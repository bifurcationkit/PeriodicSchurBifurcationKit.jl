# using Revise, Plots
using Test
	using BifurcationKit, LinearAlgebra, SparseArrays, Setfield, Parameters
	const BK = BifurcationKit

using PeriodicSchurBifurcationKit

f1(u, v) = u^2 * v
norminf(x) = norm(x, Inf)

function Fbru!(f, x, p, t = 0)
	@unpack α, β, D1, D2, l = p
	n = div(length(x), 2)
	h2 = 1.0 / n^2
	c1 = D1 / l^2 / h2
	c2 = D2 / l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]

	# Dirichlet boundary conditions
	f[1]   = c1 * (α	  - 2u[1] + u[2] ) + α - (β + 1) * u[1] + f1(u[1], v[1])
	f[end] = c2 * (v[n-1] - 2v[n] + β / α)			 + β * u[n] - f1(u[n], v[n])

	f[n]   = c1 * (u[n-1] - 2u[n] +  α   ) + α - (β + 1) * u[n] + f1(u[n], v[n])
	f[n+1] = c2 * (β / α  - 2v[1] + v[2])			 + β * u[1] - f1(u[1], v[1])

	for i=2:n-1
		  f[i] = c1 * (u[i-1] - 2u[i] + u[i+1]) + α - (β + 1) * u[i] + f1(u[i], v[i])
		f[n+i] = c2 * (v[i-1] - 2v[i] + v[i+1])			  + β * u[i] - f1(u[i], v[i])
	end
	return f
end

Fbru(x, p, t = 0) = Fbru!(similar(x), x, p, t)

function Jbru_sp(x, p)
	@unpack α, β, D1, D2, l = p
	# compute the Jacobian using a sparse representation
	n = div(length(x), 2)
	h2 = 1.0 / n^2

	c1 = D1 / p.l^2 / h2
	c2 = D2 / p.l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]

	diag   = zeros(eltype(x), 2n)
	diagp1 = zeros(eltype(x), 2n-1)
	diagm1 = zeros(eltype(x), 2n-1)

	diagpn = zeros(eltype(x), n)
	diagmn = zeros(eltype(x), n)

	@. diagmn = β - 2 * u * v
	@. diagm1[1:n-1] = c1
	@. diagm1[n+1:end] = c2

	@. diag[1:n]    = -2c1 - (β + 1) + 2 * u * v
	@. diag[n+1:2n] = -2c2 - u * u

	@. diagp1[1:n-1]   = c1
	@. diagp1[n+1:end] = c2

	@. diagpn = u * u
	J = spdiagm(0 => diag, 1 => diagp1, -1 => diagm1, n => diagpn, -n => diagmn)
	return J
end

n = 100
####################################################################################################
# different parameters to define the Brusselator model and guess for the stationary solution
par_bru = (α = 2., β = 5.45, D1 = 0.008, D2 = 0.004, l = 0.3)
	sol0 = vcat(par_bru.α * ones(n), par_bru.β/par_bru.α * ones(n))
probBif = BK.BifurcationProblem(Fbru, sol0, par_bru, (@lens _.l); J = Jbru_sp, plotSolution = (x, p; kwargs...) -> (plotsol(x; label="", kwargs... )), recordFromSolution = (x, p) -> x[div(n,2)])
####################################################################################################
eigls = EigArpack(1.1, :LM)
opts_br_eq = ContinuationPar(dsmin = 0.001, dsmax = 0.02, ds = 0.005, pMax = 1.7, detectBifurcation = 3, nev = 21, plotEveryStep = 50, newtonOptions = NewtonPar(eigsolver = eigls, tol = 1e-9), nInversion = 4)

br = @time continuation(
	probBif, PALC(),
	opts_br_eq,
	normC = norminf)

#################################################################################################### Continuation of Periodic Orbit
# Standard Shooting
using DifferentialEquations, DiffEqOperators, ForwardDiff

FOde(f, x, p, t) = Fbru!(f, x, p)
Jbru(x, dx, p) = ForwardDiff.derivative(t -> Fbru(x .+ t .* dx, p), 0.)
JacOde(J, x, p, t) = copyto!(J, Jbru_sp(x, p))

u0 = sol0 .+ 0.01 .* rand(2n)
par_hopf = (@set par_bru.l = br.specialpoint[1].param + 0.01)

jac_prototype = Jbru_sp(ones(2n), @set par_bru.β = 0)
	jac_prototype.nzval .= ones(length(jac_prototype.nzval))

using SparseDiffTools, SparseArrays
_colors = matrix_colors(jac_prototype)
vf = ODEFunction(FOde; jac_prototype = jac_prototype, colorvec = _colors)
probode = ODEProblem(vf,  sol0, (0.0, 520.), par_bru) # gives 0.22s

# sol = @time solve(probode, Rodas4P(); abstol = 1e-10, reltol = 1e-8, progress = true)
####################################################################################################
# automatic branch switching with Shooting
# linear solvers
ls = GMRESIterativeSolvers(reltol = 1e-7, maxiter = 100, verbose = false)
eig = EigKrylovKit(tol= 1e-12, x₀ = rand(2n), verbose = 0, dim = 40)
# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-9,  maxIter = 25, linsolver = ls, eigsolver = eig)
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.03, ds= 0.01, pMax = 1.5, maxSteps = 2, newtonOptions = (@set optn_po.tol = 1e-7), nev = 15, tolStability = 1e-3, detectBifurcation = 2, saveSolEveryStep = 2)

Mt = 3

probsh = ShootingProblem(Mt, probode,
		# Rodas4P();
		QNDF();
		abstol = 1e-10, reltol = 1e-8, parallel = false,
		jacobian = :FiniteDifferences,
		# jacobian = :autodiffMF,
		)

br_po0 = continuation(
	br, 1,
	# arguments for continuation
	opts_po_cont,
	probsh;
	ampfactor = 1., δp = 0.0075,
	linearAlgo = MatrixFreeBLS(@set ls.N = 2+2n*Mt),
	# eigsolver = EigKrylovKit(tol= 1e-12, x₀ = rand(2n), verbose = 0, dim = 40),
	# verbosity = 3, plot = true,
	normC = norminf)

# test passing the eig solver as an argument
br_po1 = continuation(
	br, 1,
	# arguments for continuation
	opts_po_cont, probsh;
	ampfactor = 1., δp = 0.0075,
	linearAlgo = MatrixFreeBLS(@set ls.N = 2+2n*Mt),
	eigsolver = FloquetPQZ(EigPSD_MF(tol = 1e-12, maxdim = 40, computeEigenvector = true)),
	# verbosity = 3, plot = true,
	normC = norminf)

display(br_po0.eig[end].eigenvals[1:5])
display(br_po1.eig[end].eigenvals[1:5])

@test norm(br_po0.eig[end].eigenvals[1:5] - br_po1.eig[end].eigenvals[1:5], Inf) < 1e-2

# test passing the eig solver in newtonOptions
br_po2 = continuation(
	br, 1,
	# arguments for continuation
	(@set opts_po_cont.newtonOptions.eigsolver = FloquetPQZ(EigPSD_MF(tol = 1e-12, maxdim = 40, computeEigenvector = true))),
	# opts_po_cont,
	probsh;
	ampfactor = 1., δp = 0.0075,
	linearAlgo = MatrixFreeBLS(@set ls.N = 2+2n*Mt),
	finaliseSolution = (z, tau, step, contResult; k...) -> begin
		BK.haseigenvalues(contResult) && Base.display(contResult.eig[end].eigenvals)
		return true
	end,
	normC = norminf)

@test norm(br_po0.eig[end].eigenvals[1:5] - br_po2.eig[end].eigenvals[1:5], Inf) < 1e-2
####################################################################################################
# automatic branch switching from Hopf point with Poincare Shooting
# linear solver
ls = GMRESIterativeSolvers(reltol = 1e-7, maxiter = 100, verbose = false)
# newton parameters
eig = EigKrylovKit(tol= 1e-12, x₀ = rand(2n-1), verbose = 0, dim = 40)
optn_po = NewtonPar(verbose = true, tol = 1e-7,  maxIter = 25, linsolver = ls, eigsolver = eig)
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.03, ds= 0.005, pMax = 1.5, maxSteps = 2, newtonOptions = optn_po, nev = 10, tolStability = 1e-5, detectBifurcation = 2, plotEveryStep = 2)

Mt = 2
probpsh = PoincareShootingProblem(Mt, probode,
		# Rodas4P();
		QNDF();
		abstol = 1e-10, reltol = 1e-8, parallel = false, jacobian = :FiniteDifferences)

br_po0 = continuation(
	br, 1,
	# arguments for continuation
	opts_po_cont, probpsh;
	linearAlgo = MatrixFreeBLS(@set ls.N = (2n-1)*Mt+1),
	ampfactor = 1.0, δp = 0.005,
	verbosity = 1,	plot = false,
	finaliseSolution = (z, tau, step, contResult; k...) -> begin
		BK.haseigenvalues(contResult) && Base.display(contResult.eig[end].eigenvals)
		return true
	end,
	normC = norminf)

# br_po0.eig[end].eigenvals |> display

br_po1 = continuation(
	br, 1,
	# arguments for continuation
	opts_po_cont, probpsh;
	linearAlgo = MatrixFreeBLS(@set ls.N = (2n-1)*Mt+1),
	ampfactor = 1.0, δp = 0.005,
	verbosity = 1,	plot = false,
	eigsolver = FloquetPQZ(EigPSD_MF(tol = 1e-12, maxdim = 40)),
	finaliseSolution = (z, tau, step, contResult; k...) -> begin
		BK.haseigenvalues(contResult) && Base.display(contResult.eig[end].eigenvals)
		return true
	end,
	callbackN = BK.cbMaxNorm(10),
	normC = norminf)

# br_po1.eig[end].eigenvals |> display

@test norm(br_po0.eig[end].eigenvals[1:5] - br_po1.eig[end].eigenvals[1:5], Inf) < 1e-1
