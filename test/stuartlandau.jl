# using Revise, Plots
using OrdinaryDiffEq, ForwardDiff, Test
	using BifurcationKit, LinearAlgebra, Parameters, Setfield
	const BK = BifurcationKit
	const FD = ForwardDiff

norminf(x) = norm(x, Inf)
using PeriodicSchurBifurcationKit

function Fsl!(f, u, p, t = 0)
	@unpack r, μ, ν, c3, c5 = p
	u1 = u[1]
	u2 = u[2]

	ua = u1^2 + u2^2

	f[1] = r * u1 - ν * u2 - ua * (c3 * u1 - μ * u2) - c5 * ua^2 * u1
	f[2] = r * u2 + ν * u1 - ua * (c3 * u2 + μ * u1) - c5 * ua^2 * u2
	return f
end

Fsl(x, p) = Fsl!(similar(x), x, p)
dFsl(x, dx, p) = FD.derivative(t -> Fsl(x .+ t .* dx, p), 0.)

par_sl = (r = 0.5, μ = 0., ν = 1.0, c3 = 1.0, c5 = 0.0,)
par_hopf = (@set par_sl.r = 0.1)
u0 = [.001, .001]

prob_vf = BifurcationKit.BifurcationProblem(Fsl, u0, par_hopf, (@lens _.r))

function FslMono!(f, x, p, t)
	u = x[1:2]
	du = x[3:4]
	Fsl!(f[1:2], u, p, t)
	f[3:4] .= dFsl(u, du, p)
end
####################################################################################################
# continuation
optconteq = ContinuationPar(ds = -0.01, detectBifurcation = 3, pMin = -0.5, nInversion = 4)
br = continuation(prob_vf, PALC(), optconteq)
####################################################################################################
# trapezoid
optcontpo = setproperties(optconteq; detectBifurcation = 2, tolStability = 1e-7)
@set! optcontpo.ds = -0.01
@set! optcontpo.newtonOptions.verbose = false
br_po_trap = continuation(br, 1, (@set ContinuationPar(optcontpo; ds = 0.01, saveSolEveryStep=1).newtonOptions.verbose = false),
	PeriodicOrbitTrapProblem(M = 40; jacobian = :Dense, updateSectionEveryStep = 1);
	δp = 0.1,
	usedeflation = true,
	eigsolver = FloquetPQZ(),
	)

# orthogonal collocation
optcontpo = setproperties(optconteq; detectBifurcation = 2, tolStability = 1e-7)
@set! optcontpo.ds = 0.01
@set! optcontpo.newtonOptions.verbose = false
br_po_oc = continuation(br, 1, (@set ContinuationPar(optcontpo; ds = 0.01, saveSolEveryStep=1).newtonOptions.verbose = false),
	PeriodicOrbitOCollProblem(40, 3; jacobian = :autodiffDense, updateSectionEveryStep = 1);
	δp = 0.1, verbosity = 1,
	usedeflation = true,
	eigsolver = FloquetPQZ(),
	)

@test norm(br_po_trap.eig[end].eigenvals - br_po_oc.eig[end].eigenvals, Inf) < 1e-1

####################################################################################################
prob = ODEProblem(Fsl!, u0, (0., 100.), par_hopf)
probMono = ODEProblem(FslMono!, vcat(u0, u0), (0., 100.), par_hopf)
####################################################################################################
# test automatic branch switching
br_pok2 = continuation(br, 1, optcontpo,
	ShootingProblem(1, prob, KenCarp4();  abstol = 1e-10, reltol = 1e-9, lens = (@lens _.r)); normC = norminf, verbosity = 0, eigsolver = FloquetPQZ(),)

br_pok2 = continuation(br, 1, optcontpo,
	ShootingProblem(3, prob, KenCarp4();  abstol = 1e-10, reltol = 1e-9, lens = (@lens _.r)); normC = norminf, verbosity = 0, eigsolver = FloquetPQZ(),)
####################################################################################################
# test automatic branch switching with most possible options
# calls with analytical jacobians
br_psh0 = continuation(br, 1, (@set optcontpo.ds = 0.005), PoincareShootingProblem(1, prob, KenCarp4(); abstol=1e-10, reltol=1e-9, lens = @lens _.r); normC = norminf,)

br_psh02 = continuation(br, 1, (@set optcontpo.ds = 0.005), PoincareShootingProblem(3, prob, KenCarp4(); abstol=1e-10, reltol=1e-9, lens = @lens _.r); normC = norminf,)

br_psh1 = continuation(br, 1, (@set optcontpo.ds = 0.005), PoincareShootingProblem(1, prob, KenCarp4(); abstol=1e-10, reltol=1e-9, lens = @lens _.r); normC = norminf, eigsolver = FloquetPQZ(),)

br_psh2 = continuation(br, 1, (@set optcontpo.ds = 0.005), PoincareShootingProblem(2, prob, KenCarp4(); abstol=1e-10, reltol=1e-9, lens = @lens _.r); normC = norminf, eigsolver = FloquetPQZ(),)

br_psh3 = continuation(br, 1, (@set optcontpo.ds = 0.005), PoincareShootingProblem(3, prob, KenCarp4(); abstol=1e-10, reltol=1e-9, lens = @lens _.r); normC = norminf, eigsolver = FloquetPQZ(),)

for _b in (br_psh0, br_psh02, br_psh1, br_psh2, br_psh3)
	@test real(_b.eig[1].eigenvals[1]) ≈ -0.05827903580366282 + 0.0im rtol = 1e-4
end
