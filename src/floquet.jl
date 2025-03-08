using LinearMaps
####################################################################################################
# Computation of Floquet Coefficients for periodic orbit problems

"""
	floquet = FloquetPQZ()

This composite type implements the computation of the eigenvalues of the monodromy matrix in the case of periodic orbits problems, also called the Floquet multipliers. The method, based on a Periodic Schur PeriodicSchurDecomposition (PQZ), is numerically precise for large / small Floquet exponents when the number of time sections is large.

## The arguments are as follows:
- `eigsolver::AbstractEigenSolver` solver used to compute the eigenvalues.
If `eigsolver == EigPSD()`, then the monodromy matrix is formed and its eigenvalues are computed. Otherwise, a Matrix-Free version of the monodromy is used.

!!! warning "State space"
    This solver requires the state space to be an `AbstractVector`.
"""
struct FloquetPQZ{E <: BK.AbstractEigenSolver } <: BK.AbstractFloquetSolver
	eigsolver::E
	function FloquetPQZ(eigls2::AbstractEigenSolver = EigPSD())
		return new{typeof(eigls2)}(eigls2)
	end
	FloquetPQZ(eigls::FloquetPQZ) = eigls
end
geteigenvector(eig::FloquetPQZ, vecs, n::Union{Int, Array{Int64,1}}) = geteigenvector(eig.eigsolver, vecs, n)

function (fl::FloquetPQZ)(J, nev; kwargs...)
	# this holds information on the monodromy matrix
	monodromy = MonodromyPQZ(J)
	ev = fl.eigsolver(monodromy, nev)
	vp0 = minimum(abs, ev[1])
	if (J isa BK.FloquetWrapper{ShootingProblem}) && vp0 > 1e-8
		@warn "The precision on the Floquet multipliers is $vp0. Either decrease `tolStability` in the option ContinuationPar or use a different method than `FloquetPQZ`"
	end
	return ev
end
####################################################################################################
# ShootingProblem

# Matrix free monodromy operator
# Compute the monodromy matrix at `x` explicitely, not suitable for large systems
function MonodromyPQZ(JacSH::BK.FloquetWrapper{Tpb, Tjacpb, Torbitguess, Tp}) where {Tpb <: ShootingProblem, Tjacpb, Torbitguess, Tp}
	@info "PQZ-MF"
	sh = JacSH.pb
	x = JacSH.x
	p = JacSH.par

	# period of the cycle
	T = getperiod(sh, x)

	# extract parameters
	M = BK.get_mesh_size(sh)
	N = div(length(x) - 1, M)

	# extract the time slices
	xv = @view x[1:end-1]
	xc = reshape(xv, N, M)
	TY = typeof(T)

	if M == 1
		return [LinearMaps.LinearMap{TY}(dx -> BK.evolve(sh.flow, Val(:SerialdFlow), xc[:, 1], p, dx, sh.ds[1] * T).du, N, ismutating = false)]
	else
		return [LinearMaps.LinearMap{TY}(dx -> BK.evolve(sh.flow, Val(:SerialdFlow), xc[:, ii], p, dx, sh.ds[ii] * T).du, N, ismutating = false) for ii in 1:M ]
	end
end

# Compute the monodromy matrix at `x` expression of the Jacobian of the shooting functional. We thus
# just extract the blocks needed to compute the monodromy
function MonodromyPQZ(JacSH::BK.FloquetWrapper{Tpb, Tjacpb, Torbitguess, Tp}) where {Tpb <: BK.ShootingProblem, Tjacpb <: AbstractMatrix, Torbitguess, Tp}
	J = JacSH.jacpb
	sh = JacSH.pb
	M = BK.get_mesh_size(sh)
	N = div(length(JacSH.x) - 1, M)
	As = [copy(J[1:N, 1:N])]
	if M == 1
		return [As[1] + I]
	end
	r = N
	for ii = 1:M-1
		@views push!(As, J[r+1:r+N, r+1:r+N])
		r += N
	end
	return As
end
####################################################################################################
# PoincareShooting

# matrix free evaluation of monodromy operator
function MonodromyPQZ(JacSH::BK.FloquetWrapper{Tpb, Tjacpb, Torbitguess, Tp}) where {Tpb <: BK.PoincareShootingProblem, Tjacpb, Torbitguess, Tp}
	psh = JacSH.pb
	x_bar = JacSH.x
	p = JacSH.par

	M = BK.get_mesh_size(psh)
	Nm1 = div(length(x_bar), M)

	# reshape the period orbit guess into a Matrix
	x_barc = reshape(x_bar, Nm1, M)
	TY = eltype(x_bar)

	xc = similar(x_bar, Nm1 + 1)
	As = []

	for ii in 1:M
		xci = BK.E(psh.section, view(x_barc, :, ii), ii)
		Πi = ξ -> begin
			outc = BK.dE(psh.section, ξ, ii)
			outc .= BK.diffPoincareMap(psh, xci, p, outc, ii)
			return BK.dR(psh.section, outc, ii)
		end
		push!(As, LinearMaps.LinearMap{TY}(Πi, Nm1, ismutating = false))
	end
	return As

end

# matrix based formulation of monodromy operator, not suitable for large systems
# it is based on a matrix expression of the Jacobian of the shooting functional. We thus
# just extract the blocks needed to compute the monodromy
function MonodromyPQZ(JacSH::BK.FloquetWrapper{Tpb, Tjacpb, Torbitguess, Tp}) where {Tpb <: BK.PoincareShootingProblem, Tjacpb <: AbstractMatrix, Torbitguess, Tp}
	J = JacSH.jacpb
	sh = JacSH.pb
	T = eltype(J)

	M = BK.get_mesh_size(sh)
	Nj = length(JacSH.x)
	N = div(Nj, M)

	if M == 1
		return [I - J]
	end

	r1 = mod(2N, Nj)
	r2 = N
	As = @views Matrix{T}[J[N+1:2N, 1:N]]
	for ii = 1:M-1
		# @views mul!(tmp, J[r1+1:r1+N, r2+1:r2+N], mono)
		@views push!(As, J[r1+1:r1+N, r2+1:r2+N])
		# mono .= tmp
		r1 = mod(r1 + N, Nj)
		r2 += N
	end
	return As
end

####################################################################################################
# Compute the monodromy matrix at `u0` explicitely, not suitable for large systems
function MonodromyPQZ(JacFW::BK.FloquetWrapper{Tpb, Tjacpb, Torbitguess, Tp})  where {Tpb <: BK.PeriodicOrbitTrapProblem, Tjacpb, Torbitguess, Tp}

	poPb = JacFW.pb
	u0 = JacFW.x
	par = JacFW.par

	# extraction of various constants
	M, N = size(poPb)

	# period of the cycle
	T = BK.getperiod(poPb, u0)

	# time step
	h =  T * BK.get_time_step(poPb, 1)

	u0c = BK.get_time_slices(u0, N, M)

	@views As = [(I + h/2 * BK.jacobian(poPb.prob_vf, u0c[:, 1], par))]
	@views Bs = [(I - h/2 * BK.jacobian(poPb.prob_vf, u0c[:, 1], par))]

	for ii in 2:M-1
		# also I - h/2 .* J seems to hurt (a little) the performances
		h =  T * BK.get_time_step(poPb, ii)
		@views push!(As, (I + h/2 * BK.jacobian(poPb.prob_vf, u0c[:, ii], par)))
		@views push!(Bs, (I - h/2 * BK.jacobian(poPb.prob_vf, u0c[:, ii], par)))
	end
	return (;As, Bs)
end
####################################################################################################
# orthogonal collocation

using SparseArrays

function MonodromyPQZ(JacColl::BK.FloquetWrapper{Tpb, Tjacpb, Torbitguess, Tp})  where {Tpb <: BK.PeriodicOrbitOCollProblem, Tjacpb, Torbitguess, Tp}
	prob = JacColl.pb
	𝒯 = eltype(prob)

	J = JacColl.jacpb
	n, m, Ntst = size(prob)
	nbcoll = n * m

	# condensation of parameters
	# this removes the internal unknowns of each mesh interval
	# this matrix is diagonal by blocks and each block is the L Matrix
	# which makes the corresponding J block upper triangular

	P = Matrix{𝒯}(LinearAlgebra.I(size(J, 1)))
	rg = 1:nbcoll # range
	for k = 1:Ntst
		F = lu(J[rg, rg .+ n])
		# if true
		# 	_Per = sparse(I(nbcoll))[F.p, :]
		# 	P[rg, rg] .= (_Per \ Array(F.L))
		# else
			P[rg, rg] .= (F.P \ F.L)
		# end
		rg = rg .+ m * n
	end

	Jcond = P \ J

	r1 = 1:n
	r2 = n*(m-1)+1:(m * n)

	As = Matrix{𝒯}[]
	Bs = Matrix{𝒯}[]

	for _ in 1:Ntst
		push!(As, Jcond[r2, r1])
		push!(Bs, -Jcond[r2, r1 .+ n * m])
		r1  = r1 .+ m * n
		r2  = r2 .+ m * n
	end

	return (;As, Bs)
end
