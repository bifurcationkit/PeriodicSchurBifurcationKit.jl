using Parameters

getEigenvecs(ps) = eigvecs(ps, trues(length(ps.values)))

"""
Eigensolver for cyclic matrices based on PeriodicSchurDecompositions.jl.

!!! warning "State space"
	This solver requires the state space to be an `AbstractVector`.
"""
@with_kw struct EigPSD <: BK.AbstractEigenSolver
	computeEigenvector::Bool = false
	wantZ::Bool = false
	wantT::Bool = false
end
BK.geteigenvector(::EigPSD, ev, Is) = nothing

function (eig::EigPSD)(J::NamedTuple{(:As, :Bs), 𝒯}, nev) where 𝒯
	n = size(J.As[1], 1)
	nev = min(nev, n)

	pS = gpschur(J.As, J.Bs; wantZ = eig.wantZ, wantT = eig.wantT)
	@debug "PQZ-g"
	logvals = log.(complex.(pS.values))
	I = sortperm(logvals, by = real, rev = true)[1:nev]
	# eigvec = ~eig.computeEigenvector ? nothing : getEigenvecs(pS)[1][:,I]
	return logvals[I], nothing, true, 1
end

# this solver is specific to the case where Bs = I
function (eig::EigPSD)(J, nev)
	n = size(J[1], 1)
	nev = min(nev, n)

	pS = pschur(reverse(J); wantZ = eig.wantZ, wantT = eig.wantT)
	@debug "PQZ-p"
	logvals = log.(complex.(pS.values))
	I = sortperm(logvals, by = real, rev = true)[1:nev]
	return logvals[I], nothing, true, 1
end

"""
Eigensolver for cyclic matrices based on PeriodicSchurDecompositions.jl. This is best used with a Matrix-Free formulation.

!!! warning "State space"
	This solver requires the state space to be an `AbstractVector`.
"""
@with_kw struct EigPSD_MF{T} <: BK.AbstractMFEigenSolver
	computeEigenvector::Bool = false
	tol::T = 1e-8
	maxdim::Int = 100
	restarts::Int = 100
end
BK.geteigenvector(::EigPSD_MF, ev, Is) = nothing

function (eig::EigPSD_MF)(J::NamedTuple{(:As, :Bs), 𝒯}, nev) where 𝒯
	@assert 1==0 "Not yet implemented in PeriodicSchurDecompositions.jl"
	n = size(J.As[1], 1)
	nev = min(nev, n)

	pS = gpschur(J.As, J.Bs; wantZ = eig.wantZ, wantT = eig.wantT)
	logvals = log.(complex.(pS.values))
	I = sortperm(logvals, by = abs, rev = true)[1:nev]
	eigvec = ~eig.computeEigenvector ? nothing : getEigenvecs(pS)[1][:,I]
	return logvals[I], eigvec, true, 1
end

# this solver is specific to the case where Bs = I
function (eig::EigPSD_MF)(J, nev)
	maxdim = max(eig.maxdim, max(20, 2nev))
	@info length(J)
	pS, history = partial_pschur(J, nev, PSD.LM(); maxdim = maxdim, restarts = eig.restarts, tol = eig.tol)

	logvals = log.(complex.(pS.values))

	ncv = length(logvals)
	ncv < nev && @warn "$ncv eigenvalues have converged using PeriodicSchurDecompositions.partial_pschur, you requested $nev"

	I = sortperm(logvals, by = real, rev = true)[1:nev]
	eigvec = ~eig.computeEigenvector ? nothing : getEigenvecs(pS)[1][:,I]

	return logvals[I], eigvec, true, 1
end
