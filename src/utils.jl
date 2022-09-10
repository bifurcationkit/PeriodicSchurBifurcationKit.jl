# struct to wrap a function into a LinearMap
# Tv is the type of the vector space
struct LinearMap{Tv, Tf}
    a::Tf
    N::Int
end

LinearMap(A, x::AbstractVector) = LinearMap{typeof(x), typeof(A)}(A, length(x))

Base.size(A::LinearMap) = (A.N, A.N)
Base.size(A::LinearMap, ::Int) = A.N

Base.eltype(A::LinearMap{Tv, Tf}) where {Tv, Tf} = eltype(Tv)

LinearAlgebra.mul!(out, A::LinearMap, x) = copyto!(out, A.a(x))
