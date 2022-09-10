module PeriodicSchurBifurcationKit
    using Parameters, LinearAlgebra
    using BifurcationKit
    using PeriodicSchurDecompositions

    const BK = BifurcationKit
    const PSD = PeriodicSchurDecompositions

    # using Infiltrator

    include("utils.jl")
    include("floquet.jl")
    include("eigsolver.jl")

    export FloquetPQZ, EigPSD, EigPSD_MF, EigPS, size, eltype
end
