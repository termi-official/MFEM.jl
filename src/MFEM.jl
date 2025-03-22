module MFEM
    using CxxWrap
    using MFEM_seq_jll
    using Libdl

    # Check whether the wrappers have been build locally
    gendir = normpath(joinpath(@__DIR__, "../gen"))
    if isdir(joinpath(gendir, "libMFEM"))
        include(joinpath(gendir, "libMFEM", "src", "libMFEM-export.jl"))
        @wrapmodule(()->joinpath(gendir, "libMFEM", "libMFEM.$(Libdl.dlext)"))
    else
        error("Please build libMFEM first.")
    end

    function __init__()
        @initcxx
    end
end
