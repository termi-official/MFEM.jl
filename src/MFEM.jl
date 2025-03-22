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
        #---Call Wrapper init--------------------------------------------------
        G4JL_init()
        #---Setup [data] environment-------------------------------------------
        GEANT4_DATA_DIR = Base.get(ENV, "GEANT4_DATA_DIR", Geant4_jll.data_dir)
        for line in readlines(joinpath(Geant4_jll.artifact_dir, "bin/geant4.sh"))
            m = match(r"export[ ]+(G4.*)=.*/(.*)$", line)
            if !isnothing(m)
                G4JL_setenv(String(m[1]), joinpath(GEANT4_DATA_DIR, m[2]))
            end
        end
    end

    include("SystemOfUnits.jl")
    include("PhysicalConstants.jl")
    
    include("G4Utils.jl")
    include("G4Scoring.jl")
    include("G4CallBacks.jl")
    include("G4ParticleGenerators.jl")
    include("G4JLInterface.jl")

    #---G4Vis--------------------------------------------------------------------------------------
    export draw, draw!, drawDistanceToOut, G4JLEventDisplay
    drawDistanceToOut() = "Not implemented"
    draw() = "Not implemented"
    draw!() = "Not implemented"
    G4JLEventDisplay(::Int64) = "Not implemented"

    #---G4Hist-------------------------------------------------------------------------------------
    export H1D, H2D
    H1D() = "Not implemented"  # Constructors
    H2D() = "Not implemented"  # Constructors
end
