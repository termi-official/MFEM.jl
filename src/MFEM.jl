module MFEM

using CxxWrap
using Libdl: dlext
using Preferences: @load_preference
import MPI, MPIPreferences
import HYPRE_jll, METIS_jll

const variant = @load_preference("variant", "seq")
const jlmfem_prefix = @load_preference("jlmfem_prefix", "")
const libjlmfem = joinpath(jlmfem_prefix, "lib", "libjlmfem.$dlext")

isfile(libjlmfem) || error("""
    libjlmfem not found at "$libjlmfem".
    Build it with `julia --project=jlmfem jlmfem/build.jl`, which also sets the preferences of MFEM.jl.""")

"""
    MFEM.isparallel()

Whether the loaded wrapper was built with MPI support.
"""
isparallel() = variant == "par"

module LibMFEM
    import Base.getindex
    import Base.setindex!
    using CxxWrap
    import ..libjlmfem
    @wrapmodule(() -> libjlmfem)
    __init__() = @initcxx
end

# Expose mfem::X and jlmfem::X as MFEM.X, unless X clashes with another name.
let wrapped = names(LibMFEM; all = true), own = names(@__MODULE__; all = true)
    for name in wrapped
        str = string(name)
        (startswith(str, '#') || startswith(str, "__") || name in (:LibMFEM, :eval, :include)) && continue
        alias = Symbol(replace(str, r"^(jl)?mfem!" => ""))
        (alias != name && alias in wrapped || alias in own) && (alias = name)
        @eval const $alias = LibMFEM.$name
    end
end

const real_t = Float64

include("lifetimes.jl")
include("callbacks.jl")
include("arrays.jl")
isparallel() && include("parallel.jl")

function __init__()
    if isparallel() && MPIPreferences.abi != "MPICH"
        error("The parallel MFEM.jl wrapper is built against the MPICH ABI, but MPI.jl is configured for $(MPIPreferences.abi).")
    end
end

end
