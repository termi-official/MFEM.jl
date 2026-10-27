# Classes taking an MPI_Comm as first constructor argument also accept an MPI.Comm.
for name in (:ParMesh, :CGSolver, :GMRESSolver, :FGMRESSolver, :BiCGSTABSolver, :MINRESSolver,
             :SLISolver, :NewtonSolver, :HyprePCG, :HypreGMRES, :HypreFGMRES)
    isdefined(@__MODULE__, name) || continue
    @eval $name(comm::MPI.Comm, args...) = $name(comm.val, args...)
end

# MPI.jl finalizes MPI at exit, before the finalizers of the remaining objects run. Their
# destructors may call MPI (e.g. hypre freeing communicators), so they are leaked instead.
let types = filter(x -> x isa DataType, [getfield(LibMFEM, name) for name in names(LibMFEM; all = true)])
    roots = filter(unique(types)) do T
        isabstracttype(T) && parentmodule(T) === LibMFEM && supertype(T) === Any && !Base.has_free_typevars(T)
    end
    for T in [roots; Array; Memory]
        @eval function CxxWrap.CxxWrapCore.delete(x::$T)
            MPI.Finalized() ? (x.cpp_object = C_NULL) : invoke(CxxWrap.CxxWrapCore.delete, Tuple{Any}, x)
            return nothing
        end
    end
end
