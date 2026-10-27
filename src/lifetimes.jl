# MFEM objects keep raw pointers to the objects they are built from (e.g. a FiniteElementSpace
# to its Mesh), while CxxWrap lets the Julia GC delete every object it allocated. The wrappers
# below accept these objects directly instead of CxxPtr and keep them alive as long as the
# object referring to them.

"""
    keepalive!(obj, deps...)

Keep `deps` alive at least as long as `obj`. Returns `obj`.
"""
function keepalive!(obj, deps...)
    isempty(deps) || append!(get!(Base.Vector{Any}, _dependencies, obj), deps)
    return obj
end
const _dependencies = WeakKeyDict{Any, Base.Vector{Any}}()

# Classify a CxxWrap argument type as pointer (:ptr) or reference (:ref) to a wrapped class.
function _object_argument(A)
    for U in Base.uniontypes(A)
        U isa UnionAll || continue
        R = Base.unwrap_unionall(U)
        R isa DataType || continue
        kind = R.name.wrapper in (CxxPtr, ConstCxxPtr) ? :ptr :
               R.name.wrapper in (CxxRef, ConstCxxRef) ? :ref : continue
        X = U.var.ub
        X isa DataType && parentmodule(X) === LibMFEM || continue
        return kind => X
    end
    return nothing
end

# For each method of `f`, define a method taking wrapped objects for the pointer and reference
# arguments, which calls the original one and passes the objects to `retain(result, objects...)`.
function _wrap_retaining(f, retain; skip = sig -> false)
    defined = Set{Any}()
    for m in methods(f)
        m.sig isa UnionAll && continue
        params = m.sig.parameters[2:end]
        any(p -> p isa Core.TypeofVararg, params) && continue
        objargs = map(_object_argument, params)
        any(!isnothing, objargs) || continue
        sig = [isnothing(o) ? p : o.second for (p, o) in zip(params, objargs)]
        (Base.has_free_typevars(Tuple{sig...}) || skip(sig) || sig in defined) && continue
        push!(defined, sig)
        vars = [Symbol(:arg, i) for i in eachindex(sig)]
        decl = [:($v::$T) for (v, T) in zip(vars, sig)]
        call = [isnothing(o) ? v : o.first === :ptr ? :(CxxPtr($v)) : :(CxxRef($v)) for (v, o) in zip(vars, objargs)]
        objs = [v for (v, o) in zip(vars, objargs) if !isnothing(o)]
        @eval (::Core.Typeof($f))($(decl...)) = $retain($f($(call...)), $(objs...))
    end
end

# Constructors: the new object keeps its arguments alive. Copy constructors are left alone.
let types = filter(x -> x isa DataType, [getfield(LibMFEM, name) for name in names(LibMFEM; all = true)])
    for T in unique(types)
        isabstracttype(T) && parentmodule(T) === LibMFEM && !Base.has_free_typevars(T) || continue
        _wrap_retaining(T, (obj, deps...) -> keepalive!(obj, deps...); skip = sig -> sig == [T])
    end
end

# Methods storing a reference to their argument: the receiver keeps it alive.
for f in (SetOperator, SetPreconditioner)
    _wrap_retaining(f, (_, this, deps...) -> (keepalive!(this, deps...); nothing))
end

# Forms delete the integrators they own, so they are told not to own them.
_disown_integrators(form::Union{BilinearForm, NonlinearForm}) = UseExternalIntegrators(form)
_disown_integrators(form::Union{LinearForm, MixedBilinearForm}) = LibMFEM.jlmfem!DisownIntegrators(form)

for f in (AddDomainIntegrator, AddBoundaryIntegrator, AddInteriorFaceIntegrator, AddBdrFaceIntegrator)
    _wrap_retaining(f, (_, form, deps...) -> (_disown_integrators(form); keepalive!(form, deps...); nothing);
                    skip = sig -> !hasmethod(_disown_integrators, Tuple{sig[1]}))
end
