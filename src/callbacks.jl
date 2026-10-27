# Julia implementations of MFEM coefficients and integrators. The C++ side (jlmfem.hpp) calls
# a @cfunction trampoline with a pointer to a box holding the Julia callable, so each callable
# type gets its own specialized trampoline.

"""
    PointerVector <: AbstractVector{real_t}

Non-owning vector view of memory handed out by MFEM. Only valid during the callback it is passed to.
"""
struct PointerVector <: AbstractVector{real_t}
    ptr::Ptr{real_t}
    len::Int
end
Base.size(x::PointerVector) = (x.len,)
Base.IndexStyle(::Type{PointerVector}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(x::PointerVector, i::Int)
    @boundscheck checkbounds(x, i)
    return unsafe_load(x.ptr, i)
end
Base.@propagate_inbounds function Base.setindex!(x::PointerVector, v, i::Int)
    @boundscheck checkbounds(x, i)
    unsafe_store!(x.ptr, v, i)
    return x
end

function _coefficient(box::Base.RefValue, T::CxxRef, ip::ConstCxxRef)
    return convert(real_t, box[](T, ip))
end

function _function_coefficient(box::Base.RefValue, x::Ptr{real_t}, dim::Cint, t::real_t)
    return convert(real_t, box[](PointerVector(x, dim)))
end

function _vector_function_coefficient(box::Base.RefValue, v::Ptr{real_t}, vdim::Cint, x::Ptr{real_t}, dim::Cint, t::real_t)
    box[](PointerVector(v, vdim), PointerVector(x, dim))
    return nothing
end

function _bilinear_integrator(box::Base.RefValue, el::ConstCxxRef, T::CxxRef, elmat::CxxRef)
    box[](el, T, elmat)
    return nothing
end

function _linear_integrator(box::Base.RefValue, el::ConstCxxRef, T::CxxRef, elvect::CxxRef)
    box[](el, T, elvect)
    return nothing
end

"""
    JuliaCoefficient(f)

Scalar coefficient evaluating `f(T::ElementTransformation, ip::IntegrationPoint)`.
"""
function JuliaCoefficient(f::F) where {F}
    box = Ref(f)
    cb = @cfunction(_coefficient, real_t,
        (Ref{Base.RefValue{F}}, CxxRef{ElementTransformation}, ConstCxxRef{IntegrationPoint}))
    return keepalive!(LibMFEM.jlmfem!JuliaCoefficient(cb, pointer_from_objref(box)), box)
end

"""
    FunctionCoefficient(f)

Scalar coefficient evaluating `f(x)` at the physical point `x`.
"""
function FunctionCoefficient(f::F) where {F}
    box = Ref(f)
    cb = @cfunction(_function_coefficient, real_t, (Ref{Base.RefValue{F}}, Ptr{real_t}, Cint, real_t))
    return keepalive!(LibMFEM.jlmfem!JuliaFunctionCoefficient(cb, pointer_from_objref(box)), box)
end

"""
    VectorFunctionCoefficient(vdim, f!)

Vector coefficient of dimension `vdim` evaluating `f!(v, x)`, which writes into `v` at the physical point `x`.
"""
function VectorFunctionCoefficient(vdim::Integer, f!::F) where {F}
    box = Ref(f!)
    cb = @cfunction(_vector_function_coefficient, Cvoid,
        (Ref{Base.RefValue{F}}, Ptr{real_t}, Cint, Ptr{real_t}, Cint, real_t))
    return keepalive!(LibMFEM.jlmfem!JuliaVectorFunctionCoefficient(vdim, cb, pointer_from_objref(box)), box)
end

"""
    JuliaBilinearFormIntegrator(f!)

Bilinear form integrator whose element matrix is assembled by
`f!(el::FiniteElement, T::ElementTransformation, elmat::DenseMatrix)`.
"""
function JuliaBilinearFormIntegrator(f!::F) where {F}
    box = Ref(f!)
    cb = @cfunction(_bilinear_integrator, Cvoid,
        (Ref{Base.RefValue{F}}, ConstCxxRef{FiniteElement}, CxxRef{ElementTransformation}, CxxRef{DenseMatrix}))
    return keepalive!(LibMFEM.jlmfem!JuliaBilinearFormIntegrator(cb, pointer_from_objref(box)), box)
end

"""
    JuliaLinearFormIntegrator(f!)

Linear form integrator whose element vector is assembled by
`f!(el::FiniteElement, T::ElementTransformation, elvect::Vector)`.
"""
function JuliaLinearFormIntegrator(f!::F) where {F}
    box = Ref(f!)
    cb = @cfunction(_linear_integrator, Cvoid,
        (Ref{Base.RefValue{F}}, ConstCxxRef{FiniteElement}, CxxRef{ElementTransformation}, CxxRef{Vector}))
    return keepalive!(LibMFEM.jlmfem!JuliaLinearFormIntegrator(cb, pointer_from_objref(box)), box)
end
