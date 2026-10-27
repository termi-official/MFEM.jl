"""
    unsafe_array(v::Vector)
    unsafe_array(m::DenseMatrix)

Julia array sharing memory with the MFEM object. Only valid while the object is alive and not resized.
"""
function unsafe_array(v::Union{Vector, CxxWrap.CxxWrapCore.CxxBaseRef{<:Vector}})
    return unsafe_wrap(Base.Array, GetData(v).cpp_object, Int(Size(v)))
end

function unsafe_array(m::Union{DenseMatrix, CxxWrap.CxxWrapCore.CxxBaseRef{<:DenseMatrix}})
    return unsafe_wrap(Base.Array, Data(m).cpp_object, (Int(Height(m)), Int(Width(m))))
end
