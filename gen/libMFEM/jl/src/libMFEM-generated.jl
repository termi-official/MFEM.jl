module libMFEM


import Base.getindex
import Base.setindex!

using CxxWrap
import Libdl
@wrapmodule(()->"$(@__DIR__)/../../libMFEM." * Libdl.dlext)

function __init__()
    @initcxx
end

end #module
