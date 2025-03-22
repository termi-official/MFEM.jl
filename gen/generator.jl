using WrapIt
using CxxWrap
using MFEM_seq_jll


#---Build the wrapper library----------------------------------------------------------------------
builddir = joinpath(@__DIR__, "libMFEM", "build")
sourcedir = @__DIR__
cd(@__DIR__)
mkpath(builddir)

cxxwrap_prefix = CxxWrap.prefix_path()

#---Generate the wrapper code----------------------------------------------------------------------
wit = joinpath(@__DIR__, "MFEM.wit")
witin = joinpath(@__DIR__, "MFEM.wit.in")
open(wit, "w") do f
    for l in eachline(witin)
	    println(f, replace(l, "%MFEM_INC_DIR%" => joinpath(MFEM_seq_jll.artifact_dir, "include")))
    end
end

wrapit("MFEM.wit", "./jlMFEM-veto.h", force=true)

# run(`cmake -S ./ -B build ..
#            -DCMAKE_BUILD_TYPE=Release
#            -DCMAKE_CXX_STANDARD=20
#            -DCMAKE_PREFIX_PATH=$cxxwrap_prefix\;  $sourcedir`)
# run(`cmake --build build --config Release`)

# run(`make`)
