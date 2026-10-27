# Builds MFEM, generates the WrapIt bindings and builds libjlmfem for the requested variants.
#
#   julia --project=jlmfem jlmfem/build.jl [--mfem-src=PATH] [--variants=seq,par] [--steps=mfem,wrap,lib,prefs] [-jN]
#
# The default number of parallel jobs is limited by the free memory.
#
# Results are installed to jlmfem/deps/install/<variant> and MFEM.jl's LocalPreferences.toml
# is pointed to the first requested variant.

using CxxWrap, WrapIt_jll, HYPRE_jll, METIS_jll, MPICH_jll, TOML

const JLMFEM_DIR = @__DIR__
const DEPS_DIR = joinpath(JLMFEM_DIR, "deps")

function parse_args(args)
    opts = Dict(
        "mfem-src" => normpath(joinpath(JLMFEM_DIR, "..", "..", "mfem")),
        "variants" => "seq,par",
        "steps" => "mfem,wrap,lib,prefs",
        # Compiling a generated wrapper file peaks at about 2 GiB, so 3 GiB are budgeted per job.
        "j" => string(clamp(Int(Sys.free_memory() ÷ (3 * 2^30)), 1, Sys.CPU_THREADS)),
    )
    for arg in args
        if startswith(arg, "-j")
            opts["j"] = arg[3:end]
        else
            key, val = split(lstrip(arg, '-'), '='; limit = 2)
            haskey(opts, key) || error("Unknown option --$key")
            opts[key] = val
        end
    end
    return opts
end

install_dir(variant) = joinpath(DEPS_DIR, "install", variant)
build_dir(variant, what) = joinpath(DEPS_DIR, "build", "$what-$variant")
isparallel(variant) = variant == "par"

lib(jll, name) = joinpath(jll.artifact_dir, "lib", name)
include_dir(jll) = joinpath(jll.artifact_dir, "include")

# Dependencies of the JLLs (e.g. libblastrampoline for HYPRE) must be found when linking.
rpath_link() = join(("-Wl,-rpath-link,$d" for d in (lib(MPICH_jll, ""), joinpath(Sys.BINDIR, Base.LIBDIR, "julia"))), " ")

# MPICH_jll's own compiler wrappers carry the build machine prefix, so we provide one.
function mpicxx_wrapper()
    path = joinpath(DEPS_DIR, "bin", "mpicxx")
    mkpath(dirname(path))
    write(path, """
    #!/bin/sh
    exec c++ -I$(include_dir(MPICH_jll)) "\$@" -L$(lib(MPICH_jll, "")) -lmpi -Wl,-rpath,$(lib(MPICH_jll, "")) $(rpath_link())
    """)
    chmod(path, 0o755)
    return path
end

function build_mfem(variant, mfem_src, njobs)
    flags = [
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_INSTALL_PREFIX=$(install_dir(variant))",
        "-DCMAKE_INSTALL_LIBDIR=lib",
        "-DBUILD_SHARED_LIBS=ON",
        "-DMFEM_ENABLE_TESTING=OFF",
        "-DMFEM_ENABLE_EXAMPLES=OFF",
        "-DMFEM_ENABLE_MINIAPPS=OFF",
        "-DMFEM_USE_METIS=ON",
        "-DMFEM_USE_METIS_5=ON",
        "-DMETIS_DIR=$(METIS_jll.artifact_dir)",
        "-DCMAKE_SHARED_LINKER_FLAGS=$(rpath_link())",
    ]
    if isparallel(variant)
        append!(flags, [
            "-DMFEM_USE_MPI=ON",
            "-DMPICXX=$(mpicxx_wrapper())",
            "-DHYPRE_DIR=$(HYPRE_jll.artifact_dir)",
        ])
    else
        push!(flags, "-DMFEM_USE_MPI=OFF")
    end
    bdir = build_dir(variant, "mfem")
    run(`cmake -S $mfem_src -B $bdir $flags`)
    run(`cmake --build $bdir -j$njobs`)
    run(`cmake --install $bdir`)
end

const SEQ_HEADERS = [
    "general/array.hpp",
    "mesh/mesh.hpp",
    "fem/fe_coll.hpp",
    "fem/fe/fe_base.hpp",
    "fem/fespace.hpp",
    "fem/intrules.hpp",
    "fem/eltrans.hpp",
    "fem/coefficient.hpp",
    "fem/gridfunc.hpp",
    "fem/bilinearform.hpp",
    "fem/bilininteg.hpp",
    "fem/linearform.hpp",
    "fem/lininteg.hpp",
    "fem/nonlinearform.hpp",
    "fem/nonlininteg.hpp",
    "linalg/vector.hpp",
    "linalg/densemat.hpp",
    "linalg/sparsemat.hpp",
    "linalg/solvers.hpp",
    "linalg/operator.hpp",
    "linalg/matrix.hpp",
    "linalg/sparsesmoothers.hpp",
]

const PAR_HEADERS = [
    "general/communication.hpp",
    "linalg/hypre.hpp",
    "mesh/pmesh.hpp",
    "fem/pfespace.hpp",
    "fem/pgridfunc.hpp",
    "fem/pbilinearform.hpp",
    "fem/plinearform.hpp",
]

function generate_bindings(variant)
    cfg = TOML.parsefile(joinpath(JLMFEM_DIR, "wrapit", "jlmfem.wit.in"))
    headers = isparallel(variant) ? [SEQ_HEADERS; PAR_HEADERS] : SEQ_HEADERS
    cfg["input"] = [joinpath.("mfem", headers); "jlmfem.hpp"]
    cfg["include_dirs"] = [joinpath(JLMFEM_DIR, "src"), joinpath(install_dir(variant), "include")]
    if isparallel(variant)
        append!(cfg["include_dirs"], [include_dir(MPICH_jll), include_dir(HYPRE_jll)])
    end
    bdir = build_dir(variant, "wrapit")
    staging = joinpath(bdir, "out")
    cfg["out_cxx_dir"] = joinpath(staging, "src")
    cfg["out_jl_dir"] = joinpath(staging, "jl")
    cfg["out_jl_subdir"] = ""
    cfg["veto_list"] = joinpath(JLMFEM_DIR, "wrapit", "veto.h")

    rm(staging; force = true, recursive = true)
    mkpath(bdir)
    wit = joinpath(bdir, "jlmfem.wit")
    open(io -> TOML.print(io, cfg), wit, "w")
    run(Cmd(`$(WrapIt_jll.wrapit()) --force $wit`; dir = bdir))

    # Only touch changed files to avoid recompiling all of them. The Julia module code generated
    # by WrapIt is not used, MFEM.jl defines the module itself. Source locations in the comments
    # are made machine independent.
    paths = [JLMFEM_DIR * "/" => "", (jll.artifact_dir => "<$(nameof(jll))>" for jll in (HYPRE_jll, MPICH_jll))...]
    out = joinpath(JLMFEM_DIR, "gen", variant, "src")
    mkpath(out)
    generated = readdir(cfg["out_cxx_dir"])
    for file in readdir(out)
        file in generated || rm(joinpath(out, file))
    end
    for file in generated
        code = replace(read(joinpath(cfg["out_cxx_dir"], file), String), paths...)
        dst = joinpath(out, file)
        isfile(dst) && read(dst, String) == code || write(dst, code)
    end
end

function build_jlmfem(variant, njobs)
    bdir = build_dir(variant, "jlmfem")
    run(`cmake -S $JLMFEM_DIR -B $bdir
         -DCMAKE_BUILD_TYPE=Release
         -DCMAKE_INSTALL_PREFIX=$(install_dir(variant))
         -DCMAKE_PREFIX_PATH=$(CxxWrap.prefix_path())
         -DCMAKE_SHARED_LINKER_FLAGS=$(rpath_link())
         -DMFEM_DIR=$(joinpath(install_dir(variant), "lib", "cmake", "mfem"))
         -DJLMFEM_VARIANT=$variant
         -DCMAKE_CXX_COMPILER=$(isparallel(variant) ? mpicxx_wrapper() : "c++")`)
    run(`cmake --build $bdir -j$njobs`)
    run(`cmake --install $bdir`)
end

function set_preferences(variant)
    file = joinpath(JLMFEM_DIR, "..", "LocalPreferences.toml")
    prefs = isfile(file) ? TOML.parsefile(file) : Dict{String, Any}()
    prefs["MFEM"] = Dict("variant" => variant, "jlmfem_prefix" => install_dir(variant))
    open(io -> TOML.print(io, prefs), file, "w")
    @info "MFEM.jl now uses the $variant variant from $(install_dir(variant))."
end

function main(args)
    opts = parse_args(args)
    variants = split(opts["variants"], ',')
    steps = split(opts["steps"], ',')
    for variant in variants
        variant in ("seq", "par") || error("Unknown variant $variant")
        "mfem" in steps && build_mfem(variant, opts["mfem-src"], opts["j"])
        "wrap" in steps && generate_bindings(variant)
        "lib" in steps && build_jlmfem(variant, opts["j"])
    end
    "prefs" in steps && set_preferences(first(variants))
end

main(ARGS)
