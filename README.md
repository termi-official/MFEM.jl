# MFEM.jl

Julia interface to the finite element library [MFEM](https://mfem.org).

The repository contains two parts:

* `MFEM.jl` (this package): the Julia front-end users load.
* [`jlmfem/`](jlmfem/README.md): the C++ wrapper library `libjlmfem`, generated with
  [WrapIt](https://github.com/grasph/wrapit) and built against MFEM, CxxWrap, and for the
  parallel variant MPICH_jll and HYPRE_jll, so that it shares MPI and hypre with MPI.jl and HYPRE.jl.

## Setup

There is no `jlmfem_jll` yet, so the wrapper is built locally (Linux only for now):

```sh
julia --project=jlmfem -e 'using Pkg; Pkg.instantiate()'
julia --project=jlmfem jlmfem/build.jl --mfem-src=../mfem --variants=seq,par
```

This builds MFEM and `libjlmfem` for both variants into `jlmfem/deps/install` and writes
`LocalPreferences.toml`, selecting the first listed variant. Switch with e.g.
`--variants=par --steps=prefs`. The number of compile jobs defaults to what fits into the free
memory (about 3 GiB per job); override it with `-jN`.

## Usage

MFEM classes and functions keep their C++ names, e.g. `mfem::Mesh` is `MFEM.Mesh` and
`mesh.UniformRefinement()` is `MFEM.UniformRefinement(mesh)`. Objects passed to constructors,
`Add*Integrator`, `SetOperator` and `SetPreconditioner` are kept alive by the receiving object,
so there is no need to keep Julia references to them or to wrap them in `CxxPtr`. Other
methods taking a pointer need `CxxPtr(x)` where C++ uses `&x`, e.g. `MFEM.SetIntPoint(T, CxxPtr(ip))`.

See [`examples/ex0.jl`](examples/ex0.jl) and [`examples/ex0p.jl`](examples/ex0p.jl):

```sh
julia --project examples/ex0.jl
mpiexecjl -n 4 julia --project examples/ex0p.jl   # parallel variant, see MPI.install_mpiexecjl()
```

### Julia coefficients and integrators

* `FunctionCoefficient(f)`: scalar coefficient `f(x)` of the physical point `x`.
* `VectorFunctionCoefficient(vdim, f!)`: vector coefficient filled by `f!(v, x)`.
* `JuliaCoefficient(f)`: scalar coefficient `f(T, ip)` with full access to the element transformation.
* `JuliaBilinearFormIntegrator(f!)`: element matrices assembled by `f!(el, T, elmat)`.
* `JuliaLinearFormIntegrator(f!)`: element vectors assembled by `f!(el, T, elvect)`.

`unsafe_array(v)` gives a Julia array view of an `MFEM.Vector` or `MFEM.DenseMatrix`.
See `test/runtests.jl` for a mass integrator written in Julia.
