# jlmfem

C++ wrapper library of MFEM for [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl), used by MFEM.jl.

* `wrapit/`: WrapIt configuration (`jlmfem.wit.in`) and veto list (`veto.h`).
* `src/jlmfem.hpp`: handwritten additions wrapped together with MFEM: coefficients and
  integrators calling back into Julia, and non-template versions of MFEM member templates.
* `src/templates.hpp`: template instances to wrap.
* `gen/seq`, `gen/par`: generated bindings for the serial and parallel (MPI + hypre) variants.
* `build.jl`: builds MFEM, regenerates the bindings, builds `libjlmfem` and sets the preferences of MFEM.jl.

```sh
julia --project=jlmfem jlmfem/build.jl [--mfem-src=PATH] [--variants=seq,par] [--steps=mfem,wrap,lib,prefs] [-jN]
```

The parallel variant compiles against the MPICH ABI of MPICH_jll and the hypre of HYPRE_jll,
which are the libraries MPI.jl and HYPRE.jl load by default.

When adding headers or after an MFEM update, regenerate with `--steps=wrap,lib` and extend
`wrapit/veto.h` with the signatures (given in the comments of the generated code) that fail to compile.
