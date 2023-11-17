# How to regenerate MFEM.jl

Before we start we have to build (or obtain) LLVM, Clang and LibTooling to compile (and install) [wrapit](https://github.com/grasph/wrapit/).

Make sure julia is in the PATH otherwise it won't work.

0. Make sure you have commited your git state for MFEM.jl before starting :)
1. Compile mfem (serial for now) with -fPIC
2. Set include and library paths
3. Execute makefile
4. Replace "Array<" with "mfem::Array<"
5. Comment out erroring methods again (IntRules copy ctor, code invovled in FunctionCoefficient, NURBSExt, Mesquite)
6. Deduplicate generated_cxx
7. `make all`
8. `make run_demo`

# Write up issues

* Copy constructions for global variables are generated. Can we avoid it? We can make it compile by commenting them out. E.g. IntRules -> wrapit patch
* The generated_cxx is broken for some reason. We can fix it by deleting duplicate entries.
* generated calls involving Array are broken (namespace missing)
* How to veto classes?
* functions with keyword name collisions are generated without warning
* enum class is not correctly wrapped
* include path should be excluded from files (at least I need the option)
* option to ignore deprecations
