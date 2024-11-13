# How to regenerate MFEM.jl

Linux only for now.

Before we start we have to build (or obtain) LLVM, Clang and LibTooling to compile (and install) [wrapit](https://github.com/grasph/wrapit/).

Make sure julia is in the PATH otherwise it won't work.

0. Make sure you have commited your git state for MFEM.jl before starting :)
1. Compile mfem (serial for now) with -fPIC
2. Set include and library paths
3. `make all`
4. Look at code diff for the parts that need to be fixed manually (duplicated template codegen, sort template stuff manually, comment signatures with *[], add some "invisible" functions due to templating)
6. Deduplicate generated_cxx
7. `make all`
8. `make run_demo`
9. Profit

# Write up issues

* Forward declarations of templated classes lead to duplicated codegen for some reason.
* Setters for constant global variables are generated
* `*[]` breaks namespaces
* How to enforce codegen for templated methods and constructors?
* global include path should be excluded from files (at least I need the option)
* option to ignore or propagate deprecations
* wrappers should be generated after enums#
* wrappers need to be sorted to generate code for templated classes correctly
* `t.apply` needs to be sorted
* basic types are ignored in `t.apply` list codegen
