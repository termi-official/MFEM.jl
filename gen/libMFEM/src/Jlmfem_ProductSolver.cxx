// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::ProductSolver> : std::false_type { };
  template<> struct DefaultConstructible<mfem::ProductSolver> : std::false_type { };
template<> struct SuperType<mfem::ProductSolver> { typedef mfem::Solver type; };
}

// Class generating the wrapper for type mfem::ProductSolver
// signature to use in the veto file: mfem::ProductSolver
struct Jlmfem_ProductSolver: public Wrapper {

  Jlmfem_ProductSolver(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::ProductSolver (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:1114:7
    jlcxx::TypeWrapper<mfem::ProductSolver>  t = jlModule.add_type<mfem::ProductSolver>("mfem!ProductSolver",
      jlcxx::julia_base_type<mfem::Solver>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::ProductSolver>>(new jlcxx::TypeWrapper<mfem::ProductSolver>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::ProductSolver::ProductSolver(mfem::Operator *, mfem::Solver *, mfem::Solver *, bool, bool, bool) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:1120:4
    t.constructor<mfem::Operator *, mfem::Solver *, mfem::Solver *, bool, bool, bool>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::ProductSolver::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductSolver::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:1123:17
    t.method("Mult", [](mfem::ProductSolver const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::ProductSolver const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::ProductSolver::MultTranspose(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductSolver::MultTranspose(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:1124:17
    t.method("MultTranspose", [](mfem::ProductSolver const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.MultTranspose(arg0, arg1); });
    t.method("MultTranspose", [](mfem::ProductSolver const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->MultTranspose(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::ProductSolver::SetOperator(const mfem::Operator &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductSolver::SetOperator(const mfem::Operator &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:1125:17
    t.method("SetOperator", [](mfem::ProductSolver& a, const mfem::Operator & arg0)->void { a.SetOperator(arg0); });
    t.method("SetOperator", [](mfem::ProductSolver* a, const mfem::Operator & arg0)->void { a->SetOperator(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::ProductSolver>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_ProductSolver(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_ProductSolver(module));
}
