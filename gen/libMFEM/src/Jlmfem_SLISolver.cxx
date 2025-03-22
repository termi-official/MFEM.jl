// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::SLISolver> : std::false_type { };
  template<> struct DefaultConstructible<mfem::SLISolver> : std::false_type { };
template<> struct SuperType<mfem::SLISolver> { typedef mfem::IterativeSolver type; };
}

// Class generating the wrapper for type mfem::SLISolver
// signature to use in the veto file: mfem::SLISolver
struct Jlmfem_SLISolver: public Wrapper {

  Jlmfem_SLISolver(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::SLISolver (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:433:7
    jlcxx::TypeWrapper<mfem::SLISolver>  t = jlModule.add_type<mfem::SLISolver>("mfem!SLISolver",
      jlcxx::julia_base_type<mfem::IterativeSolver>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::SLISolver>>(new jlcxx::TypeWrapper<mfem::SLISolver>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::SLISolver::SetOperator(const mfem::Operator &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::SLISolver::SetOperator(const mfem::Operator &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:447:17
    t.method("SetOperator", [](mfem::SLISolver& a, const mfem::Operator & arg0)->void { a.SetOperator(arg0); });
    t.method("SetOperator", [](mfem::SLISolver* a, const mfem::Operator & arg0)->void { a->SetOperator(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::SLISolver::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::SLISolver::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:450:17
    t.method("Mult", [](mfem::SLISolver const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::SLISolver const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::SLISolver>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_SLISolver(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_SLISolver(module));
}
