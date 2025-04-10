// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::Solver> : std::false_type { };
  template<> struct DefaultConstructible<mfem::Solver> : std::false_type { };
template<> struct SuperType<mfem::Solver> { typedef mfem::Operator type; };
}

// Class generating the wrapper for type mfem::Solver
// signature to use in the veto file: mfem::Solver
struct Jlmfem_Solver: public Wrapper {

  Jlmfem_Solver(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::Solver (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:651:7
    jlcxx::TypeWrapper<mfem::Solver>  t = jlModule.add_type<mfem::Solver>("mfem!Solver",
      jlcxx::julia_base_type<mfem::Operator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::Solver>>(new jlcxx::TypeWrapper<mfem::Solver>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;



    DEBUG_MSG("Adding wrapper for void mfem::Solver::SetOperator(const mfem::Operator &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::Solver::SetOperator(const mfem::Operator &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:669:17
    t.method("SetOperator", [](mfem::Solver& a, const mfem::Operator & arg0)->void { a.SetOperator(arg0); });
    t.method("SetOperator", [](mfem::Solver* a, const mfem::Operator & arg0)->void { a->SetOperator(arg0); });

    DEBUG_MSG("Adding iterative_mode methods  to provide read access to the field iterative_mode (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:655:9
    // signature to use in the veto list: mfem::Solver::iterative_mode
    t.method("iterative_mode", [](const mfem::Solver& a) -> bool { return a.iterative_mode; });
    t.method("iterative_mode", [](mfem::Solver& a) -> bool { return a.iterative_mode; });
    t.method("iterative_mode", [](const mfem::Solver* a) -> bool { return a->iterative_mode; });
    t.method("iterative_mode", [](mfem::Solver* a) -> bool { return a->iterative_mode; });
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:655:9
    // signature to use in the veto list: mfem::Solver::iterative_mode
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding iterative_mode! methods to provide write access to the field iterative_mode (" __HERE__ ")");
    t.method("iterative_mode!", [](mfem::Solver& a, bool val) -> bool { return a.iterative_mode = val; });

    DEBUG_MSG("Adding iterative_mode! methods to provide write access to the field iterative_mode (" __HERE__ ")");
    t.method("iterative_mode!", [](mfem::Solver* a, bool val) -> bool { return a->iterative_mode = val; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::Solver>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_Solver(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_Solver(module));
}
