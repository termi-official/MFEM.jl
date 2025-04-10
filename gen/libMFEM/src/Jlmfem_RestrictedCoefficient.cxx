// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::RestrictedCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::RestrictedCoefficient> : std::false_type { };
template<> struct SuperType<mfem::RestrictedCoefficient> { typedef mfem::Coefficient type; };
}

// Class generating the wrapper for type mfem::RestrictedCoefficient
// signature to use in the veto file: mfem::RestrictedCoefficient
struct Jlmfem_RestrictedCoefficient: public Wrapper {

  Jlmfem_RestrictedCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::RestrictedCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:410:7
    jlcxx::TypeWrapper<mfem::RestrictedCoefficient>  t = jlModule.add_type<mfem::RestrictedCoefficient>("mfem!RestrictedCoefficient",
      jlcxx::julia_base_type<mfem::Coefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::RestrictedCoefficient>>(new jlcxx::TypeWrapper<mfem::RestrictedCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::RestrictedCoefficient::RestrictedCoefficient(mfem::Coefficient &, mfem::Array<int> &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:420:4
    t.constructor<mfem::Coefficient &, mfem::Array<int> &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::RestrictedCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::RestrictedCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:424:9
    t.method("SetTime", [](mfem::RestrictedCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::RestrictedCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::RestrictedCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::RestrictedCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:427:19
    t.method("Eval", [](mfem::RestrictedCoefficient& a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a.Eval(arg0, arg1); });
    t.method("Eval", [](mfem::RestrictedCoefficient* a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a->Eval(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::RestrictedCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_RestrictedCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_RestrictedCoefficient(module));
}
