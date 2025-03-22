// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::TransformedCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::TransformedCoefficient> : std::false_type { };
template<> struct SuperType<mfem::TransformedCoefficient> { typedef mfem::Coefficient type; };
}

// Class generating the wrapper for type mfem::TransformedCoefficient
// signature to use in the veto file: mfem::TransformedCoefficient
struct Jlmfem_TransformedCoefficient: public Wrapper {

  Jlmfem_TransformedCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::TransformedCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:286:7
    jlcxx::TypeWrapper<mfem::TransformedCoefficient>  t = jlModule.add_type<mfem::TransformedCoefficient>("mfem!TransformedCoefficient",
      jlcxx::julia_base_type<mfem::Coefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::TransformedCoefficient>>(new jlcxx::TypeWrapper<mfem::TransformedCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::TransformedCoefficient::TransformedCoefficient(mfem::Coefficient *, double (*)(double)) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:295:4
    t.constructor<mfem::Coefficient *, double (*)(double)>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::TransformedCoefficient::TransformedCoefficient(mfem::Coefficient *, mfem::Coefficient *, double (*)(double, double)) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:297:4
    t.constructor<mfem::Coefficient *, mfem::Coefficient *, double (*)(double, double)>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::TransformedCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TransformedCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:302:9
    t.method("SetTime", [](mfem::TransformedCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::TransformedCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::TransformedCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::TransformedCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:305:19
    t.method("Eval", [](mfem::TransformedCoefficient& a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a.Eval(arg0, arg1); });
    t.method("Eval", [](mfem::TransformedCoefficient* a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a->Eval(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::TransformedCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_TransformedCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_TransformedCoefficient(module));
}
