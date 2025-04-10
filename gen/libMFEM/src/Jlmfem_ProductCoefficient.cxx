// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::ProductCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::ProductCoefficient> : std::false_type { };
template<> struct SuperType<mfem::ProductCoefficient> { typedef mfem::Coefficient type; };
}

// Class generating the wrapper for type mfem::ProductCoefficient
// signature to use in the veto file: mfem::ProductCoefficient
struct Jlmfem_ProductCoefficient: public Wrapper {

  Jlmfem_ProductCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::ProductCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1353:7
    jlcxx::TypeWrapper<mfem::ProductCoefficient>  t = jlModule.add_type<mfem::ProductCoefficient>("mfem!ProductCoefficient",
      jlcxx::julia_base_type<mfem::Coefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::ProductCoefficient>>(new jlcxx::TypeWrapper<mfem::ProductCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::ProductCoefficient::ProductCoefficient(double, mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1362:4
    t.constructor<double, mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::ProductCoefficient::ProductCoefficient(mfem::Coefficient &, mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1366:4
    t.constructor<mfem::Coefficient &, mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::ProductCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1370:9
    t.method("SetTime", [](mfem::ProductCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::ProductCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::ProductCoefficient::SetAConst(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductCoefficient::SetAConst(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1373:9
    t.method("SetAConst", [](mfem::ProductCoefficient& a, double arg0)->void { a.SetAConst(arg0); });
    t.method("SetAConst", [](mfem::ProductCoefficient* a, double arg0)->void { a->SetAConst(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::ProductCoefficient::GetAConst() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::ProductCoefficient::GetAConst()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1375:11
    t.method("GetAConst", [](mfem::ProductCoefficient const& a)->double { return a.GetAConst(); });
    t.method("GetAConst", [](mfem::ProductCoefficient const* a)->double { return a->GetAConst(); });

    DEBUG_MSG("Adding wrapper for void mfem::ProductCoefficient::SetACoef(mfem::Coefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductCoefficient::SetACoef(mfem::Coefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1378:9
    t.method("SetACoef", [](mfem::ProductCoefficient& a, mfem::Coefficient & arg0)->void { a.SetACoef(arg0); });
    t.method("SetACoef", [](mfem::ProductCoefficient* a, mfem::Coefficient & arg0)->void { a->SetACoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::Coefficient * mfem::ProductCoefficient::GetACoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Coefficient * mfem::ProductCoefficient::GetACoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1380:18
    t.method("GetACoef", [](mfem::ProductCoefficient const& a)->mfem::Coefficient * { return a.GetACoef(); });
    t.method("GetACoef", [](mfem::ProductCoefficient const* a)->mfem::Coefficient * { return a->GetACoef(); });

    DEBUG_MSG("Adding wrapper for void mfem::ProductCoefficient::SetBCoef(mfem::Coefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ProductCoefficient::SetBCoef(mfem::Coefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1383:9
    t.method("SetBCoef", [](mfem::ProductCoefficient& a, mfem::Coefficient & arg0)->void { a.SetBCoef(arg0); });
    t.method("SetBCoef", [](mfem::ProductCoefficient* a, mfem::Coefficient & arg0)->void { a->SetBCoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::Coefficient * mfem::ProductCoefficient::GetBCoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Coefficient * mfem::ProductCoefficient::GetBCoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1385:18
    t.method("GetBCoef", [](mfem::ProductCoefficient const& a)->mfem::Coefficient * { return a.GetBCoef(); });
    t.method("GetBCoef", [](mfem::ProductCoefficient const* a)->mfem::Coefficient * { return a->GetBCoef(); });

    DEBUG_MSG("Adding wrapper for double mfem::ProductCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::ProductCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1388:19
    t.method("Eval", [](mfem::ProductCoefficient& a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a.Eval(arg0, arg1); });
    t.method("Eval", [](mfem::ProductCoefficient* a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a->Eval(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::ProductCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_ProductCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_ProductCoefficient(module));
}
