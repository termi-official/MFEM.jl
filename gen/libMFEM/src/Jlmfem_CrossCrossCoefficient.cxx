// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::CrossCrossCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::CrossCrossCoefficient> : std::false_type { };
template<> struct SuperType<mfem::CrossCrossCoefficient> { typedef mfem::MatrixCoefficient type; };
}

// Class generating the wrapper for type mfem::CrossCrossCoefficient
// signature to use in the veto file: mfem::CrossCrossCoefficient
struct Jlmfem_CrossCrossCoefficient: public Wrapper {

  Jlmfem_CrossCrossCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::CrossCrossCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2007:7
    jlcxx::TypeWrapper<mfem::CrossCrossCoefficient>  t = jlModule.add_type<mfem::CrossCrossCoefficient>("mfem!CrossCrossCoefficient",
      jlcxx::julia_base_type<mfem::MatrixCoefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::CrossCrossCoefficient>>(new jlcxx::TypeWrapper<mfem::CrossCrossCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::CrossCrossCoefficient(double, mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2017:4
    t.constructor<double, mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::CrossCrossCoefficient(mfem::Coefficient &, mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2018:4
    t.constructor<mfem::Coefficient &, mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CrossCrossCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2021:9
    t.method("SetTime", [](mfem::CrossCrossCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::CrossCrossCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::SetAConst(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CrossCrossCoefficient::SetAConst(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2024:9
    t.method("SetAConst", [](mfem::CrossCrossCoefficient& a, double arg0)->void { a.SetAConst(arg0); });
    t.method("SetAConst", [](mfem::CrossCrossCoefficient* a, double arg0)->void { a->SetAConst(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::CrossCrossCoefficient::GetAConst() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::CrossCrossCoefficient::GetAConst()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2026:11
    t.method("GetAConst", [](mfem::CrossCrossCoefficient const& a)->double { return a.GetAConst(); });
    t.method("GetAConst", [](mfem::CrossCrossCoefficient const* a)->double { return a->GetAConst(); });

    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::SetACoef(mfem::Coefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CrossCrossCoefficient::SetACoef(mfem::Coefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2029:9
    t.method("SetACoef", [](mfem::CrossCrossCoefficient& a, mfem::Coefficient & arg0)->void { a.SetACoef(arg0); });
    t.method("SetACoef", [](mfem::CrossCrossCoefficient* a, mfem::Coefficient & arg0)->void { a->SetACoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::Coefficient * mfem::CrossCrossCoefficient::GetACoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Coefficient * mfem::CrossCrossCoefficient::GetACoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2031:18
    t.method("GetACoef", [](mfem::CrossCrossCoefficient const& a)->mfem::Coefficient * { return a.GetACoef(); });
    t.method("GetACoef", [](mfem::CrossCrossCoefficient const* a)->mfem::Coefficient * { return a->GetACoef(); });

    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::SetKCoef(mfem::VectorCoefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CrossCrossCoefficient::SetKCoef(mfem::VectorCoefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2034:9
    t.method("SetKCoef", [](mfem::CrossCrossCoefficient& a, mfem::VectorCoefficient & arg0)->void { a.SetKCoef(arg0); });
    t.method("SetKCoef", [](mfem::CrossCrossCoefficient* a, mfem::VectorCoefficient & arg0)->void { a->SetKCoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::VectorCoefficient * mfem::CrossCrossCoefficient::GetKCoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::VectorCoefficient * mfem::CrossCrossCoefficient::GetKCoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2036:24
    t.method("GetKCoef", [](mfem::CrossCrossCoefficient const& a)->mfem::VectorCoefficient * { return a.GetKCoef(); });
    t.method("GetKCoef", [](mfem::CrossCrossCoefficient const* a)->mfem::VectorCoefficient * { return a->GetKCoef(); });

    DEBUG_MSG("Adding wrapper for void mfem::CrossCrossCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CrossCrossCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:2039:17
    t.method("Eval", [](mfem::CrossCrossCoefficient& a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::CrossCrossCoefficient* a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a->Eval(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::CrossCrossCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_CrossCrossCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_CrossCrossCoefficient(module));
}
