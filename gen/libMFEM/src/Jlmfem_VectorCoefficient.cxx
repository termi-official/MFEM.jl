// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorCoefficient> : std::false_type { };
}

// Class generating the wrapper for type mfem::VectorCoefficient
// signature to use in the veto file: mfem::VectorCoefficient
struct Jlmfem_VectorCoefficient: public Wrapper {

  Jlmfem_VectorCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:432:7
    jlcxx::TypeWrapper<mfem::VectorCoefficient>  t = jlModule.add_type<mfem::VectorCoefficient>("mfem!VectorCoefficient");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorCoefficient>>(new jlcxx::TypeWrapper<mfem::VectorCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::VectorCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:443:17
    t.method("SetTime", [](mfem::VectorCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::VectorCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::VectorCoefficient::GetTime() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::VectorCoefficient::GetTime()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:446:11
    t.method("GetTime", [](mfem::VectorCoefficient& a)->double { return a.GetTime(); });
    t.method("GetTime", [](mfem::VectorCoefficient* a)->double { return a->GetTime(); });

    DEBUG_MSG("Adding wrapper for int mfem::VectorCoefficient::GetVDim() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::VectorCoefficient::GetVDim()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:449:8
    t.method("GetVDim", [](mfem::VectorCoefficient& a)->int { return a.GetVDim(); });
    t.method("GetVDim", [](mfem::VectorCoefficient* a)->int { return a->GetVDim(); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorCoefficient::Eval(mfem::Vector &, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorCoefficient::Eval(mfem::Vector &, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:456:17
    t.method("Eval", [](mfem::VectorCoefficient& a, mfem::Vector & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::VectorCoefficient* a, mfem::Vector & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a->Eval(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationRule &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationRule &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:471:17
    t.method("Eval", [](mfem::VectorCoefficient& a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationRule & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::VectorCoefficient* a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationRule & arg2)->void { a->Eval(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorCoefficient(module));
}
