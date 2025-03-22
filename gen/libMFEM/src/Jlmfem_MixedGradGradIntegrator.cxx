// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedGradGradIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedGradGradIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedGradGradIntegrator> { typedef mfem::MixedVectorIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedGradGradIntegrator
// signature to use in the veto file: mfem::MixedGradGradIntegrator
struct Jlmfem_MixedGradGradIntegrator: public Wrapper {

  Jlmfem_MixedGradGradIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedGradGradIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1105:7
    jlcxx::TypeWrapper<mfem::MixedGradGradIntegrator>  t = jlModule.add_type<mfem::MixedGradGradIntegrator>("mfem!MixedGradGradIntegrator",
      jlcxx::julia_base_type<mfem::MixedVectorIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedGradGradIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedGradGradIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MixedGradGradIntegrator::MixedGradGradIntegrator(mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1109:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MixedGradGradIntegrator::MixedGradGradIntegrator(mfem::DiagonalMatrixCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1111:4
    t.constructor<mfem::DiagonalMatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MixedGradGradIntegrator::MixedGradGradIntegrator(mfem::MatrixCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1113:4
    t.constructor<mfem::MatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for bool mfem::MixedGradGradIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::MixedGradGradIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1116:24
    t.method("VerifyFiniteElementTypes", [](mfem::MixedGradGradIntegrator const& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a.VerifyFiniteElementTypes(arg0, arg1); });
    t.method("VerifyFiniteElementTypes", [](mfem::MixedGradGradIntegrator const* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a->VerifyFiniteElementTypes(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::MixedGradGradIntegrator::FiniteElementTypeFailureMessage() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::MixedGradGradIntegrator::FiniteElementTypeFailureMessage()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1126:32
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedGradGradIntegrator const& a) { return (std::string)a.FiniteElementTypeFailureMessage(); });
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedGradGradIntegrator const* a) { return (std::string)a->FiniteElementTypeFailureMessage(); });

    DEBUG_MSG("Adding wrapper for int mfem::MixedGradGradIntegrator::GetIntegrationOrder(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::MixedGradGradIntegrator::GetIntegrationOrder(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1133:23
    t.method("GetIntegrationOrder", [](mfem::MixedGradGradIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2)->int { return a.GetIntegrationOrder(arg0, arg1, arg2); });
    t.method("GetIntegrationOrder", [](mfem::MixedGradGradIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2)->int { return a->GetIntegrationOrder(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for int mfem::MixedGradGradIntegrator::GetTrialVDim(const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::MixedGradGradIntegrator::GetTrialVDim(const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1143:23
    t.method("GetTrialVDim", [](mfem::MixedGradGradIntegrator& a, const mfem::FiniteElement & arg0)->int { return a.GetTrialVDim(arg0); });
    t.method("GetTrialVDim", [](mfem::MixedGradGradIntegrator* a, const mfem::FiniteElement & arg0)->int { return a->GetTrialVDim(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedGradGradIntegrator::CalcTrialShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedGradGradIntegrator::CalcTrialShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1146:24
    t.method("CalcTrialShape", [](mfem::MixedGradGradIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.CalcTrialShape(arg0, arg1, arg2); });
    t.method("CalcTrialShape", [](mfem::MixedGradGradIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->CalcTrialShape(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for int mfem::MixedGradGradIntegrator::GetTestVDim(const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::MixedGradGradIntegrator::GetTestVDim(const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1151:23
    t.method("GetTestVDim", [](mfem::MixedGradGradIntegrator& a, const mfem::FiniteElement & arg0)->int { return a.GetTestVDim(arg0); });
    t.method("GetTestVDim", [](mfem::MixedGradGradIntegrator* a, const mfem::FiniteElement & arg0)->int { return a->GetTestVDim(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedGradGradIntegrator::CalcTestShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedGradGradIntegrator::CalcTestShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1154:24
    t.method("CalcTestShape", [](mfem::MixedGradGradIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.CalcTestShape(arg0, arg1, arg2); });
    t.method("CalcTestShape", [](mfem::MixedGradGradIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->CalcTestShape(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedGradGradIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedGradGradIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedGradGradIntegrator(module));
}
