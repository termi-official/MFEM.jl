// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedGradDivIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedGradDivIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedGradDivIntegrator> { typedef mfem::MixedScalarVectorIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedGradDivIntegrator
// signature to use in the veto file: mfem::MixedGradDivIntegrator
struct Jlmfem_MixedGradDivIntegrator: public Wrapper {

  Jlmfem_MixedGradDivIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedGradDivIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1685:7
    jlcxx::TypeWrapper<mfem::MixedGradDivIntegrator>  t = jlModule.add_type<mfem::MixedGradDivIntegrator>("mfem!MixedGradDivIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarVectorIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedGradDivIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedGradDivIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MixedGradDivIntegrator::MixedGradDivIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1688:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for bool mfem::MixedGradDivIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::MixedGradDivIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1691:24
    t.method("VerifyFiniteElementTypes", [](mfem::MixedGradDivIntegrator const& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a.VerifyFiniteElementTypes(arg0, arg1); });
    t.method("VerifyFiniteElementTypes", [](mfem::MixedGradDivIntegrator const* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a->VerifyFiniteElementTypes(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::MixedGradDivIntegrator::FiniteElementTypeFailureMessage() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::MixedGradDivIntegrator::FiniteElementTypeFailureMessage()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1701:32
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedGradDivIntegrator const& a) { return (std::string)a.FiniteElementTypeFailureMessage(); });
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedGradDivIntegrator const* a) { return (std::string)a->FiniteElementTypeFailureMessage(); });

    DEBUG_MSG("Adding wrapper for int mfem::MixedGradDivIntegrator::GetVDim(const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::MixedGradDivIntegrator::GetVDim(const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1708:23
    t.method("GetVDim", [](mfem::MixedGradDivIntegrator& a, const mfem::FiniteElement & arg0)->int { return a.GetVDim(arg0); });
    t.method("GetVDim", [](mfem::MixedGradDivIntegrator* a, const mfem::FiniteElement & arg0)->int { return a->GetVDim(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedGradDivIntegrator::CalcVShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedGradDivIntegrator::CalcVShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1711:24
    t.method("CalcVShape", [](mfem::MixedGradDivIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.CalcVShape(arg0, arg1, arg2); });
    t.method("CalcVShape", [](mfem::MixedGradDivIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->CalcVShape(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedGradDivIntegrator::CalcShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedGradDivIntegrator::CalcShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1716:24
    t.method("CalcShape", [](mfem::MixedGradDivIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a.CalcShape(arg0, arg1, arg2); });
    t.method("CalcShape", [](mfem::MixedGradDivIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a->CalcShape(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedGradDivIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedGradDivIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedGradDivIntegrator(module));
}
