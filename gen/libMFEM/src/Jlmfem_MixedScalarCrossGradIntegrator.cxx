// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedScalarCrossGradIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedScalarCrossGradIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedScalarCrossGradIntegrator> { typedef mfem::MixedScalarVectorIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedScalarCrossGradIntegrator
// signature to use in the veto file: mfem::MixedScalarCrossGradIntegrator
struct Jlmfem_MixedScalarCrossGradIntegrator: public Wrapper {

  Jlmfem_MixedScalarCrossGradIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedScalarCrossGradIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1563:7
    jlcxx::TypeWrapper<mfem::MixedScalarCrossGradIntegrator>  t = jlModule.add_type<mfem::MixedScalarCrossGradIntegrator>("mfem!MixedScalarCrossGradIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarVectorIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedScalarCrossGradIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedScalarCrossGradIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MixedScalarCrossGradIntegrator::MixedScalarCrossGradIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1566:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for bool mfem::MixedScalarCrossGradIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::MixedScalarCrossGradIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1569:24
    t.method("VerifyFiniteElementTypes", [](mfem::MixedScalarCrossGradIntegrator const& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a.VerifyFiniteElementTypes(arg0, arg1); });
    t.method("VerifyFiniteElementTypes", [](mfem::MixedScalarCrossGradIntegrator const* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a->VerifyFiniteElementTypes(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::MixedScalarCrossGradIntegrator::FiniteElementTypeFailureMessage() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::MixedScalarCrossGradIntegrator::FiniteElementTypeFailureMessage()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1579:32
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedScalarCrossGradIntegrator const& a) { return (std::string)a.FiniteElementTypeFailureMessage(); });
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedScalarCrossGradIntegrator const* a) { return (std::string)a->FiniteElementTypeFailureMessage(); });

    DEBUG_MSG("Adding wrapper for int mfem::MixedScalarCrossGradIntegrator::GetVDim(const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::MixedScalarCrossGradIntegrator::GetVDim(const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1586:15
    t.method("GetVDim", [](mfem::MixedScalarCrossGradIntegrator& a, const mfem::FiniteElement & arg0)->int { return a.GetVDim(arg0); });
    t.method("GetVDim", [](mfem::MixedScalarCrossGradIntegrator* a, const mfem::FiniteElement & arg0)->int { return a->GetVDim(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedScalarCrossGradIntegrator::CalcVShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedScalarCrossGradIntegrator::CalcVShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1589:24
    t.method("CalcVShape", [](mfem::MixedScalarCrossGradIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.CalcVShape(arg0, arg1, arg2); });
    t.method("CalcVShape", [](mfem::MixedScalarCrossGradIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->CalcVShape(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedScalarCrossGradIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedScalarCrossGradIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedScalarCrossGradIntegrator(module));
}
