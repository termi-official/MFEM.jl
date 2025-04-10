// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedDivGradIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedDivGradIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedDivGradIntegrator> { typedef mfem::MixedScalarVectorIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedDivGradIntegrator
// signature to use in the veto file: mfem::MixedDivGradIntegrator
struct Jlmfem_MixedDivGradIntegrator: public Wrapper {

  Jlmfem_MixedDivGradIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedDivGradIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1724:7
    jlcxx::TypeWrapper<mfem::MixedDivGradIntegrator>  t = jlModule.add_type<mfem::MixedDivGradIntegrator>("mfem!MixedDivGradIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarVectorIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedDivGradIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedDivGradIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MixedDivGradIntegrator::MixedDivGradIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1727:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for bool mfem::MixedDivGradIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::MixedDivGradIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement &, const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1730:24
    t.method("VerifyFiniteElementTypes", [](mfem::MixedDivGradIntegrator const& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a.VerifyFiniteElementTypes(arg0, arg1); });
    t.method("VerifyFiniteElementTypes", [](mfem::MixedDivGradIntegrator const* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1)->bool { return a->VerifyFiniteElementTypes(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::MixedDivGradIntegrator::FiniteElementTypeFailureMessage() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::MixedDivGradIntegrator::FiniteElementTypeFailureMessage()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1741:32
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedDivGradIntegrator const& a) { return (std::string)a.FiniteElementTypeFailureMessage(); });
    t.method("FiniteElementTypeFailureMessage", [](mfem::MixedDivGradIntegrator const* a) { return (std::string)a->FiniteElementTypeFailureMessage(); });

    DEBUG_MSG("Adding wrapper for int mfem::MixedDivGradIntegrator::GetVDim(const mfem::FiniteElement &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::MixedDivGradIntegrator::GetVDim(const mfem::FiniteElement &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1748:23
    t.method("GetVDim", [](mfem::MixedDivGradIntegrator& a, const mfem::FiniteElement & arg0)->int { return a.GetVDim(arg0); });
    t.method("GetVDim", [](mfem::MixedDivGradIntegrator* a, const mfem::FiniteElement & arg0)->int { return a->GetVDim(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedDivGradIntegrator::CalcVShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedDivGradIntegrator::CalcVShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1751:24
    t.method("CalcVShape", [](mfem::MixedDivGradIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.CalcVShape(arg0, arg1, arg2); });
    t.method("CalcVShape", [](mfem::MixedDivGradIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->CalcVShape(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::MixedDivGradIntegrator::CalcShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MixedDivGradIntegrator::CalcShape(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:1756:24
    t.method("CalcShape", [](mfem::MixedDivGradIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a.CalcShape(arg0, arg1, arg2); });
    t.method("CalcShape", [](mfem::MixedDivGradIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a->CalcShape(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedDivGradIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedDivGradIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedDivGradIntegrator(module));
}
