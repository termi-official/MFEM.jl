// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::LinearFormIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::LinearFormIntegrator> : std::false_type { };
}

// Class generating the wrapper for type mfem::LinearFormIntegrator
// signature to use in the veto file: mfem::LinearFormIntegrator
struct Jlmfem_LinearFormIntegrator: public Wrapper {

  Jlmfem_LinearFormIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::LinearFormIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:22:7
    jlcxx::TypeWrapper<mfem::LinearFormIntegrator>  t = jlModule.add_type<mfem::LinearFormIntegrator>("mfem!LinearFormIntegrator");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::LinearFormIntegrator>>(new jlcxx::TypeWrapper<mfem::LinearFormIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;

    DEBUG_MSG("Adding wrapper for void mfem::LinearFormIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearFormIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:32:17
    t.method("AssembleRHSElementVect", [](mfem::LinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a.AssembleRHSElementVect(arg0, arg1, arg2); });
    t.method("AssembleRHSElementVect", [](mfem::LinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a->AssembleRHSElementVect(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearFormIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearFormIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:35:17
    t.method("AssembleRHSElementVect", [](mfem::LinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::FaceElementTransformations & arg1, mfem::Vector & arg2)->void { a.AssembleRHSElementVect(arg0, arg1, arg2); });
    t.method("AssembleRHSElementVect", [](mfem::LinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::FaceElementTransformations & arg1, mfem::Vector & arg2)->void { a->AssembleRHSElementVect(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearFormIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearFormIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:38:17
    t.method("AssembleRHSElementVect", [](mfem::LinearFormIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, mfem::Vector & arg3)->void { a.AssembleRHSElementVect(arg0, arg1, arg2, arg3); });
    t.method("AssembleRHSElementVect", [](mfem::LinearFormIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, mfem::Vector & arg3)->void { a->AssembleRHSElementVect(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearFormIntegrator::SetIntRule(const mfem::IntegrationRule *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearFormIntegrator::SetIntRule(const mfem::IntegrationRule *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:43:17
    t.method("SetIntRule", [](mfem::LinearFormIntegrator& a, const mfem::IntegrationRule * arg0)->void { a.SetIntRule(arg0); });
    t.method("SetIntRule", [](mfem::LinearFormIntegrator* a, const mfem::IntegrationRule * arg0)->void { a->SetIntRule(arg0); });

    DEBUG_MSG("Adding wrapper for const mfem::IntegrationRule * mfem::LinearFormIntegrator::GetIntRule() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::IntegrationRule * mfem::LinearFormIntegrator::GetIntRule()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:44:27
    t.method("GetIntRule", [](mfem::LinearFormIntegrator& a)->const mfem::IntegrationRule * { return a.GetIntRule(); });
    t.method("GetIntRule", [](mfem::LinearFormIntegrator* a)->const mfem::IntegrationRule * { return a->GetIntRule(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::LinearFormIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_LinearFormIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_LinearFormIntegrator(module));
}
