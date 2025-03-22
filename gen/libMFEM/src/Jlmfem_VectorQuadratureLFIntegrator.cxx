// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorQuadratureLFIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorQuadratureLFIntegrator> : std::false_type { };
template<> struct SuperType<mfem::VectorQuadratureLFIntegrator> { typedef mfem::LinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::VectorQuadratureLFIntegrator
// signature to use in the veto file: mfem::VectorQuadratureLFIntegrator
struct Jlmfem_VectorQuadratureLFIntegrator: public Wrapper {

  Jlmfem_VectorQuadratureLFIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorQuadratureLFIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:523:7
    jlcxx::TypeWrapper<mfem::VectorQuadratureLFIntegrator>  t = jlModule.add_type<mfem::VectorQuadratureLFIntegrator>("mfem!VectorQuadratureLFIntegrator",
      jlcxx::julia_base_type<mfem::LinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorQuadratureLFIntegrator>>(new jlcxx::TypeWrapper<mfem::VectorQuadratureLFIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::VectorQuadratureLFIntegrator::VectorQuadratureLFIntegrator(mfem::VectorQuadratureFunctionCoefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:529:4
    t.constructor<mfem::VectorQuadratureFunctionCoefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::VectorQuadratureLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorQuadratureLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:541:17
    t.method("AssembleRHSElementVect", [](mfem::VectorQuadratureLFIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a.AssembleRHSElementVect(arg0, arg1, arg2); });
    t.method("AssembleRHSElementVect", [](mfem::VectorQuadratureLFIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a->AssembleRHSElementVect(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorQuadratureLFIntegrator::SetIntRule(const mfem::IntegrationRule *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorQuadratureLFIntegrator::SetIntRule(const mfem::IntegrationRule *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:545:17
    t.method("SetIntRule", [](mfem::VectorQuadratureLFIntegrator& a, const mfem::IntegrationRule * arg0)->void { a.SetIntRule(arg0); });
    t.method("SetIntRule", [](mfem::VectorQuadratureLFIntegrator* a, const mfem::IntegrationRule * arg0)->void { a->SetIntRule(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorQuadratureLFIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorQuadratureLFIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorQuadratureLFIntegrator(module));
}
