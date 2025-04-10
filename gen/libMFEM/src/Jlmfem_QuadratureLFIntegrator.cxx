// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::QuadratureLFIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::QuadratureLFIntegrator> : std::false_type { };
template<> struct SuperType<mfem::QuadratureLFIntegrator> { typedef mfem::LinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::QuadratureLFIntegrator
// signature to use in the veto file: mfem::QuadratureLFIntegrator
struct Jlmfem_QuadratureLFIntegrator: public Wrapper {

  Jlmfem_QuadratureLFIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::QuadratureLFIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:554:7
    jlcxx::TypeWrapper<mfem::QuadratureLFIntegrator>  t = jlModule.add_type<mfem::QuadratureLFIntegrator>("mfem!QuadratureLFIntegrator",
      jlcxx::julia_base_type<mfem::LinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::QuadratureLFIntegrator>>(new jlcxx::TypeWrapper<mfem::QuadratureLFIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::QuadratureLFIntegrator::QuadratureLFIntegrator(mfem::QuadratureFunctionCoefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:560:4
    t.constructor<mfem::QuadratureFunctionCoefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::QuadratureLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::QuadratureLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:572:17
    t.method("AssembleRHSElementVect", [](mfem::QuadratureLFIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a.AssembleRHSElementVect(arg0, arg1, arg2); });
    t.method("AssembleRHSElementVect", [](mfem::QuadratureLFIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a->AssembleRHSElementVect(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::QuadratureLFIntegrator::SetIntRule(const mfem::IntegrationRule *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::QuadratureLFIntegrator::SetIntRule(const mfem::IntegrationRule *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:576:17
    t.method("SetIntRule", [](mfem::QuadratureLFIntegrator& a, const mfem::IntegrationRule * arg0)->void { a.SetIntRule(arg0); });
    t.method("SetIntRule", [](mfem::QuadratureLFIntegrator* a, const mfem::IntegrationRule * arg0)->void { a->SetIntRule(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::QuadratureLFIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_QuadratureLFIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_QuadratureLFIntegrator(module));
}
