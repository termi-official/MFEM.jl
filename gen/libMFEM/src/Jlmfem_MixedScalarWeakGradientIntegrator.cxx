// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedScalarWeakGradientIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedScalarWeakGradientIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedScalarWeakGradientIntegrator> { typedef mfem::MixedScalarIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedScalarWeakGradientIntegrator
// signature to use in the veto file: mfem::MixedScalarWeakGradientIntegrator
struct Jlmfem_MixedScalarWeakGradientIntegrator: public Wrapper {

  Jlmfem_MixedScalarWeakGradientIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedScalarWeakGradientIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:852:7
    jlcxx::TypeWrapper<mfem::MixedScalarWeakGradientIntegrator>  t = jlModule.add_type<mfem::MixedScalarWeakGradientIntegrator>("mfem!MixedScalarWeakGradientIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedScalarWeakGradientIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedScalarWeakGradientIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MixedScalarWeakGradientIntegrator::MixedScalarWeakGradientIntegrator(mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:856:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedScalarWeakGradientIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedScalarWeakGradientIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedScalarWeakGradientIntegrator(module));
}
