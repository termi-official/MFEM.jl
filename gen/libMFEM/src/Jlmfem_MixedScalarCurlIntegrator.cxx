// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedScalarCurlIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedScalarCurlIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedScalarCurlIntegrator> { typedef mfem::MixedScalarIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedScalarCurlIntegrator
// signature to use in the veto file: mfem::MixedScalarCurlIntegrator
struct Jlmfem_MixedScalarCurlIntegrator: public Wrapper {

  Jlmfem_MixedScalarCurlIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedScalarCurlIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:892:7
    jlcxx::TypeWrapper<mfem::MixedScalarCurlIntegrator>  t = jlModule.add_type<mfem::MixedScalarCurlIntegrator>("mfem!MixedScalarCurlIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedScalarCurlIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedScalarCurlIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MixedScalarCurlIntegrator::MixedScalarCurlIntegrator(mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:896:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedScalarCurlIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedScalarCurlIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedScalarCurlIntegrator(module));
}
