// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedVectorProductIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedVectorProductIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedVectorProductIntegrator> { typedef mfem::MixedScalarVectorIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedVectorProductIntegrator
// signature to use in the veto file: mfem::MixedVectorProductIntegrator
struct Jlmfem_MixedVectorProductIntegrator: public Wrapper {

  Jlmfem_MixedVectorProductIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedVectorProductIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:695:7
    jlcxx::TypeWrapper<mfem::MixedVectorProductIntegrator>  t = jlModule.add_type<mfem::MixedVectorProductIntegrator>("mfem!MixedVectorProductIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarVectorIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedVectorProductIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedVectorProductIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MixedVectorProductIntegrator::MixedVectorProductIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:698:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedVectorProductIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedVectorProductIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedVectorProductIntegrator(module));
}
