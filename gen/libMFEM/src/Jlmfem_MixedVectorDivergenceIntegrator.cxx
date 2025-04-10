// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MixedVectorDivergenceIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MixedVectorDivergenceIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MixedVectorDivergenceIntegrator> { typedef mfem::MixedScalarVectorIntegrator type; };
}

// Class generating the wrapper for type mfem::MixedVectorDivergenceIntegrator
// signature to use in the veto file: mfem::MixedVectorDivergenceIntegrator
struct Jlmfem_MixedVectorDivergenceIntegrator: public Wrapper {

  Jlmfem_MixedVectorDivergenceIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MixedVectorDivergenceIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:814:7
    jlcxx::TypeWrapper<mfem::MixedVectorDivergenceIntegrator>  t = jlModule.add_type<mfem::MixedVectorDivergenceIntegrator>("mfem!MixedVectorDivergenceIntegrator",
      jlcxx::julia_base_type<mfem::MixedScalarVectorIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedVectorDivergenceIntegrator>>(new jlcxx::TypeWrapper<mfem::MixedVectorDivergenceIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MixedVectorDivergenceIntegrator::MixedVectorDivergenceIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:817:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MixedVectorDivergenceIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MixedVectorDivergenceIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MixedVectorDivergenceIntegrator(module));
}
