// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::IncompressibleNeoHookeanIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::IncompressibleNeoHookeanIntegrator> : std::false_type { };
template<> struct SuperType<mfem::IncompressibleNeoHookeanIntegrator> { typedef mfem::BlockNonlinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::IncompressibleNeoHookeanIntegrator
// signature to use in the veto file: mfem::IncompressibleNeoHookeanIntegrator
struct Jlmfem_IncompressibleNeoHookeanIntegrator: public Wrapper {

  Jlmfem_IncompressibleNeoHookeanIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::IncompressibleNeoHookeanIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/nonlininteg.hpp:339:7
    jlcxx::TypeWrapper<mfem::IncompressibleNeoHookeanIntegrator>  t = jlModule.add_type<mfem::IncompressibleNeoHookeanIntegrator>("mfem!IncompressibleNeoHookeanIntegrator",
      jlcxx::julia_base_type<mfem::BlockNonlinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::IncompressibleNeoHookeanIntegrator>>(new jlcxx::TypeWrapper<mfem::IncompressibleNeoHookeanIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::IncompressibleNeoHookeanIntegrator::IncompressibleNeoHookeanIntegrator(mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/nonlininteg.hpp:348:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::IncompressibleNeoHookeanIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_IncompressibleNeoHookeanIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_IncompressibleNeoHookeanIntegrator(module));
}
