// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::IntegrationRules> : std::false_type { };
  template<> struct DefaultConstructible<mfem::IntegrationRules> : std::false_type { };
}

// Class generating the wrapper for type mfem::IntegrationRules
// signature to use in the veto file: mfem::IntegrationRules
struct Jlmfem_IntegrationRules: public Wrapper {

  Jlmfem_IntegrationRules(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::IntegrationRules (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/intrules.hpp:311:7
    jlcxx::TypeWrapper<mfem::IntegrationRules>  t = jlModule.add_type<mfem::IntegrationRules>("mfem!IntegrationRules");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::IntegrationRules>>(new jlcxx::TypeWrapper<mfem::IntegrationRules>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::IntegrationRules::IntegrationRules(int, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/intrules.hpp:364:13
    t.constructor<int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for const mfem::IntegrationRule & mfem::IntegrationRules::Get(int, int) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::IntegrationRule & mfem::IntegrationRules::Get(int, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/intrules.hpp:368:27
    t.method("Get", [](mfem::IntegrationRules& a, int arg0, int arg1)->const mfem::IntegrationRule & { return a.Get(arg0, arg1); });
    t.method("Get", [](mfem::IntegrationRules* a, int arg0, int arg1)->const mfem::IntegrationRule & { return a->Get(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::IntegrationRules::Set(int, int, mfem::IntegrationRule &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::IntegrationRules::Set(int, int, mfem::IntegrationRule &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/intrules.hpp:370:9
    t.method("Set", [](mfem::IntegrationRules& a, int arg0, int arg1, mfem::IntegrationRule & arg2)->void { a.Set(arg0, arg1, arg2); });
    t.method("Set", [](mfem::IntegrationRules* a, int arg0, int arg1, mfem::IntegrationRule & arg2)->void { a->Set(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::IntegrationRules::SetOwnRules(int) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::IntegrationRules::SetOwnRules(int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/intrules.hpp:372:9
    t.method("SetOwnRules", [](mfem::IntegrationRules& a, int arg0)->void { a.SetOwnRules(arg0); });
    t.method("SetOwnRules", [](mfem::IntegrationRules* a, int arg0)->void { a->SetOwnRules(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::IntegrationRules>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_IntegrationRules(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_IntegrationRules(module));
}
