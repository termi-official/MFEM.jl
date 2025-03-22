// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VisItFieldInfo> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VisItFieldInfo> : std::false_type { };
}

// Class generating the wrapper for type mfem::VisItFieldInfo
// signature to use in the veto file: mfem::VisItFieldInfo
struct Jlmfem_VisItFieldInfo: public Wrapper {

  Jlmfem_VisItFieldInfo(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VisItFieldInfo (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:393:7
    jlcxx::TypeWrapper<mfem::VisItFieldInfo>  t = jlModule.add_type<mfem::VisItFieldInfo>("mfem!VisItFieldInfo");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VisItFieldInfo>>(new jlcxx::TypeWrapper<mfem::VisItFieldInfo>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VisItFieldInfo::VisItFieldInfo(std::string, int, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:400:4
    t.constructor<std::string, int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<std::string, int, int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding association methods  to provide read access to the field association (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:396:16
    // signature to use in the veto list: mfem::VisItFieldInfo::association
    t.method("association", [](const mfem::VisItFieldInfo& a) -> const std::string& { return a.association; });
    t.method("association", [](mfem::VisItFieldInfo& a) -> std::string& { return a.association; });
    t.method("association", [](const mfem::VisItFieldInfo* a) -> const std::string& { return a->association; });
    t.method("association", [](mfem::VisItFieldInfo* a) -> std::string& { return a->association; });
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:396:16
    // signature to use in the veto list: mfem::VisItFieldInfo::association
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding association! methods to provide write access to the field association (" __HERE__ ")");
    t.method("association!", [](mfem::VisItFieldInfo& a, const std::string& val) -> std::string& { return a.association = val; });

    DEBUG_MSG("Adding association! methods to provide write access to the field association (" __HERE__ ")");
    t.method("association!", [](mfem::VisItFieldInfo* a, const std::string& val) -> std::string& { return a->association = val; });

    DEBUG_MSG("Adding num_components methods  to provide read access to the field num_components (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:397:8
    // signature to use in the veto list: mfem::VisItFieldInfo::num_components
    t.method("num_components", [](const mfem::VisItFieldInfo& a) -> int { return a.num_components; });
    t.method("num_components", [](mfem::VisItFieldInfo& a) -> int { return a.num_components; });
    t.method("num_components", [](const mfem::VisItFieldInfo* a) -> int { return a->num_components; });
    t.method("num_components", [](mfem::VisItFieldInfo* a) -> int { return a->num_components; });
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:397:8
    // signature to use in the veto list: mfem::VisItFieldInfo::num_components
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding num_components! methods to provide write access to the field num_components (" __HERE__ ")");
    t.method("num_components!", [](mfem::VisItFieldInfo& a, int val) -> int { return a.num_components = val; });

    DEBUG_MSG("Adding num_components! methods to provide write access to the field num_components (" __HERE__ ")");
    t.method("num_components!", [](mfem::VisItFieldInfo* a, int val) -> int { return a->num_components = val; });

    DEBUG_MSG("Adding lod methods  to provide read access to the field lod (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:398:8
    // signature to use in the veto list: mfem::VisItFieldInfo::lod
    t.method("lod", [](const mfem::VisItFieldInfo& a) -> int { return a.lod; });
    t.method("lod", [](mfem::VisItFieldInfo& a) -> int { return a.lod; });
    t.method("lod", [](const mfem::VisItFieldInfo* a) -> int { return a->lod; });
    t.method("lod", [](mfem::VisItFieldInfo* a) -> int { return a->lod; });
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/fem/datacollection.hpp:398:8
    // signature to use in the veto list: mfem::VisItFieldInfo::lod
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding lod! methods to provide write access to the field lod (" __HERE__ ")");
    t.method("lod!", [](mfem::VisItFieldInfo& a, int val) -> int { return a.lod = val; });

    DEBUG_MSG("Adding lod! methods to provide write access to the field lod (" __HERE__ ")");
    t.method("lod!", [](mfem::VisItFieldInfo* a, int val) -> int { return a->lod = val; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VisItFieldInfo>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VisItFieldInfo(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VisItFieldInfo(module));
}
