// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::NURBSExtension> : std::false_type { };
  template<> struct DefaultConstructible<mfem::NURBSExtension> : std::false_type { };
}

// Class generating the wrapper for type mfem::NURBSExtension
// signature to use in the veto file: mfem::NURBSExtension
struct Jlmfem_NURBSExtension: public Wrapper {

  Jlmfem_NURBSExtension(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::NURBSExtension (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/nurbs.hpp:166:7
    jlcxx::TypeWrapper<mfem::NURBSExtension>  t = jlModule.add_type<mfem::NURBSExtension>("mfem!NURBSExtension");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::NURBSExtension>>(new jlcxx::TypeWrapper<mfem::NURBSExtension>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::NURBSExtension>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_NURBSExtension(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_NURBSExtension(module));
}
