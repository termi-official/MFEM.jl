// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::ElementTransformation> : std::false_type { };
  template<> struct DefaultConstructible<mfem::ElementTransformation> : std::false_type { };
}

// Class generating the wrapper for type mfem::ElementTransformation
// signature to use in the veto file: mfem::ElementTransformation
struct Jlmfem_ElementTransformation: public Wrapper {

  Jlmfem_ElementTransformation(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::ElementTransformation (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/eltrans.hpp:23:7
    jlcxx::TypeWrapper<mfem::ElementTransformation>  t = jlModule.add_type<mfem::ElementTransformation>("mfem!ElementTransformation");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::ElementTransformation>>(new jlcxx::TypeWrapper<mfem::ElementTransformation>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::ElementTransformation>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_ElementTransformation(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_ElementTransformation(module));
}
