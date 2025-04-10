// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::ExtrudeCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::ExtrudeCoefficient> : std::false_type { };
template<> struct SuperType<mfem::ExtrudeCoefficient> { typedef mfem::Coefficient type; };
}

// Class generating the wrapper for type mfem::ExtrudeCoefficient
// signature to use in the veto file: mfem::ExtrudeCoefficient
struct Jlmfem_ExtrudeCoefficient: public Wrapper {

  Jlmfem_ExtrudeCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::ExtrudeCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/gridfunc.hpp:942:7
    jlcxx::TypeWrapper<mfem::ExtrudeCoefficient>  t = jlModule.add_type<mfem::ExtrudeCoefficient>("mfem!ExtrudeCoefficient",
      jlcxx::julia_base_type<mfem::Coefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::ExtrudeCoefficient>>(new jlcxx::TypeWrapper<mfem::ExtrudeCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::ExtrudeCoefficient::ExtrudeCoefficient(mfem::Mesh *, mfem::Coefficient &, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/gridfunc.hpp:949:4
    t.constructor<mfem::Mesh *, mfem::Coefficient &, int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for double mfem::ExtrudeCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::ExtrudeCoefficient::Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/gridfunc.hpp:951:19
    t.method("Eval", [](mfem::ExtrudeCoefficient& a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a.Eval(arg0, arg1); });
    t.method("Eval", [](mfem::ExtrudeCoefficient* a, mfem::ElementTransformation & arg0, const mfem::IntegrationPoint & arg1)->double { return a->Eval(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::ExtrudeCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_ExtrudeCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_ExtrudeCoefficient(module));
}
