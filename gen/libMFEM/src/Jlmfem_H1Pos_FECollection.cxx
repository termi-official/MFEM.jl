// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::H1Pos_FECollection> : std::false_type { };
  template<> struct DefaultConstructible<mfem::H1Pos_FECollection> : std::false_type { };
template<> struct SuperType<mfem::H1Pos_FECollection> { typedef mfem::H1_FECollection type; };
}

// Class generating the wrapper for type mfem::H1Pos_FECollection
// signature to use in the veto file: mfem::H1Pos_FECollection
struct Jlmfem_H1Pos_FECollection: public Wrapper {

  Jlmfem_H1Pos_FECollection(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::H1Pos_FECollection (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:256:7
    jlcxx::TypeWrapper<mfem::H1Pos_FECollection>  t = jlModule.add_type<mfem::H1Pos_FECollection>("mfem!H1Pos_FECollection",
      jlcxx::julia_base_type<mfem::H1_FECollection>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::H1Pos_FECollection>>(new jlcxx::TypeWrapper<mfem::H1Pos_FECollection>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::H1Pos_FECollection::H1Pos_FECollection(const int, const int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:259:13
    t.constructor<const int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const int, const int>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::H1Pos_FECollection>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_H1Pos_FECollection(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_H1Pos_FECollection(module));
}
