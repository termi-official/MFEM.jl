// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::LinearNonConf3DFECollection> : std::false_type { };
  template<> struct DefaultConstructible<mfem::LinearNonConf3DFECollection> : std::false_type { };
template<> struct SuperType<mfem::LinearNonConf3DFECollection> { typedef mfem::FiniteElementCollection type; };
}

// Class generating the wrapper for type mfem::LinearNonConf3DFECollection
// signature to use in the veto file: mfem::LinearNonConf3DFECollection
struct Jlmfem_LinearNonConf3DFECollection: public Wrapper {

  Jlmfem_LinearNonConf3DFECollection(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::LinearNonConf3DFECollection (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:793:7
    jlcxx::TypeWrapper<mfem::LinearNonConf3DFECollection>  t = jlModule.add_type<mfem::LinearNonConf3DFECollection>("mfem!LinearNonConf3DFECollection",
      jlcxx::julia_base_type<mfem::FiniteElementCollection>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::LinearNonConf3DFECollection>>(new jlcxx::TypeWrapper<mfem::LinearNonConf3DFECollection>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for const mfem::FiniteElement * mfem::LinearNonConf3DFECollection::FiniteElementForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::FiniteElement * mfem::LinearNonConf3DFECollection::FiniteElementForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:805:4
    t.method("FiniteElementForGeometry", [](mfem::LinearNonConf3DFECollection const& a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a.FiniteElementForGeometry(arg0); });
    t.method("FiniteElementForGeometry", [](mfem::LinearNonConf3DFECollection const* a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a->FiniteElementForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for int mfem::LinearNonConf3DFECollection::DofForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::LinearNonConf3DFECollection::DofForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:807:16
    t.method("DofForGeometry", [](mfem::LinearNonConf3DFECollection const& a, mfem::Geometry::Type arg0)->int { return a.DofForGeometry(arg0); });
    t.method("DofForGeometry", [](mfem::LinearNonConf3DFECollection const* a, mfem::Geometry::Type arg0)->int { return a->DofForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for const int * mfem::LinearNonConf3DFECollection::DofOrderForOrientation(mfem::Geometry::Type, int) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::LinearNonConf3DFECollection::DofOrderForOrientation(mfem::Geometry::Type, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:809:23
    t.method("DofOrderForOrientation", [](mfem::LinearNonConf3DFECollection const& a, mfem::Geometry::Type arg0, int arg1)->const int * { return a.DofOrderForOrientation(arg0, arg1); });
    t.method("DofOrderForOrientation", [](mfem::LinearNonConf3DFECollection const* a, mfem::Geometry::Type arg0, int arg1)->const int * { return a->DofOrderForOrientation(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::LinearNonConf3DFECollection::Name() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::LinearNonConf3DFECollection::Name()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:812:25
    t.method("Name", [](mfem::LinearNonConf3DFECollection const& a) { return (std::string)a.Name(); });
    t.method("Name", [](mfem::LinearNonConf3DFECollection const* a) { return (std::string)a->Name(); });

    DEBUG_MSG("Adding wrapper for int mfem::LinearNonConf3DFECollection::GetContType() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::LinearNonConf3DFECollection::GetContType()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:814:16
    t.method("GetContType", [](mfem::LinearNonConf3DFECollection const& a)->int { return a.GetContType(); });
    t.method("GetContType", [](mfem::LinearNonConf3DFECollection const* a)->int { return a->GetContType(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::LinearNonConf3DFECollection>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_LinearNonConf3DFECollection(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_LinearNonConf3DFECollection(module));
}
