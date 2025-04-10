// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::Const2DFECollection> : std::false_type { };
  template<> struct DefaultConstructible<mfem::Const2DFECollection> : std::false_type { };
template<> struct SuperType<mfem::Const2DFECollection> { typedef mfem::FiniteElementCollection type; };
}

// Class generating the wrapper for type mfem::Const2DFECollection
// signature to use in the veto file: mfem::Const2DFECollection
struct Jlmfem_Const2DFECollection: public Wrapper {

  Jlmfem_Const2DFECollection(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::Const2DFECollection (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:893:7
    jlcxx::TypeWrapper<mfem::Const2DFECollection>  t = jlModule.add_type<mfem::Const2DFECollection>("mfem!Const2DFECollection",
      jlcxx::julia_base_type<mfem::FiniteElementCollection>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::Const2DFECollection>>(new jlcxx::TypeWrapper<mfem::Const2DFECollection>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for const mfem::FiniteElement * mfem::Const2DFECollection::FiniteElementForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::FiniteElement * mfem::Const2DFECollection::FiniteElementForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:902:4
    t.method("FiniteElementForGeometry", [](mfem::Const2DFECollection const& a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a.FiniteElementForGeometry(arg0); });
    t.method("FiniteElementForGeometry", [](mfem::Const2DFECollection const* a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a->FiniteElementForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for int mfem::Const2DFECollection::DofForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::Const2DFECollection::DofForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:904:16
    t.method("DofForGeometry", [](mfem::Const2DFECollection const& a, mfem::Geometry::Type arg0)->int { return a.DofForGeometry(arg0); });
    t.method("DofForGeometry", [](mfem::Const2DFECollection const* a, mfem::Geometry::Type arg0)->int { return a->DofForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for const int * mfem::Const2DFECollection::DofOrderForOrientation(mfem::Geometry::Type, int) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::Const2DFECollection::DofOrderForOrientation(mfem::Geometry::Type, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:906:23
    t.method("DofOrderForOrientation", [](mfem::Const2DFECollection const& a, mfem::Geometry::Type arg0, int arg1)->const int * { return a.DofOrderForOrientation(arg0, arg1); });
    t.method("DofOrderForOrientation", [](mfem::Const2DFECollection const* a, mfem::Geometry::Type arg0, int arg1)->const int * { return a->DofOrderForOrientation(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::Const2DFECollection::Name() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::Const2DFECollection::Name()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:909:25
    t.method("Name", [](mfem::Const2DFECollection const& a) { return (std::string)a.Name(); });
    t.method("Name", [](mfem::Const2DFECollection const* a) { return (std::string)a->Name(); });

    DEBUG_MSG("Adding wrapper for int mfem::Const2DFECollection::GetContType() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::Const2DFECollection::GetContType()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:911:16
    t.method("GetContType", [](mfem::Const2DFECollection const& a)->int { return a.GetContType(); });
    t.method("GetContType", [](mfem::Const2DFECollection const* a)->int { return a->GetContType(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::Const2DFECollection>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_Const2DFECollection(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_Const2DFECollection(module));
}
