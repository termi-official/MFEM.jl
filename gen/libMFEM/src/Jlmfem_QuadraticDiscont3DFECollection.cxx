// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::QuadraticDiscont3DFECollection> : std::false_type { };
  template<> struct DefaultConstructible<mfem::QuadraticDiscont3DFECollection> : std::false_type { };
template<> struct SuperType<mfem::QuadraticDiscont3DFECollection> { typedef mfem::FiniteElementCollection type; };
}

// Class generating the wrapper for type mfem::QuadraticDiscont3DFECollection
// signature to use in the veto file: mfem::QuadraticDiscont3DFECollection
struct Jlmfem_QuadraticDiscont3DFECollection: public Wrapper {

  Jlmfem_QuadraticDiscont3DFECollection(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::QuadraticDiscont3DFECollection (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:1123:7
    jlcxx::TypeWrapper<mfem::QuadraticDiscont3DFECollection>  t = jlModule.add_type<mfem::QuadraticDiscont3DFECollection>("mfem!QuadraticDiscont3DFECollection",
      jlcxx::julia_base_type<mfem::FiniteElementCollection>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::QuadraticDiscont3DFECollection>>(new jlcxx::TypeWrapper<mfem::QuadraticDiscont3DFECollection>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for const mfem::FiniteElement * mfem::QuadraticDiscont3DFECollection::FiniteElementForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::FiniteElement * mfem::QuadraticDiscont3DFECollection::FiniteElementForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:1134:4
    t.method("FiniteElementForGeometry", [](mfem::QuadraticDiscont3DFECollection const& a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a.FiniteElementForGeometry(arg0); });
    t.method("FiniteElementForGeometry", [](mfem::QuadraticDiscont3DFECollection const* a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a->FiniteElementForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for int mfem::QuadraticDiscont3DFECollection::DofForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::QuadraticDiscont3DFECollection::DofForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:1136:16
    t.method("DofForGeometry", [](mfem::QuadraticDiscont3DFECollection const& a, mfem::Geometry::Type arg0)->int { return a.DofForGeometry(arg0); });
    t.method("DofForGeometry", [](mfem::QuadraticDiscont3DFECollection const* a, mfem::Geometry::Type arg0)->int { return a->DofForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for const int * mfem::QuadraticDiscont3DFECollection::DofOrderForOrientation(mfem::Geometry::Type, int) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::QuadraticDiscont3DFECollection::DofOrderForOrientation(mfem::Geometry::Type, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:1138:23
    t.method("DofOrderForOrientation", [](mfem::QuadraticDiscont3DFECollection const& a, mfem::Geometry::Type arg0, int arg1)->const int * { return a.DofOrderForOrientation(arg0, arg1); });
    t.method("DofOrderForOrientation", [](mfem::QuadraticDiscont3DFECollection const* a, mfem::Geometry::Type arg0, int arg1)->const int * { return a->DofOrderForOrientation(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::QuadraticDiscont3DFECollection::Name() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::QuadraticDiscont3DFECollection::Name()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:1141:25
    t.method("Name", [](mfem::QuadraticDiscont3DFECollection const& a) { return (std::string)a.Name(); });
    t.method("Name", [](mfem::QuadraticDiscont3DFECollection const* a) { return (std::string)a->Name(); });

    DEBUG_MSG("Adding wrapper for int mfem::QuadraticDiscont3DFECollection::GetContType() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::QuadraticDiscont3DFECollection::GetContType()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:1142:16
    t.method("GetContType", [](mfem::QuadraticDiscont3DFECollection const& a)->int { return a.GetContType(); });
    t.method("GetContType", [](mfem::QuadraticDiscont3DFECollection const* a)->int { return a->GetContType(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::QuadraticDiscont3DFECollection>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_QuadraticDiscont3DFECollection(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_QuadraticDiscont3DFECollection(module));
}
