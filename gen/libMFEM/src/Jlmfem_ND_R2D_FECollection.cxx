// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::ND_R2D_FECollection> : std::false_type { };
  template<> struct DefaultConstructible<mfem::ND_R2D_FECollection> : std::false_type { };
template<> struct SuperType<mfem::ND_R2D_FECollection> { typedef mfem::FiniteElementCollection type; };
}

// Class generating the wrapper for type mfem::ND_R2D_FECollection
// signature to use in the veto file: mfem::ND_R2D_FECollection
struct Jlmfem_ND_R2D_FECollection: public Wrapper {

  Jlmfem_ND_R2D_FECollection(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::ND_R2D_FECollection (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:514:7
    jlcxx::TypeWrapper<mfem::ND_R2D_FECollection>  t = jlModule.add_type<mfem::ND_R2D_FECollection>("mfem!ND_R2D_FECollection",
      jlcxx::julia_base_type<mfem::FiniteElementCollection>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::ND_R2D_FECollection>>(new jlcxx::TypeWrapper<mfem::ND_R2D_FECollection>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::ND_R2D_FECollection::ND_R2D_FECollection(const int, const int, const int, const int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:523:4
    t.constructor<const int, const int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const int, const int, const int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const int, const int, const int, const int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for const mfem::FiniteElement * mfem::ND_R2D_FECollection::FiniteElementForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::FiniteElement * mfem::ND_R2D_FECollection::FiniteElementForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:527:33
    t.method("FiniteElementForGeometry", [](mfem::ND_R2D_FECollection const& a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a.FiniteElementForGeometry(arg0); });
    t.method("FiniteElementForGeometry", [](mfem::ND_R2D_FECollection const* a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a->FiniteElementForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for int mfem::ND_R2D_FECollection::DofForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::ND_R2D_FECollection::DofForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:530:16
    t.method("DofForGeometry", [](mfem::ND_R2D_FECollection const& a, mfem::Geometry::Type arg0)->int { return a.DofForGeometry(arg0); });
    t.method("DofForGeometry", [](mfem::ND_R2D_FECollection const* a, mfem::Geometry::Type arg0)->int { return a->DofForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for const int * mfem::ND_R2D_FECollection::DofOrderForOrientation(mfem::Geometry::Type, int) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::ND_R2D_FECollection::DofOrderForOrientation(mfem::Geometry::Type, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:532:23
    t.method("DofOrderForOrientation", [](mfem::ND_R2D_FECollection const& a, mfem::Geometry::Type arg0, int arg1)->const int * { return a.DofOrderForOrientation(arg0, arg1); });
    t.method("DofOrderForOrientation", [](mfem::ND_R2D_FECollection const* a, mfem::Geometry::Type arg0, int arg1)->const int * { return a->DofOrderForOrientation(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::ND_R2D_FECollection::Name() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::ND_R2D_FECollection::Name()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:534:24
    t.method("Name", [](mfem::ND_R2D_FECollection const& a) { return (std::string)a.Name(); });
    t.method("Name", [](mfem::ND_R2D_FECollection const* a) { return (std::string)a->Name(); });

    DEBUG_MSG("Adding wrapper for int mfem::ND_R2D_FECollection::GetContType() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::ND_R2D_FECollection::GetContType()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:535:16
    t.method("GetContType", [](mfem::ND_R2D_FECollection const& a)->int { return a.GetContType(); });
    t.method("GetContType", [](mfem::ND_R2D_FECollection const* a)->int { return a->GetContType(); });

    DEBUG_MSG("Adding wrapper for mfem::FiniteElementCollection * mfem::ND_R2D_FECollection::GetTraceCollection() (" __HERE__ ")");
    // signature to use in the veto list: mfem::FiniteElementCollection * mfem::ND_R2D_FECollection::GetTraceCollection()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:536:29
    t.method("GetTraceCollection", [](mfem::ND_R2D_FECollection const& a)->mfem::FiniteElementCollection * { return a.GetTraceCollection(); });
    t.method("GetTraceCollection", [](mfem::ND_R2D_FECollection const* a)->mfem::FiniteElementCollection * { return a->GetTraceCollection(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::ND_R2D_FECollection>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_ND_R2D_FECollection(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_ND_R2D_FECollection(module));
}
