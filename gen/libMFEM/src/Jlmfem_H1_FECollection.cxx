// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::H1_FECollection> : std::false_type { };
  template<> struct DefaultConstructible<mfem::H1_FECollection> : std::false_type { };
template<> struct SuperType<mfem::H1_FECollection> { typedef mfem::FiniteElementCollection type; };
}

// Class generating the wrapper for type mfem::H1_FECollection
// signature to use in the veto file: mfem::H1_FECollection
struct Jlmfem_H1_FECollection: public Wrapper {

  Jlmfem_H1_FECollection(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::H1_FECollection (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:216:7
    jlcxx::TypeWrapper<mfem::H1_FECollection>  t = jlModule.add_type<mfem::H1_FECollection>("mfem!H1_FECollection",
      jlcxx::julia_base_type<mfem::FiniteElementCollection>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::H1_FECollection>>(new jlcxx::TypeWrapper<mfem::H1_FECollection>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::H1_FECollection::H1_FECollection(const int, const int, const int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:227:13
    t.constructor<const int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const int, const int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const int, const int, const int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for const mfem::FiniteElement * mfem::H1_FECollection::FiniteElementForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::FiniteElement * mfem::H1_FECollection::FiniteElementForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:230:33
    t.method("FiniteElementForGeometry", [](mfem::H1_FECollection const& a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a.FiniteElementForGeometry(arg0); });
    t.method("FiniteElementForGeometry", [](mfem::H1_FECollection const* a, mfem::Geometry::Type arg0)->const mfem::FiniteElement * { return a->FiniteElementForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for int mfem::H1_FECollection::DofForGeometry(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::H1_FECollection::DofForGeometry(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:232:16
    t.method("DofForGeometry", [](mfem::H1_FECollection const& a, mfem::Geometry::Type arg0)->int { return a.DofForGeometry(arg0); });
    t.method("DofForGeometry", [](mfem::H1_FECollection const* a, mfem::Geometry::Type arg0)->int { return a->DofForGeometry(arg0); });

    DEBUG_MSG("Adding wrapper for const int * mfem::H1_FECollection::DofOrderForOrientation(mfem::Geometry::Type, int) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::H1_FECollection::DofOrderForOrientation(mfem::Geometry::Type, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:234:23
    t.method("DofOrderForOrientation", [](mfem::H1_FECollection const& a, mfem::Geometry::Type arg0, int arg1)->const int * { return a.DofOrderForOrientation(arg0, arg1); });
    t.method("DofOrderForOrientation", [](mfem::H1_FECollection const* a, mfem::Geometry::Type arg0, int arg1)->const int * { return a->DofOrderForOrientation(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const char * mfem::H1_FECollection::Name() (" __HERE__ ")");
    // signature to use in the veto list: const char * mfem::H1_FECollection::Name()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:237:24
    t.method("Name", [](mfem::H1_FECollection const& a) { return (std::string)a.Name(); });
    t.method("Name", [](mfem::H1_FECollection const* a) { return (std::string)a->Name(); });

    DEBUG_MSG("Adding wrapper for int mfem::H1_FECollection::GetContType() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::H1_FECollection::GetContType()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:238:16
    t.method("GetContType", [](mfem::H1_FECollection const& a)->int { return a.GetContType(); });
    t.method("GetContType", [](mfem::H1_FECollection const* a)->int { return a->GetContType(); });

    DEBUG_MSG("Adding wrapper for int mfem::H1_FECollection::GetBasisType() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::H1_FECollection::GetBasisType()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:239:8
    t.method("GetBasisType", [](mfem::H1_FECollection const& a)->int { return a.GetBasisType(); });
    t.method("GetBasisType", [](mfem::H1_FECollection const* a)->int { return a->GetBasisType(); });

    DEBUG_MSG("Adding wrapper for mfem::FiniteElementCollection * mfem::H1_FECollection::GetTraceCollection() (" __HERE__ ")");
    // signature to use in the veto list: mfem::FiniteElementCollection * mfem::H1_FECollection::GetTraceCollection()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:241:29
    t.method("GetTraceCollection", [](mfem::H1_FECollection const& a)->mfem::FiniteElementCollection * { return a.GetTraceCollection(); });
    t.method("GetTraceCollection", [](mfem::H1_FECollection const* a)->mfem::FiniteElementCollection * { return a->GetTraceCollection(); });

    DEBUG_MSG("Adding wrapper for const int * mfem::H1_FECollection::GetDofMap(mfem::Geometry::Type) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::H1_FECollection::GetDofMap(mfem::Geometry::Type)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:244:15
    t.method("GetDofMap", [](mfem::H1_FECollection const& a, mfem::Geometry::Type arg0)->const int * { return a.GetDofMap(arg0); });
    t.method("GetDofMap", [](mfem::H1_FECollection const* a, mfem::Geometry::Type arg0)->const int * { return a->GetDofMap(arg0); });

    DEBUG_MSG("Adding wrapper for const int * mfem::H1_FECollection::GetDofMap(mfem::Geometry::Type, int) (" __HERE__ ")");
    // signature to use in the veto list: const int * mfem::H1_FECollection::GetDofMap(mfem::Geometry::Type, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:246:15
    t.method("GetDofMap", [](mfem::H1_FECollection const& a, mfem::Geometry::Type arg0, int arg1)->const int * { return a.GetDofMap(arg0, arg1); });
    t.method("GetDofMap", [](mfem::H1_FECollection const* a, mfem::Geometry::Type arg0, int arg1)->const int * { return a->GetDofMap(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for mfem::FiniteElementCollection * mfem::H1_FECollection::Clone(int) (" __HERE__ ")");
    // signature to use in the veto list: mfem::FiniteElementCollection * mfem::H1_FECollection::Clone(int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/fe_coll.hpp:248:29
    t.method("Clone", [](mfem::H1_FECollection const& a, int arg0)->mfem::FiniteElementCollection * { return a.Clone(arg0); });
    t.method("Clone", [](mfem::H1_FECollection const* a, int arg0)->mfem::FiniteElementCollection * { return a->Clone(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::H1_FECollection>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_H1_FECollection(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_H1_FECollection(module));
}
