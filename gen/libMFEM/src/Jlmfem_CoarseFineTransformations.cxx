// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::CoarseFineTransformations> : std::false_type { };
  template<> struct DefaultConstructible<mfem::CoarseFineTransformations> : std::false_type { };
}

// Class generating the wrapper for type mfem::CoarseFineTransformations
// signature to use in the veto file: mfem::CoarseFineTransformations
struct Jlmfem_CoarseFineTransformations: public Wrapper {

  Jlmfem_CoarseFineTransformations(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::CoarseFineTransformations (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:70:8
    jlcxx::TypeWrapper<mfem::CoarseFineTransformations>  t = jlModule.add_type<mfem::CoarseFineTransformations>("mfem!CoarseFineTransformations");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::CoarseFineTransformations>>(new jlcxx::TypeWrapper<mfem::CoarseFineTransformations>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::CoarseFineTransformations::MakeCoarseToFineTable(mfem::Table &, bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CoarseFineTransformations::MakeCoarseToFineTable(mfem::Table &, bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:82:9
    t.method("MakeCoarseToFineTable", [](mfem::CoarseFineTransformations const& a, mfem::Table & arg0)->void { a.MakeCoarseToFineTable(arg0); });
    t.method("MakeCoarseToFineTable", [](mfem::CoarseFineTransformations const& a, mfem::Table & arg0, bool arg1)->void { a.MakeCoarseToFineTable(arg0, arg1); });
    t.method("MakeCoarseToFineTable", [](mfem::CoarseFineTransformations const* a, mfem::Table & arg0)->void { a->MakeCoarseToFineTable(arg0); });
    t.method("MakeCoarseToFineTable", [](mfem::CoarseFineTransformations const* a, mfem::Table & arg0, bool arg1)->void { a->MakeCoarseToFineTable(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::CoarseFineTransformations::Clear() (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CoarseFineTransformations::Clear()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:85:9
    t.method("Clear", [](mfem::CoarseFineTransformations& a)->void { a.Clear(); });
    t.method("Clear", [](mfem::CoarseFineTransformations* a)->void { a->Clear(); });

    DEBUG_MSG("Adding wrapper for bool mfem::CoarseFineTransformations::IsInitialized() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::CoarseFineTransformations::IsInitialized()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:86:9
    t.method("IsInitialized", [](mfem::CoarseFineTransformations const& a)->bool { return a.IsInitialized(); });
    t.method("IsInitialized", [](mfem::CoarseFineTransformations const* a)->bool { return a->IsInitialized(); });

    DEBUG_MSG("Adding wrapper for long mfem::CoarseFineTransformations::MemoryUsage() (" __HERE__ ")");
    // signature to use in the veto list: long mfem::CoarseFineTransformations::MemoryUsage()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:87:9
    t.method("MemoryUsage", [](mfem::CoarseFineTransformations const& a)->long { return a.MemoryUsage(); });
    t.method("MemoryUsage", [](mfem::CoarseFineTransformations const* a)->long { return a->MemoryUsage(); });

    DEBUG_MSG("Adding embeddings methods  to provide read access to the field embeddings (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:73:21
    // signature to use in the veto list: mfem::CoarseFineTransformations::embeddings
    t.method("embeddings", [](const mfem::CoarseFineTransformations& a) -> const mfem::Array<mfem::Embedding>& { return a.embeddings; });
    t.method("embeddings", [](mfem::CoarseFineTransformations& a) -> mfem::Array<mfem::Embedding>& { return a.embeddings; });
    t.method("embeddings", [](const mfem::CoarseFineTransformations* a) -> const mfem::Array<mfem::Embedding>& { return a->embeddings; });
    t.method("embeddings", [](mfem::CoarseFineTransformations* a) -> mfem::Array<mfem::Embedding>& { return a->embeddings; });
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/ncmesh.hpp:73:21
    // signature to use in the veto list: mfem::CoarseFineTransformations::embeddings
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding embeddings! methods to provide write access to the field embeddings (" __HERE__ ")");
    t.method("embeddings!", [](mfem::CoarseFineTransformations& a, const mfem::Array<mfem::Embedding>& val) -> mfem::Array<mfem::Embedding>& { return a.embeddings = val; });

    DEBUG_MSG("Adding embeddings! methods to provide write access to the field embeddings (" __HERE__ ")");
    t.method("embeddings!", [](mfem::CoarseFineTransformations* a, const mfem::Array<mfem::Embedding>& val) -> mfem::Array<mfem::Embedding>& { return a->embeddings = val; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::CoarseFineTransformations>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_CoarseFineTransformations(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_CoarseFineTransformations(module));
}
