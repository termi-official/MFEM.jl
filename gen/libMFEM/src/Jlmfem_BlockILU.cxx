// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::BlockILU> : std::false_type { };
  template<> struct DefaultConstructible<mfem::BlockILU> : std::false_type { };
template<> struct SuperType<mfem::BlockILU> { typedef mfem::Solver type; };
}

// Class generating the wrapper for type mfem::BlockILU
// signature to use in the veto file: mfem::BlockILU
struct Jlmfem_BlockILU: public Wrapper {

  Jlmfem_BlockILU(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::BlockILU (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:918:7
    jlcxx::TypeWrapper<mfem::BlockILU>  t = jlModule.add_type<mfem::BlockILU>("mfem!BlockILU",
      jlcxx::julia_base_type<mfem::Solver>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::BlockILU>>(new jlcxx::TypeWrapper<mfem::BlockILU>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::BlockILU::BlockILU(int, mfem::BlockILU::Reordering, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:932:4
    t.constructor<int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, mfem::BlockILU::Reordering>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, mfem::BlockILU::Reordering, int>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::BlockILU::BlockILU(const mfem::Operator &, int, mfem::BlockILU::Reordering, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:941:4
    t.constructor<const mfem::Operator &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const mfem::Operator &, int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const mfem::Operator &, int, mfem::BlockILU::Reordering>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const mfem::Operator &, int, mfem::BlockILU::Reordering, int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::BlockILU::SetOperator(const mfem::Operator &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BlockILU::SetOperator(const mfem::Operator &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:949:9
    t.method("SetOperator", [](mfem::BlockILU& a, const mfem::Operator & arg0)->void { a.SetOperator(arg0); });
    t.method("SetOperator", [](mfem::BlockILU* a, const mfem::Operator & arg0)->void { a->SetOperator(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BlockILU::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BlockILU::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:952:9
    t.method("Mult", [](mfem::BlockILU const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::BlockILU const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for int * mfem::BlockILU::GetBlockI() (" __HERE__ ")");
    // signature to use in the veto list: int * mfem::BlockILU::GetBlockI()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:957:9
    t.method("GetBlockI", [](mfem::BlockILU& a)->int * { return a.GetBlockI(); });
    t.method("GetBlockI", [](mfem::BlockILU* a)->int * { return a->GetBlockI(); });

    DEBUG_MSG("Adding wrapper for int * mfem::BlockILU::GetBlockJ() (" __HERE__ ")");
    // signature to use in the veto list: int * mfem::BlockILU::GetBlockJ()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:962:9
    t.method("GetBlockJ", [](mfem::BlockILU& a)->int * { return a.GetBlockJ(); });
    t.method("GetBlockJ", [](mfem::BlockILU* a)->int * { return a->GetBlockJ(); });

    DEBUG_MSG("Adding wrapper for double * mfem::BlockILU::GetBlockData() (" __HERE__ ")");
    // signature to use in the veto list: double * mfem::BlockILU::GetBlockData()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:967:12
    t.method("GetBlockData", [](mfem::BlockILU& a)->double * { return a.GetBlockData(); });
    t.method("GetBlockData", [](mfem::BlockILU* a)->double * { return a->GetBlockData(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::BlockILU>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_BlockILU(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_BlockILU(module));
}
