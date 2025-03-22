// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::ConstrainedOperator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::ConstrainedOperator> : std::false_type { };
template<> struct SuperType<mfem::ConstrainedOperator> { typedef mfem::Operator type; };
}

// Class generating the wrapper for type mfem::ConstrainedOperator
// signature to use in the veto file: mfem::ConstrainedOperator
struct Jlmfem_ConstrainedOperator: public Wrapper {

  Jlmfem_ConstrainedOperator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::ConstrainedOperator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:841:7
    jlcxx::TypeWrapper<mfem::ConstrainedOperator>  t = jlModule.add_type<mfem::ConstrainedOperator>("mfem!ConstrainedOperator",
      jlcxx::julia_base_type<mfem::Operator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::ConstrainedOperator>>(new jlcxx::TypeWrapper<mfem::ConstrainedOperator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::ConstrainedOperator::ConstrainedOperator(mfem::Operator *, const mfem::Array<int> &, bool, mfem::Operator::DiagonalPolicy) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:860:4
    t.constructor<mfem::Operator *, const mfem::Array<int> &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::Operator *, const mfem::Array<int> &, bool>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::Operator *, const mfem::Array<int> &, bool, mfem::Operator::DiagonalPolicy>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for mfem::MemoryClass mfem::ConstrainedOperator::GetMemoryClass() (" __HERE__ ")");
    // signature to use in the veto list: mfem::MemoryClass mfem::ConstrainedOperator::GetMemoryClass()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:864:24
    t.method("GetMemoryClass", [](mfem::ConstrainedOperator const& a)->mfem::MemoryClass { return a.GetMemoryClass(); });
    t.method("GetMemoryClass", [](mfem::ConstrainedOperator const* a)->mfem::MemoryClass { return a->GetMemoryClass(); });

    DEBUG_MSG("Adding wrapper for void mfem::ConstrainedOperator::SetDiagonalPolicy(const mfem::Operator::DiagonalPolicy) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ConstrainedOperator::SetDiagonalPolicy(const mfem::Operator::DiagonalPolicy)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:867:9
    t.method("SetDiagonalPolicy", [](mfem::ConstrainedOperator& a, const mfem::Operator::DiagonalPolicy arg0)->void { a.SetDiagonalPolicy(arg0); });
    t.method("SetDiagonalPolicy", [](mfem::ConstrainedOperator* a, const mfem::Operator::DiagonalPolicy arg0)->void { a->SetDiagonalPolicy(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::ConstrainedOperator::AssembleDiagonal(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ConstrainedOperator::AssembleDiagonal(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:871:17
    t.method("AssembleDiagonal", [](mfem::ConstrainedOperator const& a, mfem::Vector & arg0)->void { a.AssembleDiagonal(arg0); });
    t.method("AssembleDiagonal", [](mfem::ConstrainedOperator const* a, mfem::Vector & arg0)->void { a->AssembleDiagonal(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::ConstrainedOperator::EliminateRHS(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ConstrainedOperator::EliminateRHS(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:882:9
    t.method("EliminateRHS", [](mfem::ConstrainedOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.EliminateRHS(arg0, arg1); });
    t.method("EliminateRHS", [](mfem::ConstrainedOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->EliminateRHS(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::ConstrainedOperator::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::ConstrainedOperator::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:892:17
    t.method("Mult", [](mfem::ConstrainedOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::ConstrainedOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::ConstrainedOperator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_ConstrainedOperator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_ConstrainedOperator(module));
}
