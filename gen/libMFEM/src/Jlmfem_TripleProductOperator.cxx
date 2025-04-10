// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::TripleProductOperator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::TripleProductOperator> : std::false_type { };
template<> struct SuperType<mfem::TripleProductOperator> { typedef mfem::Operator type; };
}

// Class generating the wrapper for type mfem::TripleProductOperator
// signature to use in the veto file: mfem::TripleProductOperator
struct Jlmfem_TripleProductOperator: public Wrapper {

  Jlmfem_TripleProductOperator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::TripleProductOperator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:807:7
    jlcxx::TypeWrapper<mfem::TripleProductOperator>  t = jlModule.add_type<mfem::TripleProductOperator>("mfem!TripleProductOperator",
      jlcxx::julia_base_type<mfem::Operator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::TripleProductOperator>>(new jlcxx::TypeWrapper<mfem::TripleProductOperator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::TripleProductOperator::TripleProductOperator(const mfem::Operator *, const mfem::Operator *, const mfem::Operator *, bool, bool, bool) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:817:4
    t.constructor<const mfem::Operator *, const mfem::Operator *, const mfem::Operator *, bool, bool, bool>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for mfem::MemoryClass mfem::TripleProductOperator::GetMemoryClass() (" __HERE__ ")");
    // signature to use in the veto list: mfem::MemoryClass mfem::TripleProductOperator::GetMemoryClass()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:820:24
    t.method("GetMemoryClass", [](mfem::TripleProductOperator const& a)->mfem::MemoryClass { return a.GetMemoryClass(); });
    t.method("GetMemoryClass", [](mfem::TripleProductOperator const* a)->mfem::MemoryClass { return a->GetMemoryClass(); });

    DEBUG_MSG("Adding wrapper for void mfem::TripleProductOperator::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TripleProductOperator::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:822:17
    t.method("Mult", [](mfem::TripleProductOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::TripleProductOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::TripleProductOperator::MultTranspose(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TripleProductOperator::MultTranspose(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:825:17
    t.method("MultTranspose", [](mfem::TripleProductOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.MultTranspose(arg0, arg1); });
    t.method("MultTranspose", [](mfem::TripleProductOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->MultTranspose(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::TripleProductOperator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_TripleProductOperator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_TripleProductOperator(module));
}
