// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::TransposeOperator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::TransposeOperator> : std::false_type { };
template<> struct SuperType<mfem::TransposeOperator> { typedef mfem::Operator type; };
}

// Class generating the wrapper for type mfem::TransposeOperator
// signature to use in the veto file: mfem::TransposeOperator
struct Jlmfem_TransposeOperator: public Wrapper {

  Jlmfem_TransposeOperator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::TransposeOperator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:718:7
    jlcxx::TypeWrapper<mfem::TransposeOperator>  t = jlModule.add_type<mfem::TransposeOperator>("mfem!TransposeOperator",
      jlcxx::julia_base_type<mfem::Operator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::TransposeOperator>>(new jlcxx::TypeWrapper<mfem::TransposeOperator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::TransposeOperator::TransposeOperator(const mfem::Operator *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:725:4
    t.constructor<const mfem::Operator *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::TransposeOperator::TransposeOperator(const mfem::Operator &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:729:4
    t.constructor<const mfem::Operator &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::TransposeOperator::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TransposeOperator::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:733:17
    t.method("Mult", [](mfem::TransposeOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::TransposeOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::TransposeOperator::MultTranspose(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TransposeOperator::MultTranspose(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:737:17
    t.method("MultTranspose", [](mfem::TransposeOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.MultTranspose(arg0, arg1); });
    t.method("MultTranspose", [](mfem::TransposeOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->MultTranspose(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::TransposeOperator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_TransposeOperator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_TransposeOperator(module));
}
