// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorCrossProductInterpolator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorCrossProductInterpolator> : std::false_type { };
template<> struct SuperType<mfem::VectorCrossProductInterpolator> { typedef mfem::DiscreteInterpolator type; };
}

// Class generating the wrapper for type mfem::VectorCrossProductInterpolator
// signature to use in the veto file: mfem::VectorCrossProductInterpolator
struct Jlmfem_VectorCrossProductInterpolator: public Wrapper {

  Jlmfem_VectorCrossProductInterpolator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorCrossProductInterpolator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:3443:7
    jlcxx::TypeWrapper<mfem::VectorCrossProductInterpolator>  t = jlModule.add_type<mfem::VectorCrossProductInterpolator>("mfem!VectorCrossProductInterpolator",
      jlcxx::julia_base_type<mfem::DiscreteInterpolator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorCrossProductInterpolator>>(new jlcxx::TypeWrapper<mfem::VectorCrossProductInterpolator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::VectorCrossProductInterpolator::VectorCrossProductInterpolator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:3446:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::VectorCrossProductInterpolator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorCrossProductInterpolator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:3449:17
    t.method("AssembleElementMatrix2", [](mfem::VectorCrossProductInterpolator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementMatrix2", [](mfem::VectorCrossProductInterpolator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorCrossProductInterpolator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorCrossProductInterpolator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorCrossProductInterpolator(module));
}
