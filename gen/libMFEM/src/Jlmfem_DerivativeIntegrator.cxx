// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::DerivativeIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::DerivativeIntegrator> : std::false_type { };
template<> struct SuperType<mfem::DerivativeIntegrator> { typedef mfem::BilinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::DerivativeIntegrator
// signature to use in the veto file: mfem::DerivativeIntegrator
struct Jlmfem_DerivativeIntegrator: public Wrapper {

  Jlmfem_DerivativeIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::DerivativeIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2489:7
    jlcxx::TypeWrapper<mfem::DerivativeIntegrator>  t = jlModule.add_type<mfem::DerivativeIntegrator>("mfem!DerivativeIntegrator",
      jlcxx::julia_base_type<mfem::BilinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::DerivativeIntegrator>>(new jlcxx::TypeWrapper<mfem::DerivativeIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::DerivativeIntegrator::DerivativeIntegrator(mfem::Coefficient &, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2500:4
    t.constructor<mfem::Coefficient &, int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::DerivativeIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::DerivativeIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2501:17
    t.method("AssembleElementMatrix", [](mfem::DerivativeIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::DerivativeIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::DerivativeIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::DerivativeIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2505:17
    t.method("AssembleElementMatrix2", [](mfem::DerivativeIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementMatrix2", [](mfem::DerivativeIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::DerivativeIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_DerivativeIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_DerivativeIntegrator(module));
}
