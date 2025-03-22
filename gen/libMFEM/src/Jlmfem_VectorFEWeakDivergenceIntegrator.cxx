// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorFEWeakDivergenceIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorFEWeakDivergenceIntegrator> : std::false_type { };
template<> struct SuperType<mfem::VectorFEWeakDivergenceIntegrator> { typedef mfem::BilinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::VectorFEWeakDivergenceIntegrator
// signature to use in the veto file: mfem::VectorFEWeakDivergenceIntegrator
struct Jlmfem_VectorFEWeakDivergenceIntegrator: public Wrapper {

  Jlmfem_VectorFEWeakDivergenceIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorFEWeakDivergenceIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2437:7
    jlcxx::TypeWrapper<mfem::VectorFEWeakDivergenceIntegrator>  t = jlModule.add_type<mfem::VectorFEWeakDivergenceIntegrator>("mfem!VectorFEWeakDivergenceIntegrator",
      jlcxx::julia_base_type<mfem::BilinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorFEWeakDivergenceIntegrator>>(new jlcxx::TypeWrapper<mfem::VectorFEWeakDivergenceIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorFEWeakDivergenceIntegrator::VectorFEWeakDivergenceIntegrator(mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2452:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::VectorFEWeakDivergenceIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorFEWeakDivergenceIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2453:17
    t.method("AssembleElementMatrix", [](mfem::VectorFEWeakDivergenceIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::VectorFEWeakDivergenceIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorFEWeakDivergenceIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorFEWeakDivergenceIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2456:17
    t.method("AssembleElementMatrix2", [](mfem::VectorFEWeakDivergenceIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementMatrix2", [](mfem::VectorFEWeakDivergenceIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorFEWeakDivergenceIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorFEWeakDivergenceIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorFEWeakDivergenceIntegrator(module));
}
