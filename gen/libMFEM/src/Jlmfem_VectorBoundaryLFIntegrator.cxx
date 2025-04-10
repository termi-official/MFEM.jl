// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorBoundaryLFIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorBoundaryLFIntegrator> : std::false_type { };
template<> struct SuperType<mfem::VectorBoundaryLFIntegrator> { typedef mfem::LinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::VectorBoundaryLFIntegrator
// signature to use in the veto file: mfem::VectorBoundaryLFIntegrator
struct Jlmfem_VectorBoundaryLFIntegrator: public Wrapper {

  Jlmfem_VectorBoundaryLFIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorBoundaryLFIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:241:7
    jlcxx::TypeWrapper<mfem::VectorBoundaryLFIntegrator>  t = jlModule.add_type<mfem::VectorBoundaryLFIntegrator>("mfem!VectorBoundaryLFIntegrator",
      jlcxx::julia_base_type<mfem::LinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorBoundaryLFIntegrator>>(new jlcxx::TypeWrapper<mfem::VectorBoundaryLFIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::VectorBoundaryLFIntegrator::VectorBoundaryLFIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:249:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::VectorBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:253:17
    t.method("AssembleRHSElementVect", [](mfem::VectorBoundaryLFIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a.AssembleRHSElementVect(arg0, arg1, arg2); });
    t.method("AssembleRHSElementVect", [](mfem::VectorBoundaryLFIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->void { a->AssembleRHSElementVect(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorBoundaryLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/lininteg.hpp:258:17
    t.method("AssembleRHSElementVect", [](mfem::VectorBoundaryLFIntegrator& a, const mfem::FiniteElement & arg0, mfem::FaceElementTransformations & arg1, mfem::Vector & arg2)->void { a.AssembleRHSElementVect(arg0, arg1, arg2); });
    t.method("AssembleRHSElementVect", [](mfem::VectorBoundaryLFIntegrator* a, const mfem::FiniteElement & arg0, mfem::FaceElementTransformations & arg1, mfem::Vector & arg2)->void { a->AssembleRHSElementVect(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorBoundaryLFIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorBoundaryLFIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorBoundaryLFIntegrator(module));
}
