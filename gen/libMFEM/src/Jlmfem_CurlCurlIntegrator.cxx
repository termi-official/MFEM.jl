// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::CurlCurlIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::CurlCurlIntegrator> : std::false_type { };
template<> struct SuperType<mfem::CurlCurlIntegrator> { typedef mfem::BilinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::CurlCurlIntegrator
// signature to use in the veto file: mfem::CurlCurlIntegrator
struct Jlmfem_CurlCurlIntegrator: public Wrapper {

  Jlmfem_CurlCurlIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::CurlCurlIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2512:7
    jlcxx::TypeWrapper<mfem::CurlCurlIntegrator>  t = jlModule.add_type<mfem::CurlCurlIntegrator>("mfem!CurlCurlIntegrator",
      jlcxx::julia_base_type<mfem::BilinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::CurlCurlIntegrator>>(new jlcxx::TypeWrapper<mfem::CurlCurlIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::CurlCurlIntegrator(mfem::Coefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2539:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::Coefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::CurlCurlIntegrator(mfem::DiagonalMatrixCoefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2541:4
    t.constructor<mfem::DiagonalMatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::DiagonalMatrixCoefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::CurlCurlIntegrator(mfem::MatrixCoefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2544:4
    t.constructor<mfem::MatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::MatrixCoefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::CurlCurlIntegrator(mfem::SymmetricMatrixCoefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2546:4
    t.constructor<mfem::SymmetricMatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::SymmetricMatrixCoefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CurlCurlIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2552:17
    t.method("AssembleElementMatrix", [](mfem::CurlCurlIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::CurlCurlIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::ComputeElementFlux(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, const mfem::FiniteElement &, mfem::Vector &, bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CurlCurlIntegrator::ComputeElementFlux(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, const mfem::FiniteElement &, mfem::Vector &, bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2556:17
    t.method("ComputeElementFlux", [](mfem::CurlCurlIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, const mfem::FiniteElement & arg3, mfem::Vector & arg4, bool arg5)->void { a.ComputeElementFlux(arg0, arg1, arg2, arg3, arg4, arg5); });
    t.method("ComputeElementFlux", [](mfem::CurlCurlIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, const mfem::FiniteElement & arg3, mfem::Vector & arg4, bool arg5)->void { a->ComputeElementFlux(arg0, arg1, arg2, arg3, arg4, arg5); });

    DEBUG_MSG("Adding wrapper for double mfem::CurlCurlIntegrator::ComputeFluxEnergy(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, mfem::Vector *) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::CurlCurlIntegrator::ComputeFluxEnergy(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, mfem::Vector *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2561:19
    t.method("ComputeFluxEnergy", [](mfem::CurlCurlIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->double { return a.ComputeFluxEnergy(arg0, arg1, arg2); });
    t.method("ComputeFluxEnergy", [](mfem::CurlCurlIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, mfem::Vector * arg3)->double { return a.ComputeFluxEnergy(arg0, arg1, arg2, arg3); });
    t.method("ComputeFluxEnergy", [](mfem::CurlCurlIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->double { return a->ComputeFluxEnergy(arg0, arg1, arg2); });
    t.method("ComputeFluxEnergy", [](mfem::CurlCurlIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, mfem::Vector * arg3)->double { return a->ComputeFluxEnergy(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::AssemblePA(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CurlCurlIntegrator::AssemblePA(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2566:17
    t.method("AssemblePA", [](mfem::CurlCurlIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePA(arg0); });
    t.method("AssemblePA", [](mfem::CurlCurlIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CurlCurlIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2567:17
    t.method("AddMultPA", [](mfem::CurlCurlIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultPA(arg0, arg1); });
    t.method("AddMultPA", [](mfem::CurlCurlIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultPA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::CurlCurlIntegrator::AssembleDiagonalPA(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::CurlCurlIntegrator::AssembleDiagonalPA(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2568:17
    t.method("AssembleDiagonalPA", [](mfem::CurlCurlIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalPA(arg0); });
    t.method("AssembleDiagonalPA", [](mfem::CurlCurlIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalPA(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::CurlCurlIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_CurlCurlIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_CurlCurlIntegrator(module));
}
