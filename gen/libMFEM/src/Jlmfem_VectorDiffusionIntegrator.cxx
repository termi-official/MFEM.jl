// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorDiffusionIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorDiffusionIntegrator> : std::false_type { };
template<> struct SuperType<mfem::VectorDiffusionIntegrator> { typedef mfem::BilinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::VectorDiffusionIntegrator
// signature to use in the veto file: mfem::VectorDiffusionIntegrator
struct Jlmfem_VectorDiffusionIntegrator: public Wrapper {

  Jlmfem_VectorDiffusionIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorDiffusionIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2754:7
    jlcxx::TypeWrapper<mfem::VectorDiffusionIntegrator>  t = jlModule.add_type<mfem::VectorDiffusionIntegrator>("mfem!VectorDiffusionIntegrator",
      jlcxx::julia_base_type<mfem::BilinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorDiffusionIntegrator>>(new jlcxx::TypeWrapper<mfem::VectorDiffusionIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::VectorDiffusionIntegrator(int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2781:4
    t.constructor<int>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::VectorDiffusionIntegrator(mfem::Coefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2784:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::VectorDiffusionIntegrator(mfem::Coefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2787:4
    t.constructor<mfem::Coefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::VectorDiffusionIntegrator(mfem::Coefficient &, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2798:4
    t.constructor<mfem::Coefficient &, int>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::VectorDiffusionIntegrator(mfem::VectorCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2810:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::VectorDiffusionIntegrator(mfem::MatrixCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2822:4
    t.constructor<mfem::MatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2825:17
    t.method("AssembleElementMatrix", [](mfem::VectorDiffusionIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::VectorDiffusionIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AssembleElementVector(const mfem::FiniteElement &, mfem::ElementTransformation &, const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AssembleElementVector(const mfem::FiniteElement &, mfem::ElementTransformation &, const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2828:17
    t.method("AssembleElementVector", [](mfem::VectorDiffusionIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, const mfem::Vector & arg2, mfem::Vector & arg3)->void { a.AssembleElementVector(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementVector", [](mfem::VectorDiffusionIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, const mfem::Vector & arg2, mfem::Vector & arg3)->void { a->AssembleElementVector(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AssemblePA(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AssemblePA(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2832:17
    t.method("AssemblePA", [](mfem::VectorDiffusionIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePA(arg0); });
    t.method("AssemblePA", [](mfem::VectorDiffusionIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AssembleMF(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AssembleMF(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2833:17
    t.method("AssembleMF", [](mfem::VectorDiffusionIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssembleMF(arg0); });
    t.method("AssembleMF", [](mfem::VectorDiffusionIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssembleMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AssembleDiagonalPA(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AssembleDiagonalPA(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2834:17
    t.method("AssembleDiagonalPA", [](mfem::VectorDiffusionIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalPA(arg0); });
    t.method("AssembleDiagonalPA", [](mfem::VectorDiffusionIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalPA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AssembleDiagonalMF(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AssembleDiagonalMF(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2835:17
    t.method("AssembleDiagonalMF", [](mfem::VectorDiffusionIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalMF(arg0); });
    t.method("AssembleDiagonalMF", [](mfem::VectorDiffusionIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2836:17
    t.method("AddMultPA", [](mfem::VectorDiffusionIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultPA(arg0, arg1); });
    t.method("AddMultPA", [](mfem::VectorDiffusionIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultPA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorDiffusionIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorDiffusionIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2837:17
    t.method("AddMultMF", [](mfem::VectorDiffusionIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultMF(arg0, arg1); });
    t.method("AddMultMF", [](mfem::VectorDiffusionIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultMF(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for bool mfem::VectorDiffusionIntegrator::SupportsCeed() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::VectorDiffusionIntegrator::SupportsCeed()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2838:9
    t.method("SupportsCeed", [](mfem::VectorDiffusionIntegrator const& a)->bool { return a.SupportsCeed(); });
    t.method("SupportsCeed", [](mfem::VectorDiffusionIntegrator const* a)->bool { return a->SupportsCeed(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorDiffusionIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorDiffusionIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorDiffusionIntegrator(module));
}
