// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::VectorMassIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::VectorMassIntegrator> : std::false_type { };
template<> struct SuperType<mfem::VectorMassIntegrator> { typedef mfem::BilinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::VectorMassIntegrator
// signature to use in the veto file: mfem::VectorMassIntegrator
struct Jlmfem_VectorMassIntegrator: public Wrapper {

  Jlmfem_VectorMassIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::VectorMassIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2330:7
    jlcxx::TypeWrapper<mfem::VectorMassIntegrator>  t = jlModule.add_type<mfem::VectorMassIntegrator>("mfem!VectorMassIntegrator",
      jlcxx::julia_base_type<mfem::BilinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorMassIntegrator>>(new jlcxx::TypeWrapper<mfem::VectorMassIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::VectorMassIntegrator(mfem::Coefficient &, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2356:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::Coefficient &, int>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::VectorMassIntegrator(mfem::Coefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2358:4
    t.constructor<mfem::Coefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::VectorMassIntegrator(mfem::VectorCoefficient &, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2362:4
    t.constructor<mfem::VectorCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::VectorCoefficient &, int>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::VectorMassIntegrator(mfem::MatrixCoefficient &, int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2365:4
    t.constructor<mfem::MatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::MatrixCoefficient &, int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for int mfem::VectorMassIntegrator::GetVDim() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::VectorMassIntegrator::GetVDim()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2368:8
    t.method("GetVDim", [](mfem::VectorMassIntegrator const& a)->int { return a.GetVDim(); });
    t.method("GetVDim", [](mfem::VectorMassIntegrator const* a)->int { return a->GetVDim(); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::SetVDim(int) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::SetVDim(int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2369:9
    t.method("SetVDim", [](mfem::VectorMassIntegrator& a, int arg0)->void { a.SetVDim(arg0); });
    t.method("SetVDim", [](mfem::VectorMassIntegrator* a, int arg0)->void { a->SetVDim(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2371:17
    t.method("AssembleElementMatrix", [](mfem::VectorMassIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::VectorMassIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2374:17
    t.method("AssembleElementMatrix2", [](mfem::VectorMassIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementMatrix2", [](mfem::VectorMassIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementMatrix2(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AssemblePA(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AssemblePA(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2379:17
    t.method("AssemblePA", [](mfem::VectorMassIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePA(arg0); });
    t.method("AssemblePA", [](mfem::VectorMassIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AssembleMF(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AssembleMF(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2380:17
    t.method("AssembleMF", [](mfem::VectorMassIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssembleMF(arg0); });
    t.method("AssembleMF", [](mfem::VectorMassIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssembleMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AssembleDiagonalPA(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AssembleDiagonalPA(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2381:17
    t.method("AssembleDiagonalPA", [](mfem::VectorMassIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalPA(arg0); });
    t.method("AssembleDiagonalPA", [](mfem::VectorMassIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalPA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AssembleDiagonalMF(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AssembleDiagonalMF(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2382:17
    t.method("AssembleDiagonalMF", [](mfem::VectorMassIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalMF(arg0); });
    t.method("AssembleDiagonalMF", [](mfem::VectorMassIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2383:17
    t.method("AddMultPA", [](mfem::VectorMassIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultPA(arg0, arg1); });
    t.method("AddMultPA", [](mfem::VectorMassIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultPA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::VectorMassIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::VectorMassIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2384:17
    t.method("AddMultMF", [](mfem::VectorMassIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultMF(arg0, arg1); });
    t.method("AddMultMF", [](mfem::VectorMassIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultMF(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for bool mfem::VectorMassIntegrator::SupportsCeed() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::VectorMassIntegrator::SupportsCeed()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2385:9
    t.method("SupportsCeed", [](mfem::VectorMassIntegrator const& a)->bool { return a.SupportsCeed(); });
    t.method("SupportsCeed", [](mfem::VectorMassIntegrator const* a)->bool { return a->SupportsCeed(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::VectorMassIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_VectorMassIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_VectorMassIntegrator(module));
}
