// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MassIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MassIntegrator> : std::false_type { };
template<> struct SuperType<mfem::MassIntegrator> { typedef mfem::BilinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::MassIntegrator
// signature to use in the veto file: mfem::MassIntegrator
struct Jlmfem_MassIntegrator: public Wrapper {

  Jlmfem_MassIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MassIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2170:7
    jlcxx::TypeWrapper<mfem::MassIntegrator>  t = jlModule.add_type<mfem::MassIntegrator>("mfem!MassIntegrator",
      jlcxx::julia_base_type<mfem::BilinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MassIntegrator>>(new jlcxx::TypeWrapper<mfem::MassIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::MassIntegrator(const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2185:4
    t.constructor<const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::MassIntegrator(mfem::Coefficient &, const mfem::IntegrationRule *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2189:4
    t.constructor<mfem::Coefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::Coefficient &, const mfem::IntegrationRule *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2194:17
    t.method("AssembleElementMatrix", [](mfem::MassIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::MassIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2197:17
    t.method("AssembleElementMatrix2", [](mfem::MassIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementMatrix2", [](mfem::MassIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementMatrix2(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssembleMF(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssembleMF(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2204:17
    t.method("AssembleMF", [](mfem::MassIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssembleMF(arg0); });
    t.method("AssembleMF", [](mfem::MassIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssembleMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssemblePA(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssemblePA(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2206:17
    t.method("AssemblePA", [](mfem::MassIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePA(arg0); });
    t.method("AssemblePA", [](mfem::MassIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssembleEA(const mfem::FiniteElementSpace &, mfem::Vector &, const bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssembleEA(const mfem::FiniteElementSpace &, mfem::Vector &, const bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2208:17
    t.method("AssembleEA", [](mfem::MassIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, const bool arg2)->void { a.AssembleEA(arg0, arg1, arg2); });
    t.method("AssembleEA", [](mfem::MassIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, const bool arg2)->void { a->AssembleEA(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssembleDiagonalPA(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssembleDiagonalPA(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2211:17
    t.method("AssembleDiagonalPA", [](mfem::MassIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalPA(arg0); });
    t.method("AssembleDiagonalPA", [](mfem::MassIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalPA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AssembleDiagonalMF(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AssembleDiagonalMF(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2213:17
    t.method("AssembleDiagonalMF", [](mfem::MassIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalMF(arg0); });
    t.method("AssembleDiagonalMF", [](mfem::MassIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2215:17
    t.method("AddMultMF", [](mfem::MassIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultMF(arg0, arg1); });
    t.method("AddMultMF", [](mfem::MassIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultMF(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2217:17
    t.method("AddMultPA", [](mfem::MassIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultPA(arg0, arg1); });
    t.method("AddMultPA", [](mfem::MassIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultPA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::MassIntegrator::AddMultTransposePA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MassIntegrator::AddMultTransposePA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2219:17
    t.method("AddMultTransposePA", [](mfem::MassIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultTransposePA(arg0, arg1); });
    t.method("AddMultTransposePA", [](mfem::MassIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultTransposePA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const mfem::IntegrationRule & mfem::MassIntegrator::GetRule(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &) (" __HERE__ ")");
    // signature to use in the veto list: const mfem::IntegrationRule & mfem::MassIntegrator::GetRule(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2221:34
    module_.method("mfem!MassIntegrator!GetRule", [](const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2)->const mfem::IntegrationRule & { return mfem::MassIntegrator::GetRule(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for bool mfem::MassIntegrator::SupportsCeed() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::MassIntegrator::SupportsCeed()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:2225:9
    t.method("SupportsCeed", [](mfem::MassIntegrator const& a)->bool { return a.SupportsCeed(); });
    t.method("SupportsCeed", [](mfem::MassIntegrator const* a)->bool { return a->SupportsCeed(); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MassIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MassIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MassIntegrator(module));
}
