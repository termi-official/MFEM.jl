// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::BilinearFormIntegrator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::BilinearFormIntegrator> : std::false_type { };
template<> struct SuperType<mfem::BilinearFormIntegrator> { typedef mfem::NonlinearFormIntegrator type; };
}

// Class generating the wrapper for type mfem::BilinearFormIntegrator
// signature to use in the veto file: mfem::BilinearFormIntegrator
struct Jlmfem_BilinearFormIntegrator: public Wrapper {

  Jlmfem_BilinearFormIntegrator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::BilinearFormIntegrator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:34:7
    jlcxx::TypeWrapper<mfem::BilinearFormIntegrator>  t = jlModule.add_type<mfem::BilinearFormIntegrator>("mfem!BilinearFormIntegrator",
      jlcxx::julia_base_type<mfem::NonlinearFormIntegrator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::BilinearFormIntegrator>>(new jlcxx::TypeWrapper<mfem::BilinearFormIntegrator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssemblePA(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssemblePA(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:58:17
    t.method("AssemblePA", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePA(arg0); });
    t.method("AssemblePA", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssemblePA(const mfem::FiniteElementSpace &, const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssemblePA(const mfem::FiniteElementSpace &, const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:60:17
    t.method("AssemblePA", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, const mfem::FiniteElementSpace & arg1)->void { a.AssemblePA(arg0, arg1); });
    t.method("AssemblePA", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, const mfem::FiniteElementSpace & arg1)->void { a->AssemblePA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssemblePAInteriorFaces(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssemblePAInteriorFaces(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:63:17
    t.method("AssemblePAInteriorFaces", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePAInteriorFaces(arg0); });
    t.method("AssemblePAInteriorFaces", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePAInteriorFaces(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssemblePABoundaryFaces(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssemblePABoundaryFaces(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:65:17
    t.method("AssemblePABoundaryFaces", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssemblePABoundaryFaces(arg0); });
    t.method("AssemblePABoundaryFaces", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssemblePABoundaryFaces(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleDiagonalPA(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleDiagonalPA(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:68:17
    t.method("AssembleDiagonalPA", [](mfem::BilinearFormIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalPA(arg0); });
    t.method("AssembleDiagonalPA", [](mfem::BilinearFormIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalPA(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleDiagonalPA_ADAt(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleDiagonalPA_ADAt(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:71:17
    t.method("AssembleDiagonalPA_ADAt", [](mfem::BilinearFormIntegrator& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AssembleDiagonalPA_ADAt(arg0, arg1); });
    t.method("AssembleDiagonalPA_ADAt", [](mfem::BilinearFormIntegrator* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AssembleDiagonalPA_ADAt(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AddMultPA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:80:17
    t.method("AddMultPA", [](mfem::BilinearFormIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultPA(arg0, arg1); });
    t.method("AddMultPA", [](mfem::BilinearFormIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultPA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AddMultTransposePA(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AddMultTransposePA(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:89:17
    t.method("AddMultTransposePA", [](mfem::BilinearFormIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultTransposePA(arg0, arg1); });
    t.method("AddMultTransposePA", [](mfem::BilinearFormIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultTransposePA(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleEA(const mfem::FiniteElementSpace &, mfem::Vector &, const bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleEA(const mfem::FiniteElementSpace &, mfem::Vector &, const bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:94:17
    t.method("AssembleEA", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1)->void { a.AssembleEA(arg0, arg1); });
    t.method("AssembleEA", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, const bool arg2)->void { a.AssembleEA(arg0, arg1, arg2); });
    t.method("AssembleEA", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1)->void { a->AssembleEA(arg0, arg1); });
    t.method("AssembleEA", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, const bool arg2)->void { a->AssembleEA(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleMF(const mfem::FiniteElementSpace &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleMF(const mfem::FiniteElementSpace &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:104:17
    t.method("AssembleMF", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0)->void { a.AssembleMF(arg0); });
    t.method("AssembleMF", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0)->void { a->AssembleMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AddMultMF(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:112:17
    t.method("AddMultMF", [](mfem::BilinearFormIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultMF(arg0, arg1); });
    t.method("AddMultMF", [](mfem::BilinearFormIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultMF(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AddMultTransposeMF(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AddMultTransposeMF(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:120:17
    t.method("AddMultTransposeMF", [](mfem::BilinearFormIntegrator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultTransposeMF(arg0, arg1); });
    t.method("AddMultTransposeMF", [](mfem::BilinearFormIntegrator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultTransposeMF(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleDiagonalMF(mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleDiagonalMF(mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:123:17
    t.method("AssembleDiagonalMF", [](mfem::BilinearFormIntegrator& a, mfem::Vector & arg0)->void { a.AssembleDiagonalMF(arg0); });
    t.method("AssembleDiagonalMF", [](mfem::BilinearFormIntegrator* a, mfem::Vector & arg0)->void { a->AssembleDiagonalMF(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleEAInteriorFaces(const mfem::FiniteElementSpace &, mfem::Vector &, mfem::Vector &, const bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleEAInteriorFaces(const mfem::FiniteElementSpace &, mfem::Vector &, mfem::Vector &, const bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:125:17
    t.method("AssembleEAInteriorFaces", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, mfem::Vector & arg2)->void { a.AssembleEAInteriorFaces(arg0, arg1, arg2); });
    t.method("AssembleEAInteriorFaces", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, mfem::Vector & arg2, const bool arg3)->void { a.AssembleEAInteriorFaces(arg0, arg1, arg2, arg3); });
    t.method("AssembleEAInteriorFaces", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, mfem::Vector & arg2)->void { a->AssembleEAInteriorFaces(arg0, arg1, arg2); });
    t.method("AssembleEAInteriorFaces", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, mfem::Vector & arg2, const bool arg3)->void { a->AssembleEAInteriorFaces(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleEABoundaryFaces(const mfem::FiniteElementSpace &, mfem::Vector &, const bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleEABoundaryFaces(const mfem::FiniteElementSpace &, mfem::Vector &, const bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:130:17
    t.method("AssembleEABoundaryFaces", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1)->void { a.AssembleEABoundaryFaces(arg0, arg1); });
    t.method("AssembleEABoundaryFaces", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, const bool arg2)->void { a.AssembleEABoundaryFaces(arg0, arg1, arg2); });
    t.method("AssembleEABoundaryFaces", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1)->void { a->AssembleEABoundaryFaces(arg0, arg1); });
    t.method("AssembleEABoundaryFaces", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElementSpace & arg0, mfem::Vector & arg1, const bool arg2)->void { a->AssembleEABoundaryFaces(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleElementMatrix(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:135:17
    t.method("AssembleElementMatrix", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a.AssembleElementMatrix(arg0, arg1, arg2); });
    t.method("AssembleElementMatrix", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::DenseMatrix & arg2)->void { a->AssembleElementMatrix(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleElementMatrix2(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:143:17
    t.method("AssembleElementMatrix2", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementMatrix2(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementMatrix2", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::ElementTransformation & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementMatrix2(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleFaceMatrix(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleFaceMatrix(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:148:17
    t.method("AssembleFaceMatrix", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleFaceMatrix(arg0, arg1, arg2, arg3); });
    t.method("AssembleFaceMatrix", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleFaceMatrix(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleFaceMatrix(const mfem::FiniteElement &, const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleFaceMatrix(const mfem::FiniteElement &, const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:155:17
    t.method("AssembleFaceMatrix", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, const mfem::FiniteElement & arg2, mfem::FaceElementTransformations & arg3, mfem::DenseMatrix & arg4)->void { a.AssembleFaceMatrix(arg0, arg1, arg2, arg3, arg4); });
    t.method("AssembleFaceMatrix", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, const mfem::FiniteElement & arg2, mfem::FaceElementTransformations & arg3, mfem::DenseMatrix & arg4)->void { a->AssembleFaceMatrix(arg0, arg1, arg2, arg3, arg4); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleElementVector(const mfem::FiniteElement &, mfem::ElementTransformation &, const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleElementVector(const mfem::FiniteElement &, mfem::ElementTransformation &, const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:164:17
    t.method("AssembleElementVector", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, const mfem::Vector & arg2, mfem::Vector & arg3)->void { a.AssembleElementVector(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementVector", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, const mfem::Vector & arg2, mfem::Vector & arg3)->void { a->AssembleElementVector(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleFaceVector(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleFaceVector(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:172:17
    t.method("AssembleFaceVector", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, const mfem::Vector & arg3, mfem::Vector & arg4)->void { a.AssembleFaceVector(arg0, arg1, arg2, arg3, arg4); });
    t.method("AssembleFaceVector", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, const mfem::Vector & arg3, mfem::Vector & arg4)->void { a->AssembleFaceVector(arg0, arg1, arg2, arg3, arg4); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleElementGrad(const mfem::FiniteElement &, mfem::ElementTransformation &, const mfem::Vector &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleElementGrad(const mfem::FiniteElement &, mfem::ElementTransformation &, const mfem::Vector &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:177:17
    t.method("AssembleElementGrad", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, const mfem::Vector & arg2, mfem::DenseMatrix & arg3)->void { a.AssembleElementGrad(arg0, arg1, arg2, arg3); });
    t.method("AssembleElementGrad", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, const mfem::Vector & arg2, mfem::DenseMatrix & arg3)->void { a->AssembleElementGrad(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::AssembleFaceGrad(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, const mfem::Vector &, mfem::DenseMatrix &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::AssembleFaceGrad(const mfem::FiniteElement &, const mfem::FiniteElement &, mfem::FaceElementTransformations &, const mfem::Vector &, mfem::DenseMatrix &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:182:17
    t.method("AssembleFaceGrad", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, const mfem::Vector & arg3, mfem::DenseMatrix & arg4)->void { a.AssembleFaceGrad(arg0, arg1, arg2, arg3, arg4); });
    t.method("AssembleFaceGrad", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, const mfem::FiniteElement & arg1, mfem::FaceElementTransformations & arg2, const mfem::Vector & arg3, mfem::DenseMatrix & arg4)->void { a->AssembleFaceGrad(arg0, arg1, arg2, arg3, arg4); });

    DEBUG_MSG("Adding wrapper for void mfem::BilinearFormIntegrator::ComputeElementFlux(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, const mfem::FiniteElement &, mfem::Vector &, bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::BilinearFormIntegrator::ComputeElementFlux(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, const mfem::FiniteElement &, mfem::Vector &, bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:217:17
    t.method("ComputeElementFlux", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, const mfem::FiniteElement & arg3, mfem::Vector & arg4)->void { a.ComputeElementFlux(arg0, arg1, arg2, arg3, arg4); });
    t.method("ComputeElementFlux", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, const mfem::FiniteElement & arg3, mfem::Vector & arg4, bool arg5)->void { a.ComputeElementFlux(arg0, arg1, arg2, arg3, arg4, arg5); });
    t.method("ComputeElementFlux", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, const mfem::FiniteElement & arg3, mfem::Vector & arg4)->void { a->ComputeElementFlux(arg0, arg1, arg2, arg3, arg4); });
    t.method("ComputeElementFlux", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, const mfem::FiniteElement & arg3, mfem::Vector & arg4, bool arg5)->void { a->ComputeElementFlux(arg0, arg1, arg2, arg3, arg4, arg5); });

    DEBUG_MSG("Adding wrapper for double mfem::BilinearFormIntegrator::ComputeFluxEnergy(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, mfem::Vector *) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::BilinearFormIntegrator::ComputeFluxEnergy(const mfem::FiniteElement &, mfem::ElementTransformation &, mfem::Vector &, mfem::Vector *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/bilininteg.hpp:242:19
    t.method("ComputeFluxEnergy", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->double { return a.ComputeFluxEnergy(arg0, arg1, arg2); });
    t.method("ComputeFluxEnergy", [](mfem::BilinearFormIntegrator& a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, mfem::Vector * arg3)->double { return a.ComputeFluxEnergy(arg0, arg1, arg2, arg3); });
    t.method("ComputeFluxEnergy", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2)->double { return a->ComputeFluxEnergy(arg0, arg1, arg2); });
    t.method("ComputeFluxEnergy", [](mfem::BilinearFormIntegrator* a, const mfem::FiniteElement & arg0, mfem::ElementTransformation & arg1, mfem::Vector & arg2, mfem::Vector * arg3)->double { return a->ComputeFluxEnergy(arg0, arg1, arg2, arg3); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::BilinearFormIntegrator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_BilinearFormIntegrator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_BilinearFormIntegrator(module));
}
