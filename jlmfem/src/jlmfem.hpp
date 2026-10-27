#pragma once
#include "mfem.hpp"

// MFEM subclasses forwarding their virtual methods to Julia. Each class stores a C function
// pointer (created with @cfunction on the Julia side) and an opaque context pointer, which
// is passed back as first argument and refers to the Julia object implementing the method.
namespace jlmfem
{

class JuliaCoefficient : public mfem::Coefficient
{
   using Callback = mfem::real_t (*)(void *, mfem::ElementTransformation &,
                                     const mfem::IntegrationPoint &);
   Callback f;
   void *ctx;

public:
   JuliaCoefficient(void *f, void *ctx) : f(reinterpret_cast<Callback>(f)), ctx(ctx) {}

   mfem::real_t Eval(mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override
   {
      return f(ctx, T, ip);
   }
};

class JuliaFunctionCoefficient : public mfem::Coefficient
{
   using Callback = mfem::real_t (*)(void *, const mfem::real_t *, int, mfem::real_t);
   Callback f;
   void *ctx;
   mfem::Vector x;

public:
   JuliaFunctionCoefficient(void *f, void *ctx) : f(reinterpret_cast<Callback>(f)), ctx(ctx) {}

   mfem::real_t Eval(mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override
   {
      T.Transform(ip, x);
      return f(ctx, x.GetData(), x.Size(), GetTime());
   }
};

class JuliaVectorFunctionCoefficient : public mfem::VectorCoefficient
{
   using Callback = void (*)(void *, mfem::real_t *, int, const mfem::real_t *, int, mfem::real_t);
   Callback f;
   void *ctx;
   mfem::Vector x;

public:
   JuliaVectorFunctionCoefficient(int vdim, void *f, void *ctx)
      : mfem::VectorCoefficient(vdim), f(reinterpret_cast<Callback>(f)), ctx(ctx) {}

   void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
             const mfem::IntegrationPoint &ip) override
   {
      T.Transform(ip, x);
      V.SetSize(vdim);
      f(ctx, V.GetData(), vdim, x.GetData(), x.Size(), GetTime());
   }
};

class JuliaBilinearFormIntegrator : public mfem::BilinearFormIntegrator
{
   using Callback = void (*)(void *, const mfem::FiniteElement &,
                             mfem::ElementTransformation &, mfem::DenseMatrix &);
   Callback f;
   void *ctx;

public:
   JuliaBilinearFormIntegrator(void *f, void *ctx)
      : f(reinterpret_cast<Callback>(f)), ctx(ctx) {}

   void AssembleElementMatrix(const mfem::FiniteElement &el,
                              mfem::ElementTransformation &T,
                              mfem::DenseMatrix &elmat) override
   {
      f(ctx, el, T, elmat);
   }
};

class JuliaLinearFormIntegrator : public mfem::LinearFormIntegrator
{
   using Callback = void (*)(void *, const mfem::FiniteElement &,
                             mfem::ElementTransformation &, mfem::Vector &);
   Callback f;
   void *ctx;

public:
   JuliaLinearFormIntegrator(void *f, void *ctx)
      : f(reinterpret_cast<Callback>(f)), ctx(ctx) {}

   using mfem::LinearFormIntegrator::AssembleRHSElementVect;
   void AssembleRHSElementVect(const mfem::FiniteElement &el,
                               mfem::ElementTransformation &T,
                               mfem::Vector &elvect) override
   {
      f(ctx, el, T, elvect);
   }
};

// Instances of member templates, which WrapIt cannot wrap. They are exposed as methods of
// the corresponding Julia functions (see julia_names in jlmfem.wit.in).
inline void FormLinearSystem(mfem::BilinearForm &a, const mfem::Array<int> &ess_tdof_list,
                             mfem::Vector &x, mfem::Vector &b, mfem::SparseMatrix &A,
                             mfem::Vector &X, mfem::Vector &B, int copy_interior = 0)
{
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior);
}

inline void FormSystemMatrix(mfem::BilinearForm &a, const mfem::Array<int> &ess_tdof_list,
                             mfem::SparseMatrix &A)
{
   a.FormSystemMatrix(ess_tdof_list, A);
}

#ifdef MFEM_USE_MPI
inline void FormLinearSystem(mfem::ParBilinearForm &a, const mfem::Array<int> &ess_tdof_list,
                             mfem::Vector &x, mfem::Vector &b, mfem::HypreParMatrix &A,
                             mfem::Vector &X, mfem::Vector &B, int copy_interior = 0)
{
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior);
}

inline void FormSystemMatrix(mfem::ParBilinearForm &a, const mfem::Array<int> &ess_tdof_list,
                             mfem::HypreParMatrix &A)
{
   a.FormSystemMatrix(ess_tdof_list, A);
}
#endif

// Versions without the "const IntegrationRule *irs[]" argument, which WrapIt cannot wrap.
inline mfem::real_t ComputeL2Error(const mfem::GridFunction &u, mfem::Coefficient &exsol)
{
   return u.ComputeL2Error(exsol);
}

inline mfem::real_t ComputeL2Error(const mfem::GridFunction &u, mfem::VectorCoefficient &exsol)
{
   return u.ComputeL2Error(exsol);
}

// The global variables are vetoed, as their generated setters do not compile.
inline mfem::IntegrationRules &GetIntRules() { return mfem::IntRules; }
inline mfem::IntegrationRules &GetRefinedIntRules() { return mfem::RefinedIntRules; }

// Counterparts of BilinearForm::UseExternalIntegrators, which these forms lack.
inline void DisownIntegrators(mfem::LinearForm &form)
{
   struct Access : mfem::LinearForm
   {
      static int mfem::LinearForm::*flag() { return &Access::extern_lfs; }
   };
   form.*Access::flag() = 1;
}

inline void DisownIntegrators(mfem::MixedBilinearForm &form)
{
   struct Access : mfem::MixedBilinearForm
   {
      static int mfem::MixedBilinearForm::*flag() { return &Access::extern_bfs; }
   };
   form.*Access::flag() = 1;
}

} // namespace jlmfem
