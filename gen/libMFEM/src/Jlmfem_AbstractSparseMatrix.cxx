// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::AbstractSparseMatrix> : std::false_type { };
  template<> struct DefaultConstructible<mfem::AbstractSparseMatrix> : std::false_type { };
template<> struct SuperType<mfem::AbstractSparseMatrix> { typedef mfem::Matrix type; };
}

// Class generating the wrapper for type mfem::AbstractSparseMatrix
// signature to use in the veto file: mfem::AbstractSparseMatrix
struct Jlmfem_AbstractSparseMatrix: public Wrapper {

  Jlmfem_AbstractSparseMatrix(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::AbstractSparseMatrix (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:73:7
    jlcxx::TypeWrapper<mfem::AbstractSparseMatrix>  t = jlModule.add_type<mfem::AbstractSparseMatrix>("mfem!AbstractSparseMatrix",
      jlcxx::julia_base_type<mfem::Matrix>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::AbstractSparseMatrix>>(new jlcxx::TypeWrapper<mfem::AbstractSparseMatrix>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;



    DEBUG_MSG("Adding wrapper for int mfem::AbstractSparseMatrix::NumNonZeroElems() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::AbstractSparseMatrix::NumNonZeroElems()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:83:16
    t.method("NumNonZeroElems", [](mfem::AbstractSparseMatrix const& a)->int { return a.NumNonZeroElems(); });
    t.method("NumNonZeroElems", [](mfem::AbstractSparseMatrix const* a)->int { return a->NumNonZeroElems(); });

    DEBUG_MSG("Adding wrapper for int mfem::AbstractSparseMatrix::GetRow(const int, mfem::Array<int> &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::AbstractSparseMatrix::GetRow(const int, mfem::Array<int> &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:89:16
    t.method("GetRow", [](mfem::AbstractSparseMatrix const& a, const int arg0, mfem::Array<int> & arg1, mfem::Vector & arg2)->int { return a.GetRow(arg0, arg1, arg2); });
    t.method("GetRow", [](mfem::AbstractSparseMatrix const* a, const int arg0, mfem::Array<int> & arg1, mfem::Vector & arg2)->int { return a->GetRow(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::AbstractSparseMatrix::EliminateZeroRows(const double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::AbstractSparseMatrix::EliminateZeroRows(const double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:96:17
    t.method("EliminateZeroRows", [](mfem::AbstractSparseMatrix& a)->void { a.EliminateZeroRows(); });
    t.method("EliminateZeroRows", [](mfem::AbstractSparseMatrix& a, const double arg0)->void { a.EliminateZeroRows(arg0); });
    t.method("EliminateZeroRows", [](mfem::AbstractSparseMatrix* a)->void { a->EliminateZeroRows(); });
    t.method("EliminateZeroRows", [](mfem::AbstractSparseMatrix* a, const double arg0)->void { a->EliminateZeroRows(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::AbstractSparseMatrix::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::AbstractSparseMatrix::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:99:17
    t.method("Mult", [](mfem::AbstractSparseMatrix const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::AbstractSparseMatrix const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::AbstractSparseMatrix::AddMult(const mfem::Vector &, mfem::Vector &, const double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::AbstractSparseMatrix::AddMult(const mfem::Vector &, mfem::Vector &, const double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:101:17
    t.method("AddMult", [](mfem::AbstractSparseMatrix const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMult(arg0, arg1); });
    t.method("AddMult", [](mfem::AbstractSparseMatrix const& a, const mfem::Vector & arg0, mfem::Vector & arg1, const double arg2)->void { a.AddMult(arg0, arg1, arg2); });
    t.method("AddMult", [](mfem::AbstractSparseMatrix const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMult(arg0, arg1); });
    t.method("AddMult", [](mfem::AbstractSparseMatrix const* a, const mfem::Vector & arg0, mfem::Vector & arg1, const double arg2)->void { a->AddMult(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::AbstractSparseMatrix::MultTranspose(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::AbstractSparseMatrix::MultTranspose(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:104:17
    t.method("MultTranspose", [](mfem::AbstractSparseMatrix const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.MultTranspose(arg0, arg1); });
    t.method("MultTranspose", [](mfem::AbstractSparseMatrix const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->MultTranspose(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::AbstractSparseMatrix::AddMultTranspose(const mfem::Vector &, mfem::Vector &, const double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::AbstractSparseMatrix::AddMultTranspose(const mfem::Vector &, mfem::Vector &, const double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/matrix.hpp:106:17
    t.method("AddMultTranspose", [](mfem::AbstractSparseMatrix const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.AddMultTranspose(arg0, arg1); });
    t.method("AddMultTranspose", [](mfem::AbstractSparseMatrix const& a, const mfem::Vector & arg0, mfem::Vector & arg1, const double arg2)->void { a.AddMultTranspose(arg0, arg1, arg2); });
    t.method("AddMultTranspose", [](mfem::AbstractSparseMatrix const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->AddMultTranspose(arg0, arg1); });
    t.method("AddMultTranspose", [](mfem::AbstractSparseMatrix const* a, const mfem::Vector & arg0, mfem::Vector & arg1, const double arg2)->void { a->AddMultTranspose(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::AbstractSparseMatrix>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_AbstractSparseMatrix(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_AbstractSparseMatrix(module));
}
