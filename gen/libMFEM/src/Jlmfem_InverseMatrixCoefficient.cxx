// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::InverseMatrixCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::InverseMatrixCoefficient> : std::false_type { };
template<> struct SuperType<mfem::InverseMatrixCoefficient> { typedef mfem::MatrixCoefficient type; };
}

// Class generating the wrapper for type mfem::InverseMatrixCoefficient
// signature to use in the veto file: mfem::InverseMatrixCoefficient
struct Jlmfem_InverseMatrixCoefficient: public Wrapper {

  Jlmfem_InverseMatrixCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::InverseMatrixCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1946:7
    jlcxx::TypeWrapper<mfem::InverseMatrixCoefficient>  t = jlModule.add_type<mfem::InverseMatrixCoefficient>("mfem!InverseMatrixCoefficient",
      jlcxx::julia_base_type<mfem::MatrixCoefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::InverseMatrixCoefficient>>(new jlcxx::TypeWrapper<mfem::InverseMatrixCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::InverseMatrixCoefficient::InverseMatrixCoefficient(mfem::MatrixCoefficient &) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1953:4
    t.constructor<mfem::MatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::InverseMatrixCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::InverseMatrixCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1956:9
    t.method("SetTime", [](mfem::InverseMatrixCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::InverseMatrixCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::InverseMatrixCoefficient::SetACoef(mfem::MatrixCoefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::InverseMatrixCoefficient::SetACoef(mfem::MatrixCoefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1959:9
    t.method("SetACoef", [](mfem::InverseMatrixCoefficient& a, mfem::MatrixCoefficient & arg0)->void { a.SetACoef(arg0); });
    t.method("SetACoef", [](mfem::InverseMatrixCoefficient* a, mfem::MatrixCoefficient & arg0)->void { a->SetACoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::MatrixCoefficient * mfem::InverseMatrixCoefficient::GetACoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::MatrixCoefficient * mfem::InverseMatrixCoefficient::GetACoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1961:24
    t.method("GetACoef", [](mfem::InverseMatrixCoefficient const& a)->mfem::MatrixCoefficient * { return a.GetACoef(); });
    t.method("GetACoef", [](mfem::InverseMatrixCoefficient const* a)->mfem::MatrixCoefficient * { return a->GetACoef(); });

    DEBUG_MSG("Adding wrapper for void mfem::InverseMatrixCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::InverseMatrixCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1964:17
    t.method("Eval", [](mfem::InverseMatrixCoefficient& a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::InverseMatrixCoefficient* a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a->Eval(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::InverseMatrixCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_InverseMatrixCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_InverseMatrixCoefficient(module));
}
