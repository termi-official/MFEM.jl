// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MatrixSumCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MatrixSumCoefficient> : std::false_type { };
template<> struct SuperType<mfem::MatrixSumCoefficient> { typedef mfem::MatrixCoefficient type; };
}

// Class generating the wrapper for type mfem::MatrixSumCoefficient
// signature to use in the veto file: mfem::MatrixSumCoefficient
struct Jlmfem_MatrixSumCoefficient: public Wrapper {

  Jlmfem_MatrixSumCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MatrixSumCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1810:7
    jlcxx::TypeWrapper<mfem::MatrixSumCoefficient>  t = jlModule.add_type<mfem::MatrixSumCoefficient>("mfem!MatrixSumCoefficient",
      jlcxx::julia_base_type<mfem::MatrixCoefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MatrixSumCoefficient>>(new jlcxx::TypeWrapper<mfem::MatrixSumCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::MatrixSumCoefficient(mfem::MatrixCoefficient &, mfem::MatrixCoefficient &, double, double) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1823:4
    t.constructor<mfem::MatrixCoefficient &, mfem::MatrixCoefficient &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::MatrixCoefficient &, mfem::MatrixCoefficient &, double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<mfem::MatrixCoefficient &, mfem::MatrixCoefficient &, double, double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixSumCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1827:9
    t.method("SetTime", [](mfem::MatrixSumCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::MatrixSumCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::SetACoef(mfem::MatrixCoefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixSumCoefficient::SetACoef(mfem::MatrixCoefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1830:9
    t.method("SetACoef", [](mfem::MatrixSumCoefficient& a, mfem::MatrixCoefficient & arg0)->void { a.SetACoef(arg0); });
    t.method("SetACoef", [](mfem::MatrixSumCoefficient* a, mfem::MatrixCoefficient & arg0)->void { a->SetACoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::MatrixCoefficient * mfem::MatrixSumCoefficient::GetACoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::MatrixCoefficient * mfem::MatrixSumCoefficient::GetACoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1832:24
    t.method("GetACoef", [](mfem::MatrixSumCoefficient const& a)->mfem::MatrixCoefficient * { return a.GetACoef(); });
    t.method("GetACoef", [](mfem::MatrixSumCoefficient const* a)->mfem::MatrixCoefficient * { return a->GetACoef(); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::SetBCoef(mfem::MatrixCoefficient &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixSumCoefficient::SetBCoef(mfem::MatrixCoefficient &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1835:9
    t.method("SetBCoef", [](mfem::MatrixSumCoefficient& a, mfem::MatrixCoefficient & arg0)->void { a.SetBCoef(arg0); });
    t.method("SetBCoef", [](mfem::MatrixSumCoefficient* a, mfem::MatrixCoefficient & arg0)->void { a->SetBCoef(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::MatrixCoefficient * mfem::MatrixSumCoefficient::GetBCoef() (" __HERE__ ")");
    // signature to use in the veto list: mfem::MatrixCoefficient * mfem::MatrixSumCoefficient::GetBCoef()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1837:24
    t.method("GetBCoef", [](mfem::MatrixSumCoefficient const& a)->mfem::MatrixCoefficient * { return a.GetBCoef(); });
    t.method("GetBCoef", [](mfem::MatrixSumCoefficient const* a)->mfem::MatrixCoefficient * { return a->GetBCoef(); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::SetAlpha(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixSumCoefficient::SetAlpha(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1840:9
    t.method("SetAlpha", [](mfem::MatrixSumCoefficient& a, double arg0)->void { a.SetAlpha(arg0); });
    t.method("SetAlpha", [](mfem::MatrixSumCoefficient* a, double arg0)->void { a->SetAlpha(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::MatrixSumCoefficient::GetAlpha() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::MatrixSumCoefficient::GetAlpha()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1842:11
    t.method("GetAlpha", [](mfem::MatrixSumCoefficient const& a)->double { return a.GetAlpha(); });
    t.method("GetAlpha", [](mfem::MatrixSumCoefficient const* a)->double { return a->GetAlpha(); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::SetBeta(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixSumCoefficient::SetBeta(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1845:9
    t.method("SetBeta", [](mfem::MatrixSumCoefficient& a, double arg0)->void { a.SetBeta(arg0); });
    t.method("SetBeta", [](mfem::MatrixSumCoefficient* a, double arg0)->void { a->SetBeta(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::MatrixSumCoefficient::GetBeta() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::MatrixSumCoefficient::GetBeta()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1847:11
    t.method("GetBeta", [](mfem::MatrixSumCoefficient const& a)->double { return a.GetBeta(); });
    t.method("GetBeta", [](mfem::MatrixSumCoefficient const* a)->double { return a->GetBeta(); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixSumCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixSumCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1850:17
    t.method("Eval", [](mfem::MatrixSumCoefficient& a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::MatrixSumCoefficient* a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a->Eval(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MatrixSumCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MatrixSumCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MatrixSumCoefficient(module));
}
