// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::MatrixArrayCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::MatrixArrayCoefficient> : std::false_type { };
template<> struct SuperType<mfem::MatrixArrayCoefficient> { typedef mfem::MatrixCoefficient type; };
}

// Class generating the wrapper for type mfem::MatrixArrayCoefficient
// signature to use in the veto file: mfem::MatrixArrayCoefficient
struct Jlmfem_MatrixArrayCoefficient: public Wrapper {

  Jlmfem_MatrixArrayCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::MatrixArrayCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1127:7
    jlcxx::TypeWrapper<mfem::MatrixArrayCoefficient>  t = jlModule.add_type<mfem::MatrixArrayCoefficient>("mfem!MatrixArrayCoefficient",
      jlcxx::julia_base_type<mfem::MatrixCoefficient>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::MatrixArrayCoefficient>>(new jlcxx::TypeWrapper<mfem::MatrixArrayCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::MatrixArrayCoefficient::MatrixArrayCoefficient(int) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1136:13
    t.constructor<int>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::MatrixArrayCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixArrayCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1139:9
    t.method("SetTime", [](mfem::MatrixArrayCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::MatrixArrayCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::Coefficient * mfem::MatrixArrayCoefficient::GetCoeff(int, int) (" __HERE__ ")");
    // signature to use in the veto list: mfem::Coefficient * mfem::MatrixArrayCoefficient::GetCoeff(int, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1142:17
    t.method("GetCoeff", [](mfem::MatrixArrayCoefficient& a, int arg0, int arg1)->mfem::Coefficient * { return a.GetCoeff(arg0, arg1); });
    t.method("GetCoeff", [](mfem::MatrixArrayCoefficient* a, int arg0, int arg1)->mfem::Coefficient * { return a->GetCoeff(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixArrayCoefficient::Set(int, int, mfem::Coefficient *, bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixArrayCoefficient::Set(int, int, mfem::Coefficient *, bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1147:9
    t.method("Set", [](mfem::MatrixArrayCoefficient& a, int arg0, int arg1, mfem::Coefficient * arg2)->void { a.Set(arg0, arg1, arg2); });
    t.method("Set", [](mfem::MatrixArrayCoefficient& a, int arg0, int arg1, mfem::Coefficient * arg2, bool arg3)->void { a.Set(arg0, arg1, arg2, arg3); });
    t.method("Set", [](mfem::MatrixArrayCoefficient* a, int arg0, int arg1, mfem::Coefficient * arg2)->void { a->Set(arg0, arg1, arg2); });
    t.method("Set", [](mfem::MatrixArrayCoefficient* a, int arg0, int arg1, mfem::Coefficient * arg2, bool arg3)->void { a->Set(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for double mfem::MatrixArrayCoefficient::Eval(int, int, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::MatrixArrayCoefficient::Eval(int, int, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1151:11
    t.method("Eval", [](mfem::MatrixArrayCoefficient& a, int arg0, int arg1, mfem::ElementTransformation & arg2, const mfem::IntegrationPoint & arg3)->double { return a.Eval(arg0, arg1, arg2, arg3); });
    t.method("Eval", [](mfem::MatrixArrayCoefficient* a, int arg0, int arg1, mfem::ElementTransformation & arg2, const mfem::IntegrationPoint & arg3)->double { return a->Eval(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for void mfem::MatrixArrayCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::MatrixArrayCoefficient::Eval(mfem::DenseMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1155:17
    t.method("Eval", [](mfem::MatrixArrayCoefficient& a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::MatrixArrayCoefficient* a, mfem::DenseMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a->Eval(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::MatrixArrayCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_MatrixArrayCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_MatrixArrayCoefficient(module));
}
