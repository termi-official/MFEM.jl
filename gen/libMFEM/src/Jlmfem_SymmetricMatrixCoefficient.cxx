// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::SymmetricMatrixCoefficient> : std::false_type { };
  template<> struct DefaultConstructible<mfem::SymmetricMatrixCoefficient> : std::false_type { };
}

// Class generating the wrapper for type mfem::SymmetricMatrixCoefficient
// signature to use in the veto file: mfem::SymmetricMatrixCoefficient
struct Jlmfem_SymmetricMatrixCoefficient: public Wrapper {

  Jlmfem_SymmetricMatrixCoefficient(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::SymmetricMatrixCoefficient (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1250:7
    jlcxx::TypeWrapper<mfem::SymmetricMatrixCoefficient>  t = jlModule.add_type<mfem::SymmetricMatrixCoefficient>("mfem!SymmetricMatrixCoefficient");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::SymmetricMatrixCoefficient>>(new jlcxx::TypeWrapper<mfem::SymmetricMatrixCoefficient>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void mfem::SymmetricMatrixCoefficient::SetTime(double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::SymmetricMatrixCoefficient::SetTime(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1262:17
    t.method("SetTime", [](mfem::SymmetricMatrixCoefficient& a, double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::SymmetricMatrixCoefficient* a, double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for double mfem::SymmetricMatrixCoefficient::GetTime() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::SymmetricMatrixCoefficient::GetTime()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1265:11
    t.method("GetTime", [](mfem::SymmetricMatrixCoefficient& a)->double { return a.GetTime(); });
    t.method("GetTime", [](mfem::SymmetricMatrixCoefficient* a)->double { return a->GetTime(); });

    DEBUG_MSG("Adding wrapper for int mfem::SymmetricMatrixCoefficient::GetSize() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::SymmetricMatrixCoefficient::GetSize()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1268:8
    t.method("GetSize", [](mfem::SymmetricMatrixCoefficient const& a)->int { return a.GetSize(); });
    t.method("GetSize", [](mfem::SymmetricMatrixCoefficient const* a)->int { return a->GetSize(); });

    DEBUG_MSG("Adding wrapper for void mfem::SymmetricMatrixCoefficient::Eval(mfem::DenseSymmetricMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::SymmetricMatrixCoefficient::Eval(mfem::DenseSymmetricMatrix &, mfem::ElementTransformation &, const mfem::IntegrationPoint &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/coefficient.hpp:1275:17
    t.method("Eval", [](mfem::SymmetricMatrixCoefficient& a, mfem::DenseSymmetricMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a.Eval(arg0, arg1, arg2); });
    t.method("Eval", [](mfem::SymmetricMatrixCoefficient* a, mfem::DenseSymmetricMatrix & arg0, mfem::ElementTransformation & arg1, const mfem::IntegrationPoint & arg2)->void { a->Eval(arg0, arg1, arg2); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::SymmetricMatrixCoefficient>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_SymmetricMatrixCoefficient(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_SymmetricMatrixCoefficient(module));
}
