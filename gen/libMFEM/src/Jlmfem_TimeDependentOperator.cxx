// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::TimeDependentOperator> : std::false_type { };
  template<> struct DefaultConstructible<mfem::TimeDependentOperator> : std::false_type { };
template<> struct SuperType<mfem::TimeDependentOperator> { typedef mfem::Operator type; };
}

// Class generating the wrapper for type mfem::TimeDependentOperator
// signature to use in the veto file: mfem::TimeDependentOperator
struct Jlmfem_TimeDependentOperator: public Wrapper {

  Jlmfem_TimeDependentOperator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::TimeDependentOperator (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:285:7
    jlcxx::TypeWrapper<mfem::TimeDependentOperator>  t = jlModule.add_type<mfem::TimeDependentOperator>("mfem!TimeDependentOperator",
      jlcxx::julia_base_type<mfem::Operator>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::TimeDependentOperator>>(new jlcxx::TypeWrapper<mfem::TimeDependentOperator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::TimeDependentOperator(int, double, mfem::TimeDependentOperator::Type) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:316:13
    t.constructor<int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, double, mfem::TimeDependentOperator::Type>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::TimeDependentOperator(int, int, double, mfem::TimeDependentOperator::Type) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:322:4
    t.constructor<int, int>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, int, double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<int, int, double, mfem::TimeDependentOperator::Type>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for double mfem::TimeDependentOperator::GetTime() (" __HERE__ ")");
    // signature to use in the veto list: double mfem::TimeDependentOperator::GetTime()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:326:19
    t.method("GetTime", [](mfem::TimeDependentOperator const& a)->double { return a.GetTime(); });
    t.method("GetTime", [](mfem::TimeDependentOperator const* a)->double { return a->GetTime(); });

    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::SetTime(const double) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TimeDependentOperator::SetTime(const double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:329:17
    t.method("SetTime", [](mfem::TimeDependentOperator& a, const double arg0)->void { a.SetTime(arg0); });
    t.method("SetTime", [](mfem::TimeDependentOperator* a, const double arg0)->void { a->SetTime(arg0); });

    DEBUG_MSG("Adding wrapper for bool mfem::TimeDependentOperator::isExplicit() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::TimeDependentOperator::isExplicit()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:332:9
    t.method("isExplicit", [](mfem::TimeDependentOperator const& a)->bool { return a.isExplicit(); });
    t.method("isExplicit", [](mfem::TimeDependentOperator const* a)->bool { return a->isExplicit(); });

    DEBUG_MSG("Adding wrapper for bool mfem::TimeDependentOperator::isImplicit() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::TimeDependentOperator::isImplicit()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:334:9
    t.method("isImplicit", [](mfem::TimeDependentOperator const& a)->bool { return a.isImplicit(); });
    t.method("isImplicit", [](mfem::TimeDependentOperator const* a)->bool { return a->isImplicit(); });

    DEBUG_MSG("Adding wrapper for bool mfem::TimeDependentOperator::isHomogeneous() (" __HERE__ ")");
    // signature to use in the veto list: bool mfem::TimeDependentOperator::isHomogeneous()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:336:9
    t.method("isHomogeneous", [](mfem::TimeDependentOperator const& a)->bool { return a.isHomogeneous(); });
    t.method("isHomogeneous", [](mfem::TimeDependentOperator const* a)->bool { return a->isHomogeneous(); });

    DEBUG_MSG("Adding wrapper for mfem::TimeDependentOperator::EvalMode mfem::TimeDependentOperator::GetEvalMode() (" __HERE__ ")");
    // signature to use in the veto list: mfem::TimeDependentOperator::EvalMode mfem::TimeDependentOperator::GetEvalMode()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:339:13
    t.method("GetEvalMode", [](mfem::TimeDependentOperator const& a)->mfem::TimeDependentOperator::EvalMode { return a.GetEvalMode(); });
    t.method("GetEvalMode", [](mfem::TimeDependentOperator const* a)->mfem::TimeDependentOperator::EvalMode { return a->GetEvalMode(); });

    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::SetEvalMode(const mfem::TimeDependentOperator::EvalMode) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TimeDependentOperator::SetEvalMode(const mfem::TimeDependentOperator::EvalMode)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:352:17
    t.method("SetEvalMode", [](mfem::TimeDependentOperator& a, const mfem::TimeDependentOperator::EvalMode arg0)->void { a.SetEvalMode(arg0); });
    t.method("SetEvalMode", [](mfem::TimeDependentOperator* a, const mfem::TimeDependentOperator::EvalMode arg0)->void { a->SetEvalMode(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::ExplicitMult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TimeDependentOperator::ExplicitMult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:360:17
    t.method("ExplicitMult", [](mfem::TimeDependentOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.ExplicitMult(arg0, arg1); });
    t.method("ExplicitMult", [](mfem::TimeDependentOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->ExplicitMult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::ImplicitMult(const mfem::Vector &, const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TimeDependentOperator::ImplicitMult(const mfem::Vector &, const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:367:17
    t.method("ImplicitMult", [](mfem::TimeDependentOperator const& a, const mfem::Vector & arg0, const mfem::Vector & arg1, mfem::Vector & arg2)->void { a.ImplicitMult(arg0, arg1, arg2); });
    t.method("ImplicitMult", [](mfem::TimeDependentOperator const* a, const mfem::Vector & arg0, const mfem::Vector & arg1, mfem::Vector & arg2)->void { a->ImplicitMult(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TimeDependentOperator::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:372:17
    t.method("Mult", [](mfem::TimeDependentOperator const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::TimeDependentOperator const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::TimeDependentOperator::ImplicitSolve(const double, const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::TimeDependentOperator::ImplicitSolve(const double, const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:390:17
    t.method("ImplicitSolve", [](mfem::TimeDependentOperator& a, const double arg0, const mfem::Vector & arg1, mfem::Vector & arg2)->void { a.ImplicitSolve(arg0, arg1, arg2); });
    t.method("ImplicitSolve", [](mfem::TimeDependentOperator* a, const double arg0, const mfem::Vector & arg1, mfem::Vector & arg2)->void { a->ImplicitSolve(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for mfem::Operator & mfem::TimeDependentOperator::GetImplicitGradient(const mfem::Vector &, const mfem::Vector &, double) (" __HERE__ ")");
    // signature to use in the veto list: mfem::Operator & mfem::TimeDependentOperator::GetImplicitGradient(const mfem::Vector &, const mfem::Vector &, double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:397:22
    t.method("GetImplicitGradient", [](mfem::TimeDependentOperator const& a, const mfem::Vector & arg0, const mfem::Vector & arg1, double arg2)->mfem::Operator & { return a.GetImplicitGradient(arg0, arg1, arg2); });
    t.method("GetImplicitGradient", [](mfem::TimeDependentOperator const* a, const mfem::Vector & arg0, const mfem::Vector & arg1, double arg2)->mfem::Operator & { return a->GetImplicitGradient(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for mfem::Operator & mfem::TimeDependentOperator::GetExplicitGradient(const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: mfem::Operator & mfem::TimeDependentOperator::GetExplicitGradient(const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:405:22
    t.method("GetExplicitGradient", [](mfem::TimeDependentOperator const& a, const mfem::Vector & arg0)->mfem::Operator & { return a.GetExplicitGradient(arg0); });
    t.method("GetExplicitGradient", [](mfem::TimeDependentOperator const* a, const mfem::Vector & arg0)->mfem::Operator & { return a->GetExplicitGradient(arg0); });

    DEBUG_MSG("Adding wrapper for int mfem::TimeDependentOperator::SUNImplicitSetup(const mfem::Vector &, const mfem::Vector &, int, int *, double) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::TimeDependentOperator::SUNImplicitSetup(const mfem::Vector &, const mfem::Vector &, int, int *, double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:420:16
    t.method("SUNImplicitSetup", [](mfem::TimeDependentOperator& a, const mfem::Vector & arg0, const mfem::Vector & arg1, int arg2, int * arg3, double arg4)->int { return a.SUNImplicitSetup(arg0, arg1, arg2, arg3, arg4); });
    t.method("SUNImplicitSetup", [](mfem::TimeDependentOperator* a, const mfem::Vector & arg0, const mfem::Vector & arg1, int arg2, int * arg3, double arg4)->int { return a->SUNImplicitSetup(arg0, arg1, arg2, arg3, arg4); });

    DEBUG_MSG("Adding wrapper for int mfem::TimeDependentOperator::SUNImplicitSolve(const mfem::Vector &, mfem::Vector &, double) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::TimeDependentOperator::SUNImplicitSolve(const mfem::Vector &, mfem::Vector &, double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:434:16
    t.method("SUNImplicitSolve", [](mfem::TimeDependentOperator& a, const mfem::Vector & arg0, mfem::Vector & arg1, double arg2)->int { return a.SUNImplicitSolve(arg0, arg1, arg2); });
    t.method("SUNImplicitSolve", [](mfem::TimeDependentOperator* a, const mfem::Vector & arg0, mfem::Vector & arg1, double arg2)->int { return a->SUNImplicitSolve(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for int mfem::TimeDependentOperator::SUNMassSetup() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::TimeDependentOperator::SUNMassSetup()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:442:16
    t.method("SUNMassSetup", [](mfem::TimeDependentOperator& a)->int { return a.SUNMassSetup(); });
    t.method("SUNMassSetup", [](mfem::TimeDependentOperator* a)->int { return a->SUNMassSetup(); });

    DEBUG_MSG("Adding wrapper for int mfem::TimeDependentOperator::SUNMassSolve(const mfem::Vector &, mfem::Vector &, double) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::TimeDependentOperator::SUNMassSolve(const mfem::Vector &, mfem::Vector &, double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:455:16
    t.method("SUNMassSolve", [](mfem::TimeDependentOperator& a, const mfem::Vector & arg0, mfem::Vector & arg1, double arg2)->int { return a.SUNMassSolve(arg0, arg1, arg2); });
    t.method("SUNMassSolve", [](mfem::TimeDependentOperator* a, const mfem::Vector & arg0, mfem::Vector & arg1, double arg2)->int { return a->SUNMassSolve(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for int mfem::TimeDependentOperator::SUNMassMult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: int mfem::TimeDependentOperator::SUNMassMult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/operator.hpp:466:16
    t.method("SUNMassMult", [](mfem::TimeDependentOperator& a, const mfem::Vector & arg0, mfem::Vector & arg1)->int { return a.SUNMassMult(arg0, arg1); });
    t.method("SUNMassMult", [](mfem::TimeDependentOperator* a, const mfem::Vector & arg0, mfem::Vector & arg1)->int { return a->SUNMassMult(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::TimeDependentOperator>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_TimeDependentOperator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_TimeDependentOperator(module));
}
