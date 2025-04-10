// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::OptimizationProblem> : std::false_type { };
  template<> struct DefaultConstructible<mfem::OptimizationProblem> : std::false_type { };
}

// Class generating the wrapper for type mfem::OptimizationProblem
// signature to use in the veto file: mfem::OptimizationProblem
struct Jlmfem_OptimizationProblem: public Wrapper {

  Jlmfem_OptimizationProblem(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::OptimizationProblem (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:785:7
    jlcxx::TypeWrapper<mfem::OptimizationProblem>  t = jlModule.add_type<mfem::OptimizationProblem>("mfem!OptimizationProblem");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::OptimizationProblem>>(new jlcxx::TypeWrapper<mfem::OptimizationProblem>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for double mfem::OptimizationProblem::CalcObjective(const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::OptimizationProblem::CalcObjective(const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:799:19
    t.method("CalcObjective", [](mfem::OptimizationProblem const& a, const mfem::Vector & arg0)->double { return a.CalcObjective(arg0); });
    t.method("CalcObjective", [](mfem::OptimizationProblem const* a, const mfem::Vector & arg0)->double { return a->CalcObjective(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::OptimizationProblem::CalcObjectiveGrad(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OptimizationProblem::CalcObjectiveGrad(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:801:17
    t.method("CalcObjectiveGrad", [](mfem::OptimizationProblem const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.CalcObjectiveGrad(arg0, arg1); });
    t.method("CalcObjectiveGrad", [](mfem::OptimizationProblem const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->CalcObjectiveGrad(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::OptimizationProblem::SetEqualityConstraint(const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OptimizationProblem::SetEqualityConstraint(const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:804:9
    t.method("SetEqualityConstraint", [](mfem::OptimizationProblem& a, const mfem::Vector & arg0)->void { a.SetEqualityConstraint(arg0); });
    t.method("SetEqualityConstraint", [](mfem::OptimizationProblem* a, const mfem::Vector & arg0)->void { a->SetEqualityConstraint(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::OptimizationProblem::SetInequalityConstraint(const mfem::Vector &, const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OptimizationProblem::SetInequalityConstraint(const mfem::Vector &, const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:805:9
    t.method("SetInequalityConstraint", [](mfem::OptimizationProblem& a, const mfem::Vector & arg0, const mfem::Vector & arg1)->void { a.SetInequalityConstraint(arg0, arg1); });
    t.method("SetInequalityConstraint", [](mfem::OptimizationProblem* a, const mfem::Vector & arg0, const mfem::Vector & arg1)->void { a->SetInequalityConstraint(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::OptimizationProblem::SetSolutionBounds(const mfem::Vector &, const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OptimizationProblem::SetSolutionBounds(const mfem::Vector &, const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:806:9
    t.method("SetSolutionBounds", [](mfem::OptimizationProblem& a, const mfem::Vector & arg0, const mfem::Vector & arg1)->void { a.SetSolutionBounds(arg0, arg1); });
    t.method("SetSolutionBounds", [](mfem::OptimizationProblem* a, const mfem::Vector & arg0, const mfem::Vector & arg1)->void { a->SetSolutionBounds(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for const mfem::Operator * mfem::OptimizationProblem::GetC() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Operator * mfem::OptimizationProblem::GetC()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:808:20
    t.method("GetC", [](mfem::OptimizationProblem const& a)->const mfem::Operator * { return a.GetC(); });
    t.method("GetC", [](mfem::OptimizationProblem const* a)->const mfem::Operator * { return a->GetC(); });

    DEBUG_MSG("Adding wrapper for const mfem::Operator * mfem::OptimizationProblem::GetD() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Operator * mfem::OptimizationProblem::GetD()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:809:20
    t.method("GetD", [](mfem::OptimizationProblem const& a)->const mfem::Operator * { return a.GetD(); });
    t.method("GetD", [](mfem::OptimizationProblem const* a)->const mfem::Operator * { return a->GetD(); });

    DEBUG_MSG("Adding wrapper for const mfem::Vector * mfem::OptimizationProblem::GetEqualityVec() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Vector * mfem::OptimizationProblem::GetEqualityVec()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:810:18
    t.method("GetEqualityVec", [](mfem::OptimizationProblem const& a)->const mfem::Vector * { return a.GetEqualityVec(); });
    t.method("GetEqualityVec", [](mfem::OptimizationProblem const* a)->const mfem::Vector * { return a->GetEqualityVec(); });

    DEBUG_MSG("Adding wrapper for const mfem::Vector * mfem::OptimizationProblem::GetInequalityVec_Lo() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Vector * mfem::OptimizationProblem::GetInequalityVec_Lo()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:811:18
    t.method("GetInequalityVec_Lo", [](mfem::OptimizationProblem const& a)->const mfem::Vector * { return a.GetInequalityVec_Lo(); });
    t.method("GetInequalityVec_Lo", [](mfem::OptimizationProblem const* a)->const mfem::Vector * { return a->GetInequalityVec_Lo(); });

    DEBUG_MSG("Adding wrapper for const mfem::Vector * mfem::OptimizationProblem::GetInequalityVec_Hi() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Vector * mfem::OptimizationProblem::GetInequalityVec_Hi()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:812:18
    t.method("GetInequalityVec_Hi", [](mfem::OptimizationProblem const& a)->const mfem::Vector * { return a.GetInequalityVec_Hi(); });
    t.method("GetInequalityVec_Hi", [](mfem::OptimizationProblem const* a)->const mfem::Vector * { return a->GetInequalityVec_Hi(); });

    DEBUG_MSG("Adding wrapper for const mfem::Vector * mfem::OptimizationProblem::GetBoundsVec_Lo() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Vector * mfem::OptimizationProblem::GetBoundsVec_Lo()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:813:18
    t.method("GetBoundsVec_Lo", [](mfem::OptimizationProblem const& a)->const mfem::Vector * { return a.GetBoundsVec_Lo(); });
    t.method("GetBoundsVec_Lo", [](mfem::OptimizationProblem const* a)->const mfem::Vector * { return a->GetBoundsVec_Lo(); });

    DEBUG_MSG("Adding wrapper for const mfem::Vector * mfem::OptimizationProblem::GetBoundsVec_Hi() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::Vector * mfem::OptimizationProblem::GetBoundsVec_Hi()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:814:18
    t.method("GetBoundsVec_Hi", [](mfem::OptimizationProblem const& a)->const mfem::Vector * { return a.GetBoundsVec_Hi(); });
    t.method("GetBoundsVec_Hi", [](mfem::OptimizationProblem const* a)->const mfem::Vector * { return a->GetBoundsVec_Hi(); });

    DEBUG_MSG("Adding wrapper for int mfem::OptimizationProblem::GetNumConstraints() (" __HERE__ ")");
    // signature to use in the veto list: int mfem::OptimizationProblem::GetNumConstraints()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:816:8
    t.method("GetNumConstraints", [](mfem::OptimizationProblem const& a)->int { return a.GetNumConstraints(); });
    t.method("GetNumConstraints", [](mfem::OptimizationProblem const* a)->int { return a->GetNumConstraints(); });

    DEBUG_MSG("Adding input_size methods  to provide read access to the field input_size (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:793:14
    // signature to use in the veto list: mfem::OptimizationProblem::input_size
    t.method("input_size", [](const mfem::OptimizationProblem& a) -> int { return a.input_size; });
    t.method("input_size", [](const mfem::OptimizationProblem* a) -> int { return a->input_size; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::OptimizationProblem>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_OptimizationProblem(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_OptimizationProblem(module));
}
