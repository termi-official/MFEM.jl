// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::OperatorJacobiSmoother> : std::false_type { };
  template<> struct DefaultConstructible<mfem::OperatorJacobiSmoother> : std::false_type { };
template<> struct SuperType<mfem::OperatorJacobiSmoother> { typedef mfem::Solver type; };
}

// Class generating the wrapper for type mfem::OperatorJacobiSmoother
// signature to use in the veto file: mfem::OperatorJacobiSmoother
struct Jlmfem_OperatorJacobiSmoother: public Wrapper {

  Jlmfem_OperatorJacobiSmoother(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::OperatorJacobiSmoother (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:274:7
    jlcxx::TypeWrapper<mfem::OperatorJacobiSmoother>  t = jlModule.add_type<mfem::OperatorJacobiSmoother>("mfem!OperatorJacobiSmoother",
      jlcxx::julia_base_type<mfem::Solver>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::OperatorJacobiSmoother>>(new jlcxx::TypeWrapper<mfem::OperatorJacobiSmoother>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::OperatorJacobiSmoother(const double) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:280:4
    t.constructor<const double>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::OperatorJacobiSmoother(const mfem::BilinearForm &, const mfem::Array<int> &, const double) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:291:4
    t.constructor<const mfem::BilinearForm &, const mfem::Array<int> &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const mfem::BilinearForm &, const mfem::Array<int> &, const double>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::OperatorJacobiSmoother(const mfem::Vector &, const mfem::Array<int> &, const double) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:304:4
    t.constructor<const mfem::Vector &, const mfem::Array<int> &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const mfem::Vector &, const mfem::Array<int> &, const double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::SetPositiveDiagonal(bool) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OperatorJacobiSmoother::SetPositiveDiagonal(bool)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:311:9
    t.method("SetPositiveDiagonal", [](mfem::OperatorJacobiSmoother& a)->void { a.SetPositiveDiagonal(); });
    t.method("SetPositiveDiagonal", [](mfem::OperatorJacobiSmoother& a, bool arg0)->void { a.SetPositiveDiagonal(arg0); });
    t.method("SetPositiveDiagonal", [](mfem::OperatorJacobiSmoother* a)->void { a->SetPositiveDiagonal(); });
    t.method("SetPositiveDiagonal", [](mfem::OperatorJacobiSmoother* a, bool arg0)->void { a->SetPositiveDiagonal(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::Mult(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OperatorJacobiSmoother::Mult(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:313:9
    t.method("Mult", [](mfem::OperatorJacobiSmoother const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.Mult(arg0, arg1); });
    t.method("Mult", [](mfem::OperatorJacobiSmoother const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->Mult(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::MultTranspose(const mfem::Vector &, mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OperatorJacobiSmoother::MultTranspose(const mfem::Vector &, mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:314:9
    t.method("MultTranspose", [](mfem::OperatorJacobiSmoother const& a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a.MultTranspose(arg0, arg1); });
    t.method("MultTranspose", [](mfem::OperatorJacobiSmoother const* a, const mfem::Vector & arg0, mfem::Vector & arg1)->void { a->MultTranspose(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::SetOperator(const mfem::Operator &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OperatorJacobiSmoother::SetOperator(const mfem::Operator &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:327:9
    t.method("SetOperator", [](mfem::OperatorJacobiSmoother& a, const mfem::Operator & arg0)->void { a.SetOperator(arg0); });
    t.method("SetOperator", [](mfem::OperatorJacobiSmoother* a, const mfem::Operator & arg0)->void { a->SetOperator(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::OperatorJacobiSmoother::Setup(const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::OperatorJacobiSmoother::Setup(const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/solvers.hpp:344:9
    t.method("Setup", [](mfem::OperatorJacobiSmoother& a, const mfem::Vector & arg0)->void { a.Setup(arg0); });
    t.method("Setup", [](mfem::OperatorJacobiSmoother* a, const mfem::Vector & arg0)->void { a->Setup(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::OperatorJacobiSmoother>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_OperatorJacobiSmoother(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_OperatorJacobiSmoother(module));
}
