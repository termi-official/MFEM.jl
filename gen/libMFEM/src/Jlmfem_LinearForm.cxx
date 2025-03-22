// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::LinearForm> : std::false_type { };
  template<> struct DefaultConstructible<mfem::LinearForm> : std::false_type { };
template<> struct SuperType<mfem::LinearForm> { typedef mfem::Vector type; };
}

// Class generating the wrapper for type mfem::LinearForm
// signature to use in the veto file: mfem::LinearForm
struct Jlmfem_LinearForm: public Wrapper {

  Jlmfem_LinearForm(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::LinearForm (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:23:7
    jlcxx::TypeWrapper<mfem::LinearForm>  t = jlModule.add_type<mfem::LinearForm>("mfem!LinearForm",
      jlcxx::julia_base_type<mfem::Vector>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::LinearForm>>(new jlcxx::TypeWrapper<mfem::LinearForm>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::LinearForm(mfem::FiniteElementSpace *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:78:4
    t.constructor<mfem::FiniteElementSpace *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::LinearForm(mfem::FiniteElementSpace *, mfem::LinearForm *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:88:4
    t.constructor<mfem::FiniteElementSpace *, mfem::LinearForm *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::LinearForm(mfem::FiniteElementSpace *, double *) (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:101:4
    t.constructor<mfem::FiniteElementSpace *, double *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for mfem::LinearForm & mfem::LinearForm::operator=(const mfem::LinearForm &) (" __HERE__ ")");
    // signature to use in the veto list: mfem::LinearForm & mfem::LinearForm::operator=(const mfem::LinearForm &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:110:16
    t.method("assign", [](mfem::LinearForm& a, const mfem::LinearForm & arg0)->mfem::LinearForm & { return a.operator=(arg0); });
    t.method("assign", [](mfem::LinearForm* a, const mfem::LinearForm & arg0)->mfem::LinearForm & { return a->operator=(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::FiniteElementSpace * mfem::LinearForm::FESpace() (" __HERE__ ")");
    // signature to use in the veto list: mfem::FiniteElementSpace * mfem::LinearForm::FESpace()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:118:24
    t.method("FESpace", [](mfem::LinearForm& a)->mfem::FiniteElementSpace * { return a.FESpace(); });
    t.method("FESpace", [](mfem::LinearForm* a)->mfem::FiniteElementSpace * { return a->FESpace(); });

    DEBUG_MSG("Adding wrapper for const mfem::FiniteElementSpace * mfem::LinearForm::FESpace() (" __HERE__ ")");
    // signature to use in the veto list: const mfem::FiniteElementSpace * mfem::LinearForm::FESpace()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:120:30
    t.method("FESpace", [](mfem::LinearForm const& a)->const mfem::FiniteElementSpace * { return a.FESpace(); });
    t.method("FESpace", [](mfem::LinearForm const* a)->const mfem::FiniteElementSpace * { return a->FESpace(); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddDomainIntegrator(mfem::LinearFormIntegrator *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddDomainIntegrator(mfem::LinearFormIntegrator *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:123:9
    t.method("AddDomainIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0)->void { a.AddDomainIntegrator(arg0); });
    t.method("AddDomainIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0)->void { a->AddDomainIntegrator(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddDomainIntegrator(mfem::LinearFormIntegrator *, mfem::Array<int> &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddDomainIntegrator(mfem::LinearFormIntegrator *, mfem::Array<int> &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:126:9
    t.method("AddDomainIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0, mfem::Array<int> & arg1)->void { a.AddDomainIntegrator(arg0, arg1); });
    t.method("AddDomainIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0, mfem::Array<int> & arg1)->void { a->AddDomainIntegrator(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddBoundaryIntegrator(mfem::LinearFormIntegrator *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddBoundaryIntegrator(mfem::LinearFormIntegrator *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:130:9
    t.method("AddBoundaryIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0)->void { a.AddBoundaryIntegrator(arg0); });
    t.method("AddBoundaryIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0)->void { a->AddBoundaryIntegrator(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddBoundaryIntegrator(mfem::LinearFormIntegrator *, mfem::Array<int> &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddBoundaryIntegrator(mfem::LinearFormIntegrator *, mfem::Array<int> &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:137:9
    t.method("AddBoundaryIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0, mfem::Array<int> & arg1)->void { a.AddBoundaryIntegrator(arg0, arg1); });
    t.method("AddBoundaryIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0, mfem::Array<int> & arg1)->void { a->AddBoundaryIntegrator(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddBdrFaceIntegrator(mfem::LinearFormIntegrator *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddBdrFaceIntegrator(mfem::LinearFormIntegrator *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:141:9
    t.method("AddBdrFaceIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0)->void { a.AddBdrFaceIntegrator(arg0); });
    t.method("AddBdrFaceIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0)->void { a->AddBdrFaceIntegrator(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddBdrFaceIntegrator(mfem::LinearFormIntegrator *, mfem::Array<int> &) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddBdrFaceIntegrator(mfem::LinearFormIntegrator *, mfem::Array<int> &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:148:9
    t.method("AddBdrFaceIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0, mfem::Array<int> & arg1)->void { a.AddBdrFaceIntegrator(arg0, arg1); });
    t.method("AddBdrFaceIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0, mfem::Array<int> & arg1)->void { a->AddBdrFaceIntegrator(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AddInteriorFaceIntegrator(mfem::LinearFormIntegrator *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AddInteriorFaceIntegrator(mfem::LinearFormIntegrator *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:152:9
    t.method("AddInteriorFaceIntegrator", [](mfem::LinearForm& a, mfem::LinearFormIntegrator * arg0)->void { a.AddInteriorFaceIntegrator(arg0); });
    t.method("AddInteriorFaceIntegrator", [](mfem::LinearForm* a, mfem::LinearFormIntegrator * arg0)->void { a->AddInteriorFaceIntegrator(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetDLFI() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetDLFI()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:157:34
    t.method("GetDLFI", [](mfem::LinearForm& a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a.GetDLFI(); });
    t.method("GetDLFI", [](mfem::LinearForm* a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a->GetDLFI(); });

    DEBUG_MSG("Adding wrapper for mfem::Array<mfem::DeltaLFIntegrator *> * mfem::LinearForm::GetDLFI_Delta() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Array<mfem::DeltaLFIntegrator *> * mfem::LinearForm::GetDLFI_Delta()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:161:31
    t.method("GetDLFI_Delta", [](mfem::LinearForm& a)->mfem::Array<mfem::DeltaLFIntegrator *> * { return a.GetDLFI_Delta(); });
    t.method("GetDLFI_Delta", [](mfem::LinearForm* a)->mfem::Array<mfem::DeltaLFIntegrator *> * { return a->GetDLFI_Delta(); });

    DEBUG_MSG("Adding wrapper for mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetBLFI() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetBLFI()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:164:34
    t.method("GetBLFI", [](mfem::LinearForm& a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a.GetBLFI(); });
    t.method("GetBLFI", [](mfem::LinearForm* a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a->GetBLFI(); });

    DEBUG_MSG("Adding wrapper for mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetFLFI() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetFLFI()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:167:34
    t.method("GetFLFI", [](mfem::LinearForm& a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a.GetFLFI(); });
    t.method("GetFLFI", [](mfem::LinearForm* a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a->GetFLFI(); });

    DEBUG_MSG("Adding wrapper for mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetIFLFI() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Array<mfem::LinearFormIntegrator *> * mfem::LinearForm::GetIFLFI()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:170:34
    t.method("GetIFLFI", [](mfem::LinearForm& a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a.GetIFLFI(); });
    t.method("GetIFLFI", [](mfem::LinearForm* a)->mfem::Array<mfem::LinearFormIntegrator *> * { return a->GetIFLFI(); });

    DEBUG_MSG("Adding wrapper for mfem::Array<mfem::Array<int> *> * mfem::LinearForm::GetFLFI_Marker() (" __HERE__ ")");
    // signature to use in the veto list: mfem::Array<mfem::Array<int> *> * mfem::LinearForm::GetFLFI_Marker()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:175:24
    t.method("GetFLFI_Marker", [](mfem::LinearForm& a)->mfem::Array<mfem::Array<int> *> * { return a.GetFLFI_Marker(); });
    t.method("GetFLFI_Marker", [](mfem::LinearForm* a)->mfem::Array<mfem::Array<int> *> * { return a->GetFLFI_Marker(); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::Assemble() (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::Assemble()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:178:9
    t.method("Assemble", [](mfem::LinearForm& a)->void { a.Assemble(); });
    t.method("Assemble", [](mfem::LinearForm* a)->void { a->Assemble(); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::AssembleDelta() (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::AssembleDelta()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:181:9
    t.method("AssembleDelta", [](mfem::LinearForm& a)->void { a.AssembleDelta(); });
    t.method("AssembleDelta", [](mfem::LinearForm* a)->void { a->AssembleDelta(); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::Update() (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::Update()
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:188:9
    t.method("Update", [](mfem::LinearForm& a)->void { a.Update(); });
    t.method("Update", [](mfem::LinearForm* a)->void { a->Update(); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::Update(mfem::FiniteElementSpace *) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::Update(mfem::FiniteElementSpace *)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:191:9
    t.method("Update", [](mfem::LinearForm& a, mfem::FiniteElementSpace * arg0)->void { a.Update(arg0); });
    t.method("Update", [](mfem::LinearForm* a, mfem::FiniteElementSpace * arg0)->void { a->Update(arg0); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::Update(mfem::FiniteElementSpace *, mfem::Vector &, int) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::Update(mfem::FiniteElementSpace *, mfem::Vector &, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:198:9
    t.method("Update", [](mfem::LinearForm& a, mfem::FiniteElementSpace * arg0, mfem::Vector & arg1, int arg2)->void { a.Update(arg0, arg1, arg2); });
    t.method("Update", [](mfem::LinearForm* a, mfem::FiniteElementSpace * arg0, mfem::Vector & arg1, int arg2)->void { a->Update(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void mfem::LinearForm::MakeRef(mfem::FiniteElementSpace *, mfem::Vector &, int) (" __HERE__ ")");
    // signature to use in the veto list: void mfem::LinearForm::MakeRef(mfem::FiniteElementSpace *, mfem::Vector &, int)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:208:17
    t.method("MakeRef", [](mfem::LinearForm& a, mfem::FiniteElementSpace * arg0, mfem::Vector & arg1, int arg2)->void { a.MakeRef(arg0, arg1, arg2); });
    t.method("MakeRef", [](mfem::LinearForm* a, mfem::FiniteElementSpace * arg0, mfem::Vector & arg1, int arg2)->void { a->MakeRef(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for double mfem::LinearForm::operator()(const mfem::GridFunction &) (" __HERE__ ")");
    // signature to use in the veto list: double mfem::LinearForm::operator()(const mfem::GridFunction &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:215:11
    t.method("paren", [](mfem::LinearForm const& a, const mfem::GridFunction & arg0)->double { return a.operator()(arg0); });
    t.method("paren", [](mfem::LinearForm const* a, const mfem::GridFunction & arg0)->double { return a->operator()(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::LinearForm & mfem::LinearForm::operator=(double) (" __HERE__ ")");
    // signature to use in the veto list: mfem::LinearForm & mfem::LinearForm::operator=(double)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:218:16
    t.method("assign", [](mfem::LinearForm& a, double arg0)->mfem::LinearForm & { return a.operator=(arg0); });
    t.method("assign", [](mfem::LinearForm* a, double arg0)->mfem::LinearForm & { return a->operator=(arg0); });

    DEBUG_MSG("Adding wrapper for mfem::LinearForm & mfem::LinearForm::operator=(const mfem::Vector &) (" __HERE__ ")");
    // signature to use in the veto list: mfem::LinearForm & mfem::LinearForm::operator=(const mfem::Vector &)
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/mesh/../fem/linearform.hpp:223:16
    t.method("assign", [](mfem::LinearForm& a, const mfem::Vector & arg0)->mfem::LinearForm & { return a.operator=(arg0); });
    t.method("assign", [](mfem::LinearForm* a, const mfem::Vector & arg0)->mfem::LinearForm & { return a->operator=(arg0); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::LinearForm>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_LinearForm(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_LinearForm(module));
}
