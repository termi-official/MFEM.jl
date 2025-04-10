// this file was auto-generated by wrapit v1.5.0
#include "Wrapper.h"

#include "jllibMFEM.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<mfem::DenseSymmetricMatrix> : std::false_type { };
  template<> struct DefaultConstructible<mfem::DenseSymmetricMatrix> : std::false_type { };
template<> struct SuperType<mfem::DenseSymmetricMatrix> { typedef mfem::Matrix type; };
}

// Class generating the wrapper for type mfem::DenseSymmetricMatrix
// signature to use in the veto file: mfem::DenseSymmetricMatrix
struct Jlmfem_DenseSymmetricMatrix: public Wrapper {

  Jlmfem_DenseSymmetricMatrix(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type mfem::DenseSymmetricMatrix (" __HERE__ ")");
    // defined in /home/dogiermann/.julia/artifacts/820df874853553756f46ac6dc23173c05d8db01a/include/mfem/linalg/symmat.hpp:24:7
    jlcxx::TypeWrapper<mfem::DenseSymmetricMatrix>  t = jlModule.add_type<mfem::DenseSymmetricMatrix>("mfem!DenseSymmetricMatrix",
      jlcxx::julia_base_type<mfem::Matrix>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<mfem::DenseSymmetricMatrix>>(new jlcxx::TypeWrapper<mfem::DenseSymmetricMatrix>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<mfem::DenseSymmetricMatrix>> type_;
};
std::shared_ptr<Wrapper> newJlmfem_DenseSymmetricMatrix(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlmfem_DenseSymmetricMatrix(module));
}
