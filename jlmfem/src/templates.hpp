#pragma once
#include "mfem.hpp"

// Template specializations to be wrapped. They are already instantiated in libmfem.
#ifdef WRAPIT
template class mfem::Array<int>;
template class mfem::Array<mfem::real_t>;
#endif
