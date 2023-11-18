// Here we instantiate the templates used in MFEM which should be wrapped to Julia
#include "mfem/general/array.hpp"
#include "mfem/mesh/mesh.hpp"
template class mfem::Array<int>;
// template class mfem::Array<mfem::GeometricFactors *>;
// template class mfem::Array<mfem::FaceGeometricFactors *>;