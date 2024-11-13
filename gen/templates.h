#pragma once
// Here we instantiate the templates used in MFEM which should be wrapped to Julia
// #include "mfem/general/array.hpp"
// #include "mfem/fem/intrules.hpp"
// #include "mfem/mesh/mesh.hpp"
#include "mfem/mfem.hpp"
namespace mfem 
{
    // template class Array<int>;
    // template class Array<IntegrationPoint>;
}

// template class mfem::Array<mfem::GeometricFactors *>;
// template class mfem::Array<mfem::FaceGeometricFactors *>;

// template<class T>
// using Array = mfem::Array<T>;
