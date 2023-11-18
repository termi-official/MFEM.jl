// These must be excluded because they use std::function . We provide custom wrappers for "Julia coefficients"
std::function
mfem::FunctionCoefficient
mfem::VectorFunctionCoefficient
mfem::MatrixFunctionCoefficient
mfem::SymmetricMatrixFunctionCoefficient
// wrapper generation fails and likely should not be exposed anyway
mfem::BlockArray::iterator
mfem::BlockArray::const_iterator
mfem::BlockArray

// I think these are accidentally exposed in the C++ API
mfem::Mesh::ncmesh
// NOTE: int mfem::CheckFinite(const double *, const int) is generated twice
int mfem::CheckFinite(const double *, const int)

// Deactivate these extensions for now
void mfem::Mesh::MesquiteSmooth(const int)
mfem::Mesh::NURBSext

// I don't know how to deal with these
mfem::CoefficientStorage mfem::operator|(mfem::CoefficientStorage, mfem::CoefficientStorage)
int mfem::operator&(mfem::CoefficientStorage, mfem::CoefficientStorage)

// Potential bugs
// Not sure what happened here. Codegen looks faulty
int mfem::Ordering::Map(int, int, int, int)
// mfem!IntRules! -> Don't know how to veto. The setter should not be generated because the copy ctor for the type is deleted.
// wrapper for void mfem::Mesh::Mesh(int, int, int, int, int) -> Not even sure what fails here. Maybe a conflict with the string ctor, which is modeled as "int"?
mfem::GeometricFactors
mfem::Mesh::geom_factors
const mfem::GeometricFactors * mfem::Mesh::GetGeometricFactors(const mfem::IntegrationRule &, const int, mfem::MemoryType)
void mfem::Mesh::DeleteGeometricFactors()
mfem::FaceGeometricFactors
mfem::Mesh::face_geom_factors
const mfem::FaceGeometricFactors * mfem::Mesh::GetFaceGeometricFactors(const mfem::IntegrationRule &, const int, mfem::FaceType, mfem::MemoryType)
