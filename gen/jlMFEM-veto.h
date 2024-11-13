// Auto veto breaks something... So here a manual list
// These are for Array<int>... but something is weird
Memory<T> & WrappedType::GetMemory()
const Memory<T> & WrappedType::GetMemory()

// AMR will require mork work
mfem::NCMesh
mfem::NCMesh::Master
mfem::NCMesh::Slave
mfem::NCMesh::MeshId
mfem::NCMesh::NCList

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
// wrapit bug?
mfem::Mesh::FaceInformation::(unnamed struct at /home/dogiermann/Builds/mfem-test/install/include/mfem/mesh/mesh.hpp:1862:7)
// "*[]" leads to faulty codegen
void mfem::NURBSPatchMap::SetPatchVertexMap(int, const KnotVector *[])
void mfem::NURBSPatchMap::SetPatchDofMap(int, const KnotVector *[])
void mfem::NURBSPatchMap::SetBdrPatchVertexMap(int, const KnotVector *[], int *)
void mfem::NURBSPatchMap::SetBdrPatchDofMap(int, const KnotVector *[], int *)
void mfem::NURBSExtension::MergeGridFunctions(GridFunction *[], int, mfem::GridFunction &)
void mfem::GridFunction::ProjectCoefficient(Coefficient *[])
void mfem::GridFunction::ProjectBdrCoefficient(Coefficient *[], mfem::Array<int> &)
double mfem::GridFunction::ComputeL2Error(Coefficient *[], const IntegrationRule *[], const mfem::Array<int> *)
double mfem::GridFunction::ComputeMaxError(Coefficient *[], const IntegrationRule *[])
double mfem::GridFunction::ComputeL1Error(Coefficient *[], const IntegrationRule *[])
double mfem::GridFunction::ComputeElementGradError(int, mfem::VectorCoefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeL2Error(mfem::Coefficient &, const IntegrationRule *[], const mfem::Array<int> *)
double mfem::GridFunction::ComputeL2Error(mfem::VectorCoefficient &, const IntegrationRule *[], const mfem::Array<int> *)
double mfem::GridFunction::ComputeGradError(mfem::VectorCoefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeCurlError(mfem::VectorCoefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeDivError(mfem::Coefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeDGFaceJumpError(mfem::Coefficient *, mfem::Coefficient *, mfem::JumpScaling, const IntegrationRule *[])
double mfem::GridFunction::ComputeDGFaceJumpError(mfem::Coefficient *, mfem::Coefficient *, double, const IntegrationRule *[])
double mfem::GridFunction::ComputeH1Error(mfem::Coefficient *, mfem::VectorCoefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeHDivError(mfem::VectorCoefficient *, mfem::Coefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeHCurlError(mfem::VectorCoefficient *, mfem::VectorCoefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeMaxError(mfem::Coefficient &, const IntegrationRule *[])
double mfem::GridFunction::ComputeMaxError(mfem::VectorCoefficient &, const IntegrationRule *[])
double mfem::GridFunction::ComputeL1Error(mfem::Coefficient &, const IntegrationRule *[])
double mfem::GridFunction::ComputeW11Error(mfem::Coefficient *, mfem::VectorCoefficient *, int, const mfem::Array<int> *, const IntegrationRule *[])
double mfem::GridFunction::ComputeL1Error(mfem::VectorCoefficient &, const IntegrationRule *[])
double mfem::GridFunction::ComputeLpError(const double, mfem::Coefficient &, mfem::Coefficient *, const IntegrationRule *[], const mfem::Array<int> *)
void mfem::GridFunction::ComputeElementLpErrors(const double, mfem::Coefficient &, mfem::Vector &, mfem::Coefficient *, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementL1Errors(mfem::Coefficient &, mfem::Vector &, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementL2Errors(mfem::Coefficient &, mfem::Vector &, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementMaxErrors(mfem::Coefficient &, mfem::Vector &, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementL1Errors(mfem::VectorCoefficient &, mfem::Vector &, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementL2Errors(mfem::VectorCoefficient &, mfem::Vector &, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementMaxErrors(mfem::VectorCoefficient &, mfem::Vector &, const IntegrationRule *[])
void mfem::GridFunction::ComputeElementLpErrors(const double, mfem::VectorCoefficient &, mfem::Vector &, mfem::Coefficient *, mfem::VectorCoefficient *, const IntegrationRule *[])
double mfem::GridFunction::ComputeLpError(const double, mfem::VectorCoefficient &, mfem::Coefficient *, mfem::VectorCoefficient *, const IntegrationRule *[])
// Ctor not correctly recognized
mfem::NCMesh::NCList::MeshIdAndType

// Template issues - these algorithms do not work on all types (e.g. Array<QuadraturePoint>)
int mfem::Array::Union(const T &)
int mfem::Array::Find(const T &)
void mfem::Array::Sort()
int mfem::Array::FindSorted(const T &)
void mfem::Array::DeleteFirst(const T &)
void mfem::Array::Unique()
T * mfem::Array::begin()
const T * mfem::Array::begin()
T * mfem::Array::end()
const T * mfem::Array::end()
T mfem::Array::Min()
T mfem::Array::Max()
int mfem::Array::IsSorted()
void mfem::Array::PartialSum()
T mfem::Array::Sum()
std::basic_ostream
std::basic_istream

// C function pointers
void mfem::Mesh::Transform(void (*)(const Vector &, Vector &))
