// Arrays of pointers ("T *[]") lead to faulty code generation.
/.*[*][[].*/

// std::function based coefficients are replaced by the jlmfem::Julia* classes.
std::function
mfem::FunctionCoefficient
mfem::VectorFunctionCoefficient
mfem::MatrixFunctionCoefficient
mfem::SymmetricMatrixFunctionCoefficient
mfem::TransformedCoefficient

// Streams are not wrapped. Use the file name based methods (e.g. Save).
std::basic_ostream
std::basic_istream
std::map

// C function pointers
/.*[(][*][)][(].*/

// Non-conforming meshes (AMR) and NURBS require more work.
/mfem::NCMesh.*/
mfem::Mesh::ncmesh
mfem::Mesh::NURBSext
/mfem::BlockArray.*/

// Device/geometric factor internals
mfem::GeometricFactors
mfem::FaceGeometricFactors
mfem::Mesh::geom_factors
mfem::Mesh::face_geom_factors

// Algorithms of the Array template which do not compile for all element types.
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
Memory<T> & mfem::Array::GetMemory()
const Memory<T> & mfem::Array::GetMemory()
mfem::BlockArray::iterator
mfem::BlockArray::const_iterator
T mfem::Array2D::Min()
T mfem::Array2D::Max()

// Members whose generated setters fail (deleted assignment / unique_ptr / anonymous struct)
mfem::Mesh::attribute_sets
mfem::Mesh::bdr_attribute_sets
/mfem::MeshPart.*/
/mfem::Mesh::FaceInformation.*/
mfem::ParMesh::face_nbr_elements

// Global integration rules: use jlmfem::GetIntRules() and jlmfem::GetRefinedIntRules().
mfem::IntRules
mfem::RefinedIntRules

// Internal kernel dispatch tables (they return function pointers)
/.*Kernels.*/
/.*KernelSignature.*/
/.*::kernel_name.*/
/mfem::InverseElementTransformation.*/
/mfem::BatchInverseElementTransformation.*/
/mfem::Poly_1D.*/
mfem::poly1d

// Template instances of internal containers. They would be registered before their element types.
/.*mfem::Array<[^<>]*[*]>.*/
/.*mfem::Array<mfem::Array<int> [*]>.*/
/.*mfem::Array2D<.*/
/.*mfem::Array<mfem::Geometry::Type>.*/
/.*mfem::Array<mfem::Refinement>.*/
/.*std::unordered_map.*/

// Methods using types that are not wrapped
int mfem::Mesh::FindPoints(mfem::DenseMatrix &, mfem::Array<int> &, mfem::Array<mfem::IntegrationPoint> &, bool, mfem::InverseElementTransformation *)
int mfem::ParMesh::FindPoints(mfem::DenseMatrix &, mfem::Array<int> &, mfem::Array<mfem::IntegrationPoint> &, bool, mfem::InverseElementTransformation *)
void mfem::IterativeSolver::SetMonitor(mfem::IterativeSolverMonitor &)
std::shared_ptr<const mfem::PRefinementTransferOperator> mfem::FiniteElementSpace::GetPrefUpdateOperator()

// Enums are Integers in Julia, so these overloads clash with int ones. Only host memory is used.
/.*mfem::MemoryType.*/

// Overloads that clash in Julia, where enums are Integers. The enum versions are kept, except
// for the deprecated Cartesian Mesh constructors (use Mesh!MakeCartesian2D/3D).
mfem::GSSmoother::GSSmoother(int, int)
mfem::GSSmoother::GSSmoother(const mfem::SparseMatrix &, int, int)
mfem::DSmoother::DSmoother(int, mfem::real_t, int)
mfem::DSmoother::DSmoother(const mfem::SparseMatrix &, int, mfem::real_t, int)
void mfem::GridFunction::ProjectCoefficient(mfem::VectorCoefficient &, int)
mfem::Mesh::Mesh(int, int, mfem::Element::Type, bool, mfem::real_t, mfem::real_t, bool)
mfem::Mesh::Mesh(int, int, int, mfem::Element::Type, bool, mfem::real_t, mfem::real_t, mfem::real_t, bool)

// Operators on enums break Base functions relying on Integer semantics.
mfem::CoefficientStorage ::mfem::operator|(mfem::CoefficientStorage, mfem::CoefficientStorage)
int ::mfem::operator&(mfem::CoefficientStorage, mfem::CoefficientStorage)
