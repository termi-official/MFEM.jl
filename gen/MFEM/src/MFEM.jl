module MFEM

export AddBdrElement, AddBdrPoint, AddBdrQuad, AddBdrQuadAsTriangles, AddBdrSegment, AddBdrTriangle, AddElement, AddHex
export AddHexAsPyramids, AddHexAsTets, AddHexAsWedges, AddPyramid, AddQuad, AddQuadAs4TrisWithPoints, AddQuadAs5QuadsWithPoints
export AddSegment, AddTet, AddTri, AddTriangle, AddVertex, AddVertexAtMeanCenter, AddVertexParents, AddWedge, Append, Assign
export BuildDofToArrays, CalcD2Shape, CalcDShape, CalcDnShape, CalcShape, Capacity, CartesianPartitioning, ChangeVertexDataOwnership
export CheckBdrElementOrientation, CheckDisplacements, CheckElementOrientation, CheckPartitioning, Clear, Clone, ComputeFlux, ComputeH1Error
export Conforming, ConnectBoundaries, ConvertFromConformingVDofs, ConvertToConformingVDofs, ConvertToPatches, CountElementsPerVDof
export CreatePeriodicVertexMapping, CurlDim, D2C_GlobalRestrictionMatrix, D2Const_GlobalRestrictionMatrix, DegreeElevate, DeleteAll, DeleteLast
export DerefineByError, Derivatives, Determinants, Difference, Dimension, DisableTensorProducts, DofForGeometry, DofOrderForOrientation
export DofToVDof, DofTransformationForGeometry, DofsToVDofs, ElementToEdgeTable, ElementToElementTable, ElementToFaceTable
export EnableTensorProducts, EnsureNCMesh, EnsureNodes, EulerNumber, EulerNumber2D, Eval, FEColl, FESpace, FaceIsInterior, Finalize
export FinalizeHexMesh, FinalizeMesh, FinalizeQuadMesh, FinalizeTetMesh, FinalizeTopology, FinalizeTriMesh, FinalizeWedgeMesh
export FindFaceNeighbors, FindInterpolant, FindMaxima, FindPoints, FiniteElementForDim, FiniteElementForGeometry, Flip, FlipDirection
export GeneralRefinement, GenerateBoundaryElements, GeneratePartitioning, GetAttribute, GetBE, GetBasisType, GetBdrAttribute
export GetBdrElement, GetBdrElementAdjacentElement, GetBdrElementAdjacentElement2, GetBdrElementBaseGeometry, GetBdrElementData
export GetBdrElementDofTable, GetBdrElementDofs, GetBdrElementEdgeIndex, GetBdrElementEdges, GetBdrElementFace, GetBdrElementFaceIndex
export GetBdrElementGeometry, GetBdrElementToDofTable, GetBdrElementTopo, GetBdrElementTransformation, GetBdrElementType, GetBdrElementVDofs
export GetBdrElementVertices, GetBdrFace, GetBdrFaceTransformations, GetBdrPatchKnotVectors, GetBdrPointMatrix, GetBdrValuesFrom
export GetBoundaryTrueDofs, GetBoundingBox, GetCharacteristics, GetClosedBasisType, GetCoarseToFineMap, GetConformingProlongation
export GetConformingRestriction, GetConformingVSize, GetContType, GetCurl, GetData, GetDerivMapType, GetDerivRangeType, GetDerivType
export GetDerivative, GetDivergence, GetDofMap, GetDofOrdering, GetEdgeDofs, GetEdgeElement, GetEdgeInteriorDofs
export GetEdgeInteriorVDofs, GetEdgeOrder, GetEdgeTransformation, GetEdgeVDofs, GetEdgeVertexTable, GetEdgeVertices, GetElement
export GetElementAverages, GetElementBaseGeometry, GetElementCenter, GetElementColoring, GetElementData, GetElementDofTable
export GetElementDofValues, GetElementDofs, GetElementEdges, GetElementFaces, GetElementForDof, GetElementGeometry, GetElementIJK
export GetElementIntRule, GetElementInteriorDofs, GetElementInteriorVDofs, GetElementJacobian, GetElementLocalToGlobal, GetElementOrder
export GetElementPatch, GetElementRestriction, GetElementSize, GetElementToDofTable, GetElementToFaceOrientationTable, GetElementTopo
export GetElementTransformation, GetElementType, GetElementVDofs, GetElementVertices, GetElementVolume, GetElements, GetElementsArray
export GetEntityIndex, GetEssentialTrueDofs, GetEssentialVDofs, GetFE, GetFace, GetFaceBaseGeometry, GetFaceDofs, GetFaceEdgeTable
export GetFaceEdges, GetFaceElement, GetFaceElementTransformations, GetFaceElementType, GetFaceElements, GetFaceGeometry
export GetFaceGeometryType, GetFaceInformation, GetFaceInfos, GetFaceIntRule, GetFaceInteriorDofs, GetFaceOrder
export GetFaceQuadratureInterpolator, GetFaceRestriction, GetFaceToBdrElMap, GetFaceToDofTable, GetFaceToElementTable, GetFaceTransformation
export GetFaceType, GetFaceVDofs, GetFaceValues, GetFaceVectorValues, GetFaceVertices, GetGNBE, GetGNE, GetGNV
export GetGeckoElementOrdering, GetGeometricParametersFromJacobian, GetGeometries, GetGeometry, GetGlobalNE, GetGradient, GetGradients
export GetHessians, GetHilbertElementOrdering, GetHpConformingRestriction, GetHpRestrictionMatrix, GetIntRule
export GetInteriorFaceTransformations, GetKV, GetKnotVector, GetLaplacians, GetLastOperation, GetLocalDofForDof, GetLocalFaceTransformation
export GetMapType, GetMaster, GetMaxElementOrder, GetMesh, GetNBE, GetNBP, GetNC, GetNCP, GetNConformingDofs, GetNDof, GetNDofs
export GetNE, GetNEDofs, GetNEdges, GetNF, GetNFDofs, GetNFaces, GetNFbyType, GetNKS, GetNKV, GetNP, GetNTotalDof
export GetNURBSext, GetNV, GetNVDofs, GetNodalFESpace, GetNodalValues, GetNode, GetNodes, GetNumDof, GetNumElementInteriorDofs
export GetNumFaces, GetNumFacesWithGhost, GetNumGeometries, GetOpenBasisType, GetOrder, GetOrdering, GetOrders, GetOutputLayout
export GetPatchAttribute, GetPatchBdrAttribute, GetPatchBdrElements, GetPatchDofs, GetPatchElements, GetPatchKnotVectors, GetPatchVDofs
export GetPermutedIndex, GetPointMatrix, GetProlongationMatrix, GetQuadratureInterpolator, GetRangeDim, GetRangeType
export GetRefinementTransforms, GetRestrictionMatrix, GetRestrictionOperator, GetRestrictionTransposeOperator, GetSequence, GetSize, GetSlave
export GetSpace, GetTraceCollection, GetTraceElement, GetTransferOperator, GetTransformation, GetTrueDofs
export GetTrueTransferOperator, GetTrueVSize, GetTrueVector, GetUpdateOperator, GetVDim, GetVDofs, GetVSize, GetValue, GetValues, GetValuesFrom
export GetVectorFieldNodalValues, GetVectorFieldValues, GetVectorGradient, GetVectorGradientHat, GetVectorValue, GetVectorValues, GetVertex
export GetVertexDofs, GetVertexLocalToGlobal, GetVertexToElementTable, GetVertexToVertexTable, GetVertexVDofs, GetVertices
export GetWeights, H2L_GlobalRestrictionMatrix, HasBoundaryElements, HasFaceDofs, HasGeometry, HavePatches, HostRead
export HostReadWrite, HostWrite, ImposeBounds, IsBoundary, IsConforming, IsDGSpace, IsInitialized, IsInterior, IsLocal
export IsNonconformingCoarse, IsNonconformingFine, IsOfFaceType, IsShared, IsVariableOrder, KnotInsert, Last, LoadBE, LoadFE, LoseData
export MakeCoarseToFineTable, MakeDataOwner, MakeOwner, MakeRef, MakeTRef, MakeUniformDegree, MemoryUsage, MeshGenerator, MoveNodes
export MoveVertices, Mult, MultTranspose, Name, NewElement, NewNodes, NodesUpdated, Nonconforming, OwnFEC, OwnsData, OwnsNodes
export OwnsSpace, PhysDerivatives, Prepend, PrintBdrVTU, PrintFunctions, PrintVTU, ProjectBdrCoefficient
export ProjectBdrCoefficientNormal, ProjectBdrCoefficientTangent, ProjectCoefficient, ProjectDiscCoefficient, ProjectGridFunction
export ProjectVectorFieldOn, RandomRefinement, Read, ReadWrite, RebuildElementToDofTable, ReduceInt, RefineAtVertex, RefineByError
export RefineNURBSFromFile, RemoveInternalBoundaries, RemoveUnusedVertices, ReorderByNodes, ReorderElementToDofTable, ReorderElements
export ReorientTetMesh, Reserve, Reset, RestrictConforming, Rotate, Rotate2D, Rotate3D, Save, SaveVTU, ScaleElements, ScaleSubdomains
export SetAttribute, SetAttributes, SetBdrAttribute, SetCoordsFromPatches, SetCurvature, SetElementOrder, SetFromTrueDofs
export SetFromTrueVector, SetKnotsFromPatches, SetLayer, SetNodalFESpace, SetNodalGridFunction, SetNode, SetNodes, SetNodesOwner
export SetOrder, SetOutputLayout, SetOwnsSpace, SetPatchAttribute, SetPatchBdrAttribute, SetRelaxedHpConformity, SetSize
export SetSpace, SetTrueVector, SetUpdateOperatorOwner, SetUpdateOperatorType, SetVDim, SetVertices, Size, SpaceDimension
export StealData, StealNURBSext, SubDofOrder, Swap, SwapDirections, SwapNodes, TraceFiniteElementForGeometry, Transform
export UniformRefinement, Update, UpdatesFinished, UseDevice, UsesTensorProducts, VDofToDof, Values, VectorDim, Write, assign, attributes
export attributes!, bdr_attributes, bdr_attributes!, embeddings, embeddings!, findKnotSpan, geom, geom!, getKnotLocation, ghost
export ghost!, index, index!, isElement, matrix, matrix!, mfem!Array, mfem!BilinearFormIntegrator, mfem!BoundingBox
export mfem!CoarseFineTransformations, mfem!Coefficient, mfem!ComputeElementLpDistance, mfem!Const2DFECollection, mfem!Const3DFECollection
export mfem!CrouzeixRaviartFECollection, mfem!CubicDiscont2DFECollection, mfem!CubicFECollection, mfem!DG_Interface_FECollection, mfem!DSTable
export mfem!DenseMatrix, mfem!DofTransformation, mfem!Element, mfem!Element!HEXAHEDRON, mfem!Element!POINT, mfem!Element!PYRAMID
export mfem!Element!QUADRILATERAL, mfem!Element!SEGMENT, mfem!Element!TETRAHEDRON, mfem!Element!TRIANGLE, mfem!Element!Type, mfem!Element!WEDGE
export mfem!ElementDofOrdering, mfem!ElementDofOrdering!LEXICOGRAPHIC, mfem!ElementDofOrdering!NATIVE, mfem!ElementRestrictionOperator
export mfem!ElementTransformation, mfem!Embedding, mfem!Extrude1D, mfem!Extrude1DGridFunction, mfem!Extrude2D, mfem!ExtrudeCoefficient
export mfem!FaceElementTransformations, mfem!FaceQuadratureInterpolator, mfem!FaceQuadratureInterpolator!DERIVATIVES
export mfem!FaceQuadratureInterpolator!DETERMINANTS, mfem!FaceQuadratureInterpolator!FaceEvalFlags, mfem!FaceQuadratureInterpolator!NORMALS
export mfem!FaceQuadratureInterpolator!VALUES, mfem!FaceQuadratureSpace, mfem!FaceRestriction, mfem!FaceType, mfem!FaceType!Boundary, mfem!FaceType!Interior
export mfem!FiniteElement, mfem!FiniteElementCollection, mfem!FiniteElementCollection!CONTINUOUS
export mfem!FiniteElementCollection!DISCONTINUOUS, mfem!FiniteElementCollection!NORMAL, mfem!FiniteElementCollection!New, mfem!FiniteElementCollection!TANGENTIAL
export mfem!FiniteElementSpace, mfem!FiniteElementSpace!AdjustVDofs, mfem!FiniteElementSpace!DecodeDof, mfem!FiniteElementSpace!EncodeDof
export mfem!FiniteElementSpace!ListToMarker, mfem!FiniteElementSpace!MarkerToList, mfem!GaussLinearDiscont2DFECollection
export mfem!GaussQuadraticDiscont2DFECollection, mfem!Geometry!CUBE, mfem!Geometry!INVALID, mfem!Geometry!NUM_GEOMETRIES, mfem!Geometry!POINT
export mfem!Geometry!PRISM, mfem!Geometry!PYRAMID, mfem!Geometry!SEGMENT, mfem!Geometry!SQUARE, mfem!Geometry!TETRAHEDRON
export mfem!Geometry!TRIANGLE, mfem!Geometry!Type, mfem!GetEVectorOrdering, mfem!GridFunction, mfem!GridFunction!ARITHMETIC
export mfem!GridFunction!AvgType, mfem!GridFunction!HARMONIC, mfem!H1Pos_FECollection, mfem!H1Ser_FECollection, mfem!H1_FECollection
export mfem!H1_Trace_FECollection, mfem!IntegrationPoint, mfem!IntegrationRule, mfem!InverseElementTransformation
export mfem!IsoparametricTransformation, mfem!JumpScaling, mfem!JumpScaling!CONSTANT, mfem!JumpScaling!JumpScalingType, mfem!JumpScaling!ONE_OVER_H
export mfem!JumpScaling!P_SQUARED_OVER_H, mfem!KnotVector, mfem!L2FaceValues, mfem!L2FaceValues!DoubleValued, mfem!L2FaceValues!SingleValued
export mfem!L2_FECollection, mfem!LSZZErrorEstimator, mfem!LinearDiscont2DFECollection, mfem!LinearDiscont3DFECollection
export mfem!LinearFECollection, mfem!LinearNonConf3DFECollection, mfem!Local_FECollection, mfem!MemoryType, mfem!MemoryType!DEFAULT
export mfem!MemoryType!DEVICE, mfem!MemoryType!DEVICE_DEBUG, mfem!MemoryType!DEVICE_UMPIRE, mfem!MemoryType!DEVICE_UMPIRE_2
export mfem!MemoryType!HOST, mfem!MemoryType!HOST_32, mfem!MemoryType!HOST_64, mfem!MemoryType!HOST_DEBUG, mfem!MemoryType!HOST_PINNED
export mfem!MemoryType!HOST_UMPIRE, mfem!MemoryType!MANAGED, mfem!MemoryType!PRESERVE, mfem!MemoryType!SIZE, mfem!Mesh, mfem!Mesh!DEREFINE
export mfem!Mesh!ElementConformity, mfem!Mesh!ElementConformity!Coincident, mfem!Mesh!ElementConformity!NA, mfem!Mesh!ElementConformity!Subset
export mfem!Mesh!ElementConformity!Superset, mfem!Mesh!ElementLocation, mfem!Mesh!ElementLocation!FaceNbr, mfem!Mesh!ElementLocation!Local
export mfem!Mesh!ElementLocation!NA, mfem!Mesh!FaceInfoTag, mfem!Mesh!FaceInfoTag!Boundary, mfem!Mesh!FaceInfoTag!GhostMaster
export mfem!Mesh!FaceInfoTag!GhostSlave, mfem!Mesh!FaceInfoTag!LocalConforming, mfem!Mesh!FaceInfoTag!LocalSlaveNonconforming
export mfem!Mesh!FaceInfoTag!MasterNonconforming, mfem!Mesh!FaceInfoTag!SharedConforming, mfem!Mesh!FaceInfoTag!SharedSlaveNonconforming
export mfem!Mesh!FaceInformation, mfem!Mesh!FaceTopology, mfem!Mesh!FaceTopology!Boundary, mfem!Mesh!FaceTopology!Conforming
export mfem!Mesh!FaceTopology!NA, mfem!Mesh!FaceTopology!Nonconforming, mfem!Mesh!GeometryList, mfem!Mesh!GetTransformationFEforElementType
export mfem!Mesh!LoadFromFile, mfem!Mesh!MakeCartesian1D, mfem!Mesh!MakeCartesian2D, mfem!Mesh!MakeCartesian2DWith4TrisPerQuad
export mfem!Mesh!MakeCartesian2DWith5QuadsPerQuad, mfem!Mesh!MakeCartesian3D, mfem!Mesh!MakeCartesian3DWith24TetsPerHex, mfem!Mesh!MakePeriodic
export mfem!Mesh!MakeRefined, mfem!Mesh!MakeSimplicial, mfem!Mesh!NONE, mfem!Mesh!Operation, mfem!Mesh!REBALANCE, mfem!Mesh!REFINE
export mfem!Mesh!TransformBdrElementToFace, mfem!Mesh!remove_unused_vertices, mfem!Mesh!remove_unused_vertices!, mfem!ND1_3DFECollection
export mfem!ND_FECollection, mfem!ND_R1D_FECollection, mfem!ND_R2D_FECollection, mfem!ND_R2D_Trace_FECollection, mfem!ND_Trace_FECollection
export mfem!NURBSExtension, mfem!NURBSFECollection, mfem!NURBSFECollection!VariableOrder, mfem!NURBSPatch
export mfem!NURBSPatch!Get2DRotationMatrix, mfem!NURBSPatch!Get3DRotationMatrix, mfem!NURBSPatchMap, mfem!NodeExtrudeCoefficient, mfem!Operator
export mfem!Operator!ANY_TYPE, mfem!Operator!Complex_DenseMat, mfem!Operator!Complex_Hypre_ParCSR, mfem!Operator!Complex_Operator
export mfem!Operator!Hypre_ParCSR, mfem!Operator!MFEM_Block_Matrix, mfem!Operator!MFEM_Block_Operator, mfem!Operator!MFEM_ComplexSparseMat
export mfem!Operator!MFEM_SPARSEMAT, mfem!Operator!PETSC_MATAIJ, mfem!Operator!PETSC_MATGENERIC, mfem!Operator!PETSC_MATHYPRE
export mfem!Operator!PETSC_MATIS, mfem!Operator!PETSC_MATNEST, mfem!Operator!PETSC_MATSHELL, mfem!Operator!Type, mfem!OperatorHandle
export mfem!Ordering, mfem!Ordering!Type, mfem!Ordering!byNODES, mfem!Ordering!byVDIM, mfem!P1OnQuadFECollection, mfem!QVectorLayout
export mfem!QVectorLayout!byNODES, mfem!QVectorLayout!byVDIM, mfem!QuadraticDiscont2DFECollection, mfem!QuadraticDiscont3DFECollection
export mfem!QuadraticFECollection, mfem!QuadraticPosDiscont2DFECollection, mfem!QuadraticPosFECollection, mfem!QuadratureFunction
export mfem!QuadratureInterpolator, mfem!QuadratureInterpolator!DERIVATIVES, mfem!QuadratureInterpolator!DETERMINANTS
export mfem!QuadratureInterpolator!EvalFlags, mfem!QuadratureInterpolator!MAX_ND2D, mfem!QuadratureInterpolator!MAX_ND3D
export mfem!QuadratureInterpolator!MAX_NQ2D, mfem!QuadratureInterpolator!MAX_NQ3D, mfem!QuadratureInterpolator!MAX_VDIM2D
export mfem!QuadratureInterpolator!MAX_VDIM3D, mfem!QuadratureInterpolator!PHYSICAL_DERIVATIVES, mfem!QuadratureInterpolator!VALUES, mfem!QuadratureSpace
export mfem!QuadratureSpaceBase, mfem!RT0_2DFECollection, mfem!RT0_3DFECollection, mfem!RT1_2DFECollection, mfem!RT1_3DFECollection
export mfem!RT2_2DFECollection, mfem!RT_FECollection, mfem!RT_R1D_FECollection, mfem!RT_R2D_FECollection, mfem!RT_R2D_Trace_FECollection
export mfem!RT_Trace_FECollection, mfem!RefinedLinearFECollection, mfem!Refinement, mfem!Refinement!X, mfem!Refinement!XY, mfem!Refinement!XYZ
export mfem!Refinement!XZ, mfem!Refinement!Y, mfem!Refinement!YZ, mfem!Refinement!Z, mfem!ShiftRight, mfem!SparseMatrix
export mfem!StatelessDofTransformation, mfem!Swap, mfem!Table, mfem!TensorProductLegendre, mfem!UsesTensorBasis, mfem!VTKFormat, mfem!VTKFormat!ASCII
export mfem!VTKFormat!BINARY, mfem!VTKFormat!BINARY32, mfem!Vector, mfem!VectorCoefficient, mfem!Vertex, mfem!ZZErrorEstimator, ncface
export ncface!, nx, ny, nz, paren, parent, parent!, point_matrix, ref_type, ref_type!, tag, tag!, topology, topology!

import Base.getindex
import Base.setindex!

using CxxWrap
import Libdl
@wrapmodule(()->"$(@__DIR__)/../../libMFEM/libjlMFEM." * Libdl.dlext)

function __init__()
    @initcxx
end

end #module
