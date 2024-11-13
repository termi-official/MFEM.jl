module MFEM

export AbsMult, AbsMultTranspose, ActualWidth, Add, AddBdrElement, AddBdrFaceIntegrator, AddBdrPoint, AddBdrQuad
export AddBdrQuadAsTriangles, AddBdrSegment, AddBdrTraceFaceIntegrator, AddBdrTriangle, AddBoundaryIntegrator, AddDomainIntegrator
export AddDomainInterpolator, AddElement, AddElementRule, AddElementVector, AddHex, AddHexAsPyramids, AddHexAsTets, AddHexAsWedges
export AddIntegrator, AddInteriorFaceIntegrator, AddMult, AddMultGradPA, AddMultMF, AddMultNURBSPA, AddMultPA, AddMultPatchPA
export AddMultTranspose, AddMultTransposeMF, AddMultTransposePA, AddPyramid, AddQuad, AddQuadAs4TrisWithPoints
export AddQuadAs5QuadsWithPoints, AddRow, AddSegment, AddSubMatrix, AddSubVector, AddTet, AddTraceFaceIntegrator, AddTraceFaceInterpolator
export AddTri, AddTriangle, AddVertex, AddVertexAtMeanCenter, AddVertexParents, AddWedge, All, AllocateMatrix, Append
export ApplyToKnotIntervals, Assemble, AssembleBdrElementMatrix, AssembleDelta, AssembleDeltaElementVect, AssembleDevice, AssembleDiagonal
export AssembleDiagonalMF, AssembleDiagonalPA, AssembleDiagonalPA_ADAt, AssembleDiagonal_ADAt, AssembleEA, AssembleEABoundaryFaces
export AssembleEAInteriorFaces, AssembleElementGrad, AssembleElementMatrix, AssembleElementMatrix2, AssembleElementVector, AssembleFaceGrad
export AssembleFaceMatrix, AssembleFaceVector, AssembleGradDiagonalPA, AssembleGradPA, AssembleH, AssembleMF, AssembleNURBSPA, AssemblePA
export AssemblePABoundary, AssemblePABoundaryFaces, AssemblePAInteriorFaces, AssemblePatchMatrix, AssemblePatchPA, AssembleRHSElementVect
export AssembleTraceFaceMatrix, Assign, BooleanMult, BooleanMultTranspose, BuildDofToArrays, BuildTranspose, CalcObjective, CalcObjectiveGrad
export CalcShape, CalcTestShape, CalcTrialShape, CalcVShape, Capacity, CartesianPartitioning, Center, ChangeVertexDataOwnership
export CheckBdrElementOrientation, CheckDisplacements, CheckElementOrientation, CheckFinite, CheckPartitioning, Clear, ClearColPtr, ClearCuSparse
export ClearGPUSparse, Clone, Column, Column!, ColumnsAreSorted, ComputeBdrElementMatrix, ComputeElementFlux, ComputeElementMatrices
export ComputeElementMatrix, ComputeFlux, ComputeFluxEnergy, ComputeH1Error, ComputeScalingFactor, Conforming, ConformingAssemble
export ConvertFromConformingVDofs, ConvertToConformingVDofs, CountElementsPerVDof, CountSmallElems, CreatePeriodicVertexMapping, CurlDim
export D2C_GlobalRestrictionMatrix, D2Const_GlobalRestrictionMatrix, DegreeElevate, DeleteAll, DeleteDevice, DeleteLast, DerefineByError
export DeregisterField, DeregisterQField, Destroy, DiagScale, Dimension, DistanceSquaredTo, DistanceTo, DofForGeometry
export DofOrderForOrientation, DofToVDof, DofTransformationForGeometry, DofsToVDofs, Elem, ElementToEdgeTable, ElementToElementTable
export ElementToFaceTable, EliminateBC, EliminateCol, EliminateCols, EliminateEssentialBC, EliminateEssentialBCDiag
export EliminateEssentialBCFromDofs, EliminateEssentialBCFromDofsDiag, EliminateEssentialBCFromTrialDofs, EliminateRow, EliminateRowCol
export EliminateRowColDiag, EliminateRowColMultipleRHS, EliminateTestDofs, EliminateTrialDofs, EliminateVDofs, EliminateVDofsInRHS
export EliminateZeroRows, Empty, EnableHybridization, EnableSparseMatrixSorting, EnableStaticCondensation, EnsureMultTranspose
export EnsureNCMesh, EnsureNodes, Error, Errors, EulerNumber, EulerNumber2D, Eval, EvalDelta, EvalP, EvalSymmetric, EvalW, FEColl
export FESpace, FaceIsInterior, Finalize, FinalizeHexMesh, FinalizeMesh, FinalizeQuadMesh, FinalizeTetMesh, FinalizeTopology
export FinalizeTriMesh, FinalizeWedgeMesh, Finalized, FindFaceNeighbors, FindPoints, FiniteElementForDim, FiniteElementForGeometry
export FiniteElementTypeFailureMessage, FirstAndLast, FormLinearSystem, FormRectangularLinearSystem, FormRectangularSystemMatrix, FormSystemMatrix
export FreeElementMatrices, FullAddMult, FullAddMultTranspose, FullInnerProduct, FullMult, Gauss_Seidel_back, Gauss_Seidel_forw
export GeneralRefinement, GenerateBoundaryElements, GeneratePartitioning, Get, GetA, GetACoef, GetAConst, GetAlpha, GetAlphaCoef
export GetAssemblyLevel, GetAttribute, GetB, GetBBFI, GetBBFI_Marker, GetBCoef, GetBConst, GetBE, GetBFBFI, GetBFBFI_Marker, GetBLFI
export GetBTFBFI, GetBTFBFI_Marker, GetBasisType, GetBdrAttribute, GetBdrElement, GetBdrElementAdjacentElement
export GetBdrElementAdjacentElement2, GetBdrElementBaseGeometry, GetBdrElementData, GetBdrElementDofs, GetBdrElementEdgeIndex, GetBdrElementEdges
export GetBdrElementFace, GetBdrElementFaceIndex, GetBdrElementGeometry, GetBdrElementToDofTable, GetBdrElementTransformation
export GetBdrElementType, GetBdrElementVDofs, GetBdrElementVertices, GetBdrFace, GetBdrFaceIntegrators, GetBdrFaceTransformations
export GetBdrPointMatrix, GetBdrValuesFrom, GetBeta, GetBetaCoef, GetBlockData, GetBlockI, GetBlockJ, GetBlockOffsets
export GetBlockTrueOffsets, GetBlocks, GetBoundaryTrueDofs, GetBoundingBox, GetBoundsVec_Hi, GetBoundsVec_Lo, GetC, GetCeedOp
export GetCharacteristics, GetClosedBasisType, GetCoarseToFineMap, GetCoeff, GetCoefficient, GetCoeffs, GetCollectionName
export GetConformingProlongation, GetConformingRestriction, GetConformingVSize, GetContType, GetConverged, GetCurl, GetCycle, GetD, GetDBFI
export GetDBFI_Marker, GetDI, GetDI_Marker, GetDLFI, GetDLFI_Delta, GetDLFI_Marker, GetDNFI, GetData, GetDeltaCenter
export GetDeltaCoefficient, GetDerivMapType, GetDerivRangeType, GetDerivType, GetDerivative, GetDiag, GetDim, GetDivergence, GetDofMap
export GetDofOrdering, GetEdgeDofs, GetEdgeElement, GetEdgeInteriorDofs, GetEdgeInteriorVDofs, GetEdgeOrder, GetEdgeTransformation
export GetEdgeVDofs, GetEdgeVertexTable, GetEdgeVertices, GetElement, GetElementAverages, GetElementBaseGeometry, GetElementCenter
export GetElementColoring, GetElementData, GetElementDofValues, GetElementDofs, GetElementEdges, GetElementEnergy, GetElementFaces
export GetElementForDof, GetElementGeometry, GetElementInteriorDofs, GetElementInteriorVDofs, GetElementJacobian, GetElementOrder
export GetElementRestriction, GetElementRule, GetElementSize, GetElementToDofTable, GetElementToFaceOrientationTable
export GetElementTransformation, GetElementType, GetElementVDofs, GetElementVertices, GetElementVolume, GetElementsArray, GetEnergy
export GetEqualityVec, GetEssentialTrueDofs, GetEssentialVDofs, GetExponent, GetFBFI, GetFE, GetFES, GetFLFI, GetFLFI_Marker, GetFace
export GetFaceBaseGeometry, GetFaceDofs, GetFaceEdgeTable, GetFaceEdges, GetFaceElement, GetFaceElementTransformations, GetFaceElementType
export GetFaceElements, GetFaceGeometry, GetFaceGeometryType, GetFaceInformation, GetFaceInfos, GetFaceInteriorDofs, GetFaceOrder
export GetFaceQuadratureInterpolator, GetFaceRestriction, GetFaceToBdrElMap, GetFaceToDofTable, GetFaceToElementTable, GetFaceTransformation
export GetFaceVDofs, GetFaceValues, GetFaceVectorValues, GetFaceVertices, GetField, GetFinalNorm, GetFinalRelNorm
export GetGeckoElementOrdering, GetGeometricParametersFromJacobian, GetGeometries, GetGlobalNE, GetGradient, GetGradients, GetGridFunction
export GetGridFunctionEnergy, GetHeight, GetHessians, GetHilbertElementOrdering, GetHpConformingRestriction, GetHpRestrictionMatrix
export GetHybridization, GetI, GetIFLFI, GetInequalityVec_Hi, GetInequalityVec_Lo, GetInitialNorm, GetIntRule, GetIntegrationOrder
export GetIntegrationPointFrom1D, GetIntegrationRule, GetInteriorFaceIntegrators, GetInteriorFaceTransformations, GetJ, GetJacobiScaling
export GetKCoef, GetLaplacians, GetLastOperation, GetLocalDofForDof, GetLocalFaceTransformation, GetLocalStateEnergyPA
export GetMapType, GetMatrix, GetMaxElementOrder, GetMemory, GetMemoryClass, GetMemoryData, GetMemoryI, GetMemoryJ, GetMesh
export GetNBE, GetNConformingDofs, GetNConst, GetNDofs, GetNE, GetNEDofs, GetNEdges, GetNF, GetNFDofs, GetNFaces, GetNFbyType
export GetNPoints, GetNURBSext, GetNV, GetNVDofs, GetNodalFESpace, GetNodalValues, GetNode, GetNodes, GetNumConstraints, GetNumDof
export GetNumElementInteriorDofs, GetNumFaces, GetNumFacesWithGhost, GetNumGeometries, GetNumIterations, GetOpenBasisType, GetOrder, GetOrdering
export GetOutputProlongation, GetOutputRestriction, GetOutputRestrictionTranspose, GetPatchAttribute, GetPatchBdrAttribute, GetPatchDofs
export GetPatchRule1D, GetPatchRule1D_KnotSpan, GetPatchVDofs, GetPointElement, GetPointMatrix, GetPrefixPath, GetProlongation
export GetProlongationMatrix, GetQField, GetQuadFunction, GetQuadratureInterpolator, GetRangeDim, GetRangeType, GetRefinementTransforms
export GetRestriction, GetRestrictionMatrix, GetRestrictionOperator, GetRestrictionTransposeOperator, GetRow, GetRowColumns
export GetRowEntries, GetRowNorml1, GetRowSums, GetSequence, GetSize, GetSubMatrix, GetSubVector, GetTFBFI, GetTestVDim, GetTime
export GetTimeStep, GetTraceCollection, GetTraceElement, GetTransferOperator, GetTrialVDim, GetTrueDofs, GetTrueTransferOperator
export GetTrueVSize, GetTrueVector, GetType, GetUpdateOperator, GetVDim, GetVDofs, GetVSize, GetValue, GetValues, GetValuesFrom
export GetVec, GetVectorFieldNodalValues, GetVectorFieldValues, GetVectorGradient, GetVectorGradientHat, GetVectorValue
export GetVectorValues, GetVertex, GetVertexDofs, GetVertexToElementTable, GetVertexToVertexTable, GetVertexVDofs, GetVertices
export GetWeights, GetWidth, H2L_GlobalRestrictionMatrix, HasBoundaryElements, HasFaceDofs, HasField, HasGeometry
export HasNURBSPatchIntRule, HasQField, HasSpMat, HasSpMatElim, HostRead, HostReadData, HostReadI, HostReadJ, HostReadWrite
export HostReadWriteData, HostReadWriteI, HostReadWriteJ, HostWrite, HostWriteData, HostWriteI, HostWriteJ, ImposeBounds, Init
export InnerProduct, IntPoint, Inverse, IsBinaryFormat, IsBoundary, IsConforming, IsDGSpace, IsDelta, IsInitialized, IsInterior
export IsLocal, IsNonconformingCoarse, IsNonconformingFine, IsOfFaceType, IsShared, IsSymmetric, IsVariableOrder, Iterations
export Jacobi, Jacobi2, Jacobi3, KnotInsert, Last, Load, LoseData, LoseMat, MakeCoarseToFineTable, MakeDataOwner, MakeOwner
export MakeRef, MakeTRef, Max, MaxNorm, MaxRowSize, MemoryUsage, MeshGenerator, Min, MonitorResidual, MonitorSolution
export MoveDiagonalFirst, MoveNodes, MoveVertices, Mult, MultTranspose, Name, Neg, NewDataAndSize, NewElement, NewMemoryAndSize, NewNodes
export NodesUpdated, Nonconforming, None, Norml1, Norml2, Normlinf, Normlp, NumCols, NumNonZeroElems, NumRows, OverrideSize, OwnFEC
export OwnsData, OwnsGraph, OwnsNodes, PartAddMult, PartMult, Patchwise, Prepend, Prev, Prev!, PrintBdrVTU, PrintVTU
export ProcessNewState, Project, ProjectBdrCoefficient, ProjectBdrCoefficientNormal, ProjectBdrCoefficientTangent, ProjectCoefficient
export ProjectDiscCoefficient, ProjectGridFunction, ProjectSymmetric, ProjectTranspose, ProjectVectorFieldOn, RandomRefinement, Randomize
export Read, ReadData, ReadI, ReadJ, ReadWrite, ReadWriteData, ReadWriteI, ReadWriteJ, RebuildElementToDofTable, Reciprocal
export RecoverFEMSolution, ReduceInt, RefineAtVertex, RefineByError, RefineNURBSFromFile, RegisterField, RegisterQField
export RemoveInternalBoundaries, RemoveUnusedVertices, ReorderByNodes, ReorderElementToDofTable, ReorderElements, ReorientTetMesh, Reserve
export Reset, ResetError, ResetFactors, ResetTranspose, RestrictConforming, RowIsEmpty, RowSize, SCFESpace, Save, SaveFactors
export SaveField, SaveMesh, SaveQField, SaveRootFile, Scale, ScaleColumns, ScaleElements, ScaleRow, ScaleRows, ScaleSubdomains
export SearchRow, SerialRAP, Set, Set1w, Set2, Set2w, Set3, Set3w, SetA, SetACoef, SetAConst, SetAbsTol, SetAdaptiveLinRtol
export SetAlpha, SetAlphaCoef, SetAssemblyLevel, SetAttribute, SetAttributes, SetB, SetBCoef, SetBConst, SetBdrAttribute
export SetBeta, SetBetaCoef, SetBounds, SetColPtr, SetComponent, SetCompression, SetCompressionLevel, SetConstant, SetCurvature
export SetCycle, SetData, SetDataAndSize, SetDataFormat, SetDataOwner, SetDeltaCenter, SetDeltaCoefficient, SetDiagIdentity
export SetDiagonalPolicy, SetDirection, SetElementOrder, SetElementRule, SetEqualityConstraint, SetEssentialBC, SetEssentialTrueDofs
export SetEssentialVDofs, SetExponent, SetFormat, SetFromTrueDofs, SetFromTrueVector, SetFunction, SetGraphOwner, SetGridFunction
export SetHighOrderOutput, SetHistorySize, SetInequalityConstraint, SetIntRule, SetIntegrationMode, SetIntegrationRule, SetIterativeSolver
export SetKCoef, SetKDim, SetLayer, SetLevelsOfDetail, SetLinearConstraint, SetMaxIter, SetMaxLevelsOfDetail, SetMesh
export SetMonitor, SetNURBSPatchIntRule, SetNodalFESpace, SetNodalGridFunction, SetNode, SetNodes, SetNodesOwner, SetOperator
export SetOptimizationProblem, SetOrder, SetOwnData, SetOwnRules, SetPAMemoryType, SetPadDigits, SetPadDigitsCycle, SetPadDigitsRank
export SetPatchAttribute, SetPatchBdrAttribute, SetPatchRules1D, SetPointIndices, SetPositiveDiagonal, SetPrecision, SetPreconditioner
export SetPrefixPath, SetPrintLevel, SetRelTol, SetRelaxedHpConformity, SetRow, SetScale, SetSeed, SetSize, SetSolutionBounds
export SetSolver, SetSpace, SetSpaces, SetSubMatrix, SetSubMatrixTranspose, SetSubVector, SetSubVectorComplement, SetTime
export SetTimeStep, SetTol, SetTransformation, SetTrueVector, SetUpdateOperatorOwner, SetUpdateOperatorType, SetVDim, SetVector
export SetVertices, SetWeight, SetWidth, Setup, Size, SortColumnIndices, SpMat, SpMatElim, SpaceDimension
export StaticCondensationIsEnabled, StealData, StealNURBSext, SubDofOrder, Sum, Summary, SupportsCeed, SupportsDevice, Swap, SwapNodes, Symmetrize
export SyncAliasMemory, SyncMemory, TestFESpace, Threshold, ToDenseMatrix, Tol, TraceFiniteElementForGeometry, Transform, TrialFESpace
export UniformRefinement, Update, UpdateCoefficient, UpdateCoefficients, UpdateConstants, UpdatesFinished, UseCuSparse, UseDevice
export UseExternalIntegrators, UseFastAssembly, UseGPUSparse, UsePrecomputedSparsity, UseRestartMode, UseSparsity, VDofToDof, Value, Value!
export VectorDim, VerifyFiniteElementTypes, Warnings, Weight, Write, WriteData, WriteI, WriteJ, ZeroCoefficient, _Add_, _Get_
export _Set_, add!, assign, association, association!, attributes, attributes!, bdr_attributes, bdr_attributes!, constant
export constant!, cross3D, embeddings, embeddings!, errors, errors!, fdiv!, first_and_last, first_and_last!, geom, geom!, ghost
export ghost!, index, index!, input_size, iterations, iterations!, lod, lod!, matrix, matrix!, median
export mfem!AbstractSparseMatrix, mfem!Add, mfem!Array, mfem!Array2D, mfem!AssemblyLevel, mfem!AssemblyLevel!ELEMENT, mfem!AssemblyLevel!FULL
export mfem!AssemblyLevel!LEGACY, mfem!AssemblyLevel!LEGACYFULL, mfem!AssemblyLevel!NONE, mfem!AssemblyLevel!PARTIAL, mfem!BiCGSTAB
export mfem!BiCGSTABSolver, mfem!BilinearForm, mfem!BilinearFormIntegrator, mfem!BlockILU, mfem!BlockILU!Reordering
export mfem!BlockILU!Reordering!MINIMUM_DISCARDED_FILL, mfem!BlockILU!Reordering!NONE, mfem!BlockNonlinearForm, mfem!BlockNonlinearFormIntegrator
export mfem!BoundaryFlowIntegrator, mfem!BoundaryLFIntegrator, mfem!BoundaryMassIntegrator, mfem!BoundaryNormalLFIntegrator
export mfem!BoundaryTangentialLFIntegrator, mfem!BoundingBox, mfem!CG, mfem!CGSolver, mfem!CartesianCoefficient, mfem!CartesianXCoefficient
export mfem!CartesianYCoefficient, mfem!CartesianZCoefficient, mfem!CoarseFineTransformations, mfem!Coefficient, mfem!CoefficientStorage
export mfem!CoefficientStorage!COMPRESSED, mfem!CoefficientStorage!CONSTANTS, mfem!CoefficientStorage!FULL, mfem!CoefficientStorage!SYMMETRIC
export mfem!CoefficientVector, mfem!ComputeElementLpDistance, mfem!ConservativeConvectionIntegrator, mfem!Const2DFECollection
export mfem!Const3DFECollection, mfem!ConstantCoefficient, mfem!ConvectionIntegrator, mfem!ConvectionIntegrator!GetRule
export mfem!ConvectiveVectorConvectionNLFIntegrator, mfem!CrossCrossCoefficient, mfem!CrouzeixRaviartFECollection, mfem!CubicDiscont2DFECollection
export mfem!CubicFECollection, mfem!CurlCurlIntegrator, mfem!CurlGridFunctionCoefficient, mfem!CurlInterpolator
export mfem!CylindricalAzimuthalCoefficient, mfem!CylindricalRadialCoefficient, mfem!DGDiffusionBR2Integrator, mfem!DGDiffusionIntegrator
export mfem!DGDirichletLFIntegrator, mfem!DGElasticityDirichletLFIntegrator, mfem!DGElasticityIntegrator, mfem!DGTraceIntegrator
export mfem!DGTraceIntegrator!GetRule, mfem!DG_Interface_FECollection, mfem!DSTable, mfem!DataCollection, mfem!DataCollection!Format
export mfem!DataCollection!NO_ERROR, mfem!DataCollection!No_Error, mfem!DataCollection!PARALLEL_FORMAT, mfem!DataCollection!READ_ERROR
export mfem!DataCollection!SERIAL_FORMAT, mfem!DataCollection!WRITE_ERROR, mfem!DeltaCoefficient, mfem!DeltaLFIntegrator, mfem!DenseMatrix
export mfem!DenseSymmetricMatrix, mfem!DerivativeIntegrator, mfem!DeterminantCoefficient, mfem!DiffusionIntegrator
export mfem!DiffusionIntegrator!GetRule, mfem!DirectSubBlockSolver, mfem!DiscreteInterpolator, mfem!DiscreteLinearOperator, mfem!Distance
export mfem!DistanceSquared, mfem!DivDivIntegrator, mfem!DivergenceGridFunctionCoefficient, mfem!DivergenceInterpolator
export mfem!DofTransformation, mfem!DomainLFGradIntegrator, mfem!DomainLFIntegrator, mfem!ElasticityComponentIntegrator
export mfem!ElasticityIntegrator, mfem!Element, mfem!Element!HEXAHEDRON, mfem!Element!POINT, mfem!Element!PYRAMID, mfem!Element!QUADRILATERAL
export mfem!Element!SEGMENT, mfem!Element!TETRAHEDRON, mfem!Element!TRIANGLE, mfem!Element!Type, mfem!Element!WEDGE, mfem!ElementDofOrdering
export mfem!ElementDofOrdering!LEXICOGRAPHIC, mfem!ElementDofOrdering!NATIVE, mfem!ElementRestrictionOperator, mfem!ElementTransformation, mfem!Embedding
export mfem!Extrude1D, mfem!Extrude1DGridFunction, mfem!Extrude2D, mfem!ExtrudeCoefficient, mfem!FGMRESSolver
export mfem!FaceElementTransformations, mfem!FaceQuadratureInterpolator, mfem!FaceRestriction, mfem!FaceType, mfem!FaceType!Boundary
export mfem!FaceType!Interior, mfem!FiniteElement, mfem!FiniteElementCollection, mfem!FiniteElementCollection!CONTINUOUS
export mfem!FiniteElementCollection!DISCONTINUOUS, mfem!FiniteElementCollection!NORMAL, mfem!FiniteElementCollection!New, mfem!FiniteElementCollection!TANGENTIAL
export mfem!FiniteElementSpace, mfem!FiniteElementSpace!AdjustVDofs, mfem!FiniteElementSpace!DecodeDof, mfem!FiniteElementSpace!EncodeDof
export mfem!FiniteElementSpace!ListToMarker, mfem!FiniteElementSpace!MarkerToList, mfem!GMRES, mfem!GMRESSolver, mfem!GaussLinearDiscont2DFECollection
export mfem!GaussQuadraticDiscont2DFECollection, mfem!Geometry!CUBE, mfem!Geometry!INVALID, mfem!Geometry!NUM_GEOMETRIES, mfem!Geometry!POINT
export mfem!Geometry!PRISM, mfem!Geometry!PYRAMID, mfem!Geometry!SEGMENT, mfem!Geometry!SQUARE, mfem!Geometry!TETRAHEDRON
export mfem!Geometry!TRIANGLE, mfem!Geometry!Type, mfem!GetEVectorOrdering, mfem!GradientGridFunctionCoefficient, mfem!GradientIntegrator
export mfem!GradientIntegrator!GetRule, mfem!GradientInterpolator, mfem!GridFunction, mfem!GridFunction!ARITHMETIC, mfem!GridFunction!AvgType
export mfem!GridFunction!HARMONIC, mfem!GridFunctionCoefficient, mfem!GroupConvectionIntegrator, mfem!H1Pos_FECollection, mfem!H1Ser_FECollection
export mfem!H1_FECollection, mfem!H1_Trace_FECollection, mfem!Hybridization, mfem!HyperelasticModel, mfem!HyperelasticNLFIntegrator
export mfem!IdentityInterpolator, mfem!IdentityMatrixCoefficient, mfem!IncompressibleNeoHookeanIntegrator, mfem!InnerProduct
export mfem!InnerProductCoefficient, mfem!IntRules, mfem!IntRules!, mfem!IntegrationPoint, mfem!IntegrationRule, mfem!IntegrationRules
export mfem!InverseElementTransformation, mfem!InverseHarmonicModel, mfem!InverseIntegrator, mfem!InverseMatrixCoefficient, mfem!IsFinite
export mfem!IsoparametricTransformation, mfem!IterativeSolver, mfem!IterativeSolver!PrintLevel, mfem!IterativeSolverMonitor, mfem!JumpScaling
export mfem!JumpScaling!CONSTANT, mfem!JumpScaling!JumpScalingType, mfem!JumpScaling!ONE_OVER_H, mfem!JumpScaling!P_SQUARED_OVER_H
export mfem!KnotVector, mfem!L2FaceValues, mfem!L2FaceValues!DoubleValued, mfem!L2FaceValues!SingleValued, mfem!L2_FECollection
export mfem!LBFGSSolver, mfem!LSZZErrorEstimator, mfem!LinearDiscont2DFECollection, mfem!LinearDiscont3DFECollection
export mfem!LinearFECollection, mfem!LinearForm, mfem!LinearFormIntegrator, mfem!LinearNonConf3DFECollection, mfem!Local_FECollection
export mfem!LumpedIntegrator, mfem!MINRES, mfem!MINRESSolver, mfem!MassIntegrator, mfem!MassIntegrator!GetRule, mfem!MatrixArrayCoefficient
export mfem!MatrixCoefficient, mfem!MatrixConstantCoefficient, mfem!MatrixInverse, mfem!MatrixProductCoefficient
export mfem!MatrixRestrictedCoefficient, mfem!MatrixSumCoefficient, mfem!MatrixVectorProductCoefficient, mfem!Memory, mfem!MemoryClass
export mfem!MemoryClass!DEVICE, mfem!MemoryClass!HOST, mfem!MemoryClass!HOST_32, mfem!MemoryClass!HOST_64, mfem!MemoryClass!MANAGED
export mfem!MemoryType, mfem!MemoryType!DEFAULT, mfem!MemoryType!DEVICE, mfem!MemoryType!DEVICE_DEBUG, mfem!MemoryType!DEVICE_UMPIRE
export mfem!MemoryType!DEVICE_UMPIRE_2, mfem!MemoryType!HOST, mfem!MemoryType!HOST_32, mfem!MemoryType!HOST_64, mfem!MemoryType!HOST_DEBUG
export mfem!MemoryType!HOST_PINNED, mfem!MemoryType!HOST_UMPIRE, mfem!MemoryType!MANAGED, mfem!MemoryType!PRESERVE, mfem!MemoryType!SIZE, mfem!Mesh
export mfem!Mesh!DEREFINE, mfem!Mesh!ElementConformity, mfem!Mesh!ElementConformity!Coincident, mfem!Mesh!ElementConformity!NA
export mfem!Mesh!ElementConformity!Subset, mfem!Mesh!ElementConformity!Superset, mfem!Mesh!ElementLocation, mfem!Mesh!ElementLocation!FaceNbr
export mfem!Mesh!ElementLocation!Local, mfem!Mesh!ElementLocation!NA, mfem!Mesh!FaceInfoTag, mfem!Mesh!FaceInfoTag!Boundary
export mfem!Mesh!FaceInfoTag!GhostMaster, mfem!Mesh!FaceInfoTag!GhostSlave, mfem!Mesh!FaceInfoTag!LocalConforming
export mfem!Mesh!FaceInfoTag!LocalSlaveNonconforming, mfem!Mesh!FaceInfoTag!MasterNonconforming, mfem!Mesh!FaceInfoTag!SharedConforming
export mfem!Mesh!FaceInfoTag!SharedSlaveNonconforming, mfem!Mesh!FaceInformation, mfem!Mesh!FaceTopology, mfem!Mesh!FaceTopology!Boundary
export mfem!Mesh!FaceTopology!Conforming, mfem!Mesh!FaceTopology!NA, mfem!Mesh!FaceTopology!Nonconforming, mfem!Mesh!GeometryList
export mfem!Mesh!GetTransformationFEforElementType, mfem!Mesh!LoadFromFile, mfem!Mesh!MakeCartesian1D, mfem!Mesh!MakeCartesian2D
export mfem!Mesh!MakeCartesian2DWith4TrisPerQuad, mfem!Mesh!MakeCartesian2DWith5QuadsPerQuad, mfem!Mesh!MakeCartesian3D
export mfem!Mesh!MakeCartesian3DWith24TetsPerHex, mfem!Mesh!MakePeriodic, mfem!Mesh!MakeRefined, mfem!Mesh!MakeSimplicial, mfem!Mesh!NONE, mfem!Mesh!Operation
export mfem!Mesh!REBALANCE, mfem!Mesh!REFINE, mfem!Mesh!TransformBdrElementToFace, mfem!Mesh!remove_unused_vertices
export mfem!Mesh!remove_unused_vertices!, mfem!MixedBilinearForm, mfem!MixedCrossCurlCurlIntegrator, mfem!MixedCrossCurlGradIntegrator
export mfem!MixedCrossCurlIntegrator, mfem!MixedCrossGradCurlIntegrator, mfem!MixedCrossGradGradIntegrator, mfem!MixedCrossGradIntegrator
export mfem!MixedCrossProductIntegrator, mfem!MixedCurlCurlIntegrator, mfem!MixedCurlIntegrator, mfem!MixedDirectionalDerivativeIntegrator
export mfem!MixedDivGradIntegrator, mfem!MixedDotProductIntegrator, mfem!MixedGradDivIntegrator, mfem!MixedGradGradIntegrator
export mfem!MixedScalarCrossCurlIntegrator, mfem!MixedScalarCrossGradIntegrator, mfem!MixedScalarCrossProductIntegrator, mfem!MixedScalarCurlIntegrator
export mfem!MixedScalarDerivativeIntegrator, mfem!MixedScalarDivergenceIntegrator, mfem!MixedScalarIntegrator, mfem!MixedScalarMassIntegrator
export mfem!MixedScalarVectorIntegrator, mfem!MixedScalarWeakCrossProductIntegrator, mfem!MixedScalarWeakCurlCrossIntegrator
export mfem!MixedScalarWeakCurlIntegrator, mfem!MixedScalarWeakDerivativeIntegrator, mfem!MixedScalarWeakDivergenceIntegrator
export mfem!MixedScalarWeakGradientIntegrator, mfem!MixedVectorCurlIntegrator, mfem!MixedVectorDivergenceIntegrator, mfem!MixedVectorGradientIntegrator
export mfem!MixedVectorIntegrator, mfem!MixedVectorMassIntegrator, mfem!MixedVectorProductIntegrator, mfem!MixedVectorWeakCurlIntegrator
export mfem!MixedVectorWeakDivergenceIntegrator, mfem!MixedWeakCurlCrossIntegrator, mfem!MixedWeakDivCrossIntegrator, mfem!MixedWeakGradDotIntegrator, mfem!Mult
export mfem!MultAbstractSparseMatrix, mfem!Mult_AtDA, mfem!ND1_3DFECollection, mfem!ND_FECollection, mfem!ND_R1D_FECollection
export mfem!ND_R2D_FECollection, mfem!ND_R2D_Trace_FECollection, mfem!ND_Trace_FECollection, mfem!NURBSExtension, mfem!NURBSFECollection
export mfem!NURBSFECollection!VariableOrder, mfem!NURBSMeshRules, mfem!NeoHookeanModel, mfem!NewtonSolver, mfem!NodeExtrudeCoefficient
export mfem!NonconservativeDGTraceIntegrator, mfem!NonlinearForm, mfem!NonlinearFormIntegrator, mfem!NonlinearFormIntegrator!ELEMENTWISE
export mfem!NonlinearFormIntegrator!Mode, mfem!NonlinearFormIntegrator!PATCHWISE, mfem!NonlinearFormIntegrator!PATCHWISE_REDUCED, mfem!NormalInterpolator
export mfem!NormalTraceIntegrator, mfem!NormalTraceJumpIntegrator, mfem!NormalizedVectorCoefficient, mfem!Operator, mfem!Operator!ANY_TYPE
export mfem!Operator!Complex_DenseMat, mfem!Operator!Complex_Hypre_ParCSR, mfem!Operator!Complex_Operator, mfem!Operator!DIAG_KEEP
export mfem!Operator!DIAG_ONE, mfem!Operator!DIAG_ZERO, mfem!Operator!DiagonalPolicy, mfem!Operator!Hypre_ParCSR
export mfem!Operator!MFEM_Block_Matrix, mfem!Operator!MFEM_Block_Operator, mfem!Operator!MFEM_ComplexSparseMat, mfem!Operator!MFEM_SPARSEMAT
export mfem!Operator!PETSC_MATAIJ, mfem!Operator!PETSC_MATGENERIC, mfem!Operator!PETSC_MATHYPRE, mfem!Operator!PETSC_MATIS
export mfem!Operator!PETSC_MATNEST, mfem!Operator!PETSC_MATSHELL, mfem!Operator!Type, mfem!OperatorChebyshevSmoother, mfem!OperatorHandle
export mfem!OperatorJacobiSmoother, mfem!OptimizationProblem, mfem!OptimizationSolver, mfem!Ordering, mfem!Ordering!Type, mfem!Ordering!byNODES
export mfem!Ordering!byVDIM, mfem!OrthoSolver, mfem!OuterProduct, mfem!OuterProductCoefficient, mfem!P1OnQuadFECollection, mfem!PCG
export mfem!PWCoefficient, mfem!PWConstCoefficient, mfem!PWMatrixCoefficient, mfem!PWVectorCoefficient, mfem!ParaViewDataCollection
export mfem!PositionVectorCoefficient, mfem!PowerCoefficient, mfem!ProductCoefficient, mfem!ProductSolver, mfem!QVectorLayout
export mfem!QVectorLayout!byNODES, mfem!QVectorLayout!byVDIM, mfem!QuadraticDiscont2DFECollection, mfem!QuadraticDiscont3DFECollection
export mfem!QuadraticFECollection, mfem!QuadraticPosDiscont2DFECollection, mfem!QuadraticPosFECollection, mfem!Quadrature1D
export mfem!Quadrature1D!CheckClosed, mfem!Quadrature1D!CheckOpen, mfem!Quadrature1D!ClosedGL, mfem!Quadrature1D!ClosedUniform
export mfem!Quadrature1D!GaussLegendre, mfem!Quadrature1D!GaussLobatto, mfem!Quadrature1D!Invalid, mfem!Quadrature1D!OpenHalfUniform
export mfem!Quadrature1D!OpenUniform, mfem!QuadratureFunction, mfem!QuadratureFunctionCoefficient, mfem!QuadratureFunctions1D
export mfem!QuadratureFunctions1D!ClosedGL, mfem!QuadratureFunctions1D!ClosedUniform, mfem!QuadratureFunctions1D!GaussLegendre
export mfem!QuadratureFunctions1D!GaussLobatto, mfem!QuadratureFunctions1D!GivePolyPoints, mfem!QuadratureFunctions1D!OpenHalfUniform
export mfem!QuadratureFunctions1D!OpenUniform, mfem!QuadratureInterpolator, mfem!QuadratureLFIntegrator, mfem!QuadratureSpace, mfem!QuadratureSpaceBase
export mfem!RAP, mfem!RT0_2DFECollection, mfem!RT0_3DFECollection, mfem!RT1_2DFECollection, mfem!RT1_3DFECollection
export mfem!RT2_2DFECollection, mfem!RT_FECollection, mfem!RT_R1D_FECollection, mfem!RT_R2D_FECollection, mfem!RT_R2D_Trace_FECollection
export mfem!RT_Trace_FECollection, mfem!RatioCoefficient, mfem!RefinedIntRules, mfem!RefinedIntRules!, mfem!RefinedLinearFECollection
export mfem!Refinement, mfem!Refinement!X, mfem!Refinement!XY, mfem!Refinement!XYZ, mfem!Refinement!XZ, mfem!Refinement!Y
export mfem!Refinement!YZ, mfem!Refinement!Z, mfem!ResidualBCMonitor, mfem!RestrictedCoefficient, mfem!RowNode, mfem!SLBQPOptimizer
export mfem!SLI, mfem!SLISolver, mfem!ScalarCrossProductInterpolator, mfem!ScalarMatrixProductCoefficient
export mfem!ScalarProductInterpolator, mfem!ScalarVectorProductCoefficient, mfem!ScalarVectorProductInterpolator, mfem!ShiftRight
export mfem!SkewSymmetricVectorConvectionNLFIntegrator, mfem!Solver, mfem!SparseMatrix, mfem!SparseMatrixFunction, mfem!SphericalAzimuthalCoefficient
export mfem!SphericalPolarCoefficient, mfem!SphericalRadialCoefficient, mfem!StatelessDofTransformation, mfem!SumCoefficient, mfem!SumIntegrator
export mfem!Swap, mfem!SymmetricMatrixCoefficient, mfem!SymmetricMatrixConstantCoefficient, mfem!Table
export mfem!TangentTraceIntegrator, mfem!TensorProductLegendre, mfem!TraceIntegrator, mfem!TraceJumpIntegrator, mfem!TransformedCoefficient
export mfem!Transpose, mfem!TransposeAbstractSparseMatrix, mfem!TransposeIntegrator, mfem!TransposeMatrixCoefficient
export mfem!TransposeMult, mfem!UsesTensorBasis, mfem!VTKFormat, mfem!VTKFormat!ASCII, mfem!VTKFormat!BINARY, mfem!VTKFormat!BINARY32
export mfem!Vector, mfem!VectorArrayCoefficient, mfem!VectorBoundaryFluxLFIntegrator, mfem!VectorBoundaryLFIntegrator
export mfem!VectorCoefficient, mfem!VectorConstantCoefficient, mfem!VectorConvectionNLFIntegrator, mfem!VectorConvectionNLFIntegrator!GetRule
export mfem!VectorCrossProductCoefficient, mfem!VectorCrossProductInterpolator, mfem!VectorCurlCurlIntegrator, mfem!VectorDeltaCoefficient
export mfem!VectorDiffusionIntegrator, mfem!VectorDivergenceIntegrator, mfem!VectorDivergenceIntegrator!GetRule, mfem!VectorDomainLFGradIntegrator
export mfem!VectorDomainLFIntegrator, mfem!VectorFEBoundaryFluxLFIntegrator, mfem!VectorFEBoundaryTangentLFIntegrator, mfem!VectorFECurlIntegrator
export mfem!VectorFEDivergenceIntegrator, mfem!VectorFEDomainLFCurlIntegrator, mfem!VectorFEDomainLFDivIntegrator, mfem!VectorFEDomainLFIntegrator
export mfem!VectorFEMassIntegrator, mfem!VectorFEWeakDivergenceIntegrator, mfem!VectorGridFunctionCoefficient, mfem!VectorInnerProductInterpolator
export mfem!VectorMassIntegrator, mfem!VectorQuadratureFunctionCoefficient, mfem!VectorQuadratureLFIntegrator, mfem!VectorRestrictedCoefficient
export mfem!VectorRotProductCoefficient, mfem!VectorScalarProductInterpolator, mfem!VectorSumCoefficient, mfem!Vertex, mfem!VisItDataCollection
export mfem!VisItFieldInfo, mfem!WhiteGaussianNoiseDomainLFIntegrator, mfem!ZZErrorEstimator, mfem!aGMRES, mfem!ceed!Operator
export mfem!infinity, mult!, ncface, ncface!, num_components, num_components!, paren, parent, parent!, point_matrix, ref_type
export ref_type!, sub!, summary, summary!, tag, tag!, topology, topology!, warnings, warnings!, weight, weight!, x, x!, y, y!, z
export z!

import Base.getindex
import Base.setindex!

using CxxWrap
import Libdl
@wrapmodule(()->"$(@__DIR__)/../../libMFEM/libjlMFEM." * Libdl.dlext)

function __init__()
    @initcxx
end

end #module
