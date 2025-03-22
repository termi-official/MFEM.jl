
export AbsMult, AbsMultTranspose, ActualWidth, Add, AddBdrElement, AddBdrFaceIntegrator, AddBdrPoint, AddBdrQuad
export AddBdrQuadAsTriangles, AddBdrSegment, AddBdrTraceFaceIntegrator, AddBdrTriangle, AddBoundaryIntegrator, AddDomainIntegrator
export AddDomainInterpolator, AddElement, AddElementVector, AddHex, AddHexAsPyramids, AddHexAsTets, AddHexAsWedges, AddIntegrator
export AddInteriorFaceIntegrator, AddMult, AddMultGradPA, AddMultMF, AddMultPA, AddMultTranspose, AddMultTransposeMF, AddMultTransposePA
export AddPyramid, AddQuad, AddRow, AddSegment, AddSubMatrix, AddTet, AddTraceFaceIntegrator, AddTraceFaceInterpolator, AddTri
export AddTriangle, AddVertex, AddVertexParents, AddWedge, AdjointRateMult, All, AllocateMatrix, Append, Assemble
export AssembleBdrElementMatrix, AssembleDelta, AssembleDeltaElementVect, AssembleDiagonal, AssembleDiagonalMF, AssembleDiagonalPA
export AssembleDiagonalPA_ADAt, AssembleDiagonal_ADAt, AssembleEA, AssembleEABoundaryFaces, AssembleEAInteriorFaces, AssembleElementGrad
export AssembleElementMatrix, AssembleElementMatrix2, AssembleElementVector, AssembleFaceGrad, AssembleFaceMatrix, AssembleFaceVector
export AssembleGradDiagonalPA, AssembleGradPA, AssembleH, AssembleMF, AssemblePA, AssemblePABoundaryFaces, AssemblePAInteriorFaces
export AssembleRHSElementVect, Assign, BooleanMult, BooleanMultTranspose, BuildDofToArrays, BuildTranspose, CalcObjective, CalcObjectiveGrad
export CalcShape, CalcTestShape, CalcTrialShape, CalcVShape, Capacity, CartesianPartitioning, Center, ChangeVertexDataOwnership
export CheckBdrElementOrientation, CheckDisplacements, CheckElementOrientation, CheckFinite, CheckPartitioning, Clear, ClearColPtr, ClearGPUSparse
export Clone, Column, Column!, ColumnsAreSorted, ComputeBdrElementMatrix, ComputeElementFlux, ComputeElementMatrices
export ComputeElementMatrix, ComputeFlux, ComputeFluxEnergy, ComputeH1Error, ComputeScalingFactor, Conforming, ConformingAssemble
export ConvertFrom, ConvertFromConformingVDofs, ConvertToConformingVDofs, CountSmallElems, CreatePeriodicVertexMapping, CurlDim
export D2C_GlobalRestrictionMatrix, D2Const_GlobalRestrictionMatrix, DegreeElevate, DeleteAll, DeleteDevice, DeleteLast, DerefineByError
export DeregisterField, DeregisterQField, Destroy, DiagScale, Dimension, DistanceSquaredTo, DistanceTo, DofForGeometry
export DofOrderForOrientation, DofToVDof, DofsToVDofs, Elem, ElementToEdgeTable, ElementToElementTable, ElementToFaceTable, EliminateBC
export EliminateCol, EliminateCols, EliminateEssentialBC, EliminateEssentialBCDiag, EliminateEssentialBCFromDofs
export EliminateEssentialBCFromDofsDiag, EliminateEssentialBCFromTrialDofs, EliminateRHS, EliminateRow, EliminateRowCol, EliminateRowColDiag
export EliminateRowColMultipleRHS, EliminateRows, EliminateRowsCols, EliminateTestDofs, EliminateTrialDofs, EliminateVDofs, EliminateVDofsInRHS
export EliminateZeroRows, Empty, EnableHybridization, EnableStaticCondensation, EnsureMultTranspose, EnsureNCMesh, EnsureNodes, Error
export Errors, EstimateLargestEigenvalue, EulerNumber, EulerNumber2D, Eval, EvalDelta, EvalP, EvalSymmetric, EvalW
export ExplicitMult, FEColl, FESpace, FaceIsInterior, Finalize, FinalizeHexMesh, FinalizeMesh, FinalizeQuadMesh, FinalizeTetMesh
export FinalizeTopology, FinalizeTriMesh, FinalizeWedgeMesh, Finalized, FindPoints, FiniteElementForGeometry
export FiniteElementTypeFailureMessage, FirstAndLast, FormDiscreteOperator, FormLinearSystem, FormRectangularLinearSystem, FormRectangularSystemMatrix
export FormRectangularSystemOperator, FormSystemMatrix, FormSystemOperator, FreeElementMatrices, FullAddMult, FullAddMultTranspose, FullInnerProduct
export FullMult, Gauss_Seidel_back, Gauss_Seidel_forw, GeneralRefinement, GenerateBoundaryElements, GeneratePartitioning, Get
export GetA, GetACoef, GetAConst, GetAdjointHeight, GetAlpha, GetAlphaCoef, GetAssemblyLevel, GetAttribute, GetB, GetBBFI
export GetBBFI_Marker, GetBCoef, GetBConst, GetBE, GetBFBFI, GetBFBFI_Marker, GetBLFI, GetBTFBFI, GetBTFBFI_Marker, GetBasisType
export GetBdrAttribute, GetBdrElement, GetBdrElementAdjacentElement, GetBdrElementBaseGeometry, GetBdrElementData, GetBdrElementDofs
export GetBdrElementEdgeIndex, GetBdrElementEdges, GetBdrElementFace, GetBdrElementGeometry, GetBdrElementToDofTable
export GetBdrElementTransformation, GetBdrElementType, GetBdrElementVDofs, GetBdrElementVertices, GetBdrFace, GetBdrFaceIntegrators
export GetBdrFaceTransformations, GetBdrPointMatrix, GetBdrValuesFrom, GetBeta, GetBetaCoef, GetBlockData, GetBlockI, GetBlockJ, GetBlockOffsets
export GetBlockTrueOffsets, GetBlocks, GetBoundaryTrueDofs, GetBoundingBox, GetBoundsVec_Hi, GetBoundsVec_Lo, GetC, GetCeedOp
export GetCharacteristics, GetClosedBasisType, GetCoeff, GetCoeffs, GetCollectionName, GetConformingProlongation, GetConformingRestriction
export GetConformingVSize, GetContType, GetConverged, GetCurl, GetCycle, GetD, GetDBFI, GetDI, GetDLFI, GetDLFI_Delta, GetDNFI, GetData
export GetDeltaCenter, GetDeltaCoefficient, GetDerivative, GetDiag, GetDivergence, GetDofMap, GetDofOrdering, GetEdgeDofs
export GetEdgeElement, GetEdgeInteriorDofs, GetEdgeInteriorVDofs, GetEdgeOrder, GetEdgeTransformation, GetEdgeVDofs
export GetEdgeVertexTable, GetEdgeVertices, GetElement, GetElementAverages, GetElementBaseGeometry, GetElementCenter, GetElementColoring
export GetElementData, GetElementDofValues, GetElementDofs, GetElementEdges, GetElementEnergy, GetElementFaces, GetElementForDof
export GetElementGeometry, GetElementIntRule, GetElementInteriorDofs, GetElementInteriorVDofs, GetElementOrder, GetElementRestriction
export GetElementSize, GetElementToDofTable, GetElementToFaceOrientationTable, GetElementTransformation, GetElementType
export GetElementVDofs, GetElementValues, GetElementVertices, GetElementVolume, GetElementsArray, GetEnergy, GetEqualityVec
export GetEssentialTrueDofs, GetEssentialVDofs, GetEvalMode, GetExplicitGradient, GetExponent, GetFBFI, GetFE, GetFLFI, GetFLFI_Marker
export GetFace, GetFaceBaseGeometry, GetFaceDofs, GetFaceEdgeTable, GetFaceEdges, GetFaceElement, GetFaceElementTransformations
export GetFaceElementType, GetFaceElements, GetFaceGeometry, GetFaceGeometryType, GetFaceInformation, GetFaceInfos, GetFaceInteriorDofs
export GetFaceOrder, GetFaceQuadratureInterpolator, GetFaceRestriction, GetFaceToDofTable, GetFaceToElementTable
export GetFaceTransformation, GetFaceVDofs, GetFaceValues, GetFaceVectorValues, GetFaceVertices, GetField, GetFinalNorm
export GetGeckoElementOrdering, GetGeometries, GetGlobalNE, GetGradient, GetGradients, GetGridFunction, GetGridFunctionEnergy, GetHeight
export GetHessians, GetHilbertElementOrdering, GetHpConformingRestriction, GetHpRestrictionMatrix, GetHybridization, GetI, GetIFLFI
export GetImplicitGradient, GetInequalityVec_Hi, GetInequalityVec_Lo, GetIntRule, GetIntegrationOrder, GetIntegrationRule
export GetInteriorFaceIntegrators, GetInteriorFaceTransformations, GetJ, GetJacobiScaling, GetKCoef, GetLaplacians, GetLastOperation
export GetLocalDofForDof, GetLocalFaceTransformation, GetLocalStateEnergyPA, GetMaxElementOrder, GetMemory, GetMemoryClass, GetMemoryData
export GetMemoryI, GetMemoryJ, GetMesh, GetNBE, GetNConformingDofs, GetNConst, GetNDofs, GetNE, GetNEDofs, GetNEdges, GetNF
export GetNFDofs, GetNFaces, GetNFbyType, GetNPoints, GetNURBSext, GetNV, GetNVDofs, GetNodalFESpace, GetNodalValues, GetNode
export GetNodes, GetNumConstraints, GetNumDof, GetNumElementInteriorDofs, GetNumFaces, GetNumFacesWithGhost, GetNumGeometries
export GetNumIterations, GetOpenBasisType, GetOrder, GetOrdering, GetOutputProlongation, GetOutputRestriction
export GetOutputRestrictionTranspose, GetPointMatrix, GetPrefixPath, GetProlongation, GetProlongationMatrix, GetQField, GetQuadFunction
export GetQuadratureInterpolator, GetRefinementTransforms, GetRestriction, GetRestrictionMatrix, GetRestrictionOperator
export GetRestrictionTransposeOperator, GetRow, GetRowColumns, GetRowEntries, GetRowNorml1, GetRowSums, GetSequence, GetSize, GetSpace, GetSubMatrix
export GetSubVector, GetTFBFI, GetTestVDim, GetTime, GetTimeStep, GetTraceCollection, GetTraceElement, GetTransferOperator
export GetTrialVDim, GetTrueDofs, GetTrueTransferOperator, GetTrueVSize, GetTrueVector, GetType, GetUpdateOperator, GetVDim
export GetVDofs, GetVSize, GetValue, GetValues, GetValuesFrom, GetVec, GetVectorFieldNodalValues, GetVectorFieldValues
export GetVectorGradient, GetVectorValue, GetVectorValues, GetVertex, GetVertexDofs, GetVertexToElementTable, GetVertexVDofs, GetVertices
export GetWeights, GetWidth, H2L_GlobalRestrictionMatrix, HasBoundaryElements, HasFaceDofs, HasField, HasGeometry, HasQField
export Height, HostRead, HostReadData, HostReadI, HostReadJ, HostReadWrite, HostReadWriteData, HostReadWriteI, HostReadWriteJ
export HostWrite, HostWriteData, HostWriteI, HostWriteJ, ImplicitMult, ImplicitSolve, ImposeBounds, Init, InitTVectors
export InnerProduct, IntPoint, Inverse, IsBinaryFormat, IsBoundary, IsConforming, IsDGSpace, IsDelta, IsInitialized, IsInterior
export IsLocal, IsNonconformingCoarse, IsNonconformingFine, IsOfFaceType, IsShared, IsSquare, IsSymmetric, IsVariableOrder
export Iterations, Jacobi, Jacobi2, Jacobi3, KnotInsert, Last, Load, LoseData, LoseMat, MakeCoarseToFineTable, MakeDataOwner
export MakeOwner, MakePtAP, MakeRAP, MakeRef, MakeTRef, Max, MaxNorm, MaxRowSize, MemoryUsage, MeshGenerator, Min
export MonitorResidual, MonitorSolution, MoveDiagonalFirst, MoveNodes, MoveVertices, Mult, MultTranspose, Name, Neg, NewDataAndSize
export NewElement, NewMemoryAndSize, NewNodes, Nonconforming, None, Norml1, Norml2, Normlinf, Normlp, NumCols, NumNonZeroElems
export NumRows, OwnFEC, OwnsData, OwnsGraph, OwnsNodes, OwnsOperator, OwnsSpace, PartAddMult, PartMult, Prepend, Prev, Prev!
export PrintBdrVTU, PrintVTU, ProcessNewState, ProjectBdrCoefficient, ProjectBdrCoefficientNormal, ProjectBdrCoefficientTangent
export ProjectCoefficient, ProjectDiscCoefficient, ProjectGridFunction, ProjectVectorFieldOn, Ptr, QuadratureIntegration
export QuadratureSensitivityMult, RandomRefinement, Randomize, Read, ReadData, ReadI, ReadJ, ReadWrite, ReadWriteData, ReadWriteI, ReadWriteJ
export RecoverFEMSolution, ReduceInt, RefineAtVertex, RefineByError, RegisterField, RegisterQField, RemoveInternalBoundaries
export RemoveUnusedVertices, ReorderByNodes, ReorderElementToDofTable, ReorderElements, Reserve, Reset, ResetError, ResetTranspose
export RestrictConforming, RowIsEmpty, RowSize, SCFESpace, SUNImplicitSetup, SUNImplicitSetupB, SUNImplicitSolve, SUNImplicitSolveB
export SUNMassMult, SUNMassSetup, SUNMassSolve, Save, SaveField, SaveMesh, SaveQField, SaveRootFile, SaveVTU, Scale, ScaleColumns
export ScaleElements, ScaleRow, ScaleRows, ScaleSubdomains, SearchRow, Set, Set1w, Set2, Set2w, Set3, Set3w, SetA, SetACoef
export SetAConst, SetAbsTol, SetAdaptiveLinRtol, SetAlpha, SetAlphaCoef, SetAssemblyLevel, SetAttribute, SetAttributes, SetB
export SetBCoef, SetBConst, SetBdrAttribute, SetBeta, SetBetaCoef, SetBounds, SetColPtr, SetComponent, SetCompression
export SetCompressionLevel, SetCurvature, SetCycle, SetData, SetDataAndSize, SetDataFormat, SetDataOwner, SetDeltaCenter
export SetDeltaCoefficient, SetDiagIdentity, SetDiagonalPolicy, SetDirection, SetElementOrder, SetEqualityConstraint, SetEssentialBC
export SetEssentialTrueDofs, SetEssentialVDofs, SetEvalMode, SetExponent, SetFormat, SetFromTrueDofs, SetFromTrueVector, SetFunction
export SetGraphOwner, SetGridFunction, SetHighOrderOutput, SetHistorySize, SetInequalityConstraint, SetIntRule, SetIntegrationRule
export SetIterativeSolver, SetKCoef, SetKDim, SetLayer, SetLevelsOfDetail, SetLinearConstraint, SetMaxIter, SetMaxLevelsOfDetail, SetMesh
export SetMonitor, SetNodalFESpace, SetNodalGridFunction, SetNode, SetNodes, SetNodesOwner, SetOperator, SetOperatorOwner
export SetOptimizationProblem, SetOrder, SetOwnData, SetOwnRules, SetOwnsSpace, SetPAMemoryType, SetPadDigits, SetPadDigitsCycle
export SetPadDigitsRank, SetPointIndices, SetPositiveDiagonal, SetPrecision, SetPreconditioner, SetPrefixPath, SetPrintLevel, SetRelTol
export SetRelaxedHpConformity, SetRow, SetScale, SetSize, SetSolutionBounds, SetSolver, SetSpace, SetSpaces, SetSubMatrix
export SetSubMatrixTranspose, SetSubVector, SetSubVectorComplement, SetTime, SetTimeStep, SetTol, SetTransformation, SetTrueVector, SetType
export SetUpdateOperatorOwner, SetUpdateOperatorType, SetVDim, SetVector, SetVertices, SetWeight, SetWidth, Setup, Size, SortColumnIndices
export SpMat, SpMatElim, SpaceDimension, StaticCondensationIsEnabled, StealData, StealNURBSext, SubDofOrder, Sum, Summary
export SupportsCeed, Swap, SwapNodes, Symmetrize, SyncAliasMemory, SyncMemory, TestFESpace, Threshold, ToDenseMatrix, Tol
export TraceFiniteElementForGeometry, Transform, TrialFESpace, Type, UniformRefinement, Update, UpdateCoefficient, UpdateCoefficients
export UpdateConstants, UpdatesFinished, UseDevice, UseExternalIntegrators, UseGPUSparse, UsePrecomputedSparsity, UseRestartMode
export UseSparsity, VDofToDof, Value, Value!, VectorDim, VerifyFiniteElementTypes, Warnings, Weight, Width, Write, WriteData
export WriteI, WriteJ, ZeroCoefficient, _Add_, _Get_, _Set_, add!, arrow, assign, association, association!, attributes
export attributes!, bdr_attributes, bdr_attributes!, constant, constant!, embeddings, embeddings!, errors, errors!, fdiv!
export first_and_last, first_and_last!, geom, geom!, ghost, ghost!, index, index!, input_size, isExplicit, isHomogeneous, isImplicit
export iterations, iterations!, iterative_mode, iterative_mode!, lod, lod!, matrix, matrix!, median, mfem!AbstractSparseMatrix
export mfem!Add, mfem!Array, mfem!Array2D, mfem!AssemblyLevel, mfem!AssemblyLevel!ELEMENT, mfem!AssemblyLevel!FULL
export mfem!AssemblyLevel!LEGACY, mfem!AssemblyLevel!LEGACYFULL, mfem!AssemblyLevel!NONE, mfem!AssemblyLevel!PARTIAL, mfem!BiCGSTAB
export mfem!BiCGSTABSolver, mfem!BilinearForm, mfem!BilinearFormIntegrator, mfem!BlockILU, mfem!BlockILU!Reordering
export mfem!BlockILU!Reordering!MINIMUM_DISCARDED_FILL, mfem!BlockILU!Reordering!NONE, mfem!BlockNonlinearForm, mfem!BlockNonlinearFormIntegrator
export mfem!BoundaryFlowIntegrator, mfem!BoundaryLFIntegrator, mfem!BoundaryMassIntegrator, mfem!BoundaryNormalLFIntegrator
export mfem!BoundaryTangentialLFIntegrator, mfem!CG, mfem!CGSolver, mfem!CoarseFineTransformations, mfem!Coefficient, mfem!ComputeElementLpDistance
export mfem!ConservativeConvectionIntegrator, mfem!Const2DFECollection, mfem!Const3DFECollection, mfem!ConstantCoefficient, mfem!ConstrainedOperator
export mfem!ConvectionIntegrator, mfem!ConvectionIntegrator!GetRule, mfem!ConvectiveVectorConvectionNLFIntegrator, mfem!CrossCrossCoefficient
export mfem!CrouzeixRaviartFECollection, mfem!CubicDiscont2DFECollection, mfem!CubicFECollection, mfem!CurlCurlIntegrator
export mfem!CurlGridFunctionCoefficient, mfem!CurlInterpolator, mfem!DGDiffusionBR2Integrator, mfem!DGDiffusionIntegrator, mfem!DGDirichletLFIntegrator
export mfem!DGElasticityDirichletLFIntegrator, mfem!DGElasticityIntegrator, mfem!DGTraceIntegrator, mfem!DGTraceIntegrator!GetRule
export mfem!DG_Interface_FECollection, mfem!DSmoother, mfem!DataCollection, mfem!DataCollection!Format, mfem!DataCollection!NO_ERROR
export mfem!DataCollection!PARALLEL_FORMAT, mfem!DataCollection!READ_ERROR, mfem!DataCollection!SERIAL_FORMAT, mfem!DataCollection!WRITE_ERROR
export mfem!DeltaCoefficient, mfem!DeltaLFIntegrator, mfem!DenseMatrix, mfem!DenseSymmetricMatrix, mfem!DenseTensor
export mfem!DerivativeIntegrator, mfem!DeterminantCoefficient, mfem!DiffusionIntegrator, mfem!DiffusionIntegrator!GetRule
export mfem!DirectSubBlockSolver, mfem!DiscreteInterpolator, mfem!DiscreteLinearOperator, mfem!Distance, mfem!DistanceSquared
export mfem!DivDivIntegrator, mfem!DivergenceGridFunctionCoefficient, mfem!DivergenceInterpolator, mfem!DofTransformation
export mfem!DomainLFGradIntegrator, mfem!DomainLFIntegrator, mfem!ElasticityIntegrator, mfem!Element, mfem!Element!HEXAHEDRON, mfem!Element!POINT
export mfem!Element!PYRAMID, mfem!Element!QUADRILATERAL, mfem!Element!SEGMENT, mfem!Element!TETRAHEDRON, mfem!Element!TRIANGLE
export mfem!Element!Type, mfem!Element!WEDGE, mfem!ElementDofOrdering, mfem!ElementDofOrdering!LEXICOGRAPHIC
export mfem!ElementDofOrdering!NATIVE, mfem!ElementTransformation, mfem!Embedding, mfem!Extrude1D, mfem!Extrude1DGridFunction, mfem!Extrude2D
export mfem!ExtrudeCoefficient, mfem!FGMRESSolver, mfem!FaceElementTransformations, mfem!FaceQuadratureInterpolator, mfem!FaceRestriction
export mfem!FaceType, mfem!FaceType!Boundary, mfem!FaceType!Interior, mfem!FiniteElement, mfem!FiniteElementCollection
export mfem!FiniteElementCollection!CONTINUOUS, mfem!FiniteElementCollection!DISCONTINUOUS, mfem!FiniteElementCollection!NORMAL
export mfem!FiniteElementCollection!New, mfem!FiniteElementCollection!TANGENTIAL, mfem!FiniteElementSpace, mfem!FiniteElementSpace!AdjustVDofs
export mfem!FiniteElementSpace!ListToMarker, mfem!FiniteElementSpace!MarkerToList, mfem!GMRES, mfem!GMRESSolver, mfem!GSSmoother
export mfem!GaussLinearDiscont2DFECollection, mfem!GaussQuadraticDiscont2DFECollection, mfem!Geometry!CUBE, mfem!Geometry!INVALID
export mfem!Geometry!NUM_GEOMETRIES, mfem!Geometry!POINT, mfem!Geometry!PRISM, mfem!Geometry!PYRAMID, mfem!Geometry!SEGMENT, mfem!Geometry!SQUARE
export mfem!Geometry!TETRAHEDRON, mfem!Geometry!TRIANGLE, mfem!Geometry!Type, mfem!GradientGridFunctionCoefficient, mfem!GradientIntegrator
export mfem!GradientIntegrator!GetRule, mfem!GradientInterpolator, mfem!GridFunction, mfem!GridFunction!ARITHMETIC, mfem!GridFunction!AvgType
export mfem!GridFunction!HARMONIC, mfem!GridFunctionCoefficient, mfem!GroupConvectionIntegrator, mfem!H1Pos_FECollection, mfem!H1Ser_FECollection
export mfem!H1_FECollection, mfem!H1_Trace_FECollection, mfem!HCURL_MAX_D1D, mfem!HCURL_MAX_Q1D, mfem!HDIV_MAX_D1D, mfem!HDIV_MAX_Q1D
export mfem!Hybridization, mfem!HyperelasticModel, mfem!HyperelasticNLFIntegrator, mfem!IdentityInterpolator
export mfem!IdentityMatrixCoefficient, mfem!IdentityOperator, mfem!IncompressibleNeoHookeanIntegrator, mfem!InnerProduct, mfem!InnerProductCoefficient
export mfem!IntRules, mfem!IntRules!, mfem!IntegrationPoint, mfem!IntegrationRule, mfem!IntegrationRules
export mfem!InverseElementTransformation, mfem!InverseHarmonicModel, mfem!InverseIntegrator, mfem!InverseMatrixCoefficient, mfem!IsFinite
export mfem!IsIdentityProlongation, mfem!IsoparametricTransformation, mfem!IterativeSolver, mfem!IterativeSolver!PrintLevel
export mfem!IterativeSolverMonitor, mfem!JumpScaling, mfem!JumpScaling!CONSTANT, mfem!JumpScaling!JumpScalingType, mfem!JumpScaling!ONE_OVER_H
export mfem!JumpScaling!P_SQUARED_OVER_H, mfem!L2FaceValues, mfem!L2FaceValues!DoubleValued, mfem!L2FaceValues!SingleValued, mfem!L2_FECollection
export mfem!LBFGSSolver, mfem!LinearDiscont2DFECollection, mfem!LinearDiscont3DFECollection, mfem!LinearFECollection, mfem!LinearForm
export mfem!LinearFormIntegrator, mfem!LinearNonConf3DFECollection, mfem!Local_FECollection, mfem!LumpedIntegrator, mfem!MINRES
export mfem!MINRESSolver, mfem!MassIntegrator, mfem!MassIntegrator!GetRule, mfem!Matrix, mfem!MatrixArrayCoefficient
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
export mfem!Mesh!GetTransformationFEforElementType, mfem!Mesh!LoadFromFile, mfem!Mesh!MakeCartesian1D, mfem!Mesh!MakeCartesian2D, mfem!Mesh!MakeCartesian3D
export mfem!Mesh!MakePeriodic, mfem!Mesh!MakeRefined, mfem!Mesh!MakeSimplicial, mfem!Mesh!NONE, mfem!Mesh!Operation, mfem!Mesh!REBALANCE
export mfem!Mesh!REFINE, mfem!Mesh!remove_unused_vertices, mfem!Mesh!remove_unused_vertices!, mfem!MixedBilinearForm
export mfem!MixedCrossCurlCurlIntegrator, mfem!MixedCrossCurlGradIntegrator, mfem!MixedCrossCurlIntegrator, mfem!MixedCrossGradCurlIntegrator
export mfem!MixedCrossGradGradIntegrator, mfem!MixedCrossGradIntegrator, mfem!MixedCrossProductIntegrator, mfem!MixedCurlCurlIntegrator
export mfem!MixedDirectionalDerivativeIntegrator, mfem!MixedDivGradIntegrator, mfem!MixedDotProductIntegrator, mfem!MixedGradDivIntegrator
export mfem!MixedGradGradIntegrator, mfem!MixedScalarCrossCurlIntegrator, mfem!MixedScalarCrossGradIntegrator
export mfem!MixedScalarCrossProductIntegrator, mfem!MixedScalarCurlIntegrator, mfem!MixedScalarDerivativeIntegrator, mfem!MixedScalarDivergenceIntegrator
export mfem!MixedScalarIntegrator, mfem!MixedScalarMassIntegrator, mfem!MixedScalarVectorIntegrator, mfem!MixedScalarWeakCrossProductIntegrator
export mfem!MixedScalarWeakCurlCrossIntegrator, mfem!MixedScalarWeakCurlIntegrator, mfem!MixedScalarWeakDerivativeIntegrator
export mfem!MixedScalarWeakDivergenceIntegrator, mfem!MixedScalarWeakGradientIntegrator, mfem!MixedVectorCurlIntegrator, mfem!MixedVectorDivergenceIntegrator
export mfem!MixedVectorGradientIntegrator, mfem!MixedVectorIntegrator, mfem!MixedVectorMassIntegrator, mfem!MixedVectorProductIntegrator
export mfem!MixedVectorWeakCurlIntegrator, mfem!MixedVectorWeakDivergenceIntegrator, mfem!MixedWeakCurlCrossIntegrator, mfem!MixedWeakDivCrossIntegrator
export mfem!MixedWeakGradDotIntegrator, mfem!Mult, mfem!MultAbstractSparseMatrix, mfem!Mult_AtDA, mfem!ND1_3DFECollection, mfem!ND_FECollection
export mfem!ND_R1D_FECollection, mfem!ND_R2D_FECollection, mfem!ND_R2D_Trace_FECollection, mfem!ND_Trace_FECollection, mfem!NURBSExtension
export mfem!NURBSFECollection, mfem!NURBSFECollection!VariableOrder, mfem!NeoHookeanModel, mfem!NewtonSolver, mfem!NodeExtrudeCoefficient
export mfem!NonconservativeDGTraceIntegrator, mfem!NonlinearForm, mfem!NonlinearFormIntegrator, mfem!NormalInterpolator, mfem!NormalTraceJumpIntegrator
export mfem!NormalizedVectorCoefficient, mfem!Operator, mfem!Operator!ANY_TYPE, mfem!Operator!Complex_Hypre_ParCSR, mfem!Operator!Complex_Operator
export mfem!Operator!DIAG_KEEP, mfem!Operator!DIAG_ONE, mfem!Operator!DIAG_ZERO, mfem!Operator!DiagonalPolicy, mfem!Operator!Hypre_ParCSR
export mfem!Operator!MFEM_ComplexSparseMat, mfem!Operator!MFEM_SPARSEMAT, mfem!Operator!PETSC_MATAIJ, mfem!Operator!PETSC_MATGENERIC
export mfem!Operator!PETSC_MATHYPRE, mfem!Operator!PETSC_MATIS, mfem!Operator!PETSC_MATNEST, mfem!Operator!PETSC_MATSHELL, mfem!Operator!Type
export mfem!OperatorChebyshevSmoother, mfem!OperatorHandle, mfem!OperatorJacobiSmoother, mfem!OptimizationProblem, mfem!OptimizationSolver
export mfem!Ordering, mfem!Ordering!Type, mfem!Ordering!byNODES, mfem!Ordering!byVDIM, mfem!OuterProduct
export mfem!OuterProductCoefficient, mfem!P1OnQuadFECollection, mfem!PCG, mfem!PWCoefficient, mfem!PWConstCoefficient, mfem!PWMatrixCoefficient
export mfem!PWVectorCoefficient, mfem!ParaViewDataCollection, mfem!PowerCoefficient, mfem!PowerMethod, mfem!ProductCoefficient
export mfem!ProductOperator, mfem!ProductSolver, mfem!QuadraticDiscont2DFECollection, mfem!QuadraticDiscont3DFECollection
export mfem!QuadraticFECollection, mfem!QuadraticPosDiscont2DFECollection, mfem!QuadraticPosFECollection, mfem!Quadrature1D
export mfem!Quadrature1D!CheckClosed, mfem!Quadrature1D!CheckOpen, mfem!Quadrature1D!ClosedGL, mfem!Quadrature1D!ClosedUniform
export mfem!Quadrature1D!GaussLegendre, mfem!Quadrature1D!GaussLobatto, mfem!Quadrature1D!Invalid, mfem!Quadrature1D!OpenHalfUniform
export mfem!Quadrature1D!OpenUniform, mfem!QuadratureFunction, mfem!QuadratureFunctionCoefficient, mfem!QuadratureFunctions1D
export mfem!QuadratureFunctions1D!ClosedGL, mfem!QuadratureFunctions1D!ClosedUniform, mfem!QuadratureFunctions1D!GaussLegendre
export mfem!QuadratureFunctions1D!GaussLobatto, mfem!QuadratureFunctions1D!GivePolyPoints, mfem!QuadratureFunctions1D!OpenHalfUniform
export mfem!QuadratureFunctions1D!OpenUniform, mfem!QuadratureInterpolator, mfem!QuadratureLFIntegrator, mfem!QuadratureSpace, mfem!RAP, mfem!RAPOperator
export mfem!RT0_2DFECollection, mfem!RT0_3DFECollection, mfem!RT1_2DFECollection, mfem!RT1_3DFECollection, mfem!RT2_2DFECollection
export mfem!RT_FECollection, mfem!RT_R1D_FECollection, mfem!RT_R2D_FECollection, mfem!RT_R2D_Trace_FECollection, mfem!RT_Trace_FECollection
export mfem!RatioCoefficient, mfem!RectangularConstrainedOperator, mfem!RefinedIntRules, mfem!RefinedIntRules!
export mfem!RefinedLinearFECollection, mfem!Refinement, mfem!Refinement!X, mfem!Refinement!XY, mfem!Refinement!XYZ, mfem!Refinement!XZ
export mfem!Refinement!Y, mfem!Refinement!YZ, mfem!Refinement!Z, mfem!ResidualBCMonitor, mfem!RestrictedCoefficient, mfem!RowNode
export mfem!SLBQPOptimizer, mfem!SLI, mfem!SLISolver, mfem!ScalarCrossProductInterpolator, mfem!ScalarMatrixProductCoefficient
export mfem!ScalarProductInterpolator, mfem!ScalarVectorProductCoefficient, mfem!ScalarVectorProductInterpolator, mfem!ScaledOperator
export mfem!SecondOrderTimeDependentOperator, mfem!ShiftRight, mfem!SkewSymmetricVectorConvectionNLFIntegrator, mfem!Solver, mfem!SparseMatrix
export mfem!SparseMatrixFunction, mfem!SparseSmoother, mfem!SumCoefficient, mfem!SumIntegrator, mfem!Swap, mfem!SymmetricMatrixCoefficient
export mfem!SymmetricMatrixConstantCoefficient, mfem!Table, mfem!TimeDependentAdjointOperator, mfem!TimeDependentOperator
export mfem!TimeDependentOperator!ADDITIVE_TERM_1, mfem!TimeDependentOperator!ADDITIVE_TERM_2, mfem!TimeDependentOperator!EXPLICIT
export mfem!TimeDependentOperator!EvalMode, mfem!TimeDependentOperator!HOMOGENEOUS, mfem!TimeDependentOperator!IMPLICIT, mfem!TimeDependentOperator!NORMAL
export mfem!TimeDependentOperator!Type, mfem!TraceJumpIntegrator, mfem!TransformedCoefficient, mfem!Transpose, mfem!TransposeAbstractSparseMatrix
export mfem!TransposeIntegrator, mfem!TransposeMatrixCoefficient, mfem!TransposeMult, mfem!TransposeOperator, mfem!TripleProductOperator
export mfem!UsesTensorBasis, mfem!VTKFormat, mfem!VTKFormat!ASCII, mfem!VTKFormat!BINARY, mfem!VTKFormat!BINARY32, mfem!Vector
export mfem!VectorArrayCoefficient, mfem!VectorBoundaryFluxLFIntegrator, mfem!VectorBoundaryLFIntegrator, mfem!VectorCoefficient
export mfem!VectorConstantCoefficient, mfem!VectorConvectionNLFIntegrator, mfem!VectorConvectionNLFIntegrator!GetRule
export mfem!VectorCrossProductCoefficient, mfem!VectorCrossProductInterpolator, mfem!VectorCurlCurlIntegrator, mfem!VectorDeltaCoefficient
export mfem!VectorDiffusionIntegrator, mfem!VectorDivergenceIntegrator, mfem!VectorDivergenceIntegrator!GetRule, mfem!VectorDomainLFIntegrator
export mfem!VectorFEBoundaryFluxLFIntegrator, mfem!VectorFEBoundaryTangentLFIntegrator, mfem!VectorFECurlIntegrator, mfem!VectorFEDivergenceIntegrator
export mfem!VectorFEDomainLFCurlIntegrator, mfem!VectorFEDomainLFDivIntegrator, mfem!VectorFEDomainLFIntegrator, mfem!VectorFEMassIntegrator
export mfem!VectorFEWeakDivergenceIntegrator, mfem!VectorGridFunctionCoefficient, mfem!VectorInnerProductInterpolator, mfem!VectorMassIntegrator
export mfem!VectorQuadratureFunctionCoefficient, mfem!VectorQuadratureLFIntegrator, mfem!VectorRestrictedCoefficient, mfem!VectorRotProductCoefficient
export mfem!VectorScalarProductInterpolator, mfem!VectorSumCoefficient, mfem!Vertex, mfem!VisItDataCollection, mfem!VisItFieldInfo, mfem!ZZErrorEstimator
export mfem!aGMRES, mfem!ceed!Operator, mfem!infinity, mult!, ncface, ncface!, num_components, num_components!, paren, parent
export parent!, point_matrix, ref_type, ref_type!, sub!, summary, summary!, tag, tag!, topology, topology!, warnings, warnings!
export weight, weight!