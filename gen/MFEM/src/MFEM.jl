module MFEM

export AbsMult, AbsMultTranspose, ActualWidth, Add, AddBdrElement, AddBdrFaceIntegrator, AddBdrPoint, AddBdrQuad
export AddBdrQuadAsTriangles, AddBdrSegment, AddBdrTraceFaceIntegrator, AddBdrTriangle, AddBoundaryIntegrator, AddDomainIntegrator
export AddDomainInterpolator, AddElement, AddElementRule, AddElementVector, AddHex, AddHexAsPyramids, AddHexAsTets, AddHexAsWedges
export AddIntegrator, AddInteriorFaceIntegrator, AddMult, AddMultGradPA, AddMultMF, AddMultNURBSPA, AddMultPA, AddMultPatchPA
export AddMultTranspose, AddMultTransposeMF, AddMultTransposePA, AddPyramid, AddQuad, AddQuadAs4TrisWithPoints
export AddQuadAs5QuadsWithPoints, AddRow, AddSegment, AddSubMatrix, AddSubVector, AddTet, AddTraceFaceIntegrator, AddTraceFaceInterpolator
export AddTri, AddTriangle, AddVertex, AddVertexAtMeanCenter, AddVertexParents, AddWedge, AdjointRateMult, All, AllocateMatrix
export Append, ApplyToKnotIntervals, Assemble, AssembleBdrElementMatrix, AssembleDelta, AssembleDeltaElementVect
export AssembleDevice, AssembleDiagonal, AssembleDiagonalMF, AssembleDiagonalPA, AssembleDiagonalPA_ADAt, AssembleDiagonal_ADAt
export AssembleEA, AssembleEABoundaryFaces, AssembleEAInteriorFaces, AssembleElementGrad, AssembleElementMatrix
export AssembleElementMatrix2, AssembleElementVector, AssembleFaceGrad, AssembleFaceMatrix, AssembleFaceVector, AssembleGradDiagonalPA
export AssembleGradPA, AssembleH, AssembleMF, AssembleNURBSPA, AssemblePA, AssemblePABoundary, AssemblePABoundaryFaces
export AssemblePAInteriorFaces, AssemblePatchMatrix, AssemblePatchPA, AssembleRHSElementVect, AssembleTraceFaceMatrix, Assign, BooleanMult
export BooleanMultTranspose, BuildDofToArrays, BuildTranspose, CalcObjective, CalcObjectiveGrad, CalcShape, CalcTestShape, CalcTrialShape
export CalcVShape, Capacity, CartesianPartitioning, Center, ChangeVertexDataOwnership, CheckBdrElementOrientation
export CheckDisplacements, CheckElementOrientation, CheckFinite, CheckPartitioning, Clear, ClearColPtr, ClearCuSparse, ClearGPUSparse
export Clone, Column, Column!, ColumnsAreSorted, ComputeBdrElementMatrix, ComputeElementFlux, ComputeElementMatrices
export ComputeElementMatrix, ComputeFlux, ComputeFluxEnergy, ComputeH1Error, ComputeScalingFactor, Conforming, ConformingAssemble
export ConvertFrom, ConvertFromConformingVDofs, ConvertToConformingVDofs, CountElementsPerVDof, CountSmallElems
export CreatePeriodicVertexMapping, CurlDim, D2C_GlobalRestrictionMatrix, D2Const_GlobalRestrictionMatrix, DegreeElevate, DeleteAll, DeleteDevice
export DeleteLast, DerefineByError, DeregisterField, DeregisterQField, Destroy, DiagScale, Dimension, DistanceSquaredTo
export DistanceTo, DofForGeometry, DofOrderForOrientation, DofToVDof, DofTransformationForGeometry, DofsToVDofs, Elem
export ElementToEdgeTable, ElementToElementTable, ElementToFaceTable, EliminateBC, EliminateCol, EliminateCols, EliminateEssentialBC
export EliminateEssentialBCDiag, EliminateEssentialBCFromDofs, EliminateEssentialBCFromDofsDiag, EliminateEssentialBCFromTrialDofs, EliminateRHS
export EliminateRow, EliminateRowCol, EliminateRowColDiag, EliminateRowColMultipleRHS, EliminateRows, EliminateRowsCols
export EliminateTestDofs, EliminateTrialDofs, EliminateVDofs, EliminateVDofsInRHS, EliminateZeroRows, Empty, EnableHybridization
export EnableSparseMatrixSorting, EnableStaticCondensation, EnsureMultTranspose, EnsureNCMesh, EnsureNodes, Error, Errors
export EstimateLargestEigenvalue, EulerNumber, EulerNumber2D, Eval, EvalDelta, EvalP, EvalSymmetric, EvalW, ExplicitMult, FEColl, FESpace
export FaceIsInterior, Finalize, FinalizeHexMesh, FinalizeMesh, FinalizeQuadMesh, FinalizeTetMesh, FinalizeTopology, FinalizeTriMesh
export FinalizeWedgeMesh, Finalized, FindFaceNeighbors, FindPoints, FiniteElementForDim, FiniteElementForGeometry
export FiniteElementTypeFailureMessage, FirstAndLast, FormDiscreteOperator, FormLinearSystem, FormRectangularLinearSystem, FormRectangularSystemMatrix
export FormRectangularSystemOperator, FormSystemMatrix, FormSystemOperator, FreeElementMatrices, FullAddMult, FullAddMultTranspose, FullInnerProduct
export FullMult, Gauss_Seidel_back, Gauss_Seidel_forw, GeneralRefinement, GenerateBoundaryElements, GeneratePartitioning, Get
export GetA, GetACoef, GetAConst, GetAdjointHeight, GetAlpha, GetAlphaCoef, GetAssemblyLevel, GetAttribute, GetB, GetBBFI
export GetBBFI_Marker, GetBCoef, GetBConst, GetBE, GetBFBFI, GetBFBFI_Marker, GetBLFI, GetBTFBFI, GetBTFBFI_Marker, GetBasisType
export GetBdrAttribute, GetBdrElement, GetBdrElementAdjacentElement, GetBdrElementAdjacentElement2, GetBdrElementBaseGeometry
export GetBdrElementData, GetBdrElementDofs, GetBdrElementEdgeIndex, GetBdrElementEdges, GetBdrElementFace, GetBdrElementFaceIndex
export GetBdrElementGeometry, GetBdrElementToDofTable, GetBdrElementTransformation, GetBdrElementType, GetBdrElementVDofs
export GetBdrElementVertices, GetBdrFace, GetBdrFaceIntegrators, GetBdrFaceTransformations, GetBdrPointMatrix, GetBdrValuesFrom, GetBeta
export GetBetaCoef, GetBlockData, GetBlockI, GetBlockJ, GetBlockOffsets, GetBlockTrueOffsets, GetBlocks, GetBoundaryTrueDofs
export GetBoundingBox, GetBoundsVec_Hi, GetBoundsVec_Lo, GetC, GetCeedOp, GetCharacteristics, GetClosedBasisType, GetCoarseToFineMap
export GetCoeff, GetCoefficient, GetCoeffs, GetCollectionName, GetConformingProlongation, GetConformingRestriction
export GetConformingVSize, GetContType, GetConverged, GetCurl, GetCycle, GetD, GetDBFI, GetDBFI_Marker, GetDI, GetDI_Marker, GetDLFI
export GetDLFI_Delta, GetDLFI_Marker, GetDNFI, GetData, GetDeltaCenter, GetDeltaCoefficient, GetDerivMapType, GetDerivRangeType
export GetDerivType, GetDerivative, GetDiag, GetDim, GetDivergence, GetDofMap, GetDofOrdering, GetEdgeDofs, GetEdgeElement
export GetEdgeInteriorDofs, GetEdgeInteriorVDofs, GetEdgeOrder, GetEdgeTransformation, GetEdgeVDofs, GetEdgeVertexTable, GetEdgeVertices
export GetElement, GetElementAverages, GetElementBaseGeometry, GetElementCenter, GetElementColoring, GetElementData
export GetElementDofValues, GetElementDofs, GetElementEdges, GetElementEnergy, GetElementFaces, GetElementForDof, GetElementGeometry
export GetElementInteriorDofs, GetElementInteriorVDofs, GetElementJacobian, GetElementOrder, GetElementRestriction, GetElementRule
export GetElementSize, GetElementToDofTable, GetElementToFaceOrientationTable, GetElementTransformation, GetElementType
export GetElementVDofs, GetElementVertices, GetElementVolume, GetElementsArray, GetEnergy, GetEqualityVec, GetEssentialTrueDofs
export GetEssentialVDofs, GetEvalMode, GetExplicitGradient, GetExponent, GetFBFI, GetFE, GetFES, GetFLFI, GetFLFI_Marker, GetFace
export GetFaceBaseGeometry, GetFaceDofs, GetFaceEdgeTable, GetFaceEdges, GetFaceElement, GetFaceElementTransformations, GetFaceElementType
export GetFaceElements, GetFaceGeometry, GetFaceGeometryType, GetFaceInformation, GetFaceInfos, GetFaceInteriorDofs, GetFaceOrder
export GetFaceQuadratureInterpolator, GetFaceRestriction, GetFaceToBdrElMap, GetFaceToDofTable, GetFaceToElementTable, GetFaceTransformation
export GetFaceVDofs, GetFaceValues, GetFaceVectorValues, GetFaceVertices, GetField, GetFinalNorm, GetFinalRelNorm
export GetGeckoElementOrdering, GetGeometricParametersFromJacobian, GetGeometries, GetGlobalNE, GetGradient, GetGradients, GetGridFunction
export GetGridFunctionEnergy, GetHeight, GetHessians, GetHilbertElementOrdering, GetHpConformingRestriction, GetHpRestrictionMatrix
export GetHybridization, GetI, GetIFLFI, GetImplicitGradient, GetInequalityVec_Hi, GetInequalityVec_Lo, GetInitialNorm, GetIntRule
export GetIntegrationOrder, GetIntegrationPointFrom1D, GetIntegrationRule, GetInteriorFaceIntegrators, GetInteriorFaceTransformations, GetJ
export GetJacobiScaling, GetKCoef, GetLaplacians, GetLastOperation, GetLocalDofForDof, GetLocalFaceTransformation, GetLocalStateEnergyPA
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
export HasNURBSPatchIntRule, HasQField, HasSpMat, HasSpMatElim, Height, HostRead, HostReadData, HostReadI, HostReadJ, HostReadWrite
export HostReadWriteData, HostReadWriteI, HostReadWriteJ, HostWrite, HostWriteData, HostWriteI, HostWriteJ, ImplicitMult, ImplicitSolve
export ImposeBounds, Init, InitTVectors, InnerProduct, IntPoint, Inverse, IsBinaryFormat, IsBoundary, IsConforming, IsDGSpace
export IsDelta, IsInitialized, IsInterior, IsLocal, IsNonconformingCoarse, IsNonconformingFine, IsOfFaceType, IsShared
export IsSquare, IsSymmetric, IsVariableOrder, Iterations, Jacobi, Jacobi2, Jacobi3, KnotInsert, Last, Load, LoseData, LoseMat
export MakeCoarseToFineTable, MakeDataOwner, MakeOwner, MakePtAP, MakeRAP, MakeRef, MakeTRef, Max, MaxNorm, MaxRowSize, MemoryUsage
export MeshGenerator, Min, MonitorResidual, MonitorSolution, MoveDiagonalFirst, MoveNodes, MoveVertices, Mult, MultTranspose, Name
export Neg, NewDataAndSize, NewElement, NewMemoryAndSize, NewNodes, NodesUpdated, Nonconforming, None, Norml1, Norml2
export Normlinf, Normlp, NumCols, NumNonZeroElems, NumRows, OverrideSize, OwnFEC, OwnsData, OwnsGraph, OwnsNodes, OwnsOperator
export PartAddMult, PartMult, Patchwise, Prepend, Prev, Prev!, PrintBdrVTU, PrintVTU, ProcessNewState, Project
export ProjectBdrCoefficient, ProjectBdrCoefficientNormal, ProjectBdrCoefficientTangent, ProjectCoefficient, ProjectDiscCoefficient
export ProjectGridFunction, ProjectSymmetric, ProjectTranspose, ProjectVectorFieldOn, Ptr, QuadratureIntegration, QuadratureSensitivityMult
export RandomRefinement, Randomize, Read, ReadData, ReadI, ReadJ, ReadWrite, ReadWriteData, ReadWriteI, ReadWriteJ
export RebuildElementToDofTable, Reciprocal, RecoverFEMSolution, ReduceInt, RefineAtVertex, RefineByError, RefineNURBSFromFile, RegisterField
export RegisterQField, RemoveInternalBoundaries, RemoveUnusedVertices, ReorderByNodes, ReorderElementToDofTable, ReorderElements
export ReorientTetMesh, Reserve, Reset, ResetError, ResetFactors, ResetTranspose, RestrictConforming, RowIsEmpty, RowSize, SCFESpace
export SUNImplicitSetup, SUNImplicitSetupB, SUNImplicitSolve, SUNImplicitSolveB, SUNMassMult, SUNMassSetup, SUNMassSolve, Save
export SaveFactors, SaveField, SaveMesh, SaveQField, SaveRootFile, Scale, ScaleColumns, ScaleElements, ScaleRow, ScaleRows
export ScaleSubdomains, SearchRow, SerialRAP, Set, Set1w, Set2, Set2w, Set3, Set3w, SetA, SetACoef, SetAConst, SetAbsTol
export SetAdaptiveLinRtol, SetAlpha, SetAlphaCoef, SetAssemblyLevel, SetAttribute, SetAttributes, SetB, SetBCoef, SetBConst
export SetBdrAttribute, SetBeta, SetBetaCoef, SetBounds, SetColPtr, SetComponent, SetCompression, SetCompressionLevel, SetConstant
export SetCurvature, SetCycle, SetData, SetDataAndSize, SetDataFormat, SetDataOwner, SetDeltaCenter, SetDeltaCoefficient
export SetDiagIdentity, SetDiagonalPolicy, SetDirection, SetElementOrder, SetElementRule, SetEqualityConstraint, SetEssentialBC
export SetEssentialTrueDofs, SetEssentialVDofs, SetEvalMode, SetExponent, SetFormat, SetFromTrueDofs, SetFromTrueVector, SetFunction
export SetGraphOwner, SetGridFunction, SetHighOrderOutput, SetHistorySize, SetInequalityConstraint, SetIntRule, SetIntegrationMode
export SetIntegrationRule, SetIterativeSolver, SetKCoef, SetKDim, SetLayer, SetLevelsOfDetail, SetLinearConstraint, SetMaxIter
export SetMaxLevelsOfDetail, SetMesh, SetMonitor, SetNURBSPatchIntRule, SetNodalFESpace, SetNodalGridFunction, SetNode, SetNodes
export SetNodesOwner, SetOperator, SetOperatorOwner, SetOptimizationProblem, SetOrder, SetOwnData, SetOwnRules, SetPAMemoryType
export SetPadDigits, SetPadDigitsCycle, SetPadDigitsRank, SetPatchAttribute, SetPatchBdrAttribute, SetPatchRules1D, SetPointIndices
export SetPositiveDiagonal, SetPrecision, SetPreconditioner, SetPrefixPath, SetPrintLevel, SetRelTol, SetRelaxedHpConformity, SetRow
export SetScale, SetSeed, SetSize, SetSolutionBounds, SetSolver, SetSpace, SetSpaces, SetSubMatrix, SetSubMatrixTranspose
export SetSubVector, SetSubVectorComplement, SetTime, SetTimeStep, SetTol, SetTransformation, SetTrueVector, SetType
export SetUpdateOperatorOwner, SetUpdateOperatorType, SetVDim, SetVector, SetVertices, SetWeight, SetWidth, Setup, Size, SortColumnIndices
export SpMat, SpMatElim, SpaceDimension, StaticCondensationIsEnabled, StealData, StealNURBSext, SubDofOrder, Sum, Summary
export SupportsCeed, SupportsDevice, Swap, SwapNodes, Symmetrize, SyncAliasMemory, SyncMemory, TestFESpace, Threshold, ToDenseMatrix
export Tol, TraceFiniteElementForGeometry, Transform, TrialFESpace, Type, UniformRefinement, Update, UpdateCoefficient
export UpdateCoefficients, UpdateConstants, UpdatesFinished, UseCuSparse, UseDevice, UseExternalIntegrators, UseFastAssembly, UseGPUSparse
export UsePrecomputedSparsity, UseRestartMode, UseSparsity, VDofToDof, Value, Value!, VectorDim, VerifyFiniteElementTypes, Warnings, Weight
export Width, Write, WriteData, WriteI, WriteJ, ZeroCoefficient, _Add_, _Get_, _Set_, add!, arrow, assign, association
export association!, attributes, attributes!, bdr_attributes, bdr_attributes!, constant, constant!, cross3D, embeddings, embeddings!
export errors, errors!, fdiv!, first_and_last, first_and_last!, geom, geom!, ghost, ghost!, index, index!, input_size
export isExplicit, isHomogeneous, isImplicit, iterations, iterations!, iterative_mode, iterative_mode!, lod, lod!, matrix, matrix!
export median, mfem!AbstractSparseMatrix, mfem!Add, mfem!Array, mfem!Array2D, mfem!AssemblyLevel, mfem!AssemblyLevel!ELEMENT
export mfem!AssemblyLevel!FULL, mfem!AssemblyLevel!LEGACY, mfem!AssemblyLevel!LEGACYFULL, mfem!AssemblyLevel!NONE, mfem!AssemblyLevel!PARTIAL
export mfem!BiCGSTAB, mfem!BiCGSTABSolver, mfem!BilinearForm, mfem!BilinearFormIntegrator, mfem!BlockILU, mfem!BlockILU!Reordering
export mfem!BlockILU!Reordering!MINIMUM_DISCARDED_FILL, mfem!BlockILU!Reordering!NONE, mfem!BlockNonlinearForm, mfem!BlockNonlinearFormIntegrator
export mfem!BoundaryFlowIntegrator, mfem!BoundaryLFIntegrator, mfem!BoundaryMassIntegrator, mfem!BoundaryNormalLFIntegrator
export mfem!BoundaryTangentialLFIntegrator, mfem!BoundingBox, mfem!CG, mfem!CGSolver, mfem!CartesianCoefficient, mfem!CartesianXCoefficient
export mfem!CartesianYCoefficient, mfem!CartesianZCoefficient, mfem!CoarseFineTransformations, mfem!Coefficient, mfem!CoefficientStorage
export mfem!CoefficientStorage!COMPRESSED, mfem!CoefficientStorage!CONSTANTS, mfem!CoefficientStorage!FULL, mfem!CoefficientStorage!SYMMETRIC
export mfem!CoefficientVector, mfem!ComputeElementLpDistance, mfem!ConservativeConvectionIntegrator, mfem!Const2DFECollection
export mfem!Const3DFECollection, mfem!ConstantCoefficient, mfem!ConstrainedOperator, mfem!ConvectionIntegrator
export mfem!ConvectionIntegrator!GetRule, mfem!ConvectiveVectorConvectionNLFIntegrator, mfem!CrossCrossCoefficient, mfem!CrouzeixRaviartFECollection
export mfem!CubicDiscont2DFECollection, mfem!CubicFECollection, mfem!CurlCurlIntegrator, mfem!CurlGridFunctionCoefficient, mfem!CurlInterpolator
export mfem!CylindricalAzimuthalCoefficient, mfem!CylindricalRadialCoefficient, mfem!DGDiffusionBR2Integrator, mfem!DGDiffusionIntegrator
export mfem!DGDirichletLFIntegrator, mfem!DGElasticityDirichletLFIntegrator, mfem!DGElasticityIntegrator, mfem!DGTraceIntegrator
export mfem!DGTraceIntegrator!GetRule, mfem!DG_Interface_FECollection, mfem!DSTable, mfem!DSmoother, mfem!DataCollection, mfem!DataCollection!Format
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
export mfem!FiniteElementSpace!ListToMarker, mfem!FiniteElementSpace!MarkerToList, mfem!GMRES, mfem!GMRESSolver, mfem!GSSmoother
export mfem!GaussLinearDiscont2DFECollection, mfem!GaussQuadraticDiscont2DFECollection, mfem!Geometry!CUBE, mfem!Geometry!INVALID
export mfem!Geometry!NUM_GEOMETRIES, mfem!Geometry!POINT, mfem!Geometry!PRISM, mfem!Geometry!PYRAMID, mfem!Geometry!SEGMENT, mfem!Geometry!SQUARE
export mfem!Geometry!TETRAHEDRON, mfem!Geometry!TRIANGLE, mfem!Geometry!Type, mfem!GetEVectorOrdering, mfem!GradientGridFunctionCoefficient
export mfem!GradientIntegrator, mfem!GradientIntegrator!GetRule, mfem!GradientInterpolator, mfem!GridFunction, mfem!GridFunction!ARITHMETIC
export mfem!GridFunction!AvgType, mfem!GridFunction!HARMONIC, mfem!GridFunctionCoefficient, mfem!GroupConvectionIntegrator
export mfem!H1Pos_FECollection, mfem!H1Ser_FECollection, mfem!H1_FECollection, mfem!H1_Trace_FECollection, mfem!Hybridization
export mfem!HyperelasticModel, mfem!HyperelasticNLFIntegrator, mfem!IdentityInterpolator, mfem!IdentityMatrixCoefficient
export mfem!IdentityOperator, mfem!IncompressibleNeoHookeanIntegrator, mfem!InnerProduct, mfem!InnerProductCoefficient, mfem!IntRules
export mfem!IntRules!, mfem!IntegrationPoint, mfem!IntegrationRule, mfem!IntegrationRules, mfem!InverseElementTransformation
export mfem!InverseHarmonicModel, mfem!InverseIntegrator, mfem!InverseMatrixCoefficient, mfem!IsFinite, mfem!IsIdentityProlongation
export mfem!IsoparametricTransformation, mfem!IterativeSolver, mfem!IterativeSolver!PrintLevel, mfem!IterativeSolverMonitor, mfem!JumpScaling
export mfem!JumpScaling!CONSTANT, mfem!JumpScaling!JumpScalingType, mfem!JumpScaling!ONE_OVER_H, mfem!JumpScaling!P_SQUARED_OVER_H
export mfem!KnotVector, mfem!L2FaceValues, mfem!L2FaceValues!DoubleValued, mfem!L2FaceValues!SingleValued, mfem!L2_FECollection
export mfem!LBFGSSolver, mfem!LSZZErrorEstimator, mfem!LinearDiscont2DFECollection, mfem!LinearDiscont3DFECollection
export mfem!LinearFECollection, mfem!LinearForm, mfem!LinearFormIntegrator, mfem!LinearNonConf3DFECollection, mfem!Local_FECollection
export mfem!LumpedIntegrator, mfem!MINRES, mfem!MINRESSolver, mfem!MassIntegrator, mfem!MassIntegrator!GetRule, mfem!Matrix
export mfem!MatrixArrayCoefficient, mfem!MatrixCoefficient, mfem!MatrixConstantCoefficient, mfem!MatrixInverse, mfem!MatrixProductCoefficient
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
export mfem!PositionVectorCoefficient, mfem!PowerCoefficient, mfem!PowerMethod, mfem!ProductCoefficient, mfem!ProductOperator, mfem!ProductSolver
export mfem!QVectorLayout, mfem!QVectorLayout!byNODES, mfem!QVectorLayout!byVDIM, mfem!QuadraticDiscont2DFECollection
export mfem!QuadraticDiscont3DFECollection, mfem!QuadraticFECollection, mfem!QuadraticPosDiscont2DFECollection, mfem!QuadraticPosFECollection
export mfem!Quadrature1D, mfem!Quadrature1D!CheckClosed, mfem!Quadrature1D!CheckOpen, mfem!Quadrature1D!ClosedGL
export mfem!Quadrature1D!ClosedUniform, mfem!Quadrature1D!GaussLegendre, mfem!Quadrature1D!GaussLobatto, mfem!Quadrature1D!Invalid
export mfem!Quadrature1D!OpenHalfUniform, mfem!Quadrature1D!OpenUniform, mfem!QuadratureFunction, mfem!QuadratureFunctionCoefficient
export mfem!QuadratureFunctions1D, mfem!QuadratureFunctions1D!ClosedGL, mfem!QuadratureFunctions1D!ClosedUniform
export mfem!QuadratureFunctions1D!GaussLegendre, mfem!QuadratureFunctions1D!GaussLobatto, mfem!QuadratureFunctions1D!GivePolyPoints
export mfem!QuadratureFunctions1D!OpenHalfUniform, mfem!QuadratureFunctions1D!OpenUniform, mfem!QuadratureInterpolator, mfem!QuadratureLFIntegrator
export mfem!QuadratureSpace, mfem!QuadratureSpaceBase, mfem!RAP, mfem!RAPOperator, mfem!RT0_2DFECollection, mfem!RT0_3DFECollection
export mfem!RT1_2DFECollection, mfem!RT1_3DFECollection, mfem!RT2_2DFECollection, mfem!RT_FECollection, mfem!RT_R1D_FECollection
export mfem!RT_R2D_FECollection, mfem!RT_R2D_Trace_FECollection, mfem!RT_Trace_FECollection, mfem!RatioCoefficient
export mfem!RectangularConstrainedOperator, mfem!RefinedIntRules, mfem!RefinedIntRules!, mfem!RefinedLinearFECollection, mfem!Refinement, mfem!Refinement!X
export mfem!Refinement!XY, mfem!Refinement!XYZ, mfem!Refinement!XZ, mfem!Refinement!Y, mfem!Refinement!YZ, mfem!Refinement!Z
export mfem!ResidualBCMonitor, mfem!RestrictedCoefficient, mfem!RowNode, mfem!SLBQPOptimizer, mfem!SLI, mfem!SLISolver
export mfem!ScalarCrossProductInterpolator, mfem!ScalarMatrixProductCoefficient, mfem!ScalarProductInterpolator, mfem!ScalarVectorProductCoefficient
export mfem!ScalarVectorProductInterpolator, mfem!ScaledOperator, mfem!SecondOrderTimeDependentOperator, mfem!ShiftRight
export mfem!SkewSymmetricVectorConvectionNLFIntegrator, mfem!Solver, mfem!SparseMatrix, mfem!SparseMatrixFunction, mfem!SparseSmoother
export mfem!SphericalAzimuthalCoefficient, mfem!SphericalPolarCoefficient, mfem!SphericalRadialCoefficient, mfem!StatelessDofTransformation
export mfem!SumCoefficient, mfem!SumIntegrator, mfem!Swap, mfem!SymmetricMatrixCoefficient, mfem!SymmetricMatrixConstantCoefficient
export mfem!Table, mfem!TangentTraceIntegrator, mfem!TensorProductLegendre, mfem!TimeDependentAdjointOperator
export mfem!TimeDependentOperator, mfem!TimeDependentOperator!ADDITIVE_TERM_1, mfem!TimeDependentOperator!ADDITIVE_TERM_2
export mfem!TimeDependentOperator!EXPLICIT, mfem!TimeDependentOperator!EvalMode, mfem!TimeDependentOperator!HOMOGENEOUS
export mfem!TimeDependentOperator!IMPLICIT, mfem!TimeDependentOperator!NORMAL, mfem!TimeDependentOperator!Type, mfem!TraceIntegrator
export mfem!TraceJumpIntegrator, mfem!TransformedCoefficient, mfem!Transpose, mfem!TransposeAbstractSparseMatrix, mfem!TransposeIntegrator
export mfem!TransposeMatrixCoefficient, mfem!TransposeMult, mfem!TransposeOperator, mfem!TripleProductOperator, mfem!UsesTensorBasis, mfem!VTKFormat
export mfem!VTKFormat!ASCII, mfem!VTKFormat!BINARY, mfem!VTKFormat!BINARY32, mfem!Vector, mfem!VectorArrayCoefficient
export mfem!VectorBoundaryFluxLFIntegrator, mfem!VectorBoundaryLFIntegrator, mfem!VectorCoefficient, mfem!VectorConstantCoefficient
export mfem!VectorConvectionNLFIntegrator, mfem!VectorConvectionNLFIntegrator!GetRule, mfem!VectorCrossProductCoefficient
export mfem!VectorCrossProductInterpolator, mfem!VectorCurlCurlIntegrator, mfem!VectorDeltaCoefficient, mfem!VectorDiffusionIntegrator
export mfem!VectorDivergenceIntegrator, mfem!VectorDivergenceIntegrator!GetRule, mfem!VectorDomainLFGradIntegrator, mfem!VectorDomainLFIntegrator
export mfem!VectorFEBoundaryFluxLFIntegrator, mfem!VectorFEBoundaryTangentLFIntegrator, mfem!VectorFECurlIntegrator, mfem!VectorFEDivergenceIntegrator
export mfem!VectorFEDomainLFCurlIntegrator, mfem!VectorFEDomainLFDivIntegrator, mfem!VectorFEDomainLFIntegrator, mfem!VectorFEMassIntegrator
export mfem!VectorFEWeakDivergenceIntegrator, mfem!VectorGridFunctionCoefficient, mfem!VectorInnerProductInterpolator, mfem!VectorMassIntegrator
export mfem!VectorQuadratureFunctionCoefficient, mfem!VectorQuadratureLFIntegrator, mfem!VectorRestrictedCoefficient, mfem!VectorRotProductCoefficient
export mfem!VectorScalarProductInterpolator, mfem!VectorSumCoefficient, mfem!Vertex, mfem!VisItDataCollection, mfem!VisItFieldInfo
export mfem!WhiteGaussianNoiseDomainLFIntegrator, mfem!ZZErrorEstimator, mfem!aGMRES, mfem!ceed!Operator, mfem!infinity, mult!, ncface, ncface!, num_components
export num_components!, paren, parent, parent!, point_matrix, ref_type, ref_type!, sub!, summary, summary!, tag, tag!, topology
export topology!, warnings, warnings!, weight, weight!, x, x!, y, y!, z, z!

import Base.getindex
import Base.setindex!

using CxxWrap
import Libdl
@wrapmodule(()->"$(@__DIR__)/../../libMFEM/libjlMFEM." * Libdl.dlext)

function __init__()
    @initcxx
end

end #module
