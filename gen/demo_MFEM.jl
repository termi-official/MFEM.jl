using MFEM, CxxWrap

order = 1
mesh_file = "../data/star.mesh"

mesh = MFEM.mfem!Mesh(mesh_file)
UniformRefinement(mesh)

fec = MFEM.mfem!H1_FECollection(order, Dimension(mesh))
fespace = MFEM.mfem!FiniteElementSpace(CxxPtr(mesh), CxxPtr(fec))
println("Number of unknowns: $(GetTrueVSize(fespace))")

boundary_dofs = MFEM.mfem!Array{Int32}(0)
GetBoundaryTrueDofs(fespace, boundary_dofs)
println("Num true boundry dofs: $(MFEM.Size(boundary_dofs))")

x = MFEM.mfem!GridFunction(CxxPtr(fespace))
assign(x, 0.0)

one = MFEM.mfem!ConstantCoefficient(1.0)

b   = MFEM.mfem!LinearForm(CxxPtr(fespace))
lfi = MFEM.mfem!DomainLFIntegrator(one)
MFEM.AddDomainIntegrator(b, CxxPtr(lfi))
MFEM.Assemble(b)

a   = MFEM.mfem!BilinearForm(CxxPtr(fespace))
bfi = MFEM.mfem!DiffusionIntegrator()
MFEM.AddDomainIntegrator(a, CxxPtr(bfi))
MFEM.Assemble(a)

A = MFEM.mfem!SparseMatrix()
B = MFEM.mfem!Vector()
X = MFEM.mfem!Vector()
MFEM.FormLinearSystem(a, boundary_dofs, x, b, A, X, B)

M = MFEM.mfem!GSSmoother(A)
MFEM.mfem!PCG(A, M, B, X, 1, 200, 1e-12, 0.0)

MFEM.RecoverFEMSolution(a, X, b, x)
MFEM.Save(x, "sol.gf")
MFEM.Save(mesh, "mesh.mesh")
