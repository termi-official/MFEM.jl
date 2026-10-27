# MFEM Example 0: solve -Δu = 1 with homogeneous Dirichlet boundary conditions.
#
#   julia --project examples/ex0.jl [mesh_file] [order]
using MFEM

mesh_file = get(ARGS, 1, joinpath(@__DIR__, "..", "data", "star.mesh"))
order = parse(Int, get(ARGS, 2, "1"))

mesh = MFEM.Mesh(mesh_file)
MFEM.UniformRefinement(mesh)

fec = MFEM.H1_FECollection(order, MFEM.Dimension(mesh))
fespace = MFEM.FiniteElementSpace(mesh, fec)
println("Number of unknowns: ", MFEM.GetTrueVSize(fespace))

boundary_dofs = MFEM.Array{Int32}()
MFEM.GetBoundaryTrueDofs(fespace, boundary_dofs)

x = MFEM.GridFunction(fespace)
MFEM.assign(x, 0.0)

one = MFEM.ConstantCoefficient(1.0)
b = MFEM.LinearForm(fespace)
MFEM.AddDomainIntegrator(b, MFEM.DomainLFIntegrator(one))
MFEM.Assemble(b)

a = MFEM.BilinearForm(fespace)
MFEM.AddDomainIntegrator(a, MFEM.DiffusionIntegrator())
MFEM.Assemble(a)

A = MFEM.SparseMatrix()
B = MFEM.Vector()
X = MFEM.Vector()
MFEM.FormLinearSystem(a, boundary_dofs, x, b, A, X, B)

M = MFEM.GSSmoother(A)
MFEM.PCG(A, M, B, X, 1, 200, 1e-12, 0.0)

MFEM.RecoverFEMSolution(a, X, b, x)
MFEM.Save(x, "sol.gf")
MFEM.Save(mesh, "mesh.mesh")
