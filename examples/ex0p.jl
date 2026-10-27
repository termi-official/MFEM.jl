# MFEM Example 0, parallel version: solve -Δu = 1 with homogeneous Dirichlet boundary conditions.
# Requires the parallel variant of jlmfem.
#
#   mpiexecjl -n 4 julia --project examples/ex0p.jl [mesh_file] [order]
using MFEM, MPI

MPI.Init()
MFEM.Hypre!Init()
comm = MPI.COMM_WORLD
isroot = MPI.Comm_rank(comm) == 0

mesh_file = get(ARGS, 1, joinpath(@__DIR__, "..", "data", "star.mesh"))
order = parse(Int, get(ARGS, 2, "1"))

serial_mesh = MFEM.Mesh(mesh_file)
mesh = MFEM.ParMesh(comm, serial_mesh)
MFEM.Clear(serial_mesh)
MFEM.UniformRefinement(mesh)

fec = MFEM.H1_FECollection(order, MFEM.Dimension(mesh))
fespace = MFEM.ParFiniteElementSpace(mesh, fec)
total_num_dofs = MFEM.GlobalTrueVSize(fespace)
isroot && println("Number of unknowns: ", total_num_dofs)

boundary_dofs = MFEM.Array{Int32}()
MFEM.GetBoundaryTrueDofs(fespace, boundary_dofs)

x = MFEM.ParGridFunction(fespace)
MFEM.assign(x, 0.0)

one = MFEM.ConstantCoefficient(1.0)
b = MFEM.ParLinearForm(fespace)
MFEM.AddDomainIntegrator(b, MFEM.DomainLFIntegrator(one))
MFEM.Assemble(b)

a = MFEM.ParBilinearForm(fespace)
MFEM.AddDomainIntegrator(a, MFEM.DiffusionIntegrator())
MFEM.Assemble(a)

A = MFEM.HypreParMatrix()
B = MFEM.Vector()
X = MFEM.Vector()
MFEM.FormLinearSystem(a, boundary_dofs, x, b, A, X, B)

M = MFEM.HypreBoomerAMG(A)
cg = MFEM.CGSolver(comm)
MFEM.SetRelTol(cg, 1e-12)
MFEM.SetMaxIter(cg, 2000)
MFEM.SetPrintLevel(cg, 1)
MFEM.SetPreconditioner(cg, M)
MFEM.SetOperator(cg, A)
MFEM.Mult(cg, B, X)

MFEM.RecoverFEMSolution(a, X, b, x)
MFEM.Save(x, "sol")
MFEM.Save(mesh, "mesh")
