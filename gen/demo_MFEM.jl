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
