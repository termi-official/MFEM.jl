using MFEM, CxxWrap

order = 1
mesh_file = "../data/star.mesh"

mesh = MFEM.mfem!Mesh(mesh_file)
UniformRefinement(mesh)

fec = MFEM.mfem!H1_FECollection(order, Dimension(mesh))
fespace = MFEM.mfem!FiniteElementSpace(CxxPtr(mesh), CxxPtr(fec))
println("Number of unknowns: $(GetTrueVSize(fespace))")

# boundary_dofs = nothing # TODO wrap Array<int>
# GetBoundaryTrueDofs(fespace, boundary_dofs);
