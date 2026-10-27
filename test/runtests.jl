using MFEM, CxxWrap, Test
import HYPRE, Libdl

const star_mesh = joinpath(@__DIR__, "..", "data", "star.mesh")

function h1_space(order = 2)
    mesh = MFEM.Mesh(star_mesh)
    fec = MFEM.H1_FECollection(order, MFEM.Dimension(mesh))
    return MFEM.FiniteElementSpace(mesh, fec)
end

function assembled(lfi, fes)
    b = MFEM.LinearForm(fes)
    MFEM.AddDomainIntegrator(b, lfi)
    MFEM.Assemble(b)
    return copy(MFEM.unsafe_array(b))
end

function apply(bfi, fes, x)
    a = MFEM.BilinearForm(fes)
    MFEM.AddDomainIntegrator(a, bfi)
    MFEM.Assemble(a)
    MFEM.Finalize(a)
    y = MFEM.Vector(MFEM.GetVSize(fes))
    MFEM.Mult(a, x, y)
    return copy(MFEM.unsafe_array(y))
end

function mass_integrator!(el, T, elmat)
    nd = MFEM.GetDof(el)
    shape = MFEM.Vector(nd)
    MFEM.SetSize(elmat, nd)
    MFEM.assign(elmat, 0.0)
    order = 2 * MFEM.GetOrder(el) + MFEM.OrderW(T)
    ir = MFEM.Get(MFEM.GetIntRules(), MFEM.GetGeomType(el), order)
    for i in 0:(MFEM.GetNPoints(ir) - 1)
        ip = MFEM.IntPoint(ir, i)
        MFEM.SetIntPoint(T, CxxPtr(ip))
        MFEM.CalcShape(el, ip, shape)
        MFEM.AddMult_a_VVt(MFEM.Weight(T) * MFEM.weight(ip), shape, elmat)
    end
end

function source_integrator!(el, T, elvect)
    nd = MFEM.GetDof(el)
    shape = MFEM.Vector(nd)
    MFEM.SetSize(elvect, nd)
    MFEM.assign(elvect, 0.0)
    ir = MFEM.Get(MFEM.GetIntRules(), MFEM.GetGeomType(el), 2 * MFEM.GetOrder(el))
    for i in 0:(MFEM.GetNPoints(ir) - 1)
        ip = MFEM.IntPoint(ir, i)
        MFEM.SetIntPoint(T, CxxPtr(ip))
        MFEM.CalcShape(el, ip, shape)
        MFEM.Add(elvect, MFEM.Weight(T) * MFEM.weight(ip), shape)
    end
end

@testset "MFEM.jl ($(MFEM.variant))" begin
    @testset "objects outlive their Julia references" begin
        fes = h1_space()
        GC.gc(); GC.gc()
        @test MFEM.GetNE(MFEM.GetMesh(fes)) > 0
        @test MFEM.GetTrueVSize(fes) > 0
    end

    @testset "Julia coefficients" begin
        fes = h1_space()
        reference = assembled(MFEM.DomainLFIntegrator(MFEM.ConstantCoefficient(2.0)), fes)
        @test assembled(MFEM.DomainLFIntegrator(MFEM.FunctionCoefficient(x -> 2.0)), fes) ≈ reference
        @test assembled(MFEM.DomainLFIntegrator(MFEM.JuliaCoefficient((T, ip) -> 2.0)), fes) ≈ reference

        f = x -> x[1]^2 + x[2]
        gf = MFEM.GridFunction(fes)
        MFEM.ProjectCoefficient(gf, MFEM.FunctionCoefficient(f))
        @test MFEM.ComputeL2Error(gf, MFEM.FunctionCoefficient(f)) < 1e-12

        g! = (v, x) -> (v[1] = x[2]; v[2] = -x[1])
        vgf = MFEM.GridFunction(MFEM.FiniteElementSpace(MFEM.GetMesh(fes), MFEM.FEColl(fes), 2))
        MFEM.ProjectCoefficient(vgf, MFEM.VectorFunctionCoefficient(2, g!))
        @test MFEM.ComputeL2Error(vgf, MFEM.VectorFunctionCoefficient(2, g!)) < 1e-12
    end

    @testset "Julia integrators" begin
        fes = h1_space()
        x = MFEM.Vector(MFEM.GetVSize(fes))
        MFEM.Randomize(x, 1)
        @test apply(MFEM.JuliaBilinearFormIntegrator(mass_integrator!), fes, x) ≈ apply(MFEM.MassIntegrator(), fes, x)
        @test assembled(MFEM.JuliaLinearFormIntegrator(source_integrator!), fes) ≈
              assembled(MFEM.DomainLFIntegrator(MFEM.ConstantCoefficient(1.0)), fes)
    end

    @testset "ex0" begin
        mktempdir() do dir
            cd(() -> include(joinpath(@__DIR__, "..", "examples", "ex0.jl")), dir)
            @test isfile(joinpath(dir, "sol.gf"))
        end
    end

    if MFEM.isparallel()
        @testset "shares MPI and hypre with MPI.jl and HYPRE.jl" begin
            MFEM.MPI.Init()
            HYPRE.Init()
            MFEM.Hypre!Init()
            @test count(contains("libHYPRE"), Libdl.dllist()) == 1
            @test count(contains("libmpi.so"), Libdl.dllist()) == 1
        end

        @testset "ex0p" begin
            mpiexec = MFEM.MPI.mpiexec()
            project = dirname(Base.active_project())
            script = joinpath(@__DIR__, "..", "examples", "ex0p.jl")
            mktempdir() do dir
                cmd = `$mpiexec -n 2 $(Base.julia_cmd()) --project=$project $script`
                @test success(pipeline(Cmd(cmd; dir), stdout = stdout, stderr = stderr))
                @test isfile(joinpath(dir, "sol.000001"))
            end
        end
    end
end
