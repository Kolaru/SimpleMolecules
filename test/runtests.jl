using AtomicSystems
using Graphs
using SimpleMolecules
using Test
using Unitful
using UnitfulAtomic

@testset "Angle" begin
    @test angle([1, 0, 0], [0, 0, 0], [0, 1, 0]) == π/2
    @test angle([1, 0, 1], [0, 0, 1], [0, 1, 1]) == π/2
    @test angle([1, 0, 0], [0, 0, 0], [-1, 0, 0]) ≈ π

    @test dihedral([1, 0, 0], [0, 0, 0], [0, 1, 0], [1, 1, 0]) == 0
    @test dihedral([1, 0, 0], [0, 0, 0], [0, 1, 0], [-1, 1, 0]) ≈ π


    p1 = Float64[-1, 1, 0]
    p2 = Float64[1, 1, 0]
    p3 = Float64[1, -1, 0]
    p4 = Float64[-1, -1, 0]
    p5 = Float64[2, -1, 0]
    p6 = Float64[1, -1, -1]
    p7 = Float64[1, -1, 1]

    @test dihedral(p1, p2, p3, p4) == 0
    @test dihedral(p1, p2, p3, p5) ≈ π
    @test dihedral(p1, p2, p3, p6) ≈ π/2
    @test dihedral(p1, p2, p3, p7) ≈ -π/2
end

@testset "Internal coordinates" begin
    molecule = CartesianMolecule(
        AtomicSystem([:C, :C, :C, :C]),
        0.7 * [
            1 1 -1 -1
            1 -1 -1 1
            0 0 0 0
        ]u"Å"
    )

    inc = InternalCoordinateMolecule(molecule, 1)

    @test inc.nodes[1].distances[1] == austrip(1.4u"Å")
    @test inc.nodes[2].angles[1] == π/2
end

@testset "Full molecule" begin
    molecule = read("../example/formic_acid.xyz", CartesianMolecule)

    imolecule = InternalCoordinateMolecule(molecule)
    molecule2 = CartesianMolecule(imolecule)
    imolecule2 = InternalCoordinateMolecule(molecule)
    molecule3 = CartesianMolecule(imolecule)

    @test molecule2.positions ≈ molecule3.positions

    dithiane = read("../example/dithiane.xyz", CartesianMolecule)
    idithiane = InternalCoordinateMolecule(dithiane, "S1")
    dithiane2 = CartesianMolecule(idithiane)
    idithiane2 = InternalCoordinateMolecule(dithiane2, "S1")
    dithiane3 = CartesianMolecule(idithiane2)

    @test dithiane2 ≈ dithiane3
end