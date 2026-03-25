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

# TODO Fix the fact that the two Hydrogens are at the wrong position
@testset "Full molecule" begin
    molecule = read("../example/formic_acid.xyz", CartesianMolecule)

    imolecule = InternalCoordinateMolecule(molecule)
    molecule2 = CartesianMolecule(imolecule)

    imolecule2 = InternalCoordinateMolecule(molecule)
    molecule3 = CartesianMolecule(imolecule)

    @test molecule2.positions == molecule3.positions
end