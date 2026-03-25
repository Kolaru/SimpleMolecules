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

# TODO Finish writing this test and check that InternalCoordinateMolecule is self
# consistant
@testset "Full molecule" begin
    molecule = read("../example/formic_acid.xyz")

end