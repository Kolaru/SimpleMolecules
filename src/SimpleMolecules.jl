module SimpleMolecules

using AbstractTrees
using AngleBetweenVectors
using AtomicSystems
using Colors
using Graphs
using LinearAlgebra
using PeriodicTable
using Rotations
using UnicodePlots
using Unitful
using UnitfulAtomic

import Mendeleev

include("geometry.jl")
export distance, angle, dihedral

include("molecule.jl")
export AbstractMolecule, CartesianMolecule, InternalCoordinateMolecule, AtomNode
export nbonds, get_bonds

include("loading.jl")

# TODO InternalCoordinateMolecule conversion and back still broken,
# as seen in the dithiane example -- bring it here to test (need copyright free structure data)
# TODO Implement crankshaft rotation (dihedral rotation around a given bond)
# TODO Itnerface with MoleculeMakie

end # module SimpleMolecules
