module SimpleMolecules

using AngleBetweenVectors
using AtomicSystems
using Graphs
using LinearAlgebra
using PeriodicTable
using Rotations
using UnicodePlots
using Unitful
using UnitfulAtomic

import Mendeleev

include("geometry.jl")
export angle, dihedral

include("molecule.jl")
export AbstractMolecule, CartesianMolecule, InternalCoordinateMolecule
export nbonds, get_bonds

include("loading.jl")

# TODO Implement crankshaft rotation (dihedral rotation around a given bond)
# TODO Do calculations with interval starting with internal coordinates from UED experiment
# for dithiane and see how much wiggle room there is for the Cartesian positions
# TODO Itnerface with MoleculeMakie

end # module SimpleMolecules
