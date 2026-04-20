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
export angle, dihedral

include("molecule.jl")
export AbstractMolecule, CartesianMolecule, InternalCoordinateMolecule, AtomNode
export nbonds, get_bonds

include("loading.jl")

# TODO Take into account the fact that the root (or its first children) may have multiple children:
# in this case they need to be place relative to each other
#   => use the previous children as parent and grandparent to define the angle and dihedral
# TODO Check that the dihedral is correctly oriented
#   => Should go from 0 to 2π and correctly oriented to preserve chirality
# TODO Implement crankshaft rotation (dihedral rotation around a given bond)
# TODO Do calculations with interval starting with internal coordinates from UED experiment
# for dithiane and see how much wiggle room there is for the Cartesian positions
# TODO Itnerface with MoleculeMakie

end # module SimpleMolecules
