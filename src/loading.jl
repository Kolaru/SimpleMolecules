"""
    read(filename::AbstractString, ::Type{AtomicSystem} ; units = u"Å", center = true, bond_tolerance = 0.1)

Read an XYZ file and return the corresponding CartesianMolecule.

If `center` is true, the geometry is shifted so that its center of mass is at (0, 0, 0).

`bond_tolerance` control what tolerance to use to automatically determine when a bond is present.
"""
function Base.read(filename::AbstractString, ::Type{CartesianMolecule} ;
        units = u"Å", center = true, bond_tolerance = 0.1)
    lines = readlines(filename)[3:end]
    filter!(!isempty, lines)
    natoms = length(lines)
    atoms = fill("", natoms)
    data = zeros(3, natoms)

    for (k, line) in enumerate(lines)
        atoms[k], xyz... = split(line)
        data[:, k] = parse.(Float64, xyz)
    end

    system = AtomicSystem(Symbol.(atoms))
    geometry = data*units

    if center
        masses = austrip.(atom_mass.(system))
        cm = 1/sum(masses) .* sum(geometry .* reshape(masses, 1, :) ; dims = 2)
        geometry .-= cm
    end

    return CartesianMolecule(system, geometry ; bond_tolerance)
end