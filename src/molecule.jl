Graphs.add_edge!(g::SimpleGraph, A::Atom, B::Atom) = add_edge!(g, A.index, B.index)

"""
    AbstractMolecule

Represent a AbstractMolecule.

All values used internally are in atomic units,
while all output values are Untiful quantities.
"""
abstract type AbstractMolecule end

Base.length(molecule::AbstractMolecule) = nv(molecule.topology)
nbonds(molecule::AbstractMolecule) = ne(molecule.topology)

function get_bonds(molecule::AbstractMolecule)
    map(edges(molecule.topology)) do edge
        return molecule.system[edge.src] => molecule.system[edge.dst]
    end
end

struct CartesianMolecule <: AbstractMolecule
    system::AtomicSystem
    topology::SimpleGraph
    positions::Matrix{Float64}
end

function CartesianMolecule(system::AtomicSystem, positions::AbstractMatrix{<:Quantity} ; bond_tolerance = 0.0)
    @assert size(positions, 1) == 3
    @assert size(positions, 2) == length(system)

    positions = austrip.(positions)

    topology = SimpleGraph()
    add_vertices!(topology, length(system))

    for A in system
        radius_A = austrip(Mendeleev.elements[A.element.number].covalent_radius_pyykko)
        rA = positions[:, A]
        for B in system[(A.index + 1):end]
            radius_B = austrip(Mendeleev.elements[A.element.number].covalent_radius_pyykko)
            rB = positions[:, B]

            threshold = (1 + bond_tolerance) * (radius_A + radius_B)

            if norm(rB - rA) <= threshold
                add_edge!(topology, A, B)
            end
        end
    end

    return CartesianMolecule(system, topology, positions)
end

function Base.show(io::IO, molecule::CartesianMolecule)
    print(io, "Simple molecules with $(length(molecule)) atoms and $(nbonds(molecule)) bonds.\n")

    xx = molecule.positions[1, :]
    yy = molecule.positions[2, :]

    colors = map(molecule.system) do A
        c = parse(Colorant, A.element.cpk_hex)
        return (Int(red(c) * 255), Int(green(c) * 255), Int(blue(c) * 255))
    end

    plt = scatterplot(xx, yy ;
        marker = [first(string(A.element.symbol)) for A in molecule.system],
        compact = true,
        grid = false,
        color = colors,
        border = :none,
        xticks = false,
        yticks = false
    )
    for (A, B) in get_bonds(molecule)
        midx = 0.5 * (xx[A] + xx[B])
        midy = 0.5 * (yy[A] + yy[B])
        lineplot!(plt, [xx[A], midx], [yy[A], midy] ;
            color = colors[A]
        )

        lineplot!(plt, [midx, xx[B]], [midy, yy[B]] ;
            color = colors[B]
        )
    end
    display(plt)
end

"""
    AtomNode

Represent an atom in a InternalCoordinateMolecule,
as a node which contains all relevant internal coordinates for its children.

Fields
======
- parent: The parent AtomNode, or nothing if the current AtomNode is root.
- children: The children of the AtomNode
- distances: Distance to its children
- angles: Angles (parent, self, children)
- dihedrals: Dihedral angles (grand_parent, parent, self, children)
"""
struct AtomNode
    index::Int
    parent::Union{Nothing, Int}
    children::Vector{Int}
    distances::Vector{Float64}
    angles::Vector{Float64}
    dihedrals::Vector{Float64}
end

AtomNode(atom::Atom) = AtomNode(atom.index, nothing, AtomNode[], Float64[], Float64[], Float64[])

function get_parent(tree, node)
    isnothing(node) && return nothing
    parents = inneighbors(tree, node)
    isempty(parents) && return nothing
    return only(parents)
end

get_parent(node::AtomNode) = node.parent
get_parent(::Nothing) = nothing

struct InternalCoordinateMolecule <: AbstractMolecule
    system::AtomicSystem
    topology::SimpleGraph
    root::Int
    nodes::Vector{AtomNode}
end

function InternalCoordinateMolecule(mol::CartesianMolecule, root = 1)
    system = mol.system
    topology = mol.topology
    positions = mol.positions
    nodes = Vector{AtomNode}(undef, length(system))
    root = mol.system[root]

    tree = bfs_tree(topology, root.index)
    to_process = [root.index]

    while !isempty(to_process)
        current = popfirst!(to_process)
        children = outneighbors(tree, current) 
        parent = get_parent(tree, current)
        grandparent = get_parent(tree, parent)

        append!(to_process, children)

        distances = [norm(positions[:, current] - positions[:, child]) for child in children]

        if isnothing(parent)
            angles = fill(NaN, length(children))
        else
            angles = [angle(positions[:, child], positions[:, current], positions[:, parent]) for child in children]
        end

        if isnothing(grandparent)
            dihedrals = fill(NaN, length(children))
        else
            dihedrals = map(children) do child
                dihedral(
                    positions[:, grandparent],
                    positions[:, parent],
                    positions[:, current],
                    positions[:, child]
                )
            end
        end

        nodes[current] = AtomNode(current, parent, children, distances, angles, dihedrals) 
    end
    return InternalCoordinateMolecule(system, topology, root.index, nodes)
end

function CartesianMolecule(mol::InternalCoordinateMolecule)
    positions = zeros(3, length(mol.system))

    node = nodes[mol.root]

    while !isnothing(current)
        current = node.index
        children = node.children
        parent = get_parent(node)
        grandparent = get_parent(parent)

        for (child, d, α, φ) in zip(children, node.distances, node.angles, node.dihedrals)
            if isnothing(parent)
                position[:, child] = [d, 0, 0]
                continue
            end
            
            rel_parent = positions[:, parent] - positions[:, current]
            rel = d * normalize(rel_parent)

            if isnothing(grandparent)
                if izero(norm(positions[:, parent]))
                    axis = normalize(rel × [0, -1, 0])
                end

                axis = normalize(positions[:, parent] × positions[:, current])

                positions[:, child] = positions[:, current] + RotationVec(α * axis) * rel
                continue
            end

            rel_grandparent = positions[:, grandparent] - positions[:, parent]
            axis = -normalize(rel_parent × rel_grandparent)

            positions[:, :child] = RotationVec(φ * normalize(rel_parent)) * RotationVec(α * axis) * rel
        end
    end
end
