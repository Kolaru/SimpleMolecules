
"""
    AtomNode

Represent an atom in a InternalCoordinateMolecule,
as a node which contains all relevant internal coordinates for its children.

Fields
======
- atom: Atom represented by this node
- parent: The parent AtomNode, or nothing if the current AtomNode is root.
- children: The children of the AtomNode
- distances: Distance to its children
- angles: Angles (parent, self, children)
- dihedrals: Dihedral angles (grand_parent, parent, self, children)
"""
struct AtomNode
    atom::Atom
    parent::Union{Nothing, Int}
    children::Vector{AtomNode}
    distances::Vector{Float64}
    angles::Vector{Float64}
    dihedrals::Vector{Float64}
end

AbstractTrees.children(node::AtomNode) = node.children
index(node::AtomNode) = node.atom.index

function Base.show(io::IO, node::AtomNode)
    A = node.atom
    if isempty(node.children)
        print(io, "AtomNode[$(index(node))]($(A.name) ($(A.element.name)))")
    else
        print(io, "AtomNode[$(index(node))]($(A.name) ($(A.element.name))," *
            " distances = $(round.(node.distances, digits = 2))," *
            " angles = $(round.(rad2deg.(node.angles), digits = 2))," *
            " dihedrals = $(round.(rad2deg.(node.dihedrals), digits = 2)))"
        )
    end
end

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
    print(io, "Simple molecule (cartesian coordinates) with $(length(molecule)) atoms and $(nbonds(molecule)) bonds.\n")

    for A in molecule.system
        x, y, z = round.(molecule.positions[:, A] ; digits = 2)
        println("  $(A.index). $(A.name)  ($(A.element.name))  [$x, $y, $z]")
    end

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
        color = colors
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

struct InternalCoordinateMolecule <: AbstractMolecule
    system::AtomicSystem
    topology::SimpleGraph
    root::Int
    nodes::Vector{AtomNode}
end

function Base.show(io::IO, molecule::InternalCoordinateMolecule)
    print(io, "Simple molecule (internal coordinates) with $(length(molecule)) atoms and $(nbonds(molecule)) bonds.\n")

    print_tree(io, molecule.nodes[molecule.root])
end

Graphs.add_edge!(g::SimpleGraph, A::Atom, B::Atom) = add_edge!(g, A.index, B.index)

function get_parent(tree, node)
    isnothing(node) && return nothing
    parents = inneighbors(tree, node)
    isempty(parents) && return nothing
    return only(parents)
end

function get_parent(molecule::InternalCoordinateMolecule, node::AtomNode)
    parent = node.parent
    isnothing(parent) && return nothing
    return molecule.nodes[node.parent]
end
get_parent(molecule::InternalCoordinateMolecule, node::Int) = get_parent(molecule, molecule.nodes[node])
get_parent(::InternalCoordinateMolecule, ::Nothing) = nothing
get_parent(::Nothing) = nothing

function build_node!(nodes, molecule::CartesianMolecule, tree::AbstractGraph, current::Int ; verbose = false)
    positions = molecule.positions
    system = molecule.system

    children = outneighbors(tree, current) 
    parent = get_parent(tree, current)
    grandparent = get_parent(tree, parent)

    distances = [norm(positions[:, current] - positions[:, child]) for child in children]

    angle_reference = parent
    dihedral_reference = grandparent

    if verbose
        @info "Current: $(current). $(system[current])"
        if !isnothing(angle_reference)
            @info "Angle: $(angle_reference). $(system[angle_reference])"
        end

        if !isnothing(dihedral_reference)
            @info "Dihedral: $(dihedral_reference). $(system[dihedral_reference])"
        end
        @info "Children: $children"
    end

    if isnothing(angle_reference)
        angles = fill(NaN, length(children))

        for k in 2:length(children)
            child = children[k]
            angle_reference = children[k - 1]
            
            if verbose
                @warn "No angle ref"
                @info "Child: $child. $(system[child])"
                @info "Angle: $angle_reference. $(system[angle_reference])"
            end

            angles[k] = angle(positions[:, child], positions[:, current], positions[:, angle_reference])
        end
    else
        angles = [angle(positions[:, child], positions[:, current], positions[:, angle_reference]) for child in children]
    end

    if isnothing(dihedral_reference)
        dihedrals = fill(NaN, length(children))

        if isnothing(angle_reference)
            for k in 3:length(children)
                child = children[k]
                angle_reference = children[k - 1]
                dihedral_reference = children[k - 2]

                if verbose
                    @warn "No angle ref, no dihedral ref"
                    @info "Child: $child. $(system[child])"
                    @info "Angle: $angle_reference. $(system[angle_reference])"
                    @info "Dihedral: $dihedral_reference. $(system[dihedral_reference])"
                end

                dihedrals[k] = dihedral(
                    positions[:, dihedral_reference],
                    positions[:, angle_reference],
                    positions[:, current],
                    positions[:, child]
                )
            end
        else
            for k in 2:length(children)
                child = children[k]
                dihedral_reference = children[k - 1]

                if verbose
                    @warn "No dihedral ref"
                    @info "Child: $child. $(system[child])"
                    @info "Dihedral: $dihedral_reference. $(system[dihedral_reference])"
                end

                dihedrals[k] = dihedral(
                    positions[:, dihedral_reference],
                    positions[:, angle_reference],
                    positions[:, current],
                    positions[:, child]
                )
            end
        end
    else
        dihedrals = map(children) do child
            dihedral(
                positions[:, dihedral_reference],
                positions[:, angle_reference],
                positions[:, current],
                positions[:, child]
            )
        end
    end

    if verbose
        println()
    end

    children_nodes = [build_node!(nodes, molecule, tree, C ; verbose) for C in children]
    node = AtomNode(system[current], parent, children_nodes, distances, angles, dihedrals) 
    nodes[current] = node
    return node
end

# Conversion between molecule representations

function InternalCoordinateMolecule(molecule::CartesianMolecule, root_index = 1 ; verbose = false)
    topology = molecule.topology
    tree = bfs_tree(topology, root_index)
    nodes = Vector{AtomNode}(undef, length(molecule))
    build_node!(nodes, molecule, tree, root_index ; verbose)

    return InternalCoordinateMolecule(molecule.system, topology, root_index, nodes)
end

function CartesianMolecule(molecule::InternalCoordinateMolecule)
    positions = zeros(3, length(molecule.system))

    to_process = [molecule.nodes[molecule.root]]

    while !isempty(to_process)
        node = popfirst!(to_process)
        children = node.children
        parent = get_parent(molecule, node)
        grandparent = get_parent(molecule, parent)

        append!(to_process, children)

        for (k, (child, d, α, φ)) in enumerate(zip(children, node.distances, node.angles, node.dihedrals))
            use_sibling_as_parent = false

            if !isnothing(parent)
                angle_reference = parent
            else
                if k == 1
                    angle_reference = nothing
                else
                    use_sibling_as_parent = true
                    angle_reference = children[k - 1]
                end
            end

            if !isnothing(grandparent)
                dihedral_reference = grandparent
            else
                if k <= 1 + use_sibling_as_parent
                    dihedral_reference = nothing
                else
                    if use_sibling_as_parent
                        dihedral_reference = children[k - 2]
                    else
                        dihedral_reference = children[k - 1]
                    end
                end
            end

            if isnothing(angle_reference)
                positions[:, index(child)] = [d, 0, 0]
                continue
            end

            rel_angle = positions[:, index(angle_reference)] - positions[:, index(node)]
            rel = d * normalize(rel_angle)

            if isnothing(dihedral_reference)
                if iszero(norm(positions[:, index(angle_reference)]))
                    axis = normalize(rel × [0, -1, 0])
                else
                    axis = normalize(positions[:, index(angle_reference)] × positions[:, index(node)])
                end

                positions[:, index(child)] = positions[:, index(node)] + RotationVec((α * axis)...) * rel
                continue
            end

            rel_dihedral = positions[:, index(dihedral_reference)] - positions[:, index(angle_reference)]
            axis = normalize(rel_angle × rel_dihedral)

            positions[:, index(child)] = positions[:, index(node)] + RotationVec((φ * normalize(rel_angle))...) * RotationVec((α * axis)...) * rel
        end
    end

    return CartesianMolecule(molecule.system, molecule.topology, positions)
end
