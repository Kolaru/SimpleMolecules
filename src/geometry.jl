function AngleBetweenVectors.angle(r1, r2, r3)
    r1 = convert.(Float64, r1)
    r2 = convert.(Float64, r2)
    r3 = convert.(Float64, r3)
    return angle(r1 - r2, r3 - r2)
end

function dihedral(r1, r2, r3, r4)
    u1 = r2 - r1
    u2 = r3 - r2
    u3 = r4 - r3

    return atan(
        u2 ⋅ ((u1 × u2) × (u2 × u3)),
        norm(u2) * ((u1 × u2) ⋅ (u2 × u3))
    )
end