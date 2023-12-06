
module Sources

using ..Geom

function genRays(r, bx, by)
    s = []
    for y in LinRange(by[1] + 0.01 * r, by[1] + 2 * r - 0.01 * r, 100)
        push!(s, Geom.Ray(bx[1] + 0.0001, y, 0))
    end
    return s
end

function genRaysAngular(r, bx, by)
    s = []
    for y in LinRange(by[1] + 0.01 * r, by[1] + 2 * r - 0.01 * r, 5)
        for t in LinRange(-60, 60, 14)
            push!(s, Geom.Ray(bx[1] + 0.0001, y, deg2rad(t)))
        end
    end
    return s
end

"""
    genSingleBurst(r, l, bx, by)

Generates a single, radial burst, spanning 360 degrees, slightly inside the left border of the system, at a height y. 
Takes the wire radius as an argument and the boundaries of the system
"""
function genSingleBurst(r, bx, by)
    s = []
    y = by[1] + r
    for t in LinRange(0, 359, 360)
        push!(s, Geom.Ray(bx[1] + (bx[2] - bx[1]) / 40, y, deg2rad(t)))
    end
    return s
end

end