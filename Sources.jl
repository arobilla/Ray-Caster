
module Sources

using ..Geom

"""
    genRays(r, bx, by)

Generates a series of linear ray sources parallel to the axis of the wire. Number of rays given by "linearSourcesNumber" in RunParams.json.
"""
function genRays(r, bx, by)
    s = []
    for y in LinRange(by[1] + 0.01 * r, by[1] + 2 * r - 0.01 * r, Main.RunParams["linearSourcesNumber"])
        push!(s, Geom.Ray(bx[1] + 0.0001, y, 0))
    end
    return s
end

"""
    genRaysAngular(r, bx, by)

Generates a series of conical ray sources aimed along the x-axis (the length axis of the wire).
Controls in RunParams.json: 
Uses "angularSourcesNumber" sources spaced along the radius of the wire.
Fires in a cone of angle from -"angularConeHalfWidth" to +"angularConeHalfWidth"
Uses "angularConeNumber" equally spaced beams along the cone.

"""
function genRaysAngular(r, bx, by)
    s = []
    for y in LinRange(by[1] + 0.01 * r, by[1] + 2 * r - 0.01 * r, Main.RunParams["angularSourcesNumber"])
        for t in LinRange(-Main.RunParams["angularConeHalfWidth"], Main.RunParams["angularConeHalfWidth"], Main.RunParams["angularConeNumber"])
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