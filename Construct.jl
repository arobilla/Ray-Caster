module Construct

using ..Geom

export buildkinkedstandard

function buildfrompoints(v::Vector{Geom.Point})
    if length(v) < 2 
        throw(ArgumentError("Vector for building points is too small"))
    else
        segs = []
        for i in 2:length(v)
            push!(segs, Geom.LineSegment(v[i-1], v[i]))
        end
    end
    return segs
end


"""
    buildkinkedstandard(ang, r, l)

Builds a default kinked wire system. Only runs once at start of simulation. ang is in degrees
"""
function buildkinkedstandard(ang, r, l)
    ang *= pi/180
    s=l/2
    if ang == 0
        pt1 = [Geom.Point(0,0), Geom.Point(2*s+2*l, 0)]
        pt2 = [Geom.Point(0, 2*r), Geom.Point(2*s+2*l, 2*r)]
    else
        pt1 = [Geom.Point(0, 2*r),
        Geom.Point(s + r*(-tan(ang/2)), 2*r),
        Geom.Point(s + l*cos(ang), (l*cos(ang) + r*sin(ang))*tan(ang)+r*(1+ cos(ang))),
        Geom.Point(s+2*l*cos(ang) - r*(-tan(ang/2)), 2*r),
        Geom.Point(2*s+2*l*cos(ang), 2*r)]

        pt2 = [Geom.Point(0, 0),
        Geom.Point(s+r*tan(ang/2), 0),
        Geom.Point(s+l*cos(ang), (l*cos(ang) - r*sin(ang))*tan(ang) + r*(1-cos(ang))),
        Geom.Point(s+2*l*cos(ang)-r*(tan(ang/2)), 0),
        Geom.Point(2*s + 2*l*cos(ang), 0)]
    end
    line_segs = []
    append!(line_segs, buildfrompoints(pt1))
    append!(line_segs, buildfrompoints(pt2))

    bx = [first(pt1).x, last(pt1).x]
    by = [0, (l*cos(ang)+r*sin(ang))*tan(ang) + r*(1+cos(ang))]

    return line_segs, bx, by

end

end