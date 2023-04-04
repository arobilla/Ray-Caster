module Geom

using LinearAlgebra

export LineSegment, Point, Ray, intersect

"2d Point class"
struct Point
    x::Float64
    y::Float64
end



"2d line segment class, can be constructed from points"
struct LineSegment
    x1::Float64
    y1::Float64    
    x2::Float64
    y2::Float64
    vec::Vector{Float64}
end

function LineSegment(x1, y1, x2, y2)
    LineSegment(x1, y1, x2, y2, [x2-x1, y2-y1])
end
##Point-based constructor
function LineSegment(p1::Point, p2::Point)
    LineSegment(p2.x, p2.y, p1.x, p1.y)
end

struct Line
    x1::Float64
    y1::Float64    
    x2::Float64
    y2::Float64
    vec::Vector{Float64}
end

function Line(x1, y1, x2, y2)
    Line(x1, y1, x2, y2, [x2-x1, y2-y1])
end
##Point-based constructor
function Line(p1::Point, p2::Point)
    Line(p2.x, p2.y, p1.x, p1.y)
end



struct Ray
    x::Float64
    y::Float64
    theta::Float64
end

"""
    Ray(p1::Point, theta)

TBW
"""
function Ray(p1::Point, theta)
    Ray(p1.x, p1.y, theta)
end
function Ray(p1::Point, p2::Point)
    Ray(p1.x, p1.y, tan((p2.y-p1.y)/(p2.x-p1.x)))
end
function Ray(l::LineSegment)
    Ray(l.x1, l.y1, tan(l.vec[2]/l.vec[1]))
end

#String display for ray method
Base.show(io::IO, r::Ray) = print(io, "Ray(x=$(r.x), y=$(r.y), theta=$(180/pi*r.theta))")


"""
    raytosegment(r::Ray)

creates a short line segment from a ray object
"""
function raytosegment(r::Ray)
    LineSegment(r.x, r.y, r.x+cos(r.theta), r.y+sin(r.theta))
end


"Calculates the intersection of two line segments, see Line-Line Intersection on Wikipedia"
function intersect(l1::LineSegment, l2::LineSegment)
    if l1.vec[2]/l1.vec[1] == l2.vec[2]/l2.vec[1]
        return "NoIntersection", 0, 0
    else
        M = [l1.x2-l1.x1 -(l2.x2-l2.x1)
            l1.y2-l1.y1 -(l2.y2-l2.y1)]
        b = [l2.x1-l1.x1, l2.y1-l1.y1]
        x= M \ b
        if det(M) == 0 #Lines are parallel
            return "NoIntersection", 0, 0
        end
        p = Point(l1.x1+x[1]*l1.vec[1], l1.y1+x[1]*l1.vec[2])
        if checkinsegment(p, l2) && checkinsegment(p, l1)
            return p, x[1], x[2]
        else
            return "NoIntersection", x[1], x[2]
        end
    end
end


function intersect(l1::Line, l2::Line)
    if l1.vec[2]/l1.vec[1] == l2.vec[2]/l2.vec[1]
        return "NoIntersection", 0, 0
    else
        M = [l1.x2-l1.x1 -(l2.x2-l2.x1)
            l1.y2-l1.y1 -(l2.y2-l2.y1)]
        b = [l2.x1-l1.x1, l2.y1-l1.y1]
        x= M \ b
        return Point(l1.x1+x[1]*l1.vec[1], l1.y1+x[1]*l1.vec[2]), x[1], x[2]
    end
end

function intersect(r::Ray, l2::Line)
    l1 = raytosegment(r)
    if l1.vec[2]/l1.vec[1] == l2.vec[2]/l2.vec[1]
        return "NoIntersection", 0, 0
    else
        M = [l1.x2-l1.x1 -(l2.x2-l2.x1)
            l1.y2-l1.y1 -(l2.y2-l2.y1)]
        b = [l2.x1-l1.x1, l2.y1-l1.y1]
        if det(M) == 0 #Lines are parallel
            return "NoIntersection", 0, 0
        end
        x= M \ b
        return Point(l1.x1+x[1]*l1.vec[1], l1.y1+x[1]*l1.vec[2]), x[1], x[2]
    end
end

"""
    intersect(r::Ray, l::LineSegment)

TBW
"""
function intersect(r::Ray, l2::LineSegment)
    l1 = raytosegment(r)

    if l1.vec[2]/l1.vec[1] == l2.vec[2]/l2.vec[1]
        return "NoIntersection", 0, 0
    else
        M = [l1.x2-l1.x1 -(l2.x2-l2.x1)
            l1.y2-l1.y1 -(l2.y2-l2.y1)]
        b = [l2.x1-l1.x1, l2.y1-l1.y1]
        if det(M) == 0 #Lines are parallel
            return "NoIntersection", 0, 0
        end
        try #The determinant check above should handle all the error scenarios
            x= M \ b
            p = Point(l1.x1+x[1]*l1.vec[1], l1.y1+x[1]*l1.vec[2])
            if checkinsegment(p, l2) && x[1] >=0
                return p, x[1], x[2]
            else
                return "NoIntersection", x[1], x[2]
            end
        catch e
            display(e)
            display(M)
            display(b)
            return "NoIntersection", 0, 0
        end

    end
end

function checkinsegment(p::Point, l::LineSegment)
    p1 = Point(l.x1, l.y1)
    p2 = Point(l.x2, l.y2)
    d = pointToPointDist(p1, p2)

    return (pointToPointDist(p, p1) <= d && pointToPointDist(p, p2) <= d)

end

function pointToPointDist(p1::Point, p2::Point)
    return norm([p2.x-p1.x, p2.y-p1.y])
end

function pointToPointDist(p1::Ray, p2::Point)
    return norm([p2.x-p1.x, p2.y-p1.y])
end

function pointToPointDist(p1::Ray, p2::Ray)
    return norm([p2.x-p1.x, p2.y-p1.y])
end

function translate(p::Point, v::Vector)
    return Point(p.x + v[1], p.y + v[2])

end

function getangle(l1::LineSegment, l2::LineSegment)
    theta = acos((dot(l1.vec, l2.vec))/(norm(l1.vec)*norm(l2.vec)))
end

function getangle(l1::LineSegment, l2::Line)
    theta = acos((dot(l1.vec, l2.vec))/(norm(l1.vec)*norm(l2.vec)))
end

function getangle(r::Ray, l2::Line)
    getangle(raytosegment(r), l2)
end

function getangle(r::Ray, l2::LineSegment)
    getangle(raytosegment(r), l2)
end

function remapangle(theta)
    theta = theta % (2*pi)

    return theta < 0 ? theta+2*pi : theta
end

function testintersect()
    p1 = Point(0, 0)
    p2 = Point(0, 1)
    p3 = Point(1, 1)
    p4 = Point(1, 0)

    l1 = LineSegment(p1, p3)
    l2 = LineSegment(p2, p4)
    display(intersect(l1, l2))

    l3 = LineSegment(p1, p2)
    l4 = LineSegment(p3, p4)
    display(intersect(l3, l4))
end

end