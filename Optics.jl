module Optics

using ..Geom, LinearAlgebra

function moveRay()
    leaveBounds = true
    
    ins = []
    for s in Main.dom.surfaces
        inter = Geom.intersect(Main.current_ray, s)
        if inter[1] != "NoIntersection"
            push!(ins, [inter; s; Geom.pointToPointDist(Main.current_ray, inter[1])])

        end
        #display(ins)
    end

    leaveBounds = isempty(ins)
    if !leaveBounds
        sort!(ins, by=x -> x[3])
        newtheta = reflectAngle(Main.current_ray, ins[1][2])
        #display(Geom.getangle(Main.current_ray, ins[1][2])*180/pi)
        return Geom.Ray(ins[1][1][1], newtheta), leaveBounds
    else
        ins = []
        for b in Main.dom.bounds
            inter = Geom.intersect(Main.current_ray, b)
            #display(inter)
            if inter[1] != "NoIntersection"
                push!(ins, [inter; b; Geom.pointToPointDist(Main.current_ray, inter[1])])
            end
        end
        sort!(ins, by= x -> x[3])
        newtheta = reflectAngle(Main.current_ray, ins[1][2])
        return Geom.Ray(ins[1][1][1], newtheta), leaveBounds
    end

end

function bumpRay()
    return Geom.Ray(Geom.translate(Geom.Point(Main.current_ray.x, Main.current_ray.y), Main.csize/10000*[cos(Main.current_ray.theta), sin(Main.current_ray.theta)]), Main.current_ray.theta)
end

function reflectAngle(r::Geom.Ray, l::Geom.LineSegment)
    normal = normalize([l.vec[2], -l.vec[1]])
    #display(l)
    rhat = [cos(r.theta), sin(r.theta)]
    # display(r)
    if dot(rhat, normal) > 0
        normal = -normal
    end

    rhat2 = rhat - 2*dot(rhat, normal).*normal
    # display(180/pi*r.theta)
    # display(180/pi .*[atan(rhat[2], rhat[1]), atan(rhat2[2], rhat2[1]), atan(normal[2], normal[1])])
    theta = atan(rhat2[2], rhat2[1])
    #display(rhat2)
    # display(180/pi*theta)
    # if rhat2[1] < 0 && rhat2[2] < 0
    #     return theta + pi
    # elseif rhat2[2] < 0 
    #     return theta + 2*pi
    # elseif rhat2[1] < 0
    #     return theta + pi
    # else
    #     return theta
    # end
    return theta
end

end