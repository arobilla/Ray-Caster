include("Geom.jl")
include("Construct.jl")
include("Optics.jl")
include("Sources.jl")


using Plots, .Construct, Statistics

struct Domain
    bounds::Array{Geom.LineSegment}
    surfaces::Array{Geom.LineSegment}
end

function checkSystemBoundaries()
    for bd in dom.bounds
        if Geom.intersect(current_ray, bd)[1] != "NoIntersection"
            return false
        end
    end
end




function plotSurfaces()
    x=[]
    y=[]
    for s in dom.surfaces
        append!(x, [s.x1, s.x2, NaN])
        append!(y, [s.y1, s.y2, NaN])
    end
    plot!(x,y, label="Surface", linewidth = 2)
end

function plotBounds()
    x=[]
    y=[]
    for s in dom.bounds
        append!(x, [s.x1, s.x2, NaN])
        append!(y, [s.y1, s.y2, NaN])
    end
    plot!(x,y, label="Bounds", linewidth=3)
end

function plotRays(rayset)
    x = []
    y = []
    for r in rayset
        append!(x, r.x)
        append!(y, r.y)
    end
    plot!(x,y, legend=false)
end

function plotBeams()
    for b in beams
        plotRays(b)
    end
end

function plotNormals()
    for s in dom.surfaces
        x1 = s.x1+0.5*s.vec[1]
        y1 = s.y1+0.5*s.vec[2]
        x2 = x1+0.5*s.vec[2]
        y2 = y1 - 0.5*s.vec[1]
        plot!([x1, x2], [y1,y2], linewidth=0.5, linecolor=:grey)
    end
end
    

function buildBounds(x1=-5, y1=-5, x2=5, y2=5)
    p1 = Geom.Point(x1, y1)
    p2 = Geom.Point(x1, y2)
    p3 = Geom.Point(x2,y2)
    p4 = Geom.Point(x2,y1)
    push!(dom.bounds, Geom.LineSegment(p1, p2))
    push!(dom.bounds, Geom.LineSegment(p2, p3))
    push!(dom.bounds, Geom.LineSegment(p3, p4))
    push!(dom.bounds, Geom.LineSegment(p4, p1))
    global csize = 10
end

"""
Calculates the path length of a full set of beam scattering
Also returns True if the beam leaves opposite entry, or false if it backscatters
"""
function pathLength(r_arr::Array{Any})
    exits =  (last(r_arr).x - first(r_arr).x) > 0
    path_len = 0
    for i in 2:length(r_arr)
        path_len += Geom.pointToPointDist(r_arr[i], r_arr[i-1])
    end
    return path_len, exits

end

function genRays(r, l, bx, by)
    s = []
    for y in LinRange(by[1]+0.01*r, by[1]+2*r-0.01*r, 100)
        push!(s, Geom.Ray(bx[1] + 0.0001, y, 0))
    end
    return s
end

function runSim(theta, r, l)

    global dom = Domain(Geom.LineSegment[], Geom.LineSegment[])
    global beams = []
    global rays = []
    global current_ray = Geom.Ray(-1, 0.5, 0)
    leavebd = false
    println("Running theta = $theta")

    reflect_surfs, bx, by = buildkinkedstandard(theta, r, l)


    dx = bx[2] - bx[1]
    dy = by[2] - by[1]
    buildBounds(bx[1] - 0.001*dx, by[1] - 0.05*dy, bx[2] + 0.001*dx, by[2] + 0.05*dy)

    append!(dom.surfaces, reflect_surfs)

    startrays = genRays(r, l, bx, by)


    #println(startrays)
    forward_lengths = []
    back_lengths = []
    for b in startrays
        print(".")
        n_iter = 0
        current_ray = b
        push!(rays, current_ray)
        while !leavebd && n_iter < 100
            # display(current_ray)
            current_ray = Optics.bumpRay()
            #push!(rays, current_ray)
            next_ray, leavebd = Optics.moveRay()
            #display([next_ray, leavebd])
            current_ray = next_ray
            push!(rays, current_ray)
            n_iter += 1
        end
        push!(beams, rays)
        path_len, exits = pathLength(rays)
        exits ? push!(forward_lengths, path_len) : push!(back_lengths, path_len)
        rays = []
        leavebd=false
    end


    if hdBeamPlots
        pl = plot(size=(500,500), xlims=plot_xlims, ylims = plot_ylims)
        plotSurfaces()
        plotBounds()
        plotBeams()
        plotNormals()
        savefig(pl, "figs/beamPlots/hd/phononReflect$theta.png")
    end
    if beamPlots
    #Reset and make low density plots
        beams = beams[1:10:end]
        pl = plot(size=(500,500), xlims=plot_xlims, ylims = plot_ylims)
        plotSurfaces()
        plotBounds()
        plotBeams()
        plotNormals()
        savefig(pl, "figs/beamPlots/LDphononReflect$theta.png")
    end

    println("-")
    println("Done")
    return forward_lengths, back_lengths

end

function main()
    gr()
    pl=plot()
    global beamPlots = true
    global hdBeamPlots = false
    println("Creating directories")
    println(mkpath("figs/beamPlots/hd"))
    global plot_ylims = (-1, 8)
    global plot_xlims = (-1, 8)

    x = []
    m = []
    frac = []
    for i in 0:1:75
        append!(x, i)
        beamdata = runSim(i, 2/3, 2)
        append!(m, mean(beamdata[1]))
        append!(frac, length(beamdata[1])/(length(beamdata[1])+length(beamdata[2])))
    end
    pl=plot(x, m, legend=false, xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter")
    savefig(pl, "figs/pathLengthvsAngle.png")

    pl=plot(x, frac, legend=false, xlabel="Angle (degrees)", ylabel="Fraction of Non-Backscattered Beams")
    savefig(pl, "figs/backScattervsAngle.png")

    pl=plot(x, m, xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", label="Path Length", c=:blue,
    legend=false, y_guidefontcolor=:blue, y_foreground_color_axis=:blue, right_margin=10Plots.mm)
    plot!(twinx(), frac, label="Non-Backscattered Fraction", ylabel="Fraction of Non-Backscattered Beams", c=:red,
    legend=false, y_guidefontcolor=:red, y_foreground_color_axis=:red, right_margin=10Plots.mm)
    savefig(pl, "figs/dualvsAngle.png")

    # d30 = runSim(30, 2/3, 2)
    # d45 = runSim(45, 2/3, 2)
    # d60 = runSim(60, 2/3, 2)
    # d15 = runSim(15, 2/3, 2)
    # pl = plot(size=(500, 500))
    # bin_range = range(0, 10, 30)
    # through_style = (3, :dashdot)
    # back_style = (3, :dot)
    # stephist!(d15[1], label="15 Degree Through", bins=bin_range, line=through_style, color=:blue)
    # stephist!(d15[2], label="15 Degree Backscatter", bins=bin_range, line=back_style, color=:blue)
    # stephist!(d30[1], label="30 Degree Through", bins=bin_range, line=through_style, color=:red)
    # stephist!(d30[2], label="30 Degree Backscatter", bins=bin_range, line=back_style, color=:red)
    # stephist!(d45[1], label="45 Degree Through", bins=bin_range, line=through_style, color=:green)
    # stephist!(d45[2], label="45 Degree Backscatter", bins=bin_range, line=back_style, color=:green)
    # stephist!(d60[1], label="60 Degree Through", bins=bin_range, line=through_style, color=:purple)
    # stephist!(d60[2], label="60 Degree Backscatter", bins=bin_range, line=back_style, color=:purple)
    # savefig(pl, "figs/pathLengthHisto.png")
end

main()




