include("Geom.jl")
include("Construct.jl")
include("Optics.jl")
include("Sources.jl")

import JSON

using .Construct
using Plots
using Statistics

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
    x = []
    y = []
    for s in dom.surfaces
        append!(x, [s.x1, s.x2, NaN])
        append!(y, [s.y1, s.y2, NaN])
    end
    plot!(x, y, label="Surface", linewidth=2)
end

function plotBounds()
    x = []
    y = []
    for s in dom.bounds
        append!(x, [s.x1, s.x2, NaN])
        append!(y, [s.y1, s.y2, NaN])
    end
    plot!(x, y, label="Bounds", linewidth=3)
end

function plotRays(rayset)
    x = []
    y = []
    for r in rayset
        append!(x, r.x)
        append!(y, r.y)
    end
    plot!(x, y, legend=false, linecolor=:firebrick, aspect_ratio=:equal, linestyle=:dot, linealpha=0.75, tickfontsize=RunParams["plotTickFontSize"],
     labelfontsize=RunParams["plotLabelFontSize"], xlabel='x', ylabel='y')
end

function plotRaysVector(rayset)
    x = []
    y = []
    for r in rayset
        append!(x, r.x)
        append!(y, r.y)
    end
    for i = 2:length(x)
        plot!([x[i-1], x[i]], [y[i-1], y[i]], legend=false, aspect_ratio=:equal, linecolor=:firebrick, linestyle=:dot, linealpha=0.75, tickfontsize=RunParams["plotTickFontSize"],
        labelfontsize=RunParams["plotLabelFontSize"], xlabel='x', ylabel='y', arrow=1)
    end
end

function plotBeams(pvec=true)
    if pvec
        for b in beams
            plotRaysVector(b)
        end
    else
        for b in beams
            plotRays(b)
        end
    end
end

function plotNormals()
    for s in dom.surfaces
        x1 = s.x1 + 0.5 * s.vec[1]
        y1 = s.y1 + 0.5 * s.vec[2]
        x2 = x1 + 0.5 * s.vec[2]
        y2 = y1 - 0.5 * s.vec[1]
        plot!([x1, x2], [y1, y2], linewidth=0.5, linecolor=:grey)
    end
end

function plotAll(; simple=true, pvec=true)
    plotSurfaces()
    if !simple
        plotBounds()
    end
    plotBeams(pvec)
    if !simple
        plotNormals()
    end
end


function buildBounds(x1=-5, y1=-5, x2=5, y2=5)
    p1 = Geom.Point(x1, y1)
    p2 = Geom.Point(x1, y2)
    p3 = Geom.Point(x2, y2)
    p4 = Geom.Point(x2, y1)
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
    exits = (last(r_arr).x - first(r_arr).x) > 0
    path_len = 0
    for i in 2:length(r_arr)
        path_len += Geom.pointToPointDist(r_arr[i], r_arr[i-1])
    end
    return path_len, exits

end


function iterateBeams(startrays, max_iter=100)

    leavebd = false
    rays = []
    #println(startrays)
    forward_lengths = []
    back_lengths = []
    for b in startrays
        if RunParams["progPerBeam"]
            print(".")
        end
        n_iter = 0
        global current_ray = b
        push!(rays, current_ray)
        while !leavebd && n_iter < max_iter
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
        leavebd = false
    end

    return beams, forward_lengths, back_lengths

end

function runSim(theta, r, l)

    global dom = Domain(Geom.LineSegment[], Geom.LineSegment[])
    global beams = []
    global current_ray = Geom.Ray(-1, 0.5, 0)
    println("Running theta = $theta")

    reflect_surfs, bx, by = buildkinkedstandard(theta, r, l)


    dx = bx[2] - bx[1]
    dy = by[2] - by[1]
    buildBounds(bx[1] - 0.05 * dx, by[1] - 0.05 * dy, bx[2] + 0.05 * dx, by[2] + 0.05 * dy)

    append!(dom.surfaces, reflect_surfs)

    startrays = Sources.genRays(r, bx, by)

    beams, forward_lengths, back_lengths = iterateBeams(startrays)



    if RunParams["hdBeamPlots"] && !dryRun && !RunParams["summaryPlotOnly"]
        pl = plot(size=(1000, 1000), xlims=plot_xlims, ylims=plot_ylims)
        plotAll(simple=true, pvec=RunParams["plotVectors"])
        savefig(pl, "figs/beamPlots/hd/phononReflect$theta.png")
    end
    if RunParams["beamPlots"] && !dryRun && !RunParams["summaryPlotOnly"]
        #Reset and make low density plots
        beams = beams[1:10:end]
        pl = plot(size=(1000, 1000), xlims=plot_xlims, ylims=plot_ylims)
        plotAll(simple=true, pvec=RunParams["plotVectors"])
        savefig(pl, "figs/beamPlots/LDphononReflect$theta.png")
    end

    println("-")
    return forward_lengths, back_lengths

end

function runSimAngular(theta, r, l)

    global dom = Domain(Geom.LineSegment[], Geom.LineSegment[])
    global beams = []
    global current_ray = Geom.Ray(-1, 0.5, 0)
    print("Running theta = $theta")

    reflect_surfs, bx, by = buildkinkedstandard(theta, r, l)


    dx = bx[2] - bx[1]
    dy = by[2] - by[1]
    buildBounds(bx[1] - 0.05 * dx, by[1] - 0.05 * dy, bx[2] + 0.05 * dx, by[2] + 0.05 * dy)

    append!(dom.surfaces, reflect_surfs)

    startrays = Sources.genRaysAngular(r, bx, by)


    beams, forward_lengths, back_lengths = iterateBeams(startrays)

    if RunParams["hdBeamPlots"] && !dryRun && !RunParams["summaryPlotOnly"]
        pl = plot(size=(1000, 1000), xlims=plot_xlims, ylims=plot_ylims)
        plotAll(simple=true, pvec=RunParams["plotVectors"])
        savefig(pl, "figs/beamPlots/hd/phononReflectRadial$theta.png")
    end
    if RunParams["beamPlots"] && !dryRun && !RunParams["summaryPlotOnly"]
        #Reset and make low density plots
        beams = beams[1:30:end]
        pl = plot(size=(1000, 1000), xlims=plot_xlims, ylims=plot_ylims)
        plotAll(simple=true, pvec=RunParams["plotVectors"])
        savefig(pl, "figs/beamPlots/LDphononReflectRadial$theta.png")
    end

    println("-")
    return forward_lengths, back_lengths

end

function runSimBurst(theta, r, l)

    global dom = Domain(Geom.LineSegment[], Geom.LineSegment[])
    global beams = []
    global current_ray = Geom.Ray(-1, 0.5, 0)
    print("Running theta = $theta")

    reflect_surfs, bx, by = buildkinkedstandard(theta, r, l)


    dx = bx[2] - bx[1]
    dy = by[2] - by[1]
    buildBounds(bx[1] - 0.05 * dx, by[1] - 0.05 * dy, bx[2] + 0.05 * dx, by[2] + 0.05 * dy)

    append!(dom.surfaces, reflect_surfs)

    startrays = Sources.genSingleBurst(r, bx, by)


    beams, forward_lengths, back_lengths = iterateBeams(startrays)

    if RunParams["hdBeamPlots"] && !dryRun && !RunParams["summaryPlotOnly"]
        pl = plot(size=(1000, 1000), xlims=plot_xlims, ylims=plot_ylims)
        plotAll(simple=true, pvec=RunParams["plotVectors"])
        savefig(pl, "figs/beamPlots/hd/phononReflectBurst$theta.png")
    end
    if RunParams["beamPlots"] && !dryRun && !RunParams["summaryPlotOnly"]
        #Reset and make low density plots
        beams = beams[1:3:end]
        pl = plot(size=(1000, 1000), xlims=plot_xlims, ylims=plot_ylims)
        plotAll(simple=true, pvec=RunParams["plotVectors"])
        savefig(pl, "figs/beamPlots/LDphononReflectBurst$theta.png")
    end

    println("-")
    return forward_lengths, back_lengths

end

function main()
    global RunParams = JSON.parsefile("RunParams.json")
    global dryRun = RunParams["dryRun"]
    if !dryRun
        gr()
        pl = plot()
    end
    r_sim = 2 / 3
    l_sim = 2


    global beamPlots = RunParams["beamPlots"]
    global hdBeamPlots = RunParams["hdBeamPlots"]
    println("Creating directories")
    if !dryRun
        println(mkpath("figs/beamPlots/hd"))
    end
    global plot_ylims = (-0.25, 2 * r_sim + 2.1 * l_sim * sind(60) + 0.5)
    global plot_xlims = (-0.25, 6 * l_sim + 0.5)

    x = []
    m = []
    frac = []
    if RunParams["runLinearSim"]
        for i in RunParams["startTheta"]:RunParams["deltaTheta"]:RunParams["endTheta"]
            append!(x, i)
            beamdata = runSim(i, r_sim, l_sim)
            append!(m, mean(beamdata[1]))
            append!(frac, length(beamdata[1]) / (length(beamdata[1]) + length(beamdata[2])))
        end
        #println(length(x), length(m), length(frac))
        if !dryRun
            pl = plot(x, m, legend=false, xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/pathLengthvsAngle.png")

            pl = plot(x, frac, legend=false, xlabel="Angle (degrees)", ylabel="Fraction of Non-Backscattered Beams", tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/backScattervsAngle.png")

            pl = plot(x, m, size=(1000, 800), xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", label="Path Length", c=:blue,
                legend=false, y_guidefontcolor=:blue, x_ticks=0:15:75, y_foreground_color_axis=:blue, right_margin=10Plots.mm, tickfontsize=RunParams["plotTickFontSize"])
            plot!(twinx(), x, frac, label="Non-Backscattered Fraction", ylabel="Fraction of Non-Backscattered Beams", c=:red,
                legend=false, y_guidefontcolor=:red, x_ticks=0:15:75, y_foreground_color_axis=:red, right_margin=10Plots.mm, tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/dualvsAngle.png")
        end
    end

    if RunParams["runRadialSim"]
        for i in RunParams["startTheta"]:RunParams["deltaTheta"]:RunParams["endTheta"]
            append!(x, i)
            beamdata = runSimAngular(i, r_sim, l_sim)
            append!(m, mean(beamdata[1]))
            append!(frac, length(beamdata[1]) / (length(beamdata[1]) + length(beamdata[2])))
        end
        #println(length(x), length(m), length(frac))
        if !dryRun
            pl = plot(x, m, legend=false, xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/pathLengthvsAngleRadial.png")

            pl = plot(x, frac, legend=false, xlabel="Angle (degrees)", ylabel="Fraction of Non-Backscattered Beams", tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/backScattervsAngleRadial.png")

            pl = plot(x, m, size=(1000, 800), xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", label="Path Length", c=:blue,
                legend=false, y_guidefontcolor=:blue, x_ticks=0:15:75, y_foreground_color_axis=:blue, right_margin=10Plots.mm, tickfontsize=RunParams["plotTickFontSize"])
            plot!(twinx(), x, frac, label="Non-Backscattered Fraction", ylabel="Fraction of Non-Backscattered Beams", c=:red,
                legend=false, y_guidefontcolor=:red, x_ticks=0:15:75, y_foreground_color_axis=:red, right_margin=10Plots.mm, tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/dualvsAngleRadial.png")
        end
    end

    if RunParams["runBurstSim"]
        for i in RunParams["startTheta"]:RunParams["deltaTheta"]:RunParams["endTheta"]
            append!(x, i)
            beamdata = runSimBurst(i, r_sim, l_sim)
            append!(m, mean(beamdata[1]))
            append!(frac, length(beamdata[1]) / (length(beamdata[1]) + length(beamdata[2])))
        end
        #println(length(x), length(m), length(frac))
        if !dryRun
            pl = plot(x, m, legend=false, xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/pathLengthvsAngleBurst.png")

            pl = plot(x, frac, legend=false, xlabel="Angle (degrees)", ylabel="Fraction of Non-Backscattered Beams", tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/backScattervsAngleBurst.png")

            pl = plot(x, m, size=(1000, 800), xlabel="Angle (degrees)", ylabel="Mean Path Length, No Backscatter", label="Path Length", c=:blue,
                legend=false, y_guidefontcolor=:blue, x_ticks=0:15:75, y_foreground_color_axis=:blue, right_margin=10Plots.mm, tickfontsize=RunParams["plotTickFontSize"])
            plot!(twinx(), x, frac, label="Non-Backscattered Fraction", ylabel="Fraction of Non-Backscattered Beams", c=:red,
                legend=false, y_guidefontcolor=:red, x_ticks=0:15:75, y_foreground_color_axis=:red, right_margin=10Plots.mm, tickfontsize=RunParams["plotTickFontSize"])
            savefig(pl, "figs/dualvsAngleBurst.png")
        end
    end
end

main()




