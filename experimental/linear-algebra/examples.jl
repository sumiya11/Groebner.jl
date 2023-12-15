using Nemo, Plots

s1 = 20
s2 = 16
begin
    p = plot(
        xlims=(0.0, 5.0),
        ylims=(0.0, 4),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!([(1, 1), (4, 1), (2, 3)], color=:blue, markersize=8, alpha=0.7, label=nothing)
    plot!(
        [(1, 1), (4, 1), (2, 3), (1, 1)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    savefig("ex1.png")
    # p
end

begin
    p = plot(
        xlims=(0.0, 5.0),
        ylims=(0.0, 4),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!(
        [(0, 0), (3, 0), (4, 3), (2, 2)],
        color=:blue,
        markersize=8,
        alpha=0.7,
        label=nothing
    )
    plot!(
        [(0, 0), (3, 0), (4, 3), (2, 2), (0, 0)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    savefig("ex2.png")
    #    p
end

begin
    p = plot(
        xlims=(0.0, 3),
        ylims=(0.0, 3),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!(
        [(0, 0), (0, 1), (0, 2), (2, 2), (1, 1), (1, 2)],
        color=:blue,
        markersize=8,
        alpha=0.7,
        label=nothing
    )
    plot!(
        [(0, 0), (2, 2), (0, 2), (0, 0)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    plot!(0:3, x -> x, linestyle=:dash, color=:red, linewidth=4, alpha=0.3, label=nothing)
    savefig("ex3.png")
    p
end

begin
    p = plot(
        xlims=(0.0, 3),
        ylims=(0.0, 3),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!([(0, 0), (0, 1), (1, 1)], color=:blue, markersize=8, alpha=0.7, label=nothing)
    plot!(
        [(0, 0), (0, 1), (1, 1), (0, 0)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    savefig("ex3-1.png")
    p
end

begin
    p = plot(
        xlims=(0.0, 3),
        ylims=(0.0, 3),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!(
        [(0, 0), (2, 0), (1, 1), (1, 0), (0, 1), (0, 2)],
        color=:blue,
        markersize=8,
        alpha=0.7,
        label=nothing
    )
    plot!(
        [(0, 0), (2, 0), (0, 2), (0, 0)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    savefig("ex4.png")
    p
end

begin
    p = plot(
        xlims=(0.0, 5),
        ylims=(0.0, 5),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!(
        [(0, 0), (4, 2), (3, 2), (2, 3), (2, 2), (2, 1), (1, 3), (1, 1), (0, 4), (0, 2)],
        color=:blue,
        markersize=8,
        alpha=0.7,
        label=nothing
    )
    plot!(
        [(0, 0), (4, 2), (0, 4)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    savefig("ex5.png")
    p
end

begin
    p = plot(
        xlims=(0.0, 3),
        ylims=(0.0, 3),
        dpi=200,
        xlabel="x",
        ylabel="t",
        xguidefontsize=s1,
        yguidefontsize=s1,
        legendfontsize=s2,
        xtickfontsize=s2,
        ytickfontsize=s2,
        xticks=[1, 2, 3, 4, 5],
        yticks=[0, 1, 2, 3, 4]
    )
    scatter!(
        [(0, 0), (2, 1), (1, 1), (0, 2)],
        color=:blue,
        markersize=8,
        alpha=0.7,
        label=nothing
    )
    plot!(
        [(0, 0), (2, 1), (0, 2)],
        color=:blue,
        markersize=8,
        alpha=0.2,
        fill=true,
        label=nothing
    )
    savefig("ex5-1.png")
    p
end

begin
    points = [(1, 1), (4, 1), (2, 3)]

    rad = 8
    piece = pi / 80
    ranges = [
        (center=points[1], range=(asin(2 / sqrt(5)), pi)),
        (center=points[2], range=(0pi, 3pi / 4)),
        (center=points[3], range=(-pi / 4, asin(2 / sqrt(5))))
    ]

    dranges = [(center=r.center, drange=first(r.range):piece:last(r.range)) for r in ranges]

    allranges = []
    for r in dranges
        for s in r.drange
            push!(allranges, (center=r.center, phi=s))
        end
    end

    anim = @animate for i in 1:length(allranges)
        center = allranges[i].center
        phi = allranges[i].phi

        plot(
            xlims=(0.0, 5.0),
            ylims=(0.0, 5.0),
            dpi=200,
            xlabel="x",
            ylabel="y",
            xguidefontsize=16,
            yguidefontsize=16,
            legendfontsize=14,
            xtickfontsize=12,
            ytickfontsize=12,
            xticks=[1, 2, 3, 4, 5],
            yticks=[0, 1, 2, 3, 4, 5]
        )
        scatter!(points, color=:blue, markersize=8, alpha=0.7, label=nothing)
        plot!(
            vcat(points, first(points)),
            color=:blue,
            markersize=8,
            alpha=0.2,
            fill=true,
            label=nothing
        )

        # for j in 1:length(points) - 1
        #     plot!([points[j]...], [points[j+1]...], linestyle=:dash, color=:red)
        # end
        # plot!(, points[1], linestyle=:dash, color=:red)

        plot!(
            0:5,
            x -> 1,
            linestyle=:dash,
            color=:red,
            linewidth=4,
            alpha=0.3,
            label=nothing
        )
        plot!(
            0:5,
            x -> 2x - 1,
            linestyle=:dash,
            color=:red,
            linewidth=4,
            alpha=0.3,
            label=nothing
        )
        plot!(
            0:5,
            x -> -x + 5,
            linestyle=:dash,
            color=:red,
            linewidth=4,
            alpha=0.3,
            label=nothing
        )

        e1 = rad * exp(phi * im)
        e2 = e1 * exp(pi * im)

        from_e = [real(e1), imag(e1)] .+ center
        to_e = [real(e2), imag(e2)] .+ center

        # from_e = map(rationalize, from_e)
        # to_e = map(rationalize, to_e)

        scatter!((from_e...,), label=nothing)
        scatter!((to_e...,), label=nothing)
        plot!(
            [(from_e...,), (to_e...,)],
            linewidth=4,
            color=:black,
            label="\$w \\cdot (x,y) = const\$"
        )
    end
    gif(anim, "anim_fps15.gif", fps=30, loop=0)
end

y = ax + b
w2 * y - w1 * x = 0
y = w1 / w2 * x

# plot!(points, color=:blue, fill=true, alpha=0.2)

ra = range(0.0, 2pi, 35)
rad = 1.0
j = 1

anim = @animate for i in 0:(n - 1)
    plot(xlims=(0.0, 5.0), ylims=(-1.0, 5.0))
    scatter!(orig_points_for_scatter, color=:blue, markersize=4)
    plot!(orig_points_for_scatter, color=:blue, fill=true, alpha=0.2)

    p0 = points[1]
    p1 = points[2]
    p2 = points[3]

    @info "" p0 p1 p2

    hyp01 = sqrt(sum((p1 .- p0) .^ 2))
    cos01 = (p1[1] - p0[1]) / hyp01
    sin01 = (p1[2] - p0[2]) / hyp01

    hyp12 = sqrt(sum((p2 .- p1) .^ 2))
    cos12 = (p2[1] - p1[1]) / hyp12
    sin12 = (p2[2] - p1[2]) / hyp12

    d_phi = ra[j]
    j += 1

    e1 = (cos01 + sin01 * im) * exp(d_phi * im)
    e2 = e1 * exp(pi * im)
    e3 = (cos12 + sin12 * im)

    @info "" angle(e3) angle(e2)

    if angle(e3) >= 0
        if angle(e2) >= 0 && angle(e2) >= angle(e3)
            circshift!(points, +1)
            j = 1
            continue
        end
        if angle(e2) < 0
            # pass
        end
    else
        # angle(e3) < 0
        if angle(e2) <= 0 && angle(e2) >= angle(e3)
            circshift!(points, +1)
            j = 1
            continue
        end
        if angle(e2) > 0
            # pass
        end
    end

    from_e = [real(e1), imag(e1)] .+ p1
    to_e = [real(e2), imag(e2)] .+ p1

    scatter!((from_e...,))
    scatter!((to_e...,))
    plot!([(from_e...,), (to_e...,)], linewidth=3)
end
gif(anim, "anim_fps15.gif", fps=30, loop=0)

# pyplot()

#########
# plot([(1.0, 1.0), (2.0, 2.0), (2, 0), (1, 1)], linewidth=2, fill=true, alpha=0.2)

begin
    ra = range(0.0, 2pi, 35)
    rad = 1.0
    orig_points = [(1, 1), (3, 1), (2, 3)]
    orig_points_for_scatter = [(2, 3), (1, 1), (3, 1)]

    n = length(ra)
    points = copy(orig_points)
    j = 1
    anim = @animate for i in 0:(n - 1)
        global j, orig_points, points, ra, rad, n, orig_points_for_scatter

        plot(xlims=(0.0, 5.0), ylims=(-1.0, 5.0))
        scatter!(orig_points_for_scatter, color=:blue, markersize=4)
        # plot!(orig_points_for_scatter, color=:blue, fill=true, alpha=0.2)

        p0 = points[1]
        p1 = points[2]
        p2 = points[3]

        @info "" p0 p1 p2

        hyp01 = sqrt(sum((p1 .- p0) .^ 2))
        cos01 = (p1[1] - p0[1]) / hyp01
        sin01 = (p1[2] - p0[2]) / hyp01

        hyp12 = sqrt(sum((p2 .- p1) .^ 2))
        cos12 = (p2[1] - p1[1]) / hyp12
        sin12 = (p2[2] - p1[2]) / hyp12

        d_phi = ra[j]
        j += 1

        e1 = (cos01 + sin01 * im) * exp(d_phi * im)
        e2 = e1 * exp(pi * im)
        e3 = (cos12 + sin12 * im)

        @info "" angle(e3) angle(e2)

        if angle(e3) >= 0
            if angle(e2) >= 0 && angle(e2) >= angle(e3)
                circshift!(points, +1)
                j = 1
                continue
            end
            if angle(e2) < 0
                # pass
            end
        else
            # angle(e3) < 0
            if angle(e2) <= 0 && angle(e2) >= angle(e3)
                circshift!(points, +1)
                j = 1
                continue
            end
            if angle(e2) > 0
                # pass
            end
        end

        from_e = [real(e1), imag(e1)] .+ p1
        to_e = [real(e2), imag(e2)] .+ p1

        scatter!((from_e...,))
        scatter!((to_e...,))
        plot!([(from_e...,), (to_e...,)], linewidth=3)
    end
    gif(anim, "anim_fps15.gif", fps=30, loop=0)
end

#########

R, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"], ordering=:deglex)

F = x3 * (x1 + x2^2)^2 * (x1 + x2) * (x1 * x2 + 1)
@info length(F) total_degree(F) degrees(F)

F1, F2 = x1 * x2, divexact(F, x1 * x2)

factor(F2)

##########

R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], ordering=:deglex)

A = x^2 + (y + 1) * z * x + 2
A2 = evaluate(A, [x, y, 3])
A1 = evaluate(A, [x, 4, 9])

##########
