#=
Created on Sat 30 Sep 2020
Updated on Fri 16 Feb 2024
-------------------------------------------------------------------------------
This software, like the language it is written in, is published under the MIT
License, https://opensource.org/licenses/MIT.

Copyright (c) 2023:
Alan Freed and John Clayton

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
-------------------------------------------------------------------------------

This test program verifies the B-spline implemented in BSplineKit.jl by
creating a B-spline for a sine wave, whose first derivative should produce a
cosine wave, and whose second derivative should produce a negative sine wave.
=#

#------------------------------------------------------------------------------

module testBSpline

using
    BSplineKit,
    CairoMakie

export
    run

function run(knots::Integer, figurePath::String)
    xₖ   = zeros(Float64, knots)
    yₖ   = zeros(Float64, knots)

    nodes = knots - 1
    xₙ    = zeros(Float64, nodes)
    eₙ    = zeros(Float64, nodes)
    e′ₙ   = zeros(Float64, nodes)
    e′′ₙ  = zeros(Float64, nodes)
    yₙ    = zeros(Float64, nodes)
    y′ₙ   = zeros(Float64, nodes)
    y′′ₙ  = zeros(Float64, nodes)
    zₙ    = zeros(Float64, nodes)
    z′ₙ   = zeros(Float64, nodes)
    z′′ₙ  = zeros(Float64, nodes)

    dx = 0.5π / (knots - 1)

    xₖ[1] = 0.0
    yₖ[1] = sin(xₖ[1])
    for n = 2:knots
        xₖ[n] = xₖ[n-1] + dx
        yₖ[n] = sin(xₖ[n])
    end

    # Nodes are at the mid-points between neighboring knots.
    xₙ[1]   = 0.5dx
    yₙ[1]   = sin(xₙ[1])
    y′ₙ[1]  = cos(xₙ[1])
    y′′ₙ[1] = -sin(xₙ[1])
    for n = 2:nodes
        xₙ[n]   = xₙ[n-1] + dx
        yₙ[n]   = sin(xₙ[n])
        y′ₙ[n]  = cos(xₙ[n])
        y′′ₙ[n] = -sin(xₙ[n])
    end

    # Spline the knots and its first two derivatives.
    s   = interpolate(xₖ, yₖ, BSplineOrder(4))
    s′  = diff(s)
    s′′ = diff(s′)

    # Compute spline error at the nodes, i.e., mid points between knots.
    for n = 1:nodes
        eₙ[n]   = abs(s(xₙ[n]) - yₙ[n])
        e′ₙ[n]  = abs(s′(xₙ[n]) - y′ₙ[n])
        e′′ₙ[n] = abs(s′′(xₙ[n]) - y′′ₙ[n])
        zₙ[n]   = s(xₙ[n])
        z′ₙ[n]  = s′(xₙ[n])
        z′′ₙ[n] = s′′(xₙ[n])
    end
    CairoMakie.activate!(type = "png")
    fig = Figure(; size = (1000, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = Axis(fig[1, 1];
        title  = "Values at Mid Points",
        xlabel = "x",
        ylabel = "y",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax1, xₙ, zₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "y=sin(x)")
    lines!(ax1, xₙ, z′ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "y′=cos(x)")
    lines!(ax1, xₙ, z′′ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "y′′=-sin(x)")
    axislegend("Locations",
        position = :lb)

    ax2 = Axis(fig[1,2];
        title  = "Errors at Mid Points",
        xlabel = "x",
        ylabel = "error",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20,
        yscale = log10)
    lines!(ax2, xₙ, eₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "y=sin(x)")
    lines!(ax2, xₙ, e′ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "y′=cos(x)")
    lines!(ax2, xₙ, e′′ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "y′′=-sin(x)")
    axislegend("Locations",
        position = :lc)

    figPath = string(figurePath, "testBSpline.png")
    save(figPath, fig)
end # run

end # testBSpline