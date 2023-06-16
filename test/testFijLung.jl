#=
Created on Sat 19 Feb 2022
Updated on Sun 23 Apr 2023
-------------------------------------------------------------------------------
This software, like the language it is written in, is published under the MIT
License, https://opensource.org/licenses/MIT.

Copyright (c) 2022:
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
=#
"""
Module:\n
    testFijLung\n
This is a test module for testing the FijLung module. It exports:
    figures2D
    figures3D
"""
module testFijLung

using
    CairoMakie,       # A pixel based figure construction.
    PhysicalFields,
    FijLung

import
    FijLung:
        F_loc1,
        F_loc2,
        F_loc3,
        t_loc1,
        t_loc2,
        t_loc3,
        SplineF,
        newSplineF

export
    figures2D
#   figures3D
#=
-------------------------------------------------------------------------------
=#
"""
Function:\n 
    figures2D(N)\n
where\n
    `N` is the number of nodes in the B-spline.\n 
This function tests the exported functions of FijLung for a 2D case; specifically, for the 23 plane.
"""
function figures2D(N::Integer)
    CairoMakie.activate!(type = "png")
    # These deformation gradient components associate with the 23 plane.

    println("First we create the figures using the raw data.")

    Fᵢⱼs1 = F_loc1()
    Fᵢⱼs2 = F_loc2()
    Fᵢⱼs3 = F_loc3()
    time1 = t_loc1()
    time2 = t_loc2()
    time3 = t_loc3()
    N₁ = Fᵢⱼs1.e
    N₂ = Fᵢⱼs2.e
    N₃ = Fᵢⱼs3.e
    t1 = zeros(Float64, N₁)
    for n in 1:N₁
        t1[n] = get(time1[n])
    end
    t2 = zeros(Float64, N₂)
    for n in 1:N₂
        t2[n] = get(time2[n])
    end
    t3 = zeros(Float64, N₃)
    for n in 1:N₃
        t3[n] = get(time3[n])
    end
    print("For these figures, N₁ = ", string(N₁), ", N₂ = ", string(N₂))
    println(" and N₃ = ", string(N₃), ".")

    # Create a figure for F₁₁.
    println("Working on figure F₁₁ for 2D.")
    F₁₁1 = zeros(Float64, N₁)
    F₁₁2 = zeros(Float64, N₂)
    F₁₁3 = zeros(Float64, N₃)
    for n in 1:N₁
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₁₁1[n] = get(Fᵢⱼ1[2,2])
    end
    for n in 1:N₂
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₁₁2[n] = get(Fᵢⱼ2[2,2])
    end
    for n in 1:N₃
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₁₁3[n] = get(Fᵢⱼ3[2,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 11 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₁₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₁₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₁₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF11.png", fig)

    # Create a figure for F₁₂.
    println("Working on figure F₁₂ for 2D.")
    F₁₂1 = zeros(Float64, N₁)
    F₁₂2 = zeros(Float64, N₂)
    F₁₂3 = zeros(Float64, N₃)
    for n in 1:N₁
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₁₂1[n] = get(Fᵢⱼ1[2,3])
    end
    for n in 1:N₂
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₁₂2[n] = get(Fᵢⱼ2[2,3])
    end
    for n in 1:N₃
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₁₂3[n] = get(Fᵢⱼ3[2,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 12 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₁₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₁₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₁₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :cb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF12.png", fig)

    # Create a figure for F₂₁.
    println("Working on figure F₂₁ for 2D.")
    F₂₁1 = zeros(Float64, N₁)
    F₂₁2 = zeros(Float64, N₂)
    F₂₁3 = zeros(Float64, N₃)
    for n in 1:N₁
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₂₁1[n] = get(Fᵢⱼ1[3,2])
    end
    for n in 1:N₂
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₂₁2[n] = get(Fᵢⱼ2[3,2])
    end
    for n in 1:N₃
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₂₁3[n] = get(Fᵢⱼ3[3,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 21 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₂₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₂₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₂₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF21.png", fig)

    # Create a figure for F₂₂.
    println("Working on figure F₂₂ for 2D.")
    F₂₂1 = zeros(Float64, N₁)
    F₂₂2 = zeros(Float64, N₂)
    F₂₂3 = zeros(Float64, N₃)
    for n in 1:N₁
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₂₂1[n] = get(Fᵢⱼ1[3,3])
    end
    for n in 1:N₂
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₂₂2[n] = get(Fᵢⱼ2[3,3])
    end
    for n in 1:N₃
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₂₂3[n] = get(Fᵢⱼ3[3,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 22 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₂₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₂₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₂₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF22.png", fig)

    # Create a figure for det(F).
    println("Working on figure det(F) for 2D.")
    detF1 = zeros(Float64, N₁)
    detF2 = zeros(Float64, N₂)
    detF3 = zeros(Float64, N₃)
    for n in 1:N₁
        Fᵢⱼ1 = Fᵢⱼs1[n]
        determinant = Fᵢⱼ1[2,2] * Fᵢⱼ1[3,3] - Fᵢⱼ1[3,2] * Fᵢⱼ1[2,3]
        detF1[n] = get(determinant) - 1
    end
    for n in 1:N₂
        Fᵢⱼ2 = Fᵢⱼs2[n]
        determinant = Fᵢⱼ2[2,2] * Fᵢⱼ2[3,3] - Fᵢⱼ2[3,2] * Fᵢⱼ2[2,3]
        detF2[n] = get(determinant) - 1
    end
    for n in 1:N₃
        Fᵢⱼ3 = Fᵢⱼs3[n]
        determinant = Fᵢⱼ3[2,2] * Fᵢⱼ3[3,3] - Fᵢⱼ3[3,2] * Fᵢⱼ3[2,3]
        detF3[n] = get(determinant) - 1
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Determinant of 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "det(F) - 1",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, detF1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, detF2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, detF3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DdetF.png", fig)

    # Create a figure for tr(F).
    println("Working on figure tr(F) for 2D.")
    trF1 = zeros(Float64, N₁)
    trF2 = zeros(Float64, N₂)
    trF3 = zeros(Float64, N₃)
    for n in 1:N₁
        Fᵢⱼ1 = Fᵢⱼs1[n]
        trF1[n] = get(Fᵢⱼ1[2,2] + Fᵢⱼ1[3,3]) - 2
    end
    for n in 1:N₂
        Fᵢⱼ2 = Fᵢⱼs2[n]
        trF2[n] = get(Fᵢⱼ2[2,2] + Fᵢⱼ2[3,3]) - 2
    end
    for n in 1:N₃
        Fᵢⱼ3 = Fᵢⱼs3[n]
        trF3[n] = get(Fᵢⱼ3[2,2] + Fᵢⱼ3[3,3]) - 2
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Trace of 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "tr(F) - 2",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, trF1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, trF2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, trF3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DtrF.png", fig)

    println("Now we recreate these figures using B-spline data.")

    splineF1 = newSplineF(1, N)
    splineF2 = newSplineF(2, N)
    splineF3 = newSplineF(3, N)
    t1 = zeros(Float64, N)
    t2 = zeros(Float64, N)
    t3 = zeros(Float64, N)
    for n in 1:N
        t1[n] = get(splineF1.t[n])
        t2[n] = get(splineF2.t[n])
        t3[n] = get(splineF3.t[n])
    end
    println("For these figures, N = ", string(N), ".")
    # Create figures for F₁₁ and its derivatives.
    println("Working on figure F₁₁ for 2D.")
    F₁₁1 = zeros(Float64, N)
    F₁₁2 = zeros(Float64, N)
    F₁₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₁₁1[n] = get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = splineF2.F[n]
        F₁₁2[n] = get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = splineF3.F[n]
        F₁₁3[n] = get(Fᵢⱼ3[2,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 11 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₁₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₁₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₁₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF11splined.png", fig)

    # Create a figure for dF₁₁/dt.
    println("Working on figure dF₁₁/dt for 2D.")
    dF₁₁1 = zeros(Float64, N)
    dF₁₁2 = zeros(Float64, N)
    dF₁₁3 = zeros(Float64, N)
    for n in 1:N
        dFᵢⱼ1 = splineF1.F′[n]
        dF₁₁1[n] = get(dFᵢⱼ1[2,2])
        dFᵢⱼ2 = splineF2.F′[n]
        dF₁₁2[n] = get(dFᵢⱼ2[2,2])
        dFᵢⱼ3 = splineF3.F′[n]
        dF₁₁3[n] = get(dFᵢⱼ3[2,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "dF₁₁/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, dF₁₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, dF₁₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, dF₁₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DdF11splined.png", fig)

    # Create a figure for d²F₁₁/dt².
    println("Working on figure d²F₁₁/dt² for 2D.")
    d²F₁₁1 = zeros(Float64, N)
    d²F₁₁2 = zeros(Float64, N)
    d²F₁₁3 = zeros(Float64, N)
    for n in 1:N
        d²Fᵢⱼ1 = splineF1.F′′[n]
        d²F₁₁1[n] = get(d²Fᵢⱼ1[2,2])
        d²Fᵢⱼ2 = splineF2.F′′[n]
        d²F₁₁2[n] = get(d²Fᵢⱼ2[2,2])
        d²Fᵢⱼ3 = splineF3.F′′[n]
        d²F₁₁3[n] = get(d²Fᵢⱼ3[2,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "d²F₁₁/dt² (s⁻²)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, d²F₁₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, d²F₁₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, d²F₁₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2Dd2F11splined.png", fig)

    # Create figures for F₁₂ and its derivatives.
    println("Working on figure F₁₂ for 2D.")
    F₁₂1 = zeros(Float64, N)
    F₁₂2 = zeros(Float64, N)
    F₁₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₁₂1[n] = get(Fᵢⱼ1[2,3])
        Fᵢⱼ2 = splineF2.F[n]
        F₁₂2[n] = get(Fᵢⱼ2[2,3])
        Fᵢⱼ3 = splineF3.F[n]
        F₁₂3[n] = get(Fᵢⱼ3[2,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 12 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₁₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₁₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₁₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :cb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF12splined.png", fig)

    println("Working on figure dF₁₂/dt for 2D.")
    dF₁₂1 = zeros(Float64, N)
    dF₁₂2 = zeros(Float64, N)
    dF₁₂3 = zeros(Float64, N)
    for n in 1:N
        dFᵢⱼ1 = splineF1.F′[n]
        dF₁₂1[n] = get(dFᵢⱼ1[2,3])
        dFᵢⱼ2 = splineF2.F′[n]
        dF₁₂2[n] = get(dFᵢⱼ2[2,3])
        dFᵢⱼ3 = splineF3.F′[n]
        dF₁₂3[n] = get(dFᵢⱼ3[2,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "dF₁₂/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, dF₁₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, dF₁₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, dF₁₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DdF12splined.png", fig)

    println("Working on figure d²F₁₂/dt² for 2D.")
    d²F₁₂1 = zeros(Float64, N)
    d²F₁₂2 = zeros(Float64, N)
    d²F₁₂3 = zeros(Float64, N)
    for n in 1:N
        d²Fᵢⱼ1 = splineF1.F′′[n]
        d²F₁₂1[n] = get(d²Fᵢⱼ1[2,3])
        d²Fᵢⱼ2 = splineF2.F′′[n]
        d²F₁₂2[n] = get(d²Fᵢⱼ2[2,3])
        d²Fᵢⱼ3 = splineF3.F′′[n]
        d²F₁₂3[n] = get(d²Fᵢⱼ3[2,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "d²F₁₂/dt² (s⁻²)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, d²F₁₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, d²F₁₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, d²F₁₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2Dd2F12splined.png", fig)

    # Create a figure for F₂₁.
    println("Working on figure F₂₁ for 2D.")
    F₂₁1 = zeros(Float64, N)
    F₂₁2 = zeros(Float64, N)
    F₂₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₂₁1[n] = get(Fᵢⱼ1[3,2])
        Fᵢⱼ2 = splineF2.F[n]
        F₂₁2[n] = get(Fᵢⱼ2[3,2])
        Fᵢⱼ3 = splineF3.F[n]
        F₂₁3[n] = get(Fᵢⱼ3[3,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 21 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₂₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₂₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₂₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF21splined.png", fig)

    println("Working on figure dF₂₁/dt for 2D.")
    dF₂₁1 = zeros(Float64, N)
    dF₂₁2 = zeros(Float64, N)
    dF₂₁3 = zeros(Float64, N)
    for n in 1:N
        dFᵢⱼ1 = splineF1.F′[n]
        dF₂₁1[n] = get(dFᵢⱼ1[3,2])
        dFᵢⱼ2 = splineF2.F′[n]
        dF₂₁2[n] = get(dFᵢⱼ2[3,2])
        dFᵢⱼ3 = splineF3.F′[n]
        dF₂₁3[n] = get(dFᵢⱼ3[3,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "dF₂₁/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, dF₂₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, dF₂₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, dF₂₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DdF21splined.png", fig)

    println("Working on figure d²F₂₁/dt² for 2D.")
    d²F₂₁1 = zeros(Float64, N)
    d²F₂₁2 = zeros(Float64, N)
    d²F₂₁3 = zeros(Float64, N)
    for n in 1:N
        d²Fᵢⱼ1 = splineF1.F′′[n]
        d²F₂₁1[n] = get(d²Fᵢⱼ1[3,2])
        d²Fᵢⱼ2 = splineF2.F′′[n]
        d²F₂₁2[n] = get(d²Fᵢⱼ2[3,2])
        d²Fᵢⱼ3 = splineF3.F′′[n]
        d²F₂₁3[n] = get(d²Fᵢⱼ3[3,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "d²F₂₁/dt² (s⁻²)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, d²F₂₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, d²F₂₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, d²F₂₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2Dd2F21splined.png", fig)

    # Create figures for F₂₂ and its derivatives.
    println("Working on figure F₂₂ for 2D.")
    F₂₂1 = zeros(Float64, N)
    F₂₂2 = zeros(Float64, N)
    F₂₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        F₂₂1[n] = get(Fᵢⱼ1[3,3])
        Fᵢⱼ2 = splineF2.F[n]
        F₂₂2[n] = get(Fᵢⱼ2[3,3])
        Fᵢⱼ3 = splineF3.F[n]
        F₂₂3[n] = get(Fᵢⱼ3[3,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 22 from 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, F₂₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, F₂₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, F₂₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DF22splined.png", fig)

    println("Working on figure dF₂₂/dt for 2D.")
    dF₂₂1 = zeros(Float64, N)
    dF₂₂2 = zeros(Float64, N)
    dF₂₂3 = zeros(Float64, N)
    for n in 1:N
        dFᵢⱼ1 = splineF1.F′[n]
        dF₂₂1[n] = get(dFᵢⱼ1[3,3])
        dFᵢⱼ2 = splineF2.F′[n]
        dF₂₂2[n] = get(dFᵢⱼ2[3,3])
        dFᵢⱼ3 = splineF3.F′[n]
        dF₂₂3[n] = get(dFᵢⱼ3[3,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "dF₂₂/dt (s⁻¹)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, dF₂₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, dF₂₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, dF₂₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DdF22splined.png", fig)

    println("Working on figure d²F₂₂/dt for 2D.")
    d²F₂₂1 = zeros(Float64, N)
    d²F₂₂2 = zeros(Float64, N)
    d²F₂₂3 = zeros(Float64, N)
    for n in 1:N
        d²Fᵢⱼ1 = splineF1.F′′[n]
        d²F₂₂1[n] = get(d²Fᵢⱼ1[3,3])
        d²Fᵢⱼ2 = splineF2.F′′[n]
        d²F₂₂2[n] = get(d²Fᵢⱼ2[3,3])
        d²Fᵢⱼ3 = splineF3.F′′[n]
        d²F₂₂3[n] = get(d²Fᵢⱼ3[3,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time (s)",
        ylabel = "d²F₂₂/dt² (s⁻²)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, d²F₂₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, d²F₂₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, d²F₂₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2Dd2F22splined.png", fig)

    # Create a figure for det(F).
    println("Working on figure det(F) for 2D.")
    detF1 = zeros(Float64, N)
    detF2 = zeros(Float64, N)
    detF3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        detF1[n] = get(Fᵢⱼ1[2,2] * Fᵢⱼ1[3,3] - Fᵢⱼ1[3,2] * Fᵢⱼ1[2,3]) - 1
        Fᵢⱼ2 = splineF2.F[n]
        detF2[n] = get(Fᵢⱼ2[2,2] * Fᵢⱼ2[3,3] - Fᵢⱼ2[3,2] * Fᵢⱼ2[2,3]) - 1
        Fᵢⱼ3 = splineF3.F[n]
        detF3[n] = get(Fᵢⱼ3[2,2] * Fᵢⱼ3[3,3] - Fᵢⱼ3[3,2] * Fᵢⱼ3[2,3]) - 1
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Determinant of 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "det(F) - 1",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, detF1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, detF2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, detF3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DdetFsplined.png", fig)

    # Create a figure for tr(F).
    println("Working on figure tr(F) for 2D.")
    trF1 = zeros(Float64, N)
    trF2 = zeros(Float64, N)
    trF3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = splineF1.F[n]
        trF1[n] = get(Fᵢⱼ1[2,2] + Fᵢⱼ1[3,3]) - 2
        Fᵢⱼ2 = splineF2.F[n]
        trF2[n] = get(Fᵢⱼ2[2,2] + Fᵢⱼ2[3,3]) - 2
        Fᵢⱼ3 = splineF3.F[n]
        trF3[n] = get(Fᵢⱼ3[2,2] + Fᵢⱼ3[3,3]) - 2
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Trace of 2D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "tr(F) - 2",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t1, trF1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: pleural")
    lines!(ax, t2, trF2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: interior")
    lines!(ax, t3, trF3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/2DtrFsplined.png", fig)
end

function figures3D()
    CairoMakie.activate!(type = "png")
    Fᵢⱼs1 = F_loc1()
    Fᵢⱼs2 = F_loc2()
    Fᵢⱼs3 = F_loc3()
    time1 = t_loc1()
    N = Fᵢⱼs1.e
    rows = Fᵢⱼs1.r
    cols = Fᵢⱼs1.c
    time = zeros(Float64, N)
    for n in 1:N
        time[n] = get(time1[n])
    end
    println("For these figures, N = ", string(N))

    # Create a figure for F₁₁.
    println("Working on figure F₁₁ for 3D.")
    F₁₁1 = zeros(Float64, N)
    F₁₁2 = zeros(Float64, N)
    F₁₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₁₁1[n] = get(Fᵢⱼ1[1,1])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₁₁2[n] = get(Fᵢⱼ2[1,1])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₁₁3[n] = get(Fᵢⱼ3[1,1])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 11 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₁₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₁₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₁₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F11.png", fig)

    # Create a figure for F₁₂.
    println("Working on figure F₁₂ for 3D.")
    F₁₂1 = zeros(Float64, N)
    F₁₂2 = zeros(Float64, N)
    F₁₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₁₂1[n] = get(Fᵢⱼ1[1,2])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₁₂2[n] = get(Fᵢⱼ2[1,2])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₁₂3[n] = get(Fᵢⱼ3[1,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 12 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₁₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₁₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₁₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F12.png", fig)

    # Create a figure for F₁₃.
    println("Working on figure F₁₃ for 3D.")
    F₁₃1 = zeros(Float64, N)
    F₁₃2 = zeros(Float64, N)
    F₁₃3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₁₃1[n] = get(Fᵢⱼ1[1,3])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₁₃2[n] = get(Fᵢⱼ2[1,3])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₁₃3[n] = get(Fᵢⱼ3[1,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 13 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₁₃",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₁₃1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₁₃2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₁₃3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F13.png", fig)

    # Create a figure for F₂₁.
    println("Working on figure F₂₁ for 3D.")
    F₂₁1 = zeros(Float64, N)
    F₂₁2 = zeros(Float64, N)
    F₂₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₂₁1[n] = get(Fᵢⱼ1[2,1])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₂₁2[n] = get(Fᵢⱼ2[2,1])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₂₁3[n] = get(Fᵢⱼ3[2,1])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 21 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₂₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₂₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₂₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F21.png", fig)

    # Create a figure for F₂₂.
    println("Working on figure F₂₂ for 3D.")
    F₂₂1 = zeros(Float64, N)
    F₂₂2 = zeros(Float64, N)
    F₂₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₂₂1[n] = get(Fᵢⱼ1[2,2])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₂₂2[n] = get(Fᵢⱼ2[2,2])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₂₂3[n] = get(Fᵢⱼ3[2,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 22 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₂₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₂₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₂₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F22.png", fig)

    # Create a figure for F₂₃.
    println("Working on figure F₂₃ for 3D.")
    F₂₃1 = zeros(Float64, N)
    F₂₃2 = zeros(Float64, N)
    F₂₃3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₂₃1[n] = get(Fᵢⱼ1[2,3])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₂₃2[n] = get(Fᵢⱼ2[2,3])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₂₃3[n] = get(Fᵢⱼ3[2,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 23 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₂₃",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₂₃1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₂₃2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₂₃3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :cb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F23.png", fig)

    # Create a figure for F₃₁.
    println("Working on figure F₃₁ for 3D.")
    F₃₁1 = zeros(Float64, N)
    F₃₁2 = zeros(Float64, N)
    F₃₁3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₃₁1[n] = get(Fᵢⱼ1[3,1])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₃₁2[n] = get(Fᵢⱼ2[3,1])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₃₁3[n] = get(Fᵢⱼ3[3,1])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 31 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₃₁",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₃₁1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₃₁2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₃₁3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F31.png", fig)

    # Create a figure for F₃₂.
    println("Working on figure F₃₂ for 3D.")
    F₃₂1 = zeros(Float64, N)
    F₃₂2 = zeros(Float64, N)
    F₃₂3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₃₂1[n] = get(Fᵢⱼ1[3,2])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₃₂2[n] = get(Fᵢⱼ2[3,2])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₃₂3[n] = get(Fᵢⱼ3[3,2])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 32 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₃₂",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₃₂1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₃₂2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₃₂3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F32.png", fig)

    # Create a figure for F₃₃.
    println("Working on figure F₃₃ for 3D.")
    F₃₃1 = zeros(Float64, N)
    F₃₃2 = zeros(Float64, N)
    F₃₃3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        F₃₃1[n] = get(Fᵢⱼ1[3,3])
        Fᵢⱼ2 = Fᵢⱼs2[n]
        F₃₃2[n] = get(Fᵢⱼ2[3,3])
        Fᵢⱼ3 = Fᵢⱼs3[n]
        F₃₃3[n] = get(Fᵢⱼ3[3,3])
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Component 33 from 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "F₃₃",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, F₃₃1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, F₃₃2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, F₃₃3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rt)
    save("/home/al/Documents/books/alveolar_mechanics/figures/F33.png", fig)

    # Create a figure for det(F).
    println("Working on figure det(F) for 3D.")
    detF1 = zeros(Float64, N)
    detF2 = zeros(Float64, N)
    detF3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        detF1[n] = get(det(Fᵢⱼ1)) - 1
        Fᵢⱼ2 = Fᵢⱼs2[n]
        detF2[n] = get(det(Fᵢⱼ2)) - 1
        Fᵢⱼ3 = Fᵢⱼs3[n]
        detF3[n] = get(det(Fᵢⱼ3)) - 1
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Determinant of 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "det(F) - 1",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, detF1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, detF2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, detF3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/detF.png", fig)

    # Create a figure for tr(F).
    println("Working on figure tr(F) for 3D.")
    trF1 = zeros(Float64, N)
    trF2 = zeros(Float64, N)
    trF3 = zeros(Float64, N)
    for n in 1:N
        Fᵢⱼ1 = Fᵢⱼs1[n]
        trF1[n] = get(tr(Fᵢⱼ1)) - 3
        Fᵢⱼ2 = Fᵢⱼs2[n]
        trF2[n] = get(tr(Fᵢⱼ2)) - 3
        Fᵢⱼ3 = Fᵢⱼs3[n]
        trF3[n] = get(tr(Fᵢⱼ3)) - 3
    end
    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        title  = "Trace of 3D Deformation Gradients",
        xlabel = "time (s)",
        ylabel = "tr(F) - 3",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, time, trF1;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1: surface")
    lines!(ax, time, trF2;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "2: center")
    lines!(ax, time, trF3;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "3: bronchiole")
    axislegend("Locations",
        position = :rb)
    save("/home/al/Documents/books/alveolar_mechanics/figures/trF.png", fig)
end

end  #testFijLung
