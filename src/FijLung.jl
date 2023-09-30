#=
Created on Fri 18 Feb 2022
Updated on Sat 19 Aug 2023
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
=#
"""
Module:\n
    FijLung\n
This module constructs deformation gradient histories at three locations taken from a finite element analysis of a human torso being subjected to a ballistic projectile. Location 1 lies just under the pleura membrane, and therefore resides close to the rib-cage and the location of impact. Location 2 lies deep within the lung. While location 3 lies adjacent to a bronchiole tube. These raw data sets are fit with B-splines to get continuous functions for the deformation gradient, thereby allowing its first and second derivatives in time to be approximated. This module exports:\n
    t_loc1,         # Raw data: times at location near visceral pleura.\n
    t_loc2,         # Raw data: times at location deep within the lung.\n
    t_loc3,         # Raw data: times at location next to a bronchiole tube.\n
    F_loc1,         # Raw data: deformation gradient at location 1.\n
    F_loc2,         # Raw data: deformation gradient at location 2.\n
    F_loc3,         # Raw data: deformation gradient at location 3.\n
    SplineF,        # Type used to hold these spline data.
"""
module FijLung

using
    BSplineKit,
    CSV,
    DataFrames,
    Downloads,
    PhysicalFields

export
    t_loc1,         # Raw data: times at location near visceral pleura.
    t_loc2,         # Raw data: times at location deep within the lung.
    t_loc3,         # Raw data: times at location next to a bronchiole tube.
    F_loc1,         # Raw data: deformation gradient at location 1.
    F_loc2,         # Raw data: deformation gradient at location 2.
    F_loc3,         # Raw data: deformation gradient at location 3.
    SplineF         # Type used to hold the b-spline data.

const DIMENSIONLESS = CGS_DIMENSIONLESS
const SECOND        = CGS_SECOND
const STRETCH       = CGS_STRETCH
const STRETCH_RATE  = CGS_STRETCH_RATE
const STRETCH_ACEL  = CGS_ACCELERATION - CGS_LENGTH
#=
-------------------------------------------------------------------------------
=#
"""
Function:\n
    arrayOfTimes = t_loc1()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 1 (just under the visceral pleura) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc1()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array of times.
    tLoc1 = ArrayOfPhysicalScalars(N, SECOND)
    # Populate this array with times that associate with its N nodes.
    for n in 1:N
        tLoc1[n] = PhysicalScalar(dataF11[n,1], SECOND)
    end
    return tLoc1
end

"""
Function:\n
    arrayOfF = F_loc1()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalTensors`, where each element within the array holds a deformation gradient tensor evaluated at location 1 (just under the visceral pleura) from an FEA of a human torso subjected to a high-energy ballistic impact. The time associated with each entry is supplied by the function `t_loc1()`.
"""
function F_loc1()::ArrayOfPhysicalTensors
    # Open existing files. Read in their data as DataFrames.DataFrame objects.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F12.csv"
    dataF12 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F13.csv"
    dataF13 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F21.csv"
    dataF21 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F22.csv"
    dataF22 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F23.csv"
    dataF23 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F31.csv"
    dataF31 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F32.csv"
    dataF32 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F33.csv"
    dataF33 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Size the arrays.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array with entries of type PhysicalTensor.
    FLoc1 = ArrayOfPhysicalTensors(N, 3, 3, STRETCH)
    for n in 1:N
        Fₙ = PhysicalTensor(3, 3, STRETCH)
        # Populate the deformation gradient with its components.
        Fₙ[1,1] = PhysicalScalar(dataF11[n,2], STRETCH)
        Fₙ[1,2] = PhysicalScalar(dataF12[n,2], STRETCH)
        Fₙ[1,3] = PhysicalScalar(dataF13[n,2], STRETCH)
        Fₙ[2,1] = PhysicalScalar(dataF21[n,2], STRETCH)
        Fₙ[2,2] = PhysicalScalar(dataF22[n,2], STRETCH)
        Fₙ[2,3] = PhysicalScalar(dataF23[n,2], STRETCH)
        Fₙ[3,1] = PhysicalScalar(dataF31[n,2], STRETCH)
        Fₙ[3,2] = PhysicalScalar(dataF32[n,2], STRETCH)
        Fₙ[3,3] = PhysicalScalar(dataF33[n,2], STRETCH)
        # Assign this deformation gradient to the array.
        FLoc1[n] = Fₙ
    end
    return FLoc1
end

"""
Function:\n
    arrayOfTimes = t_loc2()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 2 (roughly at the center of the right lung) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc2()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array of times whose initial time is t₁.
    tLoc2 = ArrayOfPhysicalScalars(N, SECOND)
    # Populate this array with times that associate with its N nodes.
    for n in 1:N
        tLoc2[n] = PhysicalScalar(dataF11[n,1], SECOND)
    end
    return tLoc2
end

"""
Function:\n
    arrayOfF = F_loc2()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalTensors`, where each element within the array holds a deformation gradient tensor evaluated at location 2 (roughly at the center of the right lung) from an FEA of a human torso subjected to a high-energy ballistic impact. The time associated with each entry is supplied by the function `t_loc2()`.
"""
function F_loc2()::ArrayOfPhysicalTensors
    # Open existing files. Read in their data as DataFrames.DataFrame objects.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F12.csv"
    dataF12 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F13.csv"
    dataF13 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F21.csv"
    dataF21 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F22.csv"
    dataF22 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F23.csv"
    dataF23 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F31.csv"
    dataF31 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F32.csv"
    dataF32 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F33.csv"
    dataF33 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Size the arrays.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array with entries of type PhysicalTensor.
    FLoc2 = ArrayOfPhysicalTensors(N, 3, 3, STRETCH)
    for n in 1:N
        Fₙ = PhysicalTensor(3, 3, STRETCH)
        # Populate the deformation gradient with its components.
        Fₙ[1,1] = PhysicalScalar(dataF11[n,2], STRETCH)
        Fₙ[1,2] = PhysicalScalar(dataF12[n,2], STRETCH)
        Fₙ[1,3] = PhysicalScalar(dataF13[n,2], STRETCH)
        Fₙ[2,1] = PhysicalScalar(dataF21[n,2], STRETCH)
        Fₙ[2,2] = PhysicalScalar(dataF22[n,2], STRETCH)
        Fₙ[2,3] = PhysicalScalar(dataF23[n,2], STRETCH)
        Fₙ[3,1] = PhysicalScalar(dataF31[n,2], STRETCH)
        Fₙ[3,2] = PhysicalScalar(dataF32[n,2], STRETCH)
        Fₙ[3,3] = PhysicalScalar(dataF33[n,2], STRETCH)
        # Assign this deformation gradient to the array.
        FLoc2[n] = Fₙ
    end
    return FLoc2
end

"""
Function:\n
    arrayOfTimes = t_loc3()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 3 (adjacent to a bronchiole tube) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc3()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array of times whose initial time is t₁.
    tLoc3 = ArrayOfPhysicalScalars(N, SECOND)
    # Populate this array with times that associate with its N nodes.
    for n in 1:N
        tLoc3[n] = PhysicalScalar(dataF11[n,1], SECOND)
    end
    return tLoc3
end

"""
Function:\n
    arrayOfF = F_loc3()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalTensors`, where each element within the array holds a deformation gradient tensor evaluated at location 3 (adjacent to a bronchiole tube) from an FEA of a human torso subjected to a high-energy ballistic impact. The time associated with each entry is supplied by the function `t_loc3()`.
"""
function F_loc3()::ArrayOfPhysicalTensors
    # Open existing files. Read in their data as DataFrames.DataFrame objects.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F12.csv"
    dataF12 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F13.csv"
    dataF13 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F21.csv"
    dataF21 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F22.csv"
    dataF22 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F23.csv"
    dataF23 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F31.csv"
    dataF31 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F32.csv"
    dataF32 = DataFrame(CSV.File(Downloads.download(dataURL)))
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F33.csv"
    dataF33 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Size the arrays.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Populate the array.
    FLoc3 = ArrayOfPhysicalTensors(N, 3, 3, STRETCH)
    for n in 1:N
        Fₙ = PhysicalTensor(3, 3, STRETCH)
        # Populate the deformation gradient with its components.
        Fₙ[1,1] = PhysicalScalar(dataF11[n,2], STRETCH)
        Fₙ[1,2] = PhysicalScalar(dataF12[n,2], STRETCH)
        Fₙ[1,3] = PhysicalScalar(dataF13[n,2], STRETCH)
        Fₙ[2,1] = PhysicalScalar(dataF21[n,2], STRETCH)
        Fₙ[2,2] = PhysicalScalar(dataF22[n,2], STRETCH)
        Fₙ[2,3] = PhysicalScalar(dataF23[n,2], STRETCH)
        Fₙ[3,1] = PhysicalScalar(dataF31[n,2], STRETCH)
        Fₙ[3,2] = PhysicalScalar(dataF32[n,2], STRETCH)
        Fₙ[3,3] = PhysicalScalar(dataF33[n,2], STRETCH)
        # Assign this deformation gradient to the array.
        FLoc3[n] = Fₙ
    end
    return FLoc3
end
#=
--------------------------------------------------------------------------------
=#
"""
Type:\n
    SplineF\n
        loc     # Location whereat this interpolation applies, loc ∈ {1,2,3}.\n
        N       # Number of uniform-spaced interpolation nodes spanning time.\n
        t       # Time associated with each node in this interpolation.\n
        F       # F         The deformation gradient.\n
        F′      # dF/dt     Its first derivative in time.\n
        F′′     # d²F/dt²   Its second derivative in time.\n
This data structure smooths the raw data supplied by F_loc1, F_loc2 or F_loc3, as specified by parameter `loc`, by applying B-splines to the raw data; specifically, F is fit with a cubic spline, and therefore, its first derivative dF/dt is described by a quadratic spline, and its second derivative d²F/dt² is described by a linear spline. Fields `loc` and `N` are Integers. Field `t` is an instance of type `PhysicalFields.ArrayOfPhysicalScalars`. While fields `F`, `F′` and `F′′` are instances of type `PhysicalFields.ArrayOfPhysicalTensors`.\n
The supplied deformation gradient, and its first two derivatives, are evaluated at one of three lung locations extracted from an FE analysis of a human torso subjected to a ballistic impact. Location loc = 1 associates with a location that is just under the visceral pleura. Location loc = 2 associates with a location that is deep within the lung. While location loc = 3 associates with a location that lies next to a bronchiole tube.
"""
struct SplineF
    loc::Integer                  # loc       Lung location, i.e., 1, 2 or 3.
    N::  Integer                  # N         Number of nodes spanning time.
    t::ArrayOfPhysicalScalars     # t         The time at these nodes.
    F::ArrayOfPhysicalTensors     # F         Their deformation gradients.
    F′::ArrayOfPhysicalTensors    # dF/dt     Their first derivatives in time.
    F′′::ArrayOfPhysicalTensors   # d²F/dt²   Their second derivatives in time.

    # Constructor
    """
    Constructor:\n
        splineF = SplineF(location, nodes)\n
    This constructor creates instances of type SplineF, which are returned as objects `splineF`. Each object contains interpolated values for the deformation gradient F, its first derivative dF/dt, and its second derivative d²F/dt² at `N` uniformly spaced instances in time `t`. These data associate with one of three locations in a lung from an FE analysis of a human torso subjected to a ballistic impact. Location 1 associates with a location that is just under the visceral pleura. Location 2 associates with a location that is deep within the right lung. While location 3 associates with a location that lies next to a bronchiole tube.
    """
    function SplineF(location::Integer, nodes::Integer)
        # Import the data to be fit with B-splines.
        if location == 1
            time = t_loc1()
            rawF = F_loc1()
        elseif location == 2
            time = t_loc2()
            rawF = F_loc2()
        elseif location == 3
            time = t_loc3()
            rawF = F_loc3()
        else
            msg = "An inadmissible location was submitted. It must be 1, 2 or 3."
            throw(ErrorException(msg))
        end
        knots = Int64(time.array.len)
        # Create the various static arrays that will hold these raw data.
        rawTime = zeros(Float64, knots)
        rawF₁₁  = zeros(Float64, knots)
        rawF₁₂  = zeros(Float64, knots)
        rawF₁₃  = zeros(Float64, knots)
        rawF₂₁  = zeros(Float64, knots)
        rawF₂₂  = zeros(Float64, knots)
        rawF₂₃  = zeros(Float64, knots)
        rawF₃₁  = zeros(Float64, knots)
        rawF₃₂  = zeros(Float64, knots)
        rawF₃₃  = zeros(Float64, knots)
        # Populate these arrays for B-spline fitting.
        for k in 1:knots
            rawTime[k] = get(time[k])
            rawFᵢⱼ = rawF[k]
            rawF₁₁[k] = get(rawFᵢⱼ[1,1])
            rawF₁₂[k] = get(rawFᵢⱼ[1,2])
            rawF₁₃[k] = get(rawFᵢⱼ[1,3])
            rawF₂₁[k] = get(rawFᵢⱼ[2,1])
            rawF₂₂[k] = get(rawFᵢⱼ[2,2])
            rawF₂₃[k] = get(rawFᵢⱼ[2,3])
            rawF₃₁[k] = get(rawFᵢⱼ[3,1])
            rawF₃₂[k] = get(rawFᵢⱼ[3,2])
            rawF₃₃[k] = get(rawFᵢⱼ[3,3])
        end

        # Create cubic B-splines for each component of the deformation gradient.
        splineOrder = BSplineOrder(4)
        S₁₁ = interpolate(rawTime, rawF₁₁, splineOrder)
        S₁₂ = interpolate(rawTime, rawF₁₂, splineOrder)
        S₁₃ = interpolate(rawTime, rawF₁₃, splineOrder)
        S₂₁ = interpolate(rawTime, rawF₂₁, splineOrder)
        S₂₂ = interpolate(rawTime, rawF₂₂, splineOrder)
        S₂₃ = interpolate(rawTime, rawF₂₃, splineOrder)
        S₃₁ = interpolate(rawTime, rawF₃₁, splineOrder)
        S₃₂ = interpolate(rawTime, rawF₃₂, splineOrder)
        S₃₃ = interpolate(rawTime, rawF₃₃, splineOrder)

        # Create the arrays that span time that hold these data.
        t   = ArrayOfPhysicalScalars(nodes, SECOND)
        F   = ArrayOfPhysicalTensors(nodes, 3, 3, STRETCH)
        F′  = ArrayOfPhysicalTensors(nodes, 3, 3, STRETCH_RATE)
        F′′ = ArrayOfPhysicalTensors(nodes, 3, 3, STRETCH_ACEL)

        # Establish the time increment between nodes.
        dτ = (rawTime[knots] - rawTime[1]) / (nodes - 1)
        τₙ = rawTime[1] - dτ

        # Populate the arrays holding these kinematic tensor fields.
        for n in 1:nodes
            # Create the scalar and tensor fields for the nᵗʰ node.
            tₙ   = PhysicalScalar(SECOND)
            Fₙ   = PhysicalTensor(3, 3, STRETCH)
            dFₙ  = PhysicalTensor(3, 3, STRETCH_RATE)
            d²Fₙ = PhysicalTensor(3, 3, STRETCH_ACEL)

            # Time tₙ associated with the nᵗʰ state.
            τₙ += dτ
            tₙ = PhysicalScalar(τₙ, SECOND)

            # Update the tensor components at the nᵗʰ state.

            # For the deformation gradient.
            Fₙ[1,1] = PhysicalScalar(S₁₁(τₙ), STRETCH)
            Fₙ[1,2] = PhysicalScalar(S₁₂(τₙ), STRETCH)
            Fₙ[1,3] = PhysicalScalar(S₁₃(τₙ), STRETCH)
            Fₙ[2,1] = PhysicalScalar(S₂₁(τₙ), STRETCH)
            Fₙ[2,2] = PhysicalScalar(S₂₂(τₙ), STRETCH)
            Fₙ[2,3] = PhysicalScalar(S₂₃(τₙ), STRETCH)
            Fₙ[3,1] = PhysicalScalar(S₃₁(τₙ), STRETCH)
            Fₙ[3,2] = PhysicalScalar(S₃₂(τₙ), STRETCH)
            Fₙ[3,3] = PhysicalScalar(S₃₃(τₙ), STRETCH)

            # For the first derivative of the deformation gradient.
            dFₙ[1,1] = PhysicalScalar((Derivative(1)*S₁₁)(τₙ), STRETCH_RATE)
            dFₙ[1,2] = PhysicalScalar((Derivative(1)*S₁₂)(τₙ), STRETCH_RATE)
            dFₙ[1,3] = PhysicalScalar((Derivative(1)*S₁₃)(τₙ), STRETCH_RATE)
            dFₙ[2,1] = PhysicalScalar((Derivative(1)*S₂₁)(τₙ), STRETCH_RATE)
            dFₙ[2,2] = PhysicalScalar((Derivative(1)*S₂₂)(τₙ), STRETCH_RATE)
            dFₙ[2,3] = PhysicalScalar((Derivative(1)*S₂₃)(τₙ), STRETCH_RATE)
            dFₙ[3,1] = PhysicalScalar((Derivative(1)*S₃₁)(τₙ), STRETCH_RATE)
            dFₙ[3,2] = PhysicalScalar((Derivative(1)*S₃₂)(τₙ), STRETCH_RATE)
            dFₙ[3,3] = PhysicalScalar((Derivative(1)*S₃₃)(τₙ), STRETCH_RATE)

            # For the second derivative of the deformation gradient.
            d²Fₙ[1,1] = PhysicalScalar((Derivative(2)*S₁₁)(τₙ), STRETCH_ACEL)
            d²Fₙ[1,2] = PhysicalScalar((Derivative(2)*S₁₂)(τₙ), STRETCH_ACEL)
            d²Fₙ[1,3] = PhysicalScalar((Derivative(2)*S₁₃)(τₙ), STRETCH_ACEL)
            d²Fₙ[2,1] = PhysicalScalar((Derivative(2)*S₂₁)(τₙ), STRETCH_ACEL)
            d²Fₙ[2,2] = PhysicalScalar((Derivative(2)*S₂₂)(τₙ), STRETCH_ACEL)
            d²Fₙ[2,3] = PhysicalScalar((Derivative(2)*S₂₃)(τₙ), STRETCH_ACEL)
            d²Fₙ[3,1] = PhysicalScalar((Derivative(2)*S₃₁)(τₙ), STRETCH_ACEL)
            d²Fₙ[3,2] = PhysicalScalar((Derivative(2)*S₃₂)(τₙ), STRETCH_ACEL)
            d²Fₙ[3,3] = PhysicalScalar((Derivative(2)*S₃₃)(τₙ), STRETCH_ACEL)

            # Assign these fields to their associated arrays.
            t[n]   = tₙ
            F[n]   = Fₙ
            F′[n]  = dFₙ
            F′′[n] = d²Fₙ
        end

        # Assign these B-spline arrays to the fields of object SplineF.
        new(location, nodes, t, F, F′, F′′)
    end
end # SplineF

end  # FijLung
