#=
Created on Fri 18 Feb 2022
Updated on Mon 12 Jun 2023
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
    SplineF,        # Type used to hold these splined data.\n
    newSplineF      # A constructor for type SplineF.\n
*Note:* Functions t_loc1, t_loc2, t_loc3, F_loc1, F_loc2 and F_loc3 have hardwired paths for reading in the data files. You will need to edit these paths so that they point to the appropriate raw CSV files that you have placed on your machine.
"""
module FijLung

using
    BSplineKit,
    CSV,
    DataFrames,
    PhysicalFields,
    StaticArrays

export
    t_loc1,         # Raw data: times at location near visceral pleura.
    t_loc2,         # Raw data: times at location deep within the lung.
    t_loc3,         # Raw data: times at location next to a bronchiole tube.
    F_loc1,         # Raw data: deformation gradient at location 1.
    F_loc2,         # Raw data: deformation gradient at location 2.
    F_loc3,         # Raw data: deformation gradient at location 3.
    SplineF,        # Type used to hold the b-spline data.
    newSplineF      # A constructor for instances of type SplineF.

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
    dataF11 = CSV.read("data/Loc1/F11-loc1.csv", DataFrame; header=1)
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij-loc1.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array of times whose initial time is t₁.
    t₁ = newPhysicalScalar(SECOND)
    set!(t₁, dataF11[1,1])
    tLoc1 = newArrayOfPhysicalScalars(N, t₁)
    # Populate this array with times that associate with its N nodes.
    tₙ = newPhysicalScalar(SECOND)
    for n in 2:N
        set!(tₙ, dataF11[n,1])
        tLoc1[n] = tₙ
    end
    return tLoc1
end

"""
Function:\n
    arrayOfF = F_loc1()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalTensors`, where each element within the array holds a deformation gradient tensor evaluated at location 1 (just under the visceral pleura) from an FEA of a human torso subjected to a high-energy ballistic impact. The time associated with each entry is supplied by the function `t_loc1()`.
"""
function F_loc1()::ArrayOfPhysicalTensors
    # Open existing files and read in their data as DataFrames.DataFrame objects.
    dataF11 = CSV.read("data/Loc1/F11-loc1.csv", DataFrame; header=1)
    dataF12 = CSV.read("data/Loc1/F12-loc1.csv", DataFrame; header=1)
    dataF13 = CSV.read("data/Loc1/F13-loc1.csv", DataFrame; header=1)
    dataF21 = CSV.read("data/Loc1/F21-loc1.csv", DataFrame; header=1)
    dataF22 = CSV.read("data/Loc1/F22-loc1.csv", DataFrame; header=1)
    dataF23 = CSV.read("data/Loc1/F23-loc1.csv", DataFrame; header=1)
    dataF31 = CSV.read("data/Loc1/F31-loc1.csv", DataFrame; header=1)
    dataF32 = CSV.read("data/Loc1/F32-loc1.csv", DataFrame; header=1)
    dataF33 = CSV.read("data/Loc1/F33-loc1.csv", DataFrame; header=1)
    # Size the arrays.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij-loc1.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create 9 scalar components for a deformation gradient.
    F₁₁ = newPhysicalScalar(STRETCH)
    F₁₂ = newPhysicalScalar(STRETCH)
    F₁₃ = newPhysicalScalar(STRETCH)
    F₂₁ = newPhysicalScalar(STRETCH)
    F₂₂ = newPhysicalScalar(STRETCH)
    F₂₃ = newPhysicalScalar(STRETCH)
    F₃₁ = newPhysicalScalar(STRETCH)
    F₃₂ = newPhysicalScalar(STRETCH)
    F₃₃ = newPhysicalScalar(STRETCH)
    # Assign to these scalars their values at node 1.
    set!(F₁₁, dataF11[1,2])
    set!(F₁₂, dataF12[1,2])
    set!(F₁₃, dataF13[1,2])
    set!(F₂₁, dataF21[1,2])
    set!(F₂₂, dataF22[1,2])
    set!(F₂₃, dataF23[1,2])
    set!(F₃₁, dataF31[1,2])
    set!(F₃₂, dataF32[1,2])
    set!(F₃₃, dataF33[1,2])
    # Create and populate the deformation gradient at node 1.
    F₁ = newPhysicalTensor(3, 3, STRETCH)
    F₁[1,1] = F₁₁
    F₁[1,2] = F₁₂
    F₁[1,3] = F₁₃
    F₁[2,1] = F₂₁
    F₁[2,2] = F₂₂
    F₁[2,3] = F₂₃
    F₁[3,1] = F₃₁
    F₁[3,2] = F₃₂
    F₁[3,3] = F₃₃
    # Create a new static array with entries of type PhysicalTensor.
    FLoc1 = newArrayOfPhysicalTensors(N, F₁)
    Fₙ = newPhysicalTensor(3, 3, STRETCH)
    for n in 2:N
        # Assign components to a deformation gradient evaluated at node n.
        set!(F₁₁, dataF11[n,2])
        set!(F₁₂, dataF12[n,2])
        set!(F₁₃, dataF13[n,2])
        set!(F₂₁, dataF21[n,2])
        set!(F₂₂, dataF22[n,2])
        set!(F₂₃, dataF23[n,2])
        set!(F₃₁, dataF31[n,2])
        set!(F₃₂, dataF32[n,2])
        set!(F₃₃, dataF33[n,2])
        # Populate the deformation gradient with these components.
        Fₙ[1,1] = F₁₁
        Fₙ[1,2] = F₁₂
        Fₙ[1,3] = F₁₃
        Fₙ[2,1] = F₂₁
        Fₙ[2,2] = F₂₂
        Fₙ[2,3] = F₂₃
        Fₙ[3,1] = F₃₁
        Fₙ[3,2] = F₃₂
        Fₙ[3,3] = F₃₃
        # Assign this deformation gradient to the array.
        FLoc1[n] = Fₙ
    end
    return FLoc1
end

"""
Function:\n
    arrayOfTimes = t_loc2()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 2 (roughly in the center of the lung) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc2()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataF11 = CSV.read("data/Loc2/F11-loc2.csv", DataFrame; header=1)
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij-loc2.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array of times whose initial time is t₁.
    t₁ = newPhysicalScalar(SECOND)
    set!(t₁, dataF11[1,1])
    tLoc2 = newArrayOfPhysicalScalars(N, t₁)
    # Populate this array with times that associate with its N nodes.
    tₙ = newPhysicalScalar(SECOND)
    for n in 2:N
        set!(tₙ, dataF11[n,1])
        tLoc2[n] = tₙ
    end
    return tLoc2
end

"""
Function:\n
    arrayOfF = F_loc2()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalTensors`, where each element within the array holds a deformation gradient tensor evaluated at location 2 (roughly in the center of the lung) from an FEA of a human torso subjected to a high-energy ballistic impact. The time associated with each entry is supplied by the function `t_loc2()`.
"""
function F_loc2()::ArrayOfPhysicalTensors
    # Open existing files and read in their data as DataFrames.DataFrame objects.
    dataF11 = CSV.read("data/Loc2/F11-loc2.csv", DataFrame; header=1)
    dataF12 = CSV.read("data/Loc2/F12-loc2.csv", DataFrame; header=1)
    dataF13 = CSV.read("data/Loc2/F13-loc2.csv", DataFrame; header=1)
    dataF21 = CSV.read("data/Loc2/F21-loc2.csv", DataFrame; header=1)
    dataF22 = CSV.read("data/Loc2/F22-loc2.csv", DataFrame; header=1)
    dataF23 = CSV.read("data/Loc2/F23-loc2.csv", DataFrame; header=1)
    dataF31 = CSV.read("data/Loc2/F31-loc2.csv", DataFrame; header=1)
    dataF32 = CSV.read("data/Loc2/F32-loc2.csv", DataFrame; header=1)
    dataF33 = CSV.read("data/Loc2/F33-loc2.csv", DataFrame; header=1)
    # Size the arrays.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij-loc2.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create 9 scalar components for a deformation gradient.
    F₁₁ = newPhysicalScalar(STRETCH)
    F₁₂ = newPhysicalScalar(STRETCH)
    F₁₃ = newPhysicalScalar(STRETCH)
    F₂₁ = newPhysicalScalar(STRETCH)
    F₂₂ = newPhysicalScalar(STRETCH)
    F₂₃ = newPhysicalScalar(STRETCH)
    F₃₁ = newPhysicalScalar(STRETCH)
    F₃₂ = newPhysicalScalar(STRETCH)
    F₃₃ = newPhysicalScalar(STRETCH)
    # Assign to these scalars their values at node 1.
    set!(F₁₁, dataF11[1,2])
    set!(F₁₂, dataF12[1,2])
    set!(F₁₃, dataF13[1,2])
    set!(F₂₁, dataF21[1,2])
    set!(F₂₂, dataF22[1,2])
    set!(F₂₃, dataF23[1,2])
    set!(F₃₁, dataF31[1,2])
    set!(F₃₂, dataF32[1,2])
    set!(F₃₃, dataF33[1,2])
    # Create and populate the deformation gradient at node 1.
    F₁ = newPhysicalTensor(3, 3, STRETCH)
    F₁[1,1] = F₁₁
    F₁[1,2] = F₁₂
    F₁[1,3] = F₁₃
    F₁[2,1] = F₂₁
    F₁[2,2] = F₂₂
    F₁[2,3] = F₂₃
    F₁[3,1] = F₃₁
    F₁[3,2] = F₃₂
    F₁[3,3] = F₃₃
    # Create a new static array with entries of type PhysicalTensor.
    FLoc2 = newArrayOfPhysicalTensors(N, F₁)
    Fₙ = newPhysicalTensor(3, 3, STRETCH)
    for n in 2:N
        # Assign components to a deformation gradient evaluated at node n.
        set!(F₁₁, dataF11[n,2])
        set!(F₁₂, dataF12[n,2])
        set!(F₁₃, dataF13[n,2])
        set!(F₂₁, dataF21[n,2])
        set!(F₂₂, dataF22[n,2])
        set!(F₂₃, dataF23[n,2])
        set!(F₃₁, dataF31[n,2])
        set!(F₃₂, dataF32[n,2])
        set!(F₃₃, dataF33[n,2])
        # Populate the deformation gradient with these components.
        Fₙ[1,1] = F₁₁
        Fₙ[1,2] = F₁₂
        Fₙ[1,3] = F₁₃
        Fₙ[2,1] = F₂₁
        Fₙ[2,2] = F₂₂
        Fₙ[2,3] = F₂₃
        Fₙ[3,1] = F₃₁
        Fₙ[3,2] = F₃₂
        Fₙ[3,3] = F₃₃
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
    dataF11 = CSV.read("data/Loc3/F11-loc3.csv", DataFrame; header=1)
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij-loc3.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create a new array of times whose initial time is t₁.
    t₁ = newPhysicalScalar(SECOND)
    set!(t₁, dataF11[1,1])
    tLoc3 = newArrayOfPhysicalScalars(N, t₁)
    # Populate this array with times that associate with its N nodes.
    tₙ = newPhysicalScalar(SECOND)
    for n in 2:N
        set!(tₙ, dataF11[n,1])
        tLoc3[n] = tₙ
    end
    return tLoc3
end

"""
Function:\n
    arrayOfF = F_loc3()\n
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalTensors`, where each element within the array holds a deformation gradient tensor evaluated at location 3 (adjacent to a bronchiole tube) from an FEA of a human torso subjected to a high-energy ballistic impact. The time associated with each entry is supplied by the function `t_loc3()`.
"""
function F_loc3()::ArrayOfPhysicalTensors
    # Open existing files and read in their data as DataFrames.DataFrame objects.
    dataF11 = CSV.read("data/Loc3/F11-loc3.csv", DataFrame; header=1)
    dataF12 = CSV.read("data/Loc3/F12-loc3.csv", DataFrame; header=1)
    dataF13 = CSV.read("data/Loc3/F13-loc3.csv", DataFrame; header=1)
    dataF21 = CSV.read("data/Loc3/F21-loc3.csv", DataFrame; header=1)
    dataF22 = CSV.read("data/Loc3/F22-loc3.csv", DataFrame; header=1)
    dataF23 = CSV.read("data/Loc3/F23-loc3.csv", DataFrame; header=1)
    dataF31 = CSV.read("data/Loc3/F31-loc3.csv", DataFrame; header=1)
    dataF32 = CSV.read("data/Loc3/F32-loc3.csv", DataFrame; header=1)
    dataF33 = CSV.read("data/Loc3/F33-loc3.csv", DataFrame; header=1)
    # Size the arrays.
    (N, M) = size(dataF11)
    if M ≠ 2
        msg = "Files Fij-loc3.csv have an unexpected format."
        throw(ErrorException(msg))
    end
    # Create 9 scalar components for a deformation gradient.
    F₁₁ = newPhysicalScalar(STRETCH)
    F₁₂ = newPhysicalScalar(STRETCH)
    F₁₃ = newPhysicalScalar(STRETCH)
    F₂₁ = newPhysicalScalar(STRETCH)
    F₂₂ = newPhysicalScalar(STRETCH)
    F₂₃ = newPhysicalScalar(STRETCH)
    F₃₁ = newPhysicalScalar(STRETCH)
    F₃₂ = newPhysicalScalar(STRETCH)
    F₃₃ = newPhysicalScalar(STRETCH)
    # Assign to these scalars their values at node 1.
    set!(F₁₁, dataF11[1,2])
    set!(F₁₂, dataF12[1,2])
    set!(F₁₃, dataF13[1,2])
    set!(F₂₁, dataF21[1,2])
    set!(F₂₂, dataF22[1,2])
    set!(F₂₃, dataF23[1,2])
    set!(F₃₁, dataF31[1,2])
    set!(F₃₂, dataF32[1,2])
    set!(F₃₃, dataF33[1,2])
    # Create and populate the deformation gradient at node 1.
    F₁ = newPhysicalTensor(3, 3, STRETCH)
    F₁[1,1] = F₁₁
    F₁[1,2] = F₁₂
    F₁[1,3] = F₁₃
    F₁[2,1] = F₂₁
    F₁[2,2] = F₂₂
    F₁[2,3] = F₂₃
    F₁[3,1] = F₃₁
    F₁[3,2] = F₃₂
    F₁[3,3] = F₃₃
    # Create a new static array with entries of type PhysicalTensor.
    FLoc3 = newArrayOfPhysicalTensors(N, F₁)
    Fₙ = newPhysicalTensor(3, 3, STRETCH)
    for n in 2:N
        # Assign components to a deformation gradient evaluated at node n.
        set!(F₁₁, dataF11[n,2])
        set!(F₁₂, dataF12[n,2])
        set!(F₁₃, dataF13[n,2])
        set!(F₂₁, dataF21[n,2])
        set!(F₂₂, dataF22[n,2])
        set!(F₂₃, dataF23[n,2])
        set!(F₃₁, dataF31[n,2])
        set!(F₃₂, dataF32[n,2])
        set!(F₃₃, dataF33[n,2])
        # Populate the deformation gradient with these components.
        Fₙ[1,1] = F₁₁
        Fₙ[1,2] = F₁₂
        Fₙ[1,3] = F₁₃
        Fₙ[2,1] = F₂₁
        Fₙ[2,2] = F₂₂
        Fₙ[2,3] = F₂₃
        Fₙ[3,1] = F₃₁
        Fₙ[3,2] = F₃₂
        Fₙ[3,3] = F₃₃
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
end

"""
Constructor:\n
    splineF = newSplineF(location, nodes)\n
This constructor creates instances of type SplineF, which are returned as objects `splineF`. Each object contains interpolated values for the deformation gradient F, its first derivative dF/dt, and its second derivative d²F/dt² at `N` uniformly spaced instances in time `t`. These data associate with one of three locations in a lung from an FE analysis of a human torso subjected to a ballistic impact.
"""
function newSplineF(location::Integer, nodes::Integer)::SplineF
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
    knots = Int64(time.e)
    # Create the various static arrays that will hold these raw data.
    rawTime = @MVector zeros(Float64, knots)
    rawF₁₁ = @MVector zeros(Float64, knots)
    rawF₁₂ = @MVector zeros(Float64, knots)
    rawF₁₃ = @MVector zeros(Float64, knots)
    rawF₂₁ = @MVector zeros(Float64, knots)
    rawF₂₂ = @MVector zeros(Float64, knots)
    rawF₂₃ = @MVector zeros(Float64, knots)
    rawF₃₁ = @MVector zeros(Float64, knots)
    rawF₃₂ = @MVector zeros(Float64, knots)
    rawF₃₃ = @MVector zeros(Float64, knots)
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
    # Create various scalar and tensor components needed for the first node:
    # for time
    t₁ = newPhysicalScalar(SECOND)
    # for the deformation gradient
    F₁₁ = newPhysicalScalar(STRETCH)
    F₁₂ = newPhysicalScalar(STRETCH)
    F₁₃ = newPhysicalScalar(STRETCH)
    F₂₁ = newPhysicalScalar(STRETCH)
    F₂₂ = newPhysicalScalar(STRETCH)
    F₂₃ = newPhysicalScalar(STRETCH)
    F₃₁ = newPhysicalScalar(STRETCH)
    F₃₂ = newPhysicalScalar(STRETCH)
    F₃₃ = newPhysicalScalar(STRETCH)
    # for the first derivative of the deformation gradient
    dF₁₁ = newPhysicalScalar(STRETCH_RATE)
    dF₁₂ = newPhysicalScalar(STRETCH_RATE)
    dF₁₃ = newPhysicalScalar(STRETCH_RATE)
    dF₂₁ = newPhysicalScalar(STRETCH_RATE)
    dF₂₂ = newPhysicalScalar(STRETCH_RATE)
    dF₂₃ = newPhysicalScalar(STRETCH_RATE)
    dF₃₁ = newPhysicalScalar(STRETCH_RATE)
    dF₃₂ = newPhysicalScalar(STRETCH_RATE)
    dF₃₃ = newPhysicalScalar(STRETCH_RATE)
    # for the second derivative of the deformation gradient
    d²F₁₁ = newPhysicalScalar(STRETCH_ACEL)
    d²F₁₂ = newPhysicalScalar(STRETCH_ACEL)
    d²F₁₃ = newPhysicalScalar(STRETCH_ACEL)
    d²F₂₁ = newPhysicalScalar(STRETCH_ACEL)
    d²F₂₂ = newPhysicalScalar(STRETCH_ACEL)
    d²F₂₃ = newPhysicalScalar(STRETCH_ACEL)
    d²F₃₁ = newPhysicalScalar(STRETCH_ACEL)
    d²F₃₂ = newPhysicalScalar(STRETCH_ACEL)
    d²F₃₃ = newPhysicalScalar(STRETCH_ACEL)
    # Assign components to the spline arrays to be returned for node 1:
    # for time
    τ₁ = rawTime[1]
    set!(t₁, τ₁)
    # for the deformation gradient
    set!(F₁₁, S₁₁(τ₁))
    set!(F₁₂, S₁₂(τ₁))
    set!(F₁₃, S₁₃(τ₁))
    set!(F₂₁, S₂₁(τ₁))
    set!(F₂₂, S₂₂(τ₁))
    set!(F₂₃, S₂₃(τ₁))
    set!(F₃₁, S₃₁(τ₁))
    set!(F₃₂, S₃₂(τ₁))
    set!(F₃₃, S₃₃(τ₁))
    # for the first derivative of the deformation gradient
    set!(dF₁₁, (Derivative(1) * S₁₁)(τ₁))
    set!(dF₁₂, (Derivative(1) * S₁₂)(τ₁))
    set!(dF₁₃, (Derivative(1) * S₁₃)(τ₁))
    set!(dF₂₁, (Derivative(1) * S₂₁)(τ₁))
    set!(dF₂₂, (Derivative(1) * S₂₂)(τ₁))
    set!(dF₂₃, (Derivative(1) * S₂₃)(τ₁))
    set!(dF₃₁, (Derivative(1) * S₃₁)(τ₁))
    set!(dF₃₂, (Derivative(1) * S₃₂)(τ₁))
    set!(dF₃₃, (Derivative(1) * S₃₃)(τ₁))
    # for the second derivative of the deformation gradient
    set!(d²F₁₁, (Derivative(2) * S₁₁)(τ₁))
    set!(d²F₁₂, (Derivative(2) * S₁₂)(τ₁))
    set!(d²F₁₃, (Derivative(2) * S₁₃)(τ₁))
    set!(d²F₂₁, (Derivative(2) * S₂₁)(τ₁))
    set!(d²F₂₂, (Derivative(2) * S₂₂)(τ₁))
    set!(d²F₂₃, (Derivative(2) * S₂₃)(τ₁))
    set!(d²F₃₁, (Derivative(2) * S₃₁)(τ₁))
    set!(d²F₃₂, (Derivative(2) * S₃₂)(τ₁))
    set!(d²F₃₃, (Derivative(2) * S₃₃)(τ₁))
    # Create and populate their associated tensor fields:
    # for the deformation gradient
    F₁ = newPhysicalTensor(3, 3, STRETCH)
    F₁[1,1] = F₁₁
    F₁[1,2] = F₁₂
    F₁[1,3] = F₁₃
    F₁[2,1] = F₂₁
    F₁[2,2] = F₂₂
    F₁[2,3] = F₂₃
    F₁[3,1] = F₃₁
    F₁[3,2] = F₃₂
    F₁[3,3] = F₃₃
    # for the first derivative of the deformation gradient
    dF₁ = newPhysicalTensor(3, 3, STRETCH_RATE)
    dF₁[1,1] = dF₁₁
    dF₁[1,2] = dF₁₂
    dF₁[1,3] = dF₁₃
    dF₁[2,1] = dF₂₁
    dF₁[2,2] = dF₂₂
    dF₁[2,3] = dF₂₃
    dF₁[3,1] = dF₃₁
    dF₁[3,2] = dF₃₂
    dF₁[3,3] = dF₃₃
    # for the second derivative of the deformation gradient
    d²F₁ = newPhysicalTensor(3, 3, STRETCH_ACEL)
    d²F₁[1,1] = d²F₁₁
    d²F₁[1,2] = d²F₁₂
    d²F₁[1,3] = d²F₁₃
    d²F₁[2,1] = d²F₂₁
    d²F₁[2,2] = d²F₂₂
    d²F₁[2,3] = d²F₂₃
    d²F₁[3,1] = d²F₃₁
    d²F₁[3,2] = d²F₃₂
    d²F₁[3,3] = d²F₃₃
    # Create the arrays that span time which hold these data.
    t = newArrayOfPhysicalScalars(nodes, t₁)
    F = newArrayOfPhysicalTensors(nodes, F₁)
    F′ = newArrayOfPhysicalTensors(nodes, dF₁)
    F′′ = newArrayOfPhysicalTensors(nodes, d²F₁)
    # Create the scalar and tensor fields for the nᵗʰ node.
    tₙ = newPhysicalScalar(SECOND)
    Fₙ = newPhysicalTensor(3, 3, STRETCH)
    dFₙ = newPhysicalTensor(3, 3, STRETCH_RATE)
    d²Fₙ = newPhysicalTensor(3, 3, STRETCH_ACEL)
    # Establish the time increment between nodes.
    dτ = (rawTime[knots] - rawTime[1]) / (nodes - 1)
    τₙ = τ₁
    for n in 2:nodes
        # Time tₙ associated with the nᵗʰ state.
        τₙ += dτ
        set!(tₙ, τₙ)
        # Update the tensor components at the nᵗʰ state:
        # for the deformation gradient
        set!(F₁₁, S₁₁(τₙ))
        set!(F₁₂, S₁₂(τₙ))
        set!(F₁₃, S₁₃(τₙ))
        set!(F₂₁, S₂₁(τₙ))
        set!(F₂₂, S₂₂(τₙ))
        set!(F₂₃, S₂₃(τₙ))
        set!(F₃₁, S₃₁(τₙ))
        set!(F₃₂, S₃₂(τₙ))
        set!(F₃₃, S₃₃(τₙ))
        # for the first derivative of the deformation gradient
        set!(dF₁₁, (Derivative(1) * S₁₁)(τₙ))
        set!(dF₁₂, (Derivative(1) * S₁₂)(τₙ))
        set!(dF₁₃, (Derivative(1) * S₁₃)(τₙ))
        set!(dF₂₁, (Derivative(1) * S₂₁)(τₙ))
        set!(dF₂₂, (Derivative(1) * S₂₂)(τₙ))
        set!(dF₂₃, (Derivative(1) * S₂₃)(τₙ))
        set!(dF₃₁, (Derivative(1) * S₃₁)(τₙ))
        set!(dF₃₂, (Derivative(1) * S₃₂)(τₙ))
        set!(dF₃₃, (Derivative(1) * S₃₃)(τₙ))
        # for the second derivative of the deformation gradient
        set!(d²F₁₁, (Derivative(2) * S₁₁)(τₙ))
        set!(d²F₁₂, (Derivative(2) * S₁₂)(τₙ))
        set!(d²F₁₃, (Derivative(2) * S₁₃)(τₙ))
        set!(d²F₂₁, (Derivative(2) * S₂₁)(τₙ))
        set!(d²F₂₂, (Derivative(2) * S₂₂)(τₙ))
        set!(d²F₂₃, (Derivative(2) * S₂₃)(τₙ))
        set!(d²F₃₁, (Derivative(2) * S₃₁)(τₙ))
        set!(d²F₃₂, (Derivative(2) * S₃₂)(τₙ))
        set!(d²F₃₃, (Derivative(2) * S₃₃)(τₙ))
        # Populate these tensor fields at the nᵗʰ state:
        # for the deformation gradient
        Fₙ[1,1] = F₁₁
        Fₙ[1,2] = F₁₂
        Fₙ[1,3] = F₁₃
        Fₙ[2,1] = F₂₁
        Fₙ[2,2] = F₂₂
        Fₙ[2,3] = F₂₃
        Fₙ[3,1] = F₃₁
        Fₙ[3,2] = F₃₂
        Fₙ[3,3] = F₃₃
        # for the first derivative of the deformation gradient
        dFₙ[1,1] = dF₁₁
        dFₙ[1,2] = dF₁₂
        dFₙ[1,3] = dF₁₃
        dFₙ[2,1] = dF₂₁
        dFₙ[2,2] = dF₂₂
        dFₙ[2,3] = dF₂₃
        dFₙ[3,1] = dF₃₁
        dFₙ[3,2] = dF₃₂
        dFₙ[3,3] = dF₃₃
        # for the second derivative of the deformation gradient
        d²Fₙ[1,1] = d²F₁₁
        d²Fₙ[1,2] = d²F₁₂
        d²Fₙ[1,3] = d²F₁₃
        d²Fₙ[2,1] = d²F₂₁
        d²Fₙ[2,2] = d²F₂₂
        d²Fₙ[2,3] = d²F₂₃
        d²Fₙ[3,1] = d²F₃₁
        d²Fₙ[3,2] = d²F₃₂
        d²Fₙ[3,3] = d²F₃₃
        # Assign these fields to their associated arrays.
        t[n] = tₙ
        F[n] = Fₙ
        F′[n] = dFₙ
        F′′[n] = d²Fₙ
    end
    # Assign these B-spline arrays to the fields of object SplineF.
    splineF = SplineF(location, nodes, t, F, F′, F′′)
    return splineF
end  # newSplineF

end  # FijLung