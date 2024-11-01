#=
Created on Fri 18 Feb 2022
Updated on Sat 18 May 2024
=#
"""
# FijLung

This module constructs deformation gradient histories at three locations taken from a finite element analysis of a human torso that has been subjected to a ballistic projectile. Location 1 lies just under the pleura membrane, and therefore resides close to the rib-cage and the location of impact. Location 2 lies deep within the lung. While location 3 lies adjacent to a bronchiole tube. These raw data sets are fit with B-splines to get continuous functions for the deformation gradient, thereby allowing its first and second derivatives in time to be approximated. This module exports:

    t_loc1,     # Raw data: times at location near visceral pleura.
    t_loc2,     # Raw data: times at location deep within the lung.
    t_loc3,     # Raw data: times at location next to bronchiole tube.
    F_loc1,     # Raw data: deformation gradient at location 1.
    F_loc2,     # Raw data: deformation gradient at location 2.
    F_loc3,     # Raw data: deformation gradient at location 3.
    SplineF,    # Type used to hold these spline data.
    
along with two external constructors:

    splineAtEndPoints   # Spline nodes are at end points of time intervals.
    splineAtMidPoints   # Spline nodes are at mid points of time intervals.
"""
module FijLung

using
    CubicSplines,
    CSV,
    DataFrames,
    Downloads,
    PhysicalFields

export
    t_loc1,             # Raw data: times at location near visceral pleura.
    t_loc2,             # Raw data: times at location deep within the lung.
    t_loc3,             # Raw data: times at location next to a bronchiole tube.
    F_loc1,             # Raw data: deformation gradient at location 1.
    F_loc2,             # Raw data: deformation gradient at location 2.
    F_loc3,             # Raw data: deformation gradient at location 3.
    SplineF,            # Type used to hold the b-spline data.
    splineAtEndPoints,  # Spline nodes are at the end points of time intervals.
    splineAtMidPoints   # Spline nodes are at the mid points of time intervals.

const DIMENSIONLESS = CGS_DIMENSIONLESS
const SECOND        = CGS_SECOND
const STRETCH       = CGS_STRETCH
const STRETCH_RATE  = CGS_STRETCH_RATE
const STRETCH_ACEL  = CGS_ACCELERATION - CGS_LENGTH
#=
-------------------------------------------------------------------------------
=#
"""
```julia
function t_loc1()::ArrayOfPhysicalScalars
```
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 1 (just under the visceral pleura) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc1()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc1/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        error("Files Fij.csv have an unexpected format.")
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
```julia
function F_loc1()::ArrayOfPhysicalTensors
```
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
        error("Files Fij.csv have an unexpected format.")
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
```julia
function t_loc2()::ArrayOfPhysicalScalars
```
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 2 (roughly at the center of the right lung) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc2()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc2/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        error("Files Fij.csv have an unexpected format.")
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
```julia
function F_loc2()::ArrayOfPhysicalTensors
```
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
        error("Files Fij.csv have an unexpected format.")
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
```julia
function t_loc3()::ArrayOfPhysicalScalars
```
This function returns an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, where each element within the array holds a time whereat a deformation gradient is evaluated at location 3 (adjacent to a bronchiole tube) from an FEA of a human torso subjected to a high-energy ballistic impact.
"""
function t_loc3()::ArrayOfPhysicalScalars
    # Open an existing file and read its data as a DataFrames.DataFrame object.
    dataURL = "https://raw.githubusercontent.com/AlanFreed/FijLung.jl/main/data/Loc3/F11.csv"
    dataF11 = DataFrame(CSV.File(Downloads.download(dataURL)))
    # Create an array holding the time for each data entry.
    (N, M) = size(dataF11)
    if M ≠ 2
        error("Files Fij.csv have an unexpected format.")
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
```julia
function F_loc3()::ArrayOfPhysicalTensors
```
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
        error("Files Fij.csv have an unexpected format.")
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
# SplineF

```julia
struct SplineF
    loc::Integer                  # loc       Lung location, i.e., 1, 2 or 3.
    N::  Integer                  # N         Number of nodes spanning time.
    t::ArrayOfPhysicalScalars     # t         The time at these nodes.
    F::ArrayOfPhysicalTensors     # F         Their deformation gradients.
    F′::ArrayOfPhysicalTensors    # dF/dt     Their first derivatives in time.
    F″::ArrayOfPhysicalTensors    # d²F/dt²   Their second derivatives in time.
end
```
This data structure smooths the raw data supplied by F_loc1, F_loc2 or F_loc3, as specified by parameter `loc`, by applying B-splines to the raw data; specifically, F is fit with a cubic spline, and therefore, its first derivative `dF/dt` is described by a quadratic spline, and its second derivative `d²F/dt²` is described by a linear spline. Fields `loc` and `N` are Integers. Field `t` is an instance of type `PhysicalFields.ArrayOfPhysicalScalars`, while fields `F,` `F′` and `F″` are instances of type `PhysicalFields.ArrayOfPhysicalTensors`. The data arrays are of length N+1 with entry [1] supplying its initial condition, while the remaining N entries are evenly spaced through time.

## Internal Constructor

This constructor is called by the two external constructors, which are the constructors intended for use.

```julia
splineF = SplineF(loc::Int, N::Int, t::ArrayOfPhysicalScalars, F::ArrayOfPhysicalTensors, F′::ArrayOfPhysicalTensors, F″::ArrayOfPhysicalTensors)
```

## External Constructors

These two constructors create instances of type `SplineF`, which are returned as an object `splineF`. Each object contains interpolated values for the deformation gradient F, its first derivative dF/dt, and its second derivative d²F/dt² at N uniformly spaced `nodes` over time. The arrays are of length N+1, with the first entry in each array associating with its corresponding initial condition, while the remaining entries associate with either mid-point or end-point values for N uniformly spaced time intervals.

Each data structure associates with one of three `location`s in a lung obtained from an FE analysis of a human torso subjected to a ballistic impact. Location 1 associates with a location that is just under the visceral pleura. Location 2 associates with a location that is deep within the right lung. While location 3 associates with a location that lies next to a bronchiole tube.

```julia
splineF = splineAtEndPoints(location::Int, nodes::Int)
```

```julia
splineF = splineAtMidPoints(location::Int, nodes::Int)
```
"""
struct SplineF
    loc::Int                      # loc       Lung location, i.e., 1, 2 or 3.
    N::  Int                      # N         Number of nodes spanning time.
    t::ArrayOfPhysicalScalars     # t         The time at these nodes.
    F::ArrayOfPhysicalTensors     # F         Their deformation gradients.
    F′::ArrayOfPhysicalTensors    # dF/dt     Their first derivatives in time.
    F″::ArrayOfPhysicalTensors    # d²F/dt²   Their second derivatives in time.

    # Internal Constructor

    function SplineF(loc::Int, N::Int, t::ArrayOfPhysicalScalars, F::ArrayOfPhysicalTensors, F′::ArrayOfPhysicalTensors, F″::ArrayOfPhysicalTensors)

        new(loc, N, t, F, F′, F″)::SplineF
    end
end # SplineF

# External Constructors

function splineAtEndPoints(location::Int, nodes::Int)::SplineF
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
        error("An inadmissible location was submitted. It must be 1, 2 or 3.")
    end
    knots = time.array.len

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
    S₁₁ = CubicSpline(rawTime, rawF₁₁)
    S₁₂ = CubicSpline(rawTime, rawF₁₂)
    S₁₃ = CubicSpline(rawTime, rawF₁₃)
    S₂₁ = CubicSpline(rawTime, rawF₂₁)
    S₂₂ = CubicSpline(rawTime, rawF₂₂)
    S₂₃ = CubicSpline(rawTime, rawF₂₃)
    S₃₁ = CubicSpline(rawTime, rawF₃₁)
    S₃₂ = CubicSpline(rawTime, rawF₃₂)
    S₃₃ = CubicSpline(rawTime, rawF₃₃)

    # Create the data arrays that span time which hold these data.
    t  = ArrayOfPhysicalScalars(nodes+1, SECOND)
    F  = ArrayOfPhysicalTensors(nodes+1, 3, 3, STRETCH)
    F′ = ArrayOfPhysicalTensors(nodes+1, 3, 3, STRETCH_RATE)
    F″ = ArrayOfPhysicalTensors(nodes+1, 3, 3, STRETCH_ACEL)

    # Establish the time increment between adjacent nodes.
    dt = (time[knots] - time[1]) / nodes

    # Populate the time array with uniformly spaced end-point times.
    t[1] = time[1]           # the initial time
    for n in 1:nodes
        t[n+1] = t[1] + n*dt # the end-point times
    end

    # Populate the kinematic arrays with end-point data.
    for n in 1:nodes+1
        # Create the scalar and tensor fields for the nᵗʰ node.
        Fₙ   = PhysicalTensor(3, 3, STRETCH)
        dFₙ  = PhysicalTensor(3, 3, STRETCH_RATE)
        d²Fₙ = PhysicalTensor(3, 3, STRETCH_ACEL)

        # Time tₙ associated with the nᵗʰ end-point state.
        tₙ = get(t[n])

        # Update the tensor components at the nᵗʰ end-point state.

        # For the deformation gradient.
        Fₙ[1,1] = PhysicalScalar(CubicSplines.Y(S₁₁, tₙ), STRETCH)
        Fₙ[1,2] = PhysicalScalar(CubicSplines.Y(S₁₂, tₙ), STRETCH)
        Fₙ[1,3] = PhysicalScalar(CubicSplines.Y(S₁₃, tₙ), STRETCH)
        Fₙ[2,1] = PhysicalScalar(CubicSplines.Y(S₂₁, tₙ), STRETCH)
        Fₙ[2,2] = PhysicalScalar(CubicSplines.Y(S₂₂, tₙ), STRETCH)
        Fₙ[2,3] = PhysicalScalar(CubicSplines.Y(S₂₃, tₙ), STRETCH)
        Fₙ[3,1] = PhysicalScalar(CubicSplines.Y(S₃₁, tₙ), STRETCH)
        Fₙ[3,2] = PhysicalScalar(CubicSplines.Y(S₃₂, tₙ), STRETCH)
        Fₙ[3,3] = PhysicalScalar(CubicSplines.Y(S₃₃, tₙ), STRETCH)

        # For the first derivative of the deformation gradient.
        dFₙ[1,1] = PhysicalScalar(CubicSplines.Y′(S₁₁, tₙ), STRETCH_RATE)
        dFₙ[1,2] = PhysicalScalar(CubicSplines.Y′(S₁₂, tₙ), STRETCH_RATE)
        dFₙ[1,3] = PhysicalScalar(CubicSplines.Y′(S₁₃, tₙ), STRETCH_RATE)
        dFₙ[2,1] = PhysicalScalar(CubicSplines.Y′(S₂₁, tₙ), STRETCH_RATE)
        dFₙ[2,2] = PhysicalScalar(CubicSplines.Y′(S₂₂, tₙ), STRETCH_RATE)
        dFₙ[2,3] = PhysicalScalar(CubicSplines.Y′(S₂₃, tₙ), STRETCH_RATE)
        dFₙ[3,1] = PhysicalScalar(CubicSplines.Y′(S₃₁, tₙ), STRETCH_RATE)
        dFₙ[3,2] = PhysicalScalar(CubicSplines.Y′(S₃₂, tₙ), STRETCH_RATE)
        dFₙ[3,3] = PhysicalScalar(CubicSplines.Y′(S₃₃, tₙ), STRETCH_RATE)

        # For the second derivative of the deformation gradient.
        d²Fₙ[1,1] = PhysicalScalar(CubicSplines.Y″(S₁₁, tₙ), STRETCH_ACEL)
        d²Fₙ[1,2] = PhysicalScalar(CubicSplines.Y″(S₁₂, tₙ), STRETCH_ACEL)
        d²Fₙ[1,3] = PhysicalScalar(CubicSplines.Y″(S₁₃, tₙ), STRETCH_ACEL)
        d²Fₙ[2,1] = PhysicalScalar(CubicSplines.Y″(S₂₁, tₙ), STRETCH_ACEL)
        d²Fₙ[2,2] = PhysicalScalar(CubicSplines.Y″(S₂₂, tₙ), STRETCH_ACEL)
        d²Fₙ[2,3] = PhysicalScalar(CubicSplines.Y″(S₂₃, tₙ), STRETCH_ACEL)
        d²Fₙ[3,1] = PhysicalScalar(CubicSplines.Y″(S₃₁, tₙ), STRETCH_ACEL)
        d²Fₙ[3,2] = PhysicalScalar(CubicSplines.Y″(S₃₂, tₙ), STRETCH_ACEL)
        d²Fₙ[3,3] = PhysicalScalar(CubicSplines.Y″(S₃₃, tₙ), STRETCH_ACEL)

        # Assign these fields to their associated arrays.
        F[n]  = Fₙ
        F′[n] = dFₙ
        F″[n] = d²Fₙ
    end

    # Assign these B-spline arrays to the fields of object SplineF.
    splineF = SplineF(location, nodes, t, F, F′, F″)

    return splineF
end # splineAtEndPoints

function splineAtMidPoints(location::Int, nodes::Int)::SplineF
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
        error("An inadmissible location was submitted. It must be 1, 2 or 3.")
    end
    knots = time.array.len

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
    S₁₁ = CubicSpline(rawTime, rawF₁₁)
    S₁₂ = CubicSpline(rawTime, rawF₁₂)
    S₁₃ = CubicSpline(rawTime, rawF₁₃)
    S₂₁ = CubicSpline(rawTime, rawF₂₁)
    S₂₂ = CubicSpline(rawTime, rawF₂₂)
    S₂₃ = CubicSpline(rawTime, rawF₂₃)
    S₃₁ = CubicSpline(rawTime, rawF₃₁)
    S₃₂ = CubicSpline(rawTime, rawF₃₂)
    S₃₃ = CubicSpline(rawTime, rawF₃₃)

    # Create the data arrays that span time which hold these data.
    t  = ArrayOfPhysicalScalars(nodes+1, SECOND)
    F  = ArrayOfPhysicalTensors(nodes+1, 3, 3, STRETCH)
    F′ = ArrayOfPhysicalTensors(nodes+1, 3, 3, STRETCH_RATE)
    F″ = ArrayOfPhysicalTensors(nodes+1, 3, 3, STRETCH_ACEL)

    # Establish the time increment between adjacent nodes.
    dt = (time[knots] - time[1]) / nodes

    # Populate the time array with uniformly spaced mid-point times.
    t[1] = time[1]                  # the initial time
    t[2] = t[1] + 0.5dt             # the first mid-point time
    for n in 2:nodes
        t[n+1] = t[2] + (n-1)*dt    # the remaining mid-point times
    end

    # Populate the kinematic arrays with mid-point data.
    for n in 1:nodes+1
        # Create the scalar and tensor fields for the nᵗʰ node.
        Fₙ   = PhysicalTensor(3, 3, STRETCH)
        dFₙ  = PhysicalTensor(3, 3, STRETCH_RATE)
        d²Fₙ = PhysicalTensor(3, 3, STRETCH_ACEL)

        # Time tₙ associated with the nᵗʰ mid-point state.
        tₙ = get(t[n])

        # Update the tensor components at the nᵗʰ mid-point state.

        # For the deformation gradient.
        Fₙ[1,1] = PhysicalScalar(CubicSplines.Y(S₁₁, tₙ), STRETCH)
        Fₙ[1,2] = PhysicalScalar(CubicSplines.Y(S₁₂, tₙ), STRETCH)
        Fₙ[1,3] = PhysicalScalar(CubicSplines.Y(S₁₃, tₙ), STRETCH)
        Fₙ[2,1] = PhysicalScalar(CubicSplines.Y(S₂₁, tₙ), STRETCH)
        Fₙ[2,2] = PhysicalScalar(CubicSplines.Y(S₂₂, tₙ), STRETCH)
        Fₙ[2,3] = PhysicalScalar(CubicSplines.Y(S₂₃, tₙ), STRETCH)
        Fₙ[3,1] = PhysicalScalar(CubicSplines.Y(S₃₁, tₙ), STRETCH)
        Fₙ[3,2] = PhysicalScalar(CubicSplines.Y(S₃₂, tₙ), STRETCH)
        Fₙ[3,3] = PhysicalScalar(CubicSplines.Y(S₃₃, tₙ), STRETCH)

        # For the first derivative of the deformation gradient.
        dFₙ[1,1] = PhysicalScalar(CubicSplines.Y′(S₁₁, tₙ), STRETCH_RATE)
        dFₙ[1,2] = PhysicalScalar(CubicSplines.Y′(S₁₂, tₙ), STRETCH_RATE)
        dFₙ[1,3] = PhysicalScalar(CubicSplines.Y′(S₁₃, tₙ), STRETCH_RATE)
        dFₙ[2,1] = PhysicalScalar(CubicSplines.Y′(S₂₁, tₙ), STRETCH_RATE)
        dFₙ[2,2] = PhysicalScalar(CubicSplines.Y′(S₂₂, tₙ), STRETCH_RATE)
        dFₙ[2,3] = PhysicalScalar(CubicSplines.Y′(S₂₃, tₙ), STRETCH_RATE)
        dFₙ[3,1] = PhysicalScalar(CubicSplines.Y′(S₃₁, tₙ), STRETCH_RATE)
        dFₙ[3,2] = PhysicalScalar(CubicSplines.Y′(S₃₂, tₙ), STRETCH_RATE)
        dFₙ[3,3] = PhysicalScalar(CubicSplines.Y′(S₃₃, tₙ), STRETCH_RATE)

        # For the second derivative of the deformation gradient.
        d²Fₙ[1,1] = PhysicalScalar(CubicSplines.Y″(S₁₁, tₙ), STRETCH_ACEL)
        d²Fₙ[1,2] = PhysicalScalar(CubicSplines.Y″(S₁₂, tₙ), STRETCH_ACEL)
        d²Fₙ[1,3] = PhysicalScalar(CubicSplines.Y″(S₁₃, tₙ), STRETCH_ACEL)
        d²Fₙ[2,1] = PhysicalScalar(CubicSplines.Y″(S₂₁, tₙ), STRETCH_ACEL)
        d²Fₙ[2,2] = PhysicalScalar(CubicSplines.Y″(S₂₂, tₙ), STRETCH_ACEL)
        d²Fₙ[2,3] = PhysicalScalar(CubicSplines.Y″(S₂₃, tₙ), STRETCH_ACEL)
        d²Fₙ[3,1] = PhysicalScalar(CubicSplines.Y″(S₃₁, tₙ), STRETCH_ACEL)
        d²Fₙ[3,2] = PhysicalScalar(CubicSplines.Y″(S₃₂, tₙ), STRETCH_ACEL)
        d²Fₙ[3,3] = PhysicalScalar(CubicSplines.Y″(S₃₃, tₙ), STRETCH_ACEL)

        # Assign these fields to their associated arrays.
        F[n]  = Fₙ
        F′[n] = dFₙ
        F″[n] = d²Fₙ
    end

    # Assign these B-spline arrays to the fields of object SplineF.
    splineF = SplineF(location, nodes, t, F, F′, F″)

    return splineF
end # splineAtMidPoints

end  # FijLung
