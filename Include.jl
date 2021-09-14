# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# includes -
include("Types.jl")
include("Kinetics.jl")
include("Control.jl")
include("Inputs.jl")
include("Data.jl")
include("SolveBalances.jl")
include("Balances.jl")
include("Utility.jl")
include("Error.jl")
include("POETs.jl/src/POETs.jl")
include("POETs.jl/src/Main.jl")

# system packages - these are required to be installed -
# check - if they are installed, all is good with the world.
# if not then install them
using LinearAlgebra # pre-installed w/Julia
using Statistics    # pre-installed w/Julia
using Pkg           # pre-installed w/Julia
installed_package_set = keys(Pkg.installed())

# Do we have DifferentialEquations?
if (in("DifferentialEquations",installed_package_set) == false)
    Pkg.add("DifferentialEquations")
end

# Do we have DelimitedFiles?
if (in("DelimitedFiles",installed_package_set) == false)
    Pkg.add("DelimitedFiles")
end

# Do we have JSON?
if (in("JSON",installed_package_set) == false)
    Pkg.add("JSON")
end

# Do we have ProgressMeter?
if (in("ProgressMeter",installed_package_set) == false)
    Pkg.add("ProgressMeter")
end

if (in("Optim",installed_package_set) == false)
    Pkg.add("Optim")
end

if (in("DataFrames",installed_package_set) == false)
    Pkg.add("DataFrames")
end

if (in("CSV",installed_package_set) == false)
    Pkg.add("CSV")
end

if (in("Interpolations",installed_package_set) == false)
    Pkg.add("Interpolations")
end

if (in("PyPlot",installed_package_set) == false)
    Pkg.add("PyPlot")
end

if (in("NumericalIntegration",installed_package_set) == false)
    Pkg.add("NumericalIntegration")
end

if (in("DiffEqSensitivity",installed_package_set) == false)
    Pkg.add("DiffEqSensitivity")
end

if (in("GlobalSensitivity",installed_package_set) == false)
    Pkg.add("GlobalSensitivity")
end

# Ok to use now ...
using DifferentialEquations
using DelimitedFiles
using JSON
using ProgressMeter
using Optim
using DataFrames
using CSV
using Interpolations
using PyPlot
using NumericalIntegration
using DiffEqSensitivity
using GlobalSensitivity
