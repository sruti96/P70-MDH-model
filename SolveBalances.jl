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

# ----------------------------------------------------------------------------------- #
# SolveAdjBalances: Solves adjoint model equations from TSTART to TSTOP given
# parameters in data_dictionary.
#
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
# parameter_index - index of parameter that we want to look at
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# TSIM - Simulation time vector
# X - Simulation state array (NTIME x NSPECIES)
# ----------------------------------------------------------------------------------- #
function SolveAdjBalances(TSTART,TSTOP,Ts,parameter_index,data_dictionary)

  # Get required stuff from DataFile struct -
  time_span = (TSTART,TSTOP)
  initial_condition_array = data_dictionary["initial_condition_array"];

  # which parameter are we looking at?
  data_dictionary["parameter_index"] = parameter_index

  # build problem object -
  problem_object = ODEProblem(AdjBalances, initial_condition_array, time_span, data_dictionary)

  # solve -
  solution = solve(problem_object, alg_hints=[:auto], reltol=1e-8,abstol=1e-8)

  # pull solution apart -
  T = solution.t

  # initialize the state array -
  number_of_times_steps = length(T)
  number_of_states = length(initial_condition_array)
  X = zeros(number_of_times_steps,number_of_states)
  for step_index=1:number_of_times_steps

      # grab the solution 0
      soln_array = solution.u[step_index]
      for state_index = 1:number_of_states
          X[step_index, state_index] = soln_array[state_index]
      end
  end

  # return -
  return (T,X)
end

# ----------------------------------------------------------------------------------- #
# SolveBalances: Solves model equations from TSTART to TSTOP given parameters in data_dictionary.
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# TSIM - Simulation time vector
# X - Simulation state array (NTIME x NSPECIES)
# ----------------------------------------------------------------------------------- #
function SolveBalances(TSTART,TSTOP,Ts,data_dictionary)

    # Get required stuff from DataFile struct -
    time_span = (TSTART,TSTOP)
    initial_condition_array = data_dictionary["initial_condition_array"];

    # build problem object -
    problem_object = ODEProblem(Balances, initial_condition_array, time_span, data_dictionary)

    # solve -
    solution = DifferentialEquations.solve(problem_object, AutoTsit5(Rosenbrock23(autodiff=false)), reltol=1e-8,abstol=1e-8,saveat=Ts)

    # pull solution apart -
    T = solution.t

    # initialize the state array -
    number_of_times_steps = length(T)
    number_of_states = length(initial_condition_array)
    X = zeros(number_of_times_steps,number_of_states)
    for step_index=1:number_of_times_steps

        # grab the solution 0
        soln_array = solution.u[step_index]
        for state_index = 1:number_of_states
            X[step_index, state_index] = soln_array[state_index]
        end
    end

    # return -
    return (T,X)
end
