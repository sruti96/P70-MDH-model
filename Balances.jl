# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
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
# estimate_steady_state: estimates the steady-state vector for the dynamic system
#
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# XSS - steady state vector
# ----------------------------------------------------------------------------------- #
function estimate_steady_state(data_dictionary)

    # Get the initial conditions -
    initial_condition_array = data_dictionary["initial_condition_array"];

    # calculate the steady-state soln -
    steady_state_prob = SteadyStateProblem(Balances, initial_condition_array, data_dictionary)
    steady_state_soln  = solve(steady_state_prob, DynamicSS(AutoTsit5(Rosenbrock23(autodiff=false));abstol=1e-8,reltol=1e-6,tspan=Inf))

    # return -
    return abs.(steady_state_soln.u)
end

# ----------------------------------------------------------------------------------- #
# AdjBalances: Evaluates adjoint model equations given time, state, parameter_index,
# and the data_dictionary.
#
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# t  - current time
# dx - derivative vector
# x  - state array
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# N/A - derivative array at current time step is filled in the function
# ----------------------------------------------------------------------------------- #
function AdjBalances(dx,x,data_dictionary,t)

    # look up the number of states -
    number_of_states = data_dictionary["number_of_states"]

    # grab the parameter index -
    parameter_index = data_dictionary["parameter_index"]

    # partition that state -
    state_array = x[1:number_of_states]
    sensitivity_array = x[(number_of_states+1):end]

    # call the balances -
    dxdt_array = calculate_balances(t,state_array,data_dictionary)

    # Calculate the sensitivity states -
    local_data_dictionary = deepcopy(data_dictionary)
    JM = calculate_jacobian(t,state_array,local_data_dictionary)
    BM = calculate_bmatrix(t,state_array,local_data_dictionary)

    # calulate the sensitivity state -
    dsdt_array = JM*sensitivity_array+BM[:,parameter_index]
    r_array = [dxdt_array ; dsdt_array]

    # package -
    number_of_total_states = length(r_array)
    for index = 1:number_of_total_states
        dx[index] = r_array[index]
    end
end

# ----------------------------------------------------------------------------------- #
# calculate_balances: Evaluates model equations given time, state and the data_dictionary.
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# t  - current time
# x  - state array
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# dxdt - derivative array at current time step
# ----------------------------------------------------------------------------------- #
function calculate_balances(t,x,data_dictionary)

    # get the structure arrays from the data_dictionary -
    AM = data_dictionary["dilution_degradation_matrix"]
    SM = data_dictionary["stoichiometric_matrix"]
    half_life_translation = data_dictionary["half_life_translation_capacity"]

    # what is my system size?
    number_of_states = data_dictionary["number_of_states"]

    # calculate the kinetics array -
    kinetics_array = calculate_txtl_kinetics_array(t,x,data_dictionary)

    # calculate the control array -
    transcription_control_array = calculate_transcription_control_array(t,x,data_dictionary)
    translation_control_array = calculate_translation_control_array(t,x,data_dictionary)

    # control_array is the transcription and translation control arrays -
    control_array = [transcription_control_array ; translation_control_array]

    # modfiy the kinetics -
    rV = kinetics_array.*control_array

    # compute the dxdt vector -
    dxdt = AM*x+SM*rV

    # return -
    return dxdt
end

# ----------------------------------------------------------------------------------- #
# Balances: Evaluates model equations given time, state and the data_dictionary.
# Type: GRN-JULIA
# Version: 1.0
#
# Input arguments:
# dx - derivative array
# t  - current time
# x  - state array
# data_dictionary  - Data dictionary instance (holds model parameters)
#
# Return arguments:
# N/A - dx is filled in the function and returned to the caller
# ----------------------------------------------------------------------------------- #
function Balances(dx,x,data_dictionary,t)

    # get the structure arrays from the data_dictionary -
    AM = data_dictionary["dilution_degradation_matrix"]
    SM = data_dictionary["stoichiometric_matrix"]
    half_life_translation = data_dictionary["half_life_translation_capacity"]

    # what is my system size?
    number_of_states = data_dictionary["number_of_states"]

    # calculate the kinetics array -
    kinetics_array = calculate_txtl_kinetics_array(t,x,data_dictionary)

    # calculate the control array -
    transcription_control_array = calculate_transcription_control_array(t,x,data_dictionary)
    translation_control_array = calculate_translation_control_array(t,x,data_dictionary)

    # control_array is the transcription and translation control arrays -
    control_array = [transcription_control_array ; translation_control_array]

    # modfiy the kinetics -
    rV = kinetics_array.*control_array

    # remove last n states -
    x_main = x[1:6]

    # compute the dxdt vector -
    dxdt = AM*x_main+SM*rV

    # package -
    for index = 1:number_of_states
        dx[index] = dxdt[index]
    end

    # extra states -
    dx[7] = -(log(2)/half_life_translation)*x[7]
end
