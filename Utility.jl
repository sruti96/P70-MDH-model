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

function contains(string,token)
    return occursin(token,string)
end

function eye(number_of_rows,number_of_cols)
    return Matrix{Float64}(I,number_of_rows, number_of_cols)
end

# load the files -
searchdir(path,key) = filter(x->contains(x,key),readdir(path))

function sort_parameter_ensemble_array(parameter_ensemble_array::Array{Float64,2})

  # The last row is the error -
  error_array = parameter_ensemble_array[end,:]

  # what is the perm index?
  idx_sort = sortperm(error_array) # default: small to big

  # shuffle cols -
  sorted_parameter_array = parameter_ensemble_array[:,idx_sort]

  # return -
  return sorted_parameter_array
end

function parameter_bounds_array_from_ensemble(path_to_ensemble_file::String)

  # setup ranges -
  sample_bounds_array = Array{Tuple,1}()
  ensemble_array = readdlm(path_to_ensemble_file)
  (number_of_parameters,number_of_trials) = size(ensemble_array)
  for parameter_index = 1:(number_of_parameters-1)

      # get row of parameters -
      parameter_row = ensemble_array[parameter_index,:]
      min_value = minimum(parameter_row)
      max_value = maximum(parameter_row)

      # create the tuple -
      tmp_tuple = (min_value,max_value)

      # cache -
      push!(sample_bounds_array,tmp_tuple)
  end

  # return -
  return sample_bounds_array
end

function load_experimental_data_dictionary(base_path::String)

  # create experimental data dictionary -
  exp_data_dictionary = Dict{String,DataFrame}()

  # load data -
  full_data_table_path = "$(base_path)/Data3.csv"
  full_data_table = CSV.read(full_data_table_path)

  # mRNA -
  T = full_data_table[!,:time_hr]
  mean_mRNA = full_data_table[!,:mean_mRNA_nM]
  stdev_mRNA = full_data_table[!,:stdev_mRNA_nM]
  mRNA_data_array = [T mean_mRNA stdev_mRNA]

  # protein -
  T = full_data_table[!,:time_hr]
  mean_MDH = full_data_table[!,:mean_MDH_uM]
  stdev_MDH = full_data_table[!,:stdev_MDH_uM]
  prot_data_array = [T mean_MDH stdev_MDH]

  # package -------------------------------------------------- #
  exp_data_dictionary["mRNA_data_array"] = mRNA_data_array
  exp_data_dictionary["prot_data_array"] = prot_data_array
  return exp_data_dictionary
  # ---------------------------------------------------------- #
end

function calculate_fisher_information_matrix(scaled_sensitivity_array,measurement_weight_array,id_parameter_index_array)

  # Get the size of the system -
  (number_of_rows,number_of_cols) = size(scaled_sensitivity_array);

  # grab the id cols -
  id_sensitivity_array = scaled_sensitivity_array[:,id_parameter_index_array]

  # Create the FIM (fisher_information_matrix) -
  FIM = transpose(id_sensitivity_array)*measurement_weight_array*id_sensitivity_array

  # return the FIM to the caller -
  return FIM
end

function estimate_identifiable_parameters(scaled_sensitivity_array,epsilon)

  # Get the size of the system -
  (number_of_rows,number_of_cols) = size(scaled_sensitivity_array);
  X = zeros(number_of_rows,1)
  pset = Int64[]

  R = scaled_sensitivity_array

  for col_index = 1:number_of_cols

    # get mag of col -
    local_m_array = zeros(number_of_cols)
    for inner_col_index = 1:number_of_cols

      value = R[:,inner_col_index]'*R[:,inner_col_index]
      local_m_array[inner_col_index] = value[1]
    end

    # what is the maximum element?
    max_value = maximum(local_m_array)
    max_index = indmax(local_m_array)

    # check tolerance -
    if max_value>epsilon

      # Grab this parameter -
      push!(pset,max_index)

      # Grab this col -
      X = [X scaled_sensitivity_array[:,max_index]]

      # remove leading col if first time through -
      if (col_index == 1)
        X = X[:,2:end]
      end

      # Create local array -
      Shat=X*inv(X'*X)*X'*scaled_sensitivity_array

      # Update R -
      R = scaled_sensitivity_array - Shat
    end
  end

  return sort!(pset)
end

function calculate_sensitivity_array(path_to_senstivity_files,file_pattern,time_skip,data_dictionary)

  # what is my system dimension?
  number_of_states = data_dictionary["number_of_states"]

  # load the files -
  searchdir(path,key) = filter(x->contains(x,key),readdir(path))

  # block_dictionary -
  block_dictionary = Dict()

  # build file list -
  list_of_files = searchdir(path_to_senstivity_files,file_pattern)
  number_of_files = length(list_of_files)
  time_array = []
  for file_index = 1:number_of_files

    # Build path -
    path_to_data_file = path_to_senstivity_files*"/"*file_pattern*string(file_index)*".dat"

    # Load file -
    local_data_array = readdlm(path_to_data_file)

    # split -
    time_array = local_data_array[:,1]
    X = local_data_array[:,2:end]
    state_array = X[:,1:number_of_states]
    sensitivity_array = X[:,(number_of_states+1):end]
    scaled_sensitivity_block = scale_sensitivity_array(time_array,state_array,sensitivity_array,file_index,data_dictionary)

    # store the transpose -
    key_symbol = file_pattern*string(file_index)*".dat"
    block_dictionary[key_symbol] = transpose(scaled_sensitivity_block)
  end

  # what is my system dimension?
  number_of_timesteps = length(time_array)
  number_of_parameters = number_of_files

  # initialize -
  sensitivity_array = zeros(number_of_states,number_of_parameters)
  sample_time_array = Float64[]
  for time_step_index = 1:time_skip:number_of_timesteps

    time_value = time_array[time_step_index]
    push!(sample_time_array,time_value)

    local_sens_block = zeros(number_of_states,number_of_parameters)
    for parameter_index = 1:number_of_parameters

      # get the block from the dictionary -
      key_symbol = file_pattern*string(parameter_index)*".dat"
      block = block_dictionary[key_symbol]

      # grab the col -
      block_col = block[:,time_step_index]

      for state_index = 1:number_of_states
        local_sens_block[state_index,parameter_index] = block_col[state_index]
      end
    end

    # add block to s array -
    sensitivity_array = [sensitivity_array ; local_sens_block]
  end

  # cutoff leading block -
  sensitivity_array = sensitivity_array[(number_of_states+1):end,:]
  return (sample_time_array,sensitivity_array)
end

function calculate_average_scaled_sensitivity_array(path_to_senstivity_files,file_pattern,data_dictionary)

  # what is my system dimension?
  number_of_states = data_dictionary["number_of_states"]

  # initialize -
  average_scaled_sensitivity_array = zeros(number_of_states,1)

  # load the files -
  searchdir(path,key) = filter(x->contains(x,key),readdir(path))

  # build file list -
  list_of_files = searchdir(path_to_senstivity_files,file_pattern)
  number_of_files = length(list_of_files)
  for file_index = 1:number_of_files

    # Build path -
    path_to_data_file = path_to_senstivity_files*"/"*file_pattern*string(file_index)*".dat"

    # Load file -
    local_data_array = readdlm(path_to_data_file)

    # split -
    time_array = local_data_array[:,1]
    X = local_data_array[:,2:end]
    state_array = X[:,1:number_of_states]
    sensitivity_array = X[:,(number_of_states+1):end]
    scaled_sensitivity_array = scale_sensitivity_array(time_array,state_array,sensitivity_array,file_index,data_dictionary)

    # time average -
    average_sensitivity_col = time_average_array(time_array,scaled_sensitivity_array)

    # grab -
    average_scaled_sensitivity_array = [average_scaled_sensitivity_array average_sensitivity_col]
  end

  # trim leading col -
  average_scaled_sensitivity_array = average_scaled_sensitivity_array[:,2:end]
  return average_scaled_sensitivity_array
end

function time_average_array(time_array,data_array)

  # what is the delta T?
  delta_time = (time_array[end] - time_array[1])

  # initialize -
  average_array = Float64[]

  # what is the size of the array?
  (number_of_timesteps,number_of_states) = size(data_array)
  for state_index = 1:number_of_states

    # grab the data column -
    data_col = data_array[:,state_index]

    # average -
    average_value = (1/delta_time)*trapz(time_array,data_col)

    # push -
    push!(average_array,average_value)
  end

  return average_array

end

function scale_sensitivity_array(time_array,state_array,sensitivity_array,parameter_index,data_dictionary)

  # what is small?
  epsilon = 1e-6

  # initialize -
  (number_of_timesteps,number_of_states) = size(state_array)
  scaled_sensitivity_array = zeros(number_of_timesteps,number_of_states)

  # What is the nominal parameter value?
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  parameter_name = parameter_name_mapping_array[parameter_index]

  # build the total parameter dictionary -
  binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  txtl_parameter_dictionary = data_dictionary["biophysical_constants_dictionary"]
  total_parameter_dictionary = merge(binding_parameter_dictionary,control_parameter_dictionary,txtl_parameter_dictionary)

  # Grab the default value -
  default_parameter_value = total_parameter_dictionary[parameter_name]

  # main loop -
  for (time_index,time_value) in enumerate(time_array)

    for state_index = 1:number_of_states

      state_value = state_array[time_index,state_index]
      if (state_value<epsilon)
        state_value = epsilon
      end

      old_value = sensitivity_array[time_index,state_index]
      new_value = old_value*(default_parameter_value/state_value)
      scaled_sensitivity_array[time_index,state_index] = new_value
    end
  end

  return scaled_sensitivity_array
end

function calculate_jacobian(time,state_array,data_dictionary)

  # what is the size of the system?
  number_of_states = length(state_array)

  # calculate each row of the jacobian -
  jacobian_array = zeros(1,number_of_states)
  for (state_index,state_value) in enumerate(state_array)

    jacobian_row = calculate_jacobian_row(time,state_array,state_index,data_dictionary)
    jacobian_array = [jacobian_array  ; transpose(jacobian_row)]
  end

  jacobian_array = jacobian_array[2:end,:]
  return jacobian_array
end

function calculate_bmatrix(time,state_array,data_dictionary)

  # what is the size of the system?
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  number_of_parameters = length(parameter_name_mapping_array)
  number_of_states = length(state_array)

  # calculate each row of the jacobian -
  b_array = zeros(number_of_parameters,1)
  for (state_index,state_value) in enumerate(state_array)

    bmatrtix_row = calculate_bmatrix_row(time,state_array,state_index,data_dictionary)
    b_array = [b_array bmatrtix_row]
  end

  b_array = b_array[:,2:end]
  return transpose(b_array)
end

function calculate_bmatrix_row(time,state_array,balance_index,data_dictionary)

  # define some constants -
  epsilon = 1e-6
  delta = 0.01

  # parameter name dictionary -
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  txtl_parameter_dictionary = data_dictionary["biophysical_constants_dictionary"]
  number_of_parameters = length(parameter_name_mapping_array)

  # how many binding parameters do we have?
  number_of_binding_parameters = length(binding_parameter_dictionary)
  number_of_control_parameters = length(control_parameter_dictionary)

  # create a mega dictionary -
  total_parameter_dictionary = merge(binding_parameter_dictionary,control_parameter_dictionary,txtl_parameter_dictionary)

  # create delta parameter array -
  parameter_delta_array = Float64[]
  for (parameter_index,parameter_name) in enumerate(parameter_name_mapping_array)

    # Grab the default value -
    default_parameter_value = total_parameter_dictionary[parameter_name]

    # state perturbation -
    peturbed_parameter_value = default_parameter_value*(delta);

    #@show (peturbed_parameter_value,default_parameter_value)

    # check -
    if (peturbed_parameter_value<epsilon)
      peturbed_parameter_value = epsilon
    end

    # capture -
    push!(parameter_delta_array,(peturbed_parameter_value))
  end

  # Create the diag array -
  diag_delta_array = Diagonal(vec(parameter_delta_array))

  # Create bVec -
  f_nominal = calculate_balances(time,vec(state_array),data_dictionary)

  # estimate the perturbed balances -
  rhs_delta_array = Float64[]
  for (parameter_index,parameter_name) in enumerate(parameter_name_mapping_array)

    # copy -
    local_data_dictionary = deepcopy(data_dictionary)

    # Grab the default value -
    default_parameter_value = total_parameter_dictionary[parameter_name]

    # update the state -
    perturbed_parameter_array = zeros(number_of_parameters)
    for local_index = 1:number_of_parameters
      if (parameter_index == local_index)
        perturbed_parameter_array[parameter_index] = default_parameter_value*(1+delta);
      else
        local_parameter_name = parameter_name_mapping_array[local_index]
        perturbed_parameter_array[local_index] = total_parameter_dictionary[local_parameter_name]
      end
    end

    if (parameter_index<=number_of_binding_parameters)

      # we are in the binding section -
      local_data_dictionary["binding_parameter_dictionary"][parameter_name] = perturbed_parameter_array[parameter_index]
    elseif (parameter_index>number_of_binding_parameters && parameter_index<=(number_of_binding_parameters+number_of_control_parameters))
      # we are in the control section -
      local_data_dictionary["control_parameter_dictionary"][parameter_name] = perturbed_parameter_array[parameter_index]
    else
      # we are in the bar parameter section -
      local_data_dictionary[parameter_name] = perturbed_parameter_array[parameter_index]
    end

    # calculate the perturbed balances -
    f_perturbed = calculate_balances(time,vec(state_array),local_data_dictionary)
    f_perturbed = f_perturbed[balance_index] - f_nominal[balance_index]

    # capture -
    push!(rhs_delta_array,f_perturbed)
  end

  # calculate the bmatrix row -
  bmatrix_row = diag_delta_array\rhs_delta_array

  # return -
  return bmatrix_row
end

function finite_diff_jacobian(time,state_array,data_dictionary)

  # define some constants -
  epsilon = 1e-6
  delta = 0.05
  number_of_states = length(state_array)

  # initialize -
  jacobian_array = zeros(number_of_states,number_of_states)

  # nominal -
  f_nominal = calculate_balances(time,vec(state_array),data_dictionary)

  #@show state_array

  for row_index = 1:number_of_states
    for col_index = 1:number_of_states

      perturbed_state_array = zeros(number_of_states)
      for perturbation_index = 1:number_of_states

        if (col_index == perturbation_index)
          perturbed_state_array[col_index] = state_array[col_index]*(1+delta)
        else
          perturbed_state_array[perturbation_index] = state_array[perturbation_index]
        end
      end

      #@show perturbed_state_array

      # calculate the balances -
      f_perturbed = calculate_balances(time,vec(perturbed_state_array),data_dictionary)

      # calculate the entry -
      perturbation_size = state_array[col_index]*delta
      if (perturbation_size<epsilon)
        perturbation_size = epsilon
      end
      jacobian_array[row_index,col_index] = (f_perturbed[row_index] - f_nominal[row_index])/(perturbation_size)
    end
  end

  return jacobian_array
end

function calculate_jacobian_row(time,state_array,balance_index,data_dictionary)

  # define some constants -
  epsilon = 1e-6
  delta = 0.001
  number_of_states = length(state_array)

  # Create the delta state array -
  state_delta_array = Float64[]
  for (state_index,state_value) in enumerate(state_array)

    # state perturbation -
    peturbed_state = state_value*(delta);

    # check -
    if (peturbed_state<epsilon)
      peturbed_state = epsilon
    end

    # capture -
    push!(state_delta_array,peturbed_state)
  end

  # Create the diag array -
  diag_delta_array = Diagonal(vec(state_delta_array))

  # Create bVec -
  f_nominal = calculate_balances(time,vec(state_array),data_dictionary)

  # estimate the perturbed balances -
  rhs_delta_array = Float64[]
  for (state_index,state_value) in enumerate(state_array)

    # update the state -
    perturbed_state_array = zeros(number_of_states)

    for local_index = 1:number_of_states
      if (state_index == local_index)
        perturbed_state_array[state_index] = state_value*(1+delta);
      else
        perturbed_state_array[local_index] = state_array[local_index]
      end
    end


    # calculate the perturbed balances -
    f_perturbed = calculate_balances(time,vec(perturbed_state_array),data_dictionary)
    f_perturbed_delta = f_perturbed[balance_index] - f_nominal[balance_index]

    # capture -
    push!(rhs_delta_array,f_perturbed_delta)
  end

  # calculate the jacobian row -
  jacobian_row = diag_delta_array\rhs_delta_array

  # return -
  return jacobian_row
end

function build_biophysical_dictionary(path_to_biophysical_constants_file::String, host_type::Symbol)::Dict{String,Any}

    # check the path - is it legit?
    check_file_existence(path_to_biophysical_constants_file)

    # Load the model dictionary -
    model_dictionary = JSON.parsefile(path_to_biophysical_constants_file)

    # check the host type in the biophysical file - need to match
    file_host_type = Symbol(model_dictionary["host_type"])
    if (file_host_type != host_type)
        throw(error("inconsistent host_type. Expected $(host_type) found $(file_host_type)"))
    end

    # initialize -
    biophysical_dictionary = Dict{String,Any}()

    # convert -
    list_of_constant_dictionaries = model_dictionary["biophysical_constants"]
    for (key_string, local_dictionary) in list_of_constant_dictionaries

        # get the value -
        value = parse(Float64, local_dictionary["value"])

        # cache -
        biophysical_dictionary[key_string] = value
    end

    # return -
    return biophysical_dictionary
end

function trapz(x::Array{Float64}, y::Array{Float64})

    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end

    r = 0
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end
