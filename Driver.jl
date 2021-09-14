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

# include -
include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 16.0
time_step_size = 0.01

# what is the host_type?
host_type = :cell_free

# path to parameters -
path_to_biophysical_constants_file = "./CellFree.json"

# Load the data dictionary (uses the default biophysical_constants file)
default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
data_dictionary = customize_data_dictionary(default_data_dictionary, host_type)

# Species_ensemble

for i in 1:size(PC)[2]

  parameter_guess_array = PC[:,i]

  R = data_dictionary["R"]
  T_K = data_dictionary["T_K"]

  # compute W -
  tmp_W_array = Float64[]
  for index = 1:2
      parameter_guess = parameter_guess_array[index]
      value = exp(-1*(parameter_guess/100)/(R*T_K))
      push!(tmp_W_array,value)
  end

    # update the control W's -
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary["W_MDH_RNAP"] = tmp_W_array[1]       # 3
  control_parameter_dictionary["W_MDH_s70"] = tmp_W_array[2]   # 4
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

  binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
  binding_parameter_dictionary["n_MDH_s70"] = parameter_guess_array[3]
  binding_parameter_dictionary["K_MDH_s70"] = parameter_guess_array[4]
  data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    # time constant modifiers -
  time_constant_modifier_array = [
    0.0	                        ;	# 1	MDH
    0.0	                        ;	# 2	s70
    parameter_guess_array[5] 	;	# 3	mRNA_MDH
    1.0	                        ;	# 4	mRNA_s70
    parameter_guess_array[6] 	;	# 5	protein_MDH
    1.0	                        ;	# 6	protein_s70
    ]
    data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # setup degradation_modifier_array -
    degradation_modifier_array = [
    0.0	;	# 1	MDH
    0.0	;	# 2	s70
    parameter_guess_array[7]	;	# 3	mRNA_MDH
    1.0	;	# 4	mRNA_s70
    parameter_guess_array[8]	;	# 5	protein_MDH
    parameter_guess_array[9]	;	# 6	protein_s70
    ]

  data_dictionary["degradation_modifier_array"] = degradation_modifier_array

  # update the translation time -
  data_dictionary["half_life_translation_capacity"] = parameter_guess_array[10]

  # lastly, update KL -
  biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
  biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[11]

  # grab defaults -
  species_symbol_type_array = data_dictionary["species_symbol_type_array"]
  protein_coding_length_array = data_dictionary["protein_coding_length_array"]
  gene_coding_length_array = data_dictionary["gene_coding_length_array"]
  time_constant_modifier_array = data_dictionary["time_constant_modifier_array"]
  initial_condition_array = data_dictionary["initial_condition_array"]

  # # get gene IC -
  idx_gene = findall(x->x==:gene,species_symbol_type_array)
  gene_abundance_array = initial_condition_array[idx_gene]

  # Precompute the translation parameters -
  translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
  data_dictionary["translation_parameter_array"] = translation_parameter_array

  # Precompute the kinetic limit of transcription -
  transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
  data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

  # Dilution degrdation matrix -
  dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
  data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

  # Solve the model equations -
  (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)

  # Species_ensemble[i,:] = X

  # T_exp = [0, 2, 4, 6, 8, 16]
  # X_exp = [0, 6.284153543, 11.46068615, 15.45104172, 18.37365078, 21.29625984]

  # plot protein
  T_exp = [0, 4, 6]
  X_exp = [0, 11.46068615, 15.45104172]

  plot(T_exp, X_exp, "o", color="black")
  plot(T,X[:,5])
  ylabel("MDH Protein (μm)")
  xlabel("Time (hr)")

  # # plot mRNA
  # plot(T,X[:,3]*1e3)
  # ylabel("MDH mRNA (nm)")
  # xlabel("Time (hr)")


end

# Species_mean = mean(Species_ensemble, dims = 3)
# Species_err = std(Species_ensemble, dims = 3)
# Species_pos = Species_mean + 1.96*Species_err
# Species_neg = Species_mean - 1.96*Species_err
