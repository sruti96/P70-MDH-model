# ----------------------------------------------------------------------------------- #
# Copyright (c) 2020 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
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
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2020-02-05T11:38:24.145
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{String,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = "./Default.json", host_type::Symbol = :bacteria)::Dict{String,Any}

	# load the biophysical_constants dictionary
	biophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of species types -
	species_symbol_type_array = [
		:gene	;	# 1	MDH
		:gene	;	# 2	s70
		:mrna	;	# 3	mRNA_MDH
		:mrna	;	# 4	mRNA_s70
		:protein	;	# 5	protein_MDH
		:protein	;	# 6	protein_s70
	]

	# we need to store the species symbol array for later -
	biophysical_constants_dictionary["species_symbol_type_array"] = species_symbol_type_array

	# array of gene lengths -
	gene_coding_length_array = [
		1000.0	;	# 1	MDH
		1000.0	;	# 2	s70
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 3	1	mRNA_MDH
		gene_coding_length_array[2]	;	# 4	2	mRNA_s70
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 5	1	protein_MDH
		round((0.33)*mRNA_coding_length_array[2])	;	# 6	2	protein_s70
	]

	# array of gene concentrations -
	gene_abundance_array = [
		5.0	;	# (nM) 1	MDH
		5.0	;	# (nM) 2	s70
	]

	# initial condition array -
	initial_condition_array = [
		gene_abundance_array[1]	;	# 1	MDH
		gene_abundance_array[2]	;	# 2	s70
		0.0	;	# 3	mRNA_MDH
		0.0	;	# 4	mRNA_s70
		0.0	;	# 5	protein_MDH
		0.0	;	# 6	protein_s70

		# translation capacity -
		100.0	;	# translation capacity
	]

	binding_parameter_dictionary = Dict{String,Float64}()
	binding_parameter_dictionary["n_MDH_s70"] = 1.0
	binding_parameter_dictionary["K_MDH_s70"] = 0.05

	# Alias the control function parameters -
	control_parameter_dictionary = Dict{String,Float64}()
	control_parameter_dictionary["W_MDH_RNAP"] = 0.001
	control_parameter_dictionary["W_MDH_s70"] = 1.0
	control_parameter_dictionary["W_s70_RNAP"] = 0.001

	# degradation modifiers -
	degradation_modifier_array = [
		0.0	;	# 1	MDH
		0.0	;	# 2	s70
		1.0	;	# 3	mRNA_MDH
		1.0	;	# 4	mRNA_s70
		1.0	;	# 5	protein_MDH
		1.0	;	# 6	protein_s70
	]

	# time constant modifiers -
	time_constant_modifier_array = [
		0.0	;	# 1	MDH
		0.0	;	# 2	s70
		1.0	;	# 3	mRNA_MDH
		1.0	;	# 4	mRNA_s70
		1.0	;	# 5	protein_MDH
		1.0	;	# 6	protein_s70
	]

	# Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)

	# Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array, host_type)

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)

	# Parameter name index array -
	parameter_name_mapping_array = [
		"n_MDH_s70"	;	# 1
		"K_MDH_s70"	;	# 2
		"W_MDH_RNAP"	;	# 3
		"W_MDH_s70"	;	# 4
		"W_s70_RNAP"	;	# 5
		"rnapII_concentration"	;	# 6
		"ribosome_concentration"	;	# 7
		"degradation_constant_mRNA"	;	# 8
		"degradation_constant_protein"	;	# 9
		"kcat_transcription"	;	# 10
		"kcat_translation"	;	# 11
		"maximum_specific_growth_rate"	;	# 12
		"saturation_constant_transcription"	;	# 13
		"saturation_constant_translation"	;	# 14
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["species_symbol_type_array"] = species_symbol_type_array
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
	data_dictionary["translation_parameter_array"] = translation_parameter_array
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

	# extra stuff -
	data_dictionary["R"] = 8.314 			# J mol^-1 K^-1
	data_dictionary["T_K"] = 273.15 + 29.0 	# K
	data_dictionary["half_life_translation_capacity"] = 5.17	 # hr
	data_dictionary["base_line_weight"] = (1.0/1.0)
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end

function customize_data_dictionary(default_data_dictionary::Dict{String,Any}, host_type::Symbol)::Dict{String,Any}

	# make a deepcopy -
	customized_data_dictionary = deepcopy(default_data_dictionary)

	# set the gene lengths -
	gene_length_array = customized_data_dictionary["gene_coding_length_array"]
	gene_length_array[1] = 906.0    # MDH
	gene_length_array[2] = 720.0    # S70 => doesn't matter, gene expression = 0
	customized_data_dictionary["gene_coding_length_array"] = gene_length_array

	# set the protein lengths -
	protein_coding_length_array = customized_data_dictionary["protein_coding_length_array"]
	protein_coding_length_array[1] = 299.0  # MDH
	protein_coding_length_array[2] = 240.0  # S70 => doesn't matter, translation = 0
	customized_data_dictionary["protein_coding_length_array"] = protein_coding_length_array

	# set ICs -
	initial_condition_array = customized_data_dictionary["initial_condition_array"]
	initial_condition_array[1] = 0.01  # MDH gene
	initial_condition_array[2] = 0.0    # no S70 gene, just S70 protein in the extract
	initial_condition_array[6] = 0.035 	# S70 protein in muM
	customized_data_dictionary["initial_condition_array"] = initial_condition_array

	# setup the W's -
	R = customized_data_dictionary["R"]
	T_K = customized_data_dictionary["T_K"]
	control_parameter_dictionary = customized_data_dictionary["control_parameter_dictionary"]
  	control_parameter_dictionary["W_MDH_RNAP"] = exp(-1*(47309/100)/(R*T_K))
  	control_parameter_dictionary["W_MDH_s70"] = exp(-1*(-39458/100)/(R*T_K))
	control_parameter_dictionary["W_s70_RNAP"] = 0.0
	customized_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

	# setup binding parameters -
	binding_parameter_dictionary = customized_data_dictionary["binding_parameter_dictionary"]
	binding_parameter_dictionary["n_MDH_s70"] = 0.5
	binding_parameter_dictionary["K_MDH_s70"] = 20
	customized_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

	# setup degradation_modifier_array -
	degradation_modifier_array = [
		0.0		;	# 1	MDH
		0.0		;	# 2	s70
		1.34	;	# 3	mRNA_MDH
		1.0	;	# 4	mRNA_s70
		3.16		;	# 5	protein_MDH
		4.77	;	# 6	protein_s70
	]
	customized_data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	biophysical_constants_dictionary = customized_data_dictionary["biophysical_constants_dictionary"]
	species_symbol_type_array = customized_data_dictionary["species_symbol_type_array"]
	protein_coding_length_array = customized_data_dictionary["protein_coding_length_array"]
	gene_coding_length_array = customized_data_dictionary["gene_coding_length_array"]

	# get gene IC -
	idx_gene = findall(x->x==:gene,species_symbol_type_array)
	gene_abundance_array = initial_condition_array[idx_gene]

	# no s70 gene -
	gene_abundance_array[2] = 0.0

	# Dilution degrdation matrix -
  	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)
	customized_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

	# time constant modifiers -
	time_constant_modifier_array = [
		0.0	;	# 1	MDH
		0.0	;	# 2	s70
		0.43	;	# 3	mRNA_MDH
		1.0	;	# 4	mRNA_s70
		1.52	;	# 5	protein_MDH
		1.0 ;	# 6	protein_s70
	]
	customized_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

	# Precompute the translation parameters -
  	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
	customized_data_dictionary["translation_parameter_array"] = translation_parameter_array

  	# Precompute the kinetic limit of transcription -
  	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
	customized_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

	# return -
	return customized_data_dictionary
  end
