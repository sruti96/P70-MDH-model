# load my files -
include("Include.jl")

function update_model_dictionary(parameter_array,default_model_dictionary)

    # what is the host_type?
    host_type = :cell_free

    # make a deepcopy -
    model_data_dictionary = deepcopy(default_model_dictionary)

    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

    # compute W -
    tmp_W_array = Float64[]
    for index = 1:2
        parameter_guess = parameter_array[index]
        value = exp(-1*(parameter_guess/100)/(R*T_K))
        push!(tmp_W_array,value)
    end

    # update the control W's -
    control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
	control_parameter_dictionary["W_MDH_RNAP"] = tmp_W_array[1]       # 3
	control_parameter_dictionary["W_MDH_s70"] = tmp_W_array[2]   # 4
    model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
	binding_parameter_dictionary["n_MDH_s70"] = parameter_array[3]
	binding_parameter_dictionary["K_MDH_s70"] = parameter_array[4]
    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    # time constant modifier -
    # time constant modifiers -
	time_constant_modifier_array = [
		0.0	                        ;	# 1	MDH
		0.0	                        ;	# 2	s70
		parameter_array[5] 	        ;	# 3	mRNA_MDH
		1.0	                        ;	# 4	mRNA_s70
		parameter_array[6] 	        ;	# 5	protein_MDH
		1.0	                        ;	# 6	protein_s70
    ]
    model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # setup degradation_modifier_array -
    degradation_modifier_array = [
		0.0	;	# 1	MDH
		0.0	;	# 2	s70
		parameter_array[7]	;	# 3	mRNA_MDH
		1.0	;	# 4	mRNA_s70
		parameter_array[8]	;	# 5	protein_MDH
		parameter_array[9]	;	# 6	protein_s70
    ]

    # update the translation time -
    model_data_dictionary["half_life_translation_capacity"] = parameter_array[10]

    # lastly, update KL -
    biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_array[11]
    model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

    # grab defaults -
    species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
    protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
    gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
    time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
    initial_condition_array = model_data_dictionary["initial_condition_array"]

    # # get gene IC -
    idx_gene = findall(x->x==:gene,species_symbol_type_array)
    gene_abundance_array = initial_condition_array[idx_gene]

    # Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
    model_data_dictionary["translation_parameter_array"] = translation_parameter_array

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
    model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

    # Dilution degrdation matrix -
    dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
    model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

    # return -
    return model_data_dictionary
end


function main(path_to_ensemble_file::String, path_to_sim_dir::String)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free
    path_to_biophysical_constants_file = "./CellFree.json"

    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
    customized_data_dictionary = customize_data_dictionary(default_data_dictionary, host_type)

    # load the ensemble file -
    parameter_ensemble_array = readdlm(path_to_ensemble_file)
    (number_of_parameters, number_of_samples) = size(parameter_ensemble_array)

    for i = 1:100

        # grab the parameter array -
        local_parameter_array = parameter_ensemble_array[:,i]

		# sample TXTL parameters
		biophysical_constants_dictionary = customized_data_dictionary["biophysical_constants_dictionary"]
		biophysical_constants_dictionary["RNAPII_concentration"] = 0.060 + rand(1)[1]*(0.075-0.060); #RNAP_concentration  #uM # 60-75nM (ACS SynBio Garamella 2016)
	    biophysical_constants_dictionary["transcription_elongation_rate"] = 20.0 + rand(1)[1]*(30-20);  #max_transcription_rate  # >5 NT/s (ACS SynBio Garamella 2016)
	    biophysical_constants_dictionary["ribosome_concentration"] = 2.0 + rand(1)[1]*(2.5-2.0); #RIBOSOME_concentration #uM  <0.0023mM (ACS SynBio Garamella 2016)
	    biophysical_constants_dictionary["translation_elongation_rate"] = 12.0 + rand(1)[1]*(18-10); #max_translation_rate #>1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)
		customized_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

        model_data_dictionary = update_model_dictionary(local_parameter_array, customized_data_dictionary)

		# solve the model equations -
        (T,X) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)

        # dump -
        data_array = [T X]
        filename = "$(path_to_sim_dir)/simulation-P$(i).dat"
        writedlm(filename, data_array)

        # give user some notification -
        @show i
    end
end


path_to_sim_file = "$(pwd())/simulations"
path_to_ensemble_file = "$(pwd())/Poets-T10.dat"
main(path_to_ensemble_file, path_to_sim_file)
