# include packages -
include("Include.jl")

function objective_function(parameter_guess_array,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # what is the host_type?
    host_type = :cell_free

    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

    # compute W -
    tmp_W_array = Float64[]
    for index = 1:2
        parameter_guess = parameter_guess_array[index]
        value = exp(-1*parameter_guess/(R*T_K))
        push!(tmp_W_array,value)
    end

    # update the control W's -
    control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
	control_parameter_dictionary["W_MDH_RNAP"] = tmp_W_array[1]       # 3
	control_parameter_dictionary["W_MDH_s70"] = tmp_W_array[2]   # 4
    model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
	binding_parameter_dictionary["n_MDH_s70"] = parameter_guess_array[3]
	binding_parameter_dictionary["K_MDH_s70"] = parameter_guess_array[4]
    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    # time constant modifier -
    # time constant modifiers -
	time_constant_modifier_array = [
		0.0	                        ;	# 1	MDH
		0.0	                        ;	# 2	s70
		parameter_guess_array[5] 	;	# 3	mRNA_MDH
		1.0	                        ;	# 4	mRNA_s70
		parameter_guess_array[6] 	;	# 5	protein_MDH
		1.0	                        ;	# 6	protein_s70
    ]
    model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # setup degradation_modifier_array -
    degradation_modifier_array = [
		0.0	;	# 1	MDH
		0.0	;	# 2	s70
		parameter_guess_array[7]	;	# 3	mRNA_MDH
		1.0	;	# 4	mRNA_s70
		parameter_guess_array[8]	;	# 5	protein_MDH
		parameter_guess_array[9]	;	# 6	protein_s70
    ]

    # update the translation time -
    model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[10]

    # lastly, update KL -
    biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[11]
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
    # ===================================================================================================== #

    # Phase 2:  solve model equations ===================================================================== #
    # solve the balance equations -
    (TSIM,XSIM) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)
    # ===================================================================================================== #

    # Phase 3:  compute simulation error ================================================================== #

	tsim_exp = exp_data_dictionary["mRNA_data_array"][:,1]

    # GFP protein -
    itp_gfp_protein =  LinearInterpolation(TSIM, XSIM[:,5]);
    protein_GFP_sim = itp_gfp_protein[tsim_exp]

    # get experimental data -
    mRNA_GFP_exp = exp_data_dictionary["mRNA_data_array"][:,2]          # mean is col 2 nM
    mRNA_GFP_std_exp = exp_data_dictionary["mRNA_data_array"][:,3]      # stdev is col 3 nM

    protein_GFP_exp = exp_data_dictionary["prot_data_array"][:,2]       # mean is col 2 muM
    protein_GFP_std_exp = exp_data_dictionary["prot_data_array"][:,3]   # stdev is col 3 muM

	# print(protein_GFP_exp, protein_GFP_sim)

	# compute error terms -
    # protein MDH -
    error_vector = (protein_GFP_exp .- protein_GFP_sim)
    error_term_array = transpose(error_vector)*error_vector
    # ===================================================================================================== #

    # error total -
    error_total = error_term_array

    # return -
    return error_total
end

function main(path_to_data_dir::String; algorithm::Symbol = :LBFGS)

    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # load the default data_dictionary -
    host_type = :cell_free
    path_to_biophysical_constants_file = "./CellFree.json"
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the objective function -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # setup paramter bounds -
    pvec_bounds = [

        # dG's -
        20000.0  80000.0    ;   # 1     W_MDH_RNAP
        -50000.0  -100.0    ;   # 2     W_MDH_s70

        # binding parameters -
        0.5 10.0            ;   # 3     n_MDH_s70
        0.001 100.0         ;   # 4     K_MDH_s70

        # time constants -
		0.001 100.0         ;	# 5	    mRNA_MDH
		0.001 100.0         ;	# 6	    protein_MDH

        # degradation mods -
		0.01 100.0 	    ;	# 7	    mRNA_MDH
		0.01 100.0 	    ;	# 8	    protein_MDH
        0.01 50.0 	    ;	# 9	    protein_s70

         # w -
         4.0 10.0           ;   # 10    translation capacity half-life

        # KL value -
        100.0 500.0         ;   # 11    KL in muM
    ];

    # setup initial condition vector -
    pvec_initial = [

        # dG's -
        40000.0     ;   # 1     W_MDH_RNAP
        -25000.0    ;   # 2     W_MDH_s70

        # binding parameters -
        2.0         ;   # 3     n_MDH_s70
        0.15        ;   # 4     K_MDH_s70

        # time constants -
		1.0         ;	# 5	    mRNA_MDH
		3.0         ;	# 6	    protein_MDH

        # degradation mods -
		0.1 	    ;	# 7	    mRNA_MDH
		3.0 	    ;	# 8	    protein_MDH
        1.0 	    ;	# 9	    protein_s70

        # w -
        8.0         ;   # 10    translation capacity half-life

        # KL -
        250.0       ;   # 11    KL in muM
    ];

    # check bounds -
    number_of_parameters = length(pvec_initial)
    for parameter_index = 1:number_of_parameters

        # what is the parameter value?
        p_i = pvec_initial[parameter_index]

        # is p_i outside of the bounds?
        lb_value = pvec_bounds[parameter_index,1]
        ub_value = pvec_bounds[parameter_index,2]

        if (p_i<lb_value)
            pvec_initial[parameter_index,1] = lb_value
        end

        if (p_i>ub_value)
            pvec_initial[parameter_index,1] = ub_value
        end
    end

    # search -
    options = Optim.Options(show_trace=true,show_every=100,iterations=8000)
    result = optimize(OF,pvec_bounds[:,1],pvec_bounds[:,2],pvec_initial,SAMIN(),options)

    return result
end

# setup paths -
path_to_data_dir = "$(pwd())/data"

# generate an ensemble of parameters -
number_of_p_sets = 5
number_of_parameters = 11+1 # end + 1 hold the error -
ensemble_array = zeros(number_of_parameters,number_of_p_sets)
algorithm_flag = :SimulatedAnnealing
for pset_index = 1:number_of_p_sets

    # call main -
    opt_result = main(path_to_data_dir; algorithm = algorithm_flag)

    # grab the best pset, and then dump to disk -
    pV = Optim.minimizer(opt_result)
    ensemble_array[end,pset_index] = Optim.minimum(opt_result)
    for (p_index,p_value) in enumerate(pV)
        ensemble_array[p_index,pset_index] = p_value
    end
end

# dump -
writedlm("Ensemble-T4.dat",ensemble_array)
