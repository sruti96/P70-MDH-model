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
        value = exp(-1*(parameter_guess/100)/(R*T_K))
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

	model_data_dictionary["degradation_modifier_array"] = degradation_modifier_array

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
    # compute the error - we need to do a bunch of interpolation -
    tsim_exp = exp_data_dictionary["mRNA_data_array"][:,1]

    # GFP mRNA -
    itp_gfp_mRNA =  LinearInterpolation(TSIM, (1000)*XSIM[:,3]);
    mRNA_GFP_sim = itp_gfp_mRNA[tsim_exp]  # convert to muM from nM

    # GFP protein -
    itp_gfp_protein =  LinearInterpolation(TSIM, XSIM[:,5]);
    protein_GFP_sim = itp_gfp_protein[tsim_exp]

    # get experimental data -
    mRNA_GFP_exp = exp_data_dictionary["mRNA_data_array"][:,2]          # mean is col 2 nM
    mRNA_GFP_std_exp = exp_data_dictionary["mRNA_data_array"][:,3]      # stdev is col 3 nM

    protein_GFP_exp = exp_data_dictionary["prot_data_array"][:,2]       # mean is col 2 muM
    protein_GFP_std_exp = exp_data_dictionary["prot_data_array"][:,3]   # stdev is col 3 muM

    # compute error terms -
    error_term_array = zeros(1,1)

    # mRNA GFP -
    #mRNA_GFP_std_exp[1] = 1.0   # we have 0 ic
    #tmp_arr = 1.0./((mRNA_GFP_std_exp).^2)
    #W_mRNA = diagm(tmp_arr)
    # error_vector_1 = (mRNA_GFP_exp .- mRNA_GFP_sim)
    # error_term_array[1] = transpose(error_vector_1)*error_vector_1

    # protein MDH -
    #protein_GFP_std_exp[1] = 1.0
    #tmp_arr = 1.0./((protein_GFP_std_exp).^2)
    #W_prot = diagm(tmp_arr)
    error_vector_2 = (protein_GFP_exp .- protein_GFP_sim)
    error_term_array[1] = transpose(error_vector_2)*error_vector_2
    # ===================================================================================================== #

	# print(protein_GFP_exp, protein_GFP_sim)
    # return -
    return error_term_array
end

# Evaluates the objective function values -
function local_refienment_step(path_to_data_dir, parameter_array; sigma=0.05, iteration_max=100)

    # inner functions -
    function _compute_error_total(objective_array,W)
        value = transpose(objective_array)*W*objective_array
        return value[1]
    end

    # initialize -
    number_of_parameters = length(parameter_array)
    BIG = 1e10

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # wght array -
    W = diagm(ones(1))

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)

    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # calculate the starting error -
    parameter_array_best = parameter_array
    error_array = BIG*ones(4)
    error_array[1] = _compute_error_total(OF(parameter_array_best), W)

    # main refinement loop -
    iteration_counter = 1
    while (iteration_counter<iteration_max)

        # take a step up -
        parameter_up = parameter_array_best.*(1 .+ sigma*rand(number_of_parameters))
        parameter_up = check_parameter_bounds(parameter_up)

        # take a step down -
        parameter_down = parameter_array_best.*(1 .- sigma*rand(number_of_parameters))
        parameter_down = check_parameter_bounds(parameter_down)

        # Evaluate the obj function -
        error_array[2] = _compute_error_total(OF(parameter_up),W)
        error_array[3] = _compute_error_total(OF(parameter_down),W)

        # Calculate a correction factor -
        a = error_array[2]+error_array[3] - 2.0*error_array[1]
        parameter_corrected = parameter_array_best
        if (a>0.0)
            amda = -0.5*(error_array[3] - error_array[2])/a
            parameter_corrected = parameter_array_best .+ amda*rand(number_of_parameters)
            parameter_corrected = check_parameter_bounds(parameter_corrected)
            error_array[4] = _compute_error_total(OF(parameter_corrected), W)
        end

        # Which step has the min error?
        min_index = argmin(error_array)
        if (min_index == 1)
            parameter_array_best = parameter_array_best
        elseif (min_index == 2)
            parameter_array_best = parameter_up
        elseif (min_index == 3)
            parameter_array_best = parameter_down
        elseif (min_index == 4)
            parameter_array_best = parameter_corrected
        end

        # Update the local error
        error_array[1] = error_array[min_index]

        @show iteration_counter,error_array[min_index]

        # update local counter -
        iteration_counter = iteration_counter + 1
    end

    return parameter_array_best
end

function check_parameter_bounds(parameter_array)

    # setup paramter bounds -
    pvec_bounds = [

        # dG's -
        20000.0  80000.0    ;   # 1     W_MDH_RNAP
        -50000.0  -100.0    ;   # 2     W_MDH_s70

        # binding parameters -
        0.5 10.0            ;   # 3     n_MDH_s70
        0.001 10.0         ;   # 4     K_MDH_s70

        # time constants -
		0.001 10.0         ;	# 5	    mRNA_MDH
		0.001 10.0         ;	# 6	    protein_MDH

        # degradation mods -
		0.001 10.0 	    ;	# 7	    mRNA_MDH
		0.001 10.0 	    ;	# 8	    protein_MDH
        0.001 10.0 	    ;	# 9	    protein_s70

         # w -
         4.0 10.0           ;   # 10    translation capacity half-life

        # KL value -
        100.0 500.0         ;   # 11    KL in muM
    ];

    # tmp -
    pvec_initial = parameter_array

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

    # return -
    return pvec_initial
end

function neighbor_function(parameter_array; sigma=0.05)

    # setup -
    number_of_parameters = length(parameter_array)

    # calculate new parameter array -
    new_parameter_array = parameter_array.*(1 .+ sigma*randn(number_of_parameters))

    # check the bounds and return -
    return check_parameter_bounds(new_parameter_array)
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end

function acceptance_probability_function(rank_array,temperature)
    return (exp(-rank_array[end]/temperature))
end

function main(path_to_data_dir::String, initial_parameter_array::Array{Float64,1}; rank_cutoff::Int64=4, maximum_number_of_iterations::Int64=100)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 6.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free
    number_of_parameters = length(initial_parameter_array)

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)

    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)
    NF(P) = neighbor_function(P;sigma=0.01)

    # make call to POETs -
    (EC,PC,RA) = estimate_ensemble(OF,NF,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=rank_cutoff,maximum_number_of_iterations=maximum_number_of_iterations)

    # return -
    return (EC,PC,RA)
end

 # setup initial condition vector -
 pvec_initial = [

    # dG's -
    40000.0     ;   # 1     W_MDH_RNAP
    -25000.0    ;   # 2     W_MDH_s70

    # binding parameters -
    1.0         ;   # 3     n_MDH_s70
    30.0        ;   # 4     K_MDH_s70

    # time constants -
    1.0         ;	# 5	    mRNA_MDH
    1.0         ;	# 6	    protein_MDH

    # degradation mods -
    1.0 	    ;	# 7	    mRNA_MDH
    1.0 	    ;	# 8	    protein_MDH
    1.0 	    ;	# 9	    protein_s70

    # w -
    8.0         ;   # 10    translation capacity half-life

    # KL -
    250.0       ;   # 11    KL in muM
];

# setup -
path_to_data_dir = "$(pwd())/data"
pV = neighbor_function(pvec_initial; sigma=0.25)
EC = 0
PC = 0
RA = 0

# execute -
number_of_trials = 10
for trial_index = 1:number_of_trials

    global pV
    global EC
    global PC
    global RA

   # do a local step -
    if (mod(trial_index,2) == 0)

        # find the lowest score pV -
        sum_error_array = sum(EC,dims=1)
        best_p_index = argmin(vec(sum_error_array))
        pV_best = PC[:,best_p_index]

        # local refine -
        pV = local_refienment_step(path_to_data_dir, pV_best; iteration_max=500)
    end

    # main -
    (EC,PC,RA) = main(path_to_data_dir, vec(pV); rank_cutoff=3,maximum_number_of_iterations=100)

    # dump results to disk -
    fname = "./poets_ensemble/RA_T$(trial_index).dat"
    writedlm(fname,RA)
    fname = "./poets_ensemble/EC_T$(trial_index).dat"
    writedlm(fname,EC)
    fname = "./poets_ensemble/PC_T$(trial_index).dat"
    writedlm(fname,PC)
end

# include("Driver1.jl")
