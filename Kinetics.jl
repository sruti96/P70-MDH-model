# -- INTERNAL METHODS BELOW ---------------------------------------------------------- #
function precompute_translation_parameter_array_bacteria(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})::Array{TranslationParameters,1}

    # get some stuff from the biophysical_dictionary -
    mass_of_single_cell = biophysical_dictionary["mass_of_single_cell"]
    fraction_of_water_per_cell = biophysical_dictionary["fraction_of_water_per_cell"]
    ribosome_copy_number = biophysical_dictionary["copies_of_ribosome_per_cell"]
    translation_elongation_rate = biophysical_dictionary["translation_elongation_rate"]
    characteristic_initiation_time_translation = biophysical_dictionary["characteristic_initiation_time_translation"]
    KL = biophysical_dictionary["translation_saturation_constant"]
    species_symbol_type_array = biophysical_dictionary["species_symbol_type_array"]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the ribosome concentration -
    ribosome_concentration = ribosome_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # which indices are protein -
    idx_protein = findall(x->x==:protein,species_symbol_type_array)

    # initialize -
    translation_parameters_array = Array{TranslationParameters,1}()

    # compute the translation parameters -
    for (index,protein_length) in enumerate(protein_coding_length_array)

        # compute kE -
        kE = translation_elongation_rate*(1/protein_length)

        # compute kI -
        kI = (1/characteristic_initiation_time_translation)

        # what is the correct index?
        correction_index = idx_protein[index]
        time_scale_parameter = time_constant_modifier_array[correction_index]

        # compute the tau factor -
        tau_factor = (kE/kI)*time_scale_parameter

        # compute the vmax -
        translation_vmax = kE*ribosome_concentration

        # build -
        translationParameters = TranslationParameters()
        translationParameters.vmax = translation_vmax*(3600)    # convert to hr
        translationParameters.tau_factor = tau_factor
        translationParameters.KL = KL
        push!(translation_parameters_array, translationParameters)
    end

    return translation_parameters_array
end

function precompute_translation_parameter_array_cell_free(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})::Array{TranslationParameters,1}

    # get some stuff from the biophysical_dictionary -
    ribosome_concentration = biophysical_dictionary["ribosome_concentration"]
    translation_elongation_rate = biophysical_dictionary["translation_elongation_rate"]
    characteristic_initiation_time_translation = biophysical_dictionary["characteristic_initiation_time_translation"]
    KL = biophysical_dictionary["translation_saturation_constant"]
    species_symbol_type_array = biophysical_dictionary["species_symbol_type_array"]

    # initialize -
    translation_parameters_array = Array{TranslationParameters,1}()

    # which indices are protein -
    idx_protein = findall(x->x==:protein,species_symbol_type_array)

    # compute the translation parameters -
    for (index,protein_length) in enumerate(protein_coding_length_array)

        # compute kE -
        kE = translation_elongation_rate*(1/protein_length)

        # compute kI -
        kI = (1/characteristic_initiation_time_translation)

        # what is the correct index?
        correction_index = idx_protein[index]
        time_scale_parameter = time_constant_modifier_array[correction_index]

        # compute the tau factor -
        tau_factor = (kE/kI)*time_scale_parameter

        # compute the vmax -
        # use a polysome factor of 6 -
        translation_vmax = kE*ribosome_concentration

        # build -
        translationParameters = TranslationParameters()
        translationParameters.vmax = translation_vmax*(3600)    # convert to hr
        translationParameters.tau_factor = tau_factor
        translationParameters.KL = KL
        push!(translation_parameters_array, translationParameters)
    end

    return translation_parameters_array
end

function precompute_translation_parameter_array_mammalian(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})::Array{TranslationParameters,1}

    # get some stuff from the biophysical_dictionary -
    mass_of_single_cell = biophysical_dictionary["mass_of_single_cell"]
    fraction_of_water_per_cell = biophysical_dictionary["fraction_of_water_per_cell"]
    ribosome_copy_number = biophysical_dictionary["copies_of_ribosome_per_cell"]
    translation_elongation_rate = biophysical_dictionary["translation_elongation_rate"]
    characteristic_initiation_time_translation = biophysical_dictionary["characteristic_initiation_time_translation"]
    KL = biophysical_dictionary["translation_saturation_constant"]
    species_symbol_type_array = biophysical_dictionary["species_symbol_type_array"]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the ribosome concentration -
    ribosome_concentration = ribosome_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # which indices are protein -
    idx_protein = findall(x->x==:protein,species_symbol_type_array)

    # initialize -
    translation_parameters_array = Array{TranslationParameters,1}()

    # compute the translation parameters -
    for (index,protein_length) in enumerate(protein_coding_length_array)

        # compute kE -
        kE = translation_elongation_rate*(1/protein_length)

        # compute kI -
        kI = (1/characteristic_initiation_time_translation)

        # what is the correct index?
        correction_index = idx_protein[index]
        time_scale_parameter = time_constant_modifier_array[correction_index]

        # compute the tau factor -
        tau_factor = (kE/kI)*time_scale_parameter

        # compute the vmax -
        translation_vmax = kE*ribosome_concentration

        # build -
        translationParameters = TranslationParameters()
        translationParameters.vmax = translation_vmax*(3600)    # convert to hr
        translationParameters.tau_factor = tau_factor
        translationParameters.KL = KL
        push!(translation_parameters_array, translationParameters)
    end

    return translation_parameters_array
end


function precompute_translation_parameter_array(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1}, host_type::Symbol = :bacteria)::Array{TranslationParameters,1}

    # we need to check - which host_type are we?
    if host_type == :bacteria
        return precompute_translation_parameter_array_bacteria(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})
    elseif host_type == :mammalian
        return precompute_translation_parameter_array_mammalian(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})
    elseif host_type == :cell_free
        return precompute_translation_parameter_array_cell_free(biophysical_dictionary::Dict{String,Any}, protein_coding_length_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})
    else
        throw(error("unsupported host_type value. Expected {:bacteria,:mammalian,:cell_free}. Got $(host_type)"))
    end
end

function precompute_transcription_kinetic_limit_array_cell_free(biophysical_dictionary::Dict{String,Any}, gene_coding_length_array::Array{Float64,1}, gene_abundance_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})::Array{Float64,1}

    # get stuff from the biophysical dictionary -
    characteristic_length = biophysical_dictionary["characteristic_transcript_length"]
    transcription_elongation_rate = biophysical_dictionary["transcription_elongation_rate"]
    characteristic_initiation_time = biophysical_dictionary["characteristic_initiation_time_transcription"]
    KX = biophysical_dictionary["transcription_saturation_constant"]
    RNAPII_concentration = biophysical_dictionary["RNAPII_concentration"]
    species_symbol_type_array = biophysical_dictionary["species_symbol_type_array"]

    # initialize -
    transcription_kinetics_array = Float64[]

    # which indices are mRNA (the production of transcription)?
    idx_mrna = findall(x->x==:mrna,species_symbol_type_array)

    # compute the kinetic limit -
    for (gene_index, gene_length) in enumerate(gene_coding_length_array)

        # how much gene do we have?
        # for host_type == :cell_free this is in concentration units -
        G = gene_abundance_array[gene_index]

        # what is the correction index?
        correction_index = idx_mrna[gene_index]
        time_scale_parameter = time_constant_modifier_array[correction_index]

        # compute kE -
        kE = transcription_elongation_rate*(1/gene_length)

        # compute kI -
        kI = (1/characteristic_initiation_time)

        # compute the tau factor -
        tau_factor = (kE/kI)*time_scale_parameter

        # Compute the rate -
        sat_term = (G/(KX*tau_factor+(1+tau_factor)*G))
        value = kE*RNAPII_concentration*sat_term*(3600) # nmol/gDW-hr

        # push -
        push!(transcription_kinetics_array, value)
    end

    # return -
    return transcription_kinetics_array
end

function precompute_transcription_kinetic_limit_array_bacteria(biophysical_dictionary::Dict{String,Any}, gene_coding_length_array::Array{Float64,1}, gene_copy_number_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})::Array{Float64,1}

    # get stuff from the biophysical dictionary -
    RNAPII_copy_number = biophysical_dictionary["copies_of_rnapII_per_cell"]
    characteristic_length = biophysical_dictionary["characteristic_transcript_length"]
    transcription_elongation_rate = biophysical_dictionary["transcription_elongation_rate"]
    characteristic_initiation_time = biophysical_dictionary["characteristic_initiation_time_transcription"]
    KX = biophysical_dictionary["transcription_saturation_constant"]
    mass_of_single_cell = biophysical_dictionary["mass_of_single_cell"]
    fraction_of_water_per_cell = biophysical_dictionary["fraction_of_water_per_cell"]
    species_symbol_type_array = biophysical_dictionary["species_symbol_type_array"]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the RNAP concentration -
    RNAPII_concentration = RNAPII_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # which indices are mRNA (the production of transcription)?
    idx_mrna = findall(x->x==:mrna,species_symbol_type_array)

    # initialize -
    transcription_kinetics_array = Float64[]

    # compute the kinetic limit -
    for (gene_index, gene_length) in enumerate(gene_coding_length_array)

        # compute elongation constant -
        gene_copy_number = gene_copy_number_array[gene_index]

        # compute kE -
        kE = transcription_elongation_rate*(1/gene_length)

        # compute kI -
        kI = (1/characteristic_initiation_time)

        # what is the correction index?
        correction_index = idx_mrna[gene_index]
        time_scale_parameter = time_constant_modifier_array[correction_index]

        # compute the tau factor -
        tau_factor = (kE/kI)*time_scale_parameter

        # compute the gene concentration -
        G = gene_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell)   # nmol/gDW

        # Compute the rate -
        value = kE*RNAPII_concentration*(G/(KX*tau_factor+(1+tau_factor)*G))*(3600) # nmol/gDW-hr

        # push -
        push!(transcription_kinetics_array, value)
    end

    # return -
    return transcription_kinetics_array
end

function precompute_transcription_kinetic_limit_array_mammalian(biophysical_dictionary::Dict{String,Any}, gene_coding_length_array::Array{Float64,1}, gene_copy_number_array::Array{Float64,1},  time_constant_modifier_array::Array{Float64,1})::Array{Float64,1}

    # get stuff from the biophysical dictionary -
    RNAPII_copy_number = biophysical_dictionary["copies_of_rnapII_per_cell"]
    characteristic_length = biophysical_dictionary["characteristic_transcript_length"]
    transcription_elongation_rate = biophysical_dictionary["transcription_elongation_rate"]
    characteristic_initiation_time = biophysical_dictionary["characteristic_initiation_time_transcription"]
    KX = biophysical_dictionary["transcription_saturation_constant"]
    mass_of_single_cell = biophysical_dictionary["mass_of_single_cell"]
    fraction_of_water_per_cell = biophysical_dictionary["fraction_of_water_per_cell"]
    species_symbol_type_array = biophysical_dictionary["species_symbol_type_array"]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the RNAP concentration -
    RNAPII_concentration = RNAPII_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # initialize -
    transcription_kinetics_array = Float64[]

    # which indices are mRNA (the production of transcription)?
    idx_mrna = findall(x->x==:mrna,species_symbol_type_array)

    # compute the kinetic limit -
    for (gene_index, gene_length) in enumerate(gene_coding_length_array)

        # compute elongation constant -
        gene_copy_number = gene_copy_number_array[gene_index]

        # compute kE -
        kE = transcription_elongation_rate*(1/gene_length)

        # compute kI -
        kI = (1/characteristic_initiation_time)

        # what is the correction index?
        correction_index = idx_mrna[gene_index]
        time_scale_parameter = time_constant_modifier_array[correction_index]

        # compute the tau factor -
        tau_factor = (kE/kI)*time_scale_parameter

        # compute the gene concentration -
        G = gene_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell)   # nmol/gDW

        # Compute the rate -
        value = kE*RNAPII_concentration*(G/(KX*tau_factor+(1+tau_factor)*G))*(3600) # nmol/gDW-hr

        # push -
        push!(transcription_kinetics_array, value)
    end

    # return -
    return transcription_kinetics_array
end

function precompute_transcription_kinetic_limit_array(biophysical_dictionary::Dict{String,Any}, gene_coding_length_array::Array{Float64,1}, gene_abundance_array::Array{Float64,1}, time_constant_modifier_array::Array{Float64,1}, host_type::Symbol)::Array{Float64,1}

    # swicth on the host_type flag -
    if host_type == :bacteria
        return precompute_transcription_kinetic_limit_array_bacteria(biophysical_dictionary,gene_coding_length_array,gene_abundance_array, time_constant_modifier_array)
    elseif host_type == :mammalian
        return precompute_transcription_kinetic_limit_array_mammalian(biophysical_dictionary,gene_coding_length_array,gene_abundance_array, time_constant_modifier_array)
    elseif host_type == :cell_free
        return precompute_transcription_kinetic_limit_array_cell_free(biophysical_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array)
    else
        throw(error("unsupported host_type value. Expected {:bacteria,:mammalian,:cell_free}. Got $(host_type)"))
    end
end

function calculate_translation_kinetics_array(t,x,mRNA_index_array,translation_parameters_array, data_dictionary)

    # we need to calculate the translation rate on the fly
    # since it depends upon mRNA -
    translation_kinetics_array = Float64[]

    # grab the translation parameters from the data dictionary -
    number_of_translation_rates = length(translation_parameters_array)
    counter = 1
    for index = 1:number_of_translation_rates

        # use the mRNA_index_array to index into the state vector -
        mRNA_concentration = x[mRNA_index_array[index]]

        # grab the parameters struct -
        parameters_struct = translation_parameters_array[index]
        vmax = parameters_struct.vmax
        tau = parameters_struct.tau_factor
        KL = parameters_struct.KL

        # build the rate -
        value = vmax*(mRNA_concentration/(KL*tau+(1+tau)*mRNA_concentration))

        # package -
        push!(translation_kinetics_array, value)
    end

    # return -
    return translation_kinetics_array
end

function build_dilution_degradation_matrix(biophysical_dictionary::Dict{String,Any}, species_type_array::Array{Symbol,1}, degradation_modifier_array::Array{Float64,1})::Array{Float64,2}

    # how many species do we have?
    total_number_of_species = length(species_type_array)

    # estimate *max* growth -
    # get the cell doubling time -
    cell_doubling_time = biophysical_dictionary["cell_doubling_time"]
    mugmax = 0.0
    if cell_doubling_time != 0.0

        # compute the mumax -
        mugmax = log(2)/cell_doubling_time
    end

    # calculate the mRNA and protein half life -
    mRNA_half_life_in_hr = biophysical_dictionary["mRNA_half_life_in_hr"]
    protein_half_life_in_hr = biophysical_dictionary["protein_half_life_in_hr"]

    # compute the degradation rate constant -
    kDX = -(1/mRNA_half_life_in_hr)*log(0.5)
    kDL = -(1/protein_half_life_in_hr)*log(0.5)

    # initialize array which holds the (mu+d) terms -
    dilution_degradation_term_array = Array{Float64,1}()
    for (species_index,species_type) in enumerate(species_type_array)

        # for species, compute a degradation term -
        if species_type == :gene || species_type == :constant
            push!(dilution_degradation_term_array,0.0)
        elseif species_type == :mrna

            # get degradation modifier -
            degradation_modifier_value = degradation_modifier_array[species_index]

            # compute -
            value = -1*(mugmax + degradation_modifier_value*kDX)

            # grab -
            push!(dilution_degradation_term_array,value)

        elseif species_type == :protein

            # get degradation modifier -
            degradation_modifier_value = degradation_modifier_array[species_index]

            # compute -
            value = -1*(mugmax + degradation_modifier_value*kDL)

            # grab -
            push!(dilution_degradation_term_array,value)
        end
    end

    # return the dilution array -
    return Diagonal(dilution_degradation_term_array)
end
# -- INTERNAL METHODS ABOVE ---------------------------------------------------------- #

# -- PUBLIC METHODS BELOW ------------------------------------------------------------ #
function calculate_txtl_kinetics_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})::Array{Float64,1}

    # initialize -
    calculated_txtl_kinetics_array = Array{Float64,1}()

    # we have precomputed the translation kinetic parameters, and the kinetic limit of the transcription rate -
    transcription_kinetic_limit_array = data_dictionary["transcription_kinetic_limit_array"]
    translation_parameters_array = data_dictionary["translation_parameter_array"]
    species_symbol_type_array = data_dictionary["species_symbol_type_array"]

    # find mRNA and protein indices -
    mRNA_index_array = findall(x->x==:mrna,species_symbol_type_array)
    protein_index_array = findall(x->x==:protein,species_symbol_type_array)

    # compute the rate of translation -
    translation_kinetics_array = calculate_translation_kinetics_array(t,x,mRNA_index_array,translation_parameters_array, data_dictionary)

    # package -
    for (index,value) in enumerate(transcription_kinetic_limit_array)
        push!(calculated_txtl_kinetics_array,value)
    end

    for (index,value) in enumerate(translation_kinetics_array)
        push!(calculated_txtl_kinetics_array,value)
    end

    # return - TX rates, then TL rates -
    return calculated_txtl_kinetics_array
end
# -- PUBLIC METHODS ABOVE ------------------------------------------------------------ #
