include("Include.jl")

function plot_prot_MDH_data(time_array, simulation_data_array, experimental_data_dictionary)
    #
    # (NR,NC) = size(time_array)
    # for index = 1:NC
    #     plot(time_array[:,index],simulation_data_array[:,index],color="tan",alpha=0.80,lw=0.5)
    # end

    figure(figsize=(4,3.5))

    μ = mean(simulation_data_array,dims=2)
    σ = std(simulation_data_array,dims=2)
    LB = μ .- (1.96/sqrt(1))*σ
    UB = μ .+ (1.96/sqrt(1))*σ
    fill_between(time_array[:,1], vec(UB), vec(LB), color="powderblue", alpha=0.60, label="95% CI of Ensemble")

    # Plot mean -
    plot(time_array[:,1],μ,"-",color="black",lw=2.0, label="Mean of Ensemble")

    # plot the experimemtal data -
    TEXP = experimental_data_dictionary["prot_data_array"][:,1]
    DATA = experimental_data_dictionary["prot_data_array"][:,2]
    STD = (1.96/sqrt(3))*experimental_data_dictionary["prot_data_array"][:,3]
    yerr_array = transpose([STD STD])
    errorbar(TEXP, DATA,yerr=yerr_array,fmt="o",color="black", markersize = 5, capsize=5, elinewidth=2, label="Experimental Data")

    # labels -
    xlabel("Time (hr)", fontsize=11)
    ylabel(L"MDH Protein Concentration ($\mu$M)", fontsize=11)

    axis([-0.5,16.5,-1,30])
    xticks([0,4,8,12,16], fontsize=11)
    yticks([0,10,20,30], fontsize=11)

    tight_layout()
    legend(loc=2, frameon=false, fontsize=9)

    # wwrite -
    savefig("./plots/Prot-MDH-Ensemble.pdf")
end


function main(path_to_simulation_dir::String, path_to_plot_file::String)

    # what index is prot MDH?
    state_index = 5

    # load the experimemtal data -
    exp_dd = load_experimental_data_dictionary("./data")

    # how many files?
    file_name_array = searchdir(path_to_simulation_dir, ".dat")
    number_of_trials = length(file_name_array)

    # initialize -> hardcode the dimension for now
    data_array = zeros(1601,number_of_trials)
    time_array = zeros(1601,number_of_trials)

    # read the simulation dir -
    for (file_index,file_name) in enumerate(file_name_array)

        # load -
        file_path = "$(path_to_simulation_dir)/$(file_name)"
        sim_data_array = readdlm(file_path)

        # what size?
        (NR,NC) = size(sim_data_array)
        for step_index = 1:NR
            data_array[step_index,file_index] = sim_data_array[step_index,(state_index+1)]
            time_array[step_index,file_index] = sim_data_array[step_index,1]
        end
    end

    # plot -
    plot_prot_MDH_data(time_array, data_array, exp_dd)
end

path_to_simulation_dir = "./simulations"
path_to_plot_file = "./plots"
main(path_to_simulation_dir, path_to_plot_file)
