include("Include.jl")

function plot_mRNA_MDH_data(time_array, simulation_data_array)

    figure(figsize=(4,3.5))

    μ = mean(simulation_data_array,dims=2)
    σ = std(simulation_data_array,dims=2)
    LB = μ .- (1.96/sqrt(1))*σ
    UB = μ .+ (1.96/sqrt(1))*σ
    fill_between(time_array[:,1], vec(UB), vec(LB), color="powderblue", alpha=0.60)

    # Plot mean -
    plot(time_array[:,1],μ,"-",color="black",lw=2.0)

    # labels -
    xlabel("Time (hr)", fontsize=11)
    ylabel("MDH mRNA concentration (nM)", fontsize=11)

    axis([-0.5,16.5,-50,1600])
    xticks([0,4,8,12,16], fontsize=11)
    yticks([0,500,1000,1500], fontsize=11)

    tight_layout()

    # wwrite -
    savefig("./plots/mRNA-MDH-Ensemble.pdf")
end


function main(path_to_simulation_dir::String, path_to_plot_file::String)

    # what index is prot MDH?
    state_index = 3

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
            data_array[step_index,file_index] = sim_data_array[step_index,(state_index+1)].*10e3
            time_array[step_index,file_index] = sim_data_array[step_index,1]
        end
    end

    # plot -
    plot_mRNA_MDH_data(time_array, data_array)
end

path_to_simulation_dir = "./simulations"
path_to_plot_file = "./plots"
main(path_to_simulation_dir, path_to_plot_file)
