using PyPlot
using CSV
using DataFrames
using Statistics
using DelimitedFiles

#===================Find mean and CI for Ensemble==================#
t_sim = collect(0.0:0.1:16)

S1 = readdlm("simulations/simulation-P1.dat")
no_time, no_species = size(S1)
no_samples = 700
species_ensemble = zeros(no_time,no_species,no_samples)
for i = 1:no_samples
    species_ensemble[:,:,i] = readdlm("/simulations/simulation-P$i.dat")
end
species_mean = mean(species_ensemble, dims = 3)
species_err = std(species_ensemble, dims = 3)
species_pos = species_mean + 1.96*species_err
species_neg = species_mean - 1.96*species_err

t = [0;2;4;6;8;16]

markercolor = "darkred"
shade = "lightpink"
lcolor = "darkred"

#0 = no
#1 = yes
plotline = 1
plotsave =  0
filedir = "Figs/all"
savecase = "all"
#--------------------mRNA & GFP----------------------------#
gfp = CSV.read("data/Data3.csv")
mean_gfp = gfp[1:5,:]
error_gfp = gfp[6:10,:]
gfp_sim_idx = [3;5] #species index for mRNA and GFP

row = 1
col = 2

figure(1,figsize=(6.5,2.5))
for idx = 1:2
    subplot(row,col,idx)
    ylabel(names(gfp)[idx])
    errorbar(t, mean_gfp[:,idx], yerr=error_gfp[:,idx],fmt="o",markersize = 3,capsize=2,elinewidth=1,color=markercolor)
    if idx == 1
        if plotline == 1
             plot(t_sim,species_mean[:,gfp_sim_idx[idx]]*1e6,color=lcolor)
             axis([0,16.0,0,30])
             yticks([0,300,600,900])

        end
        fill_between(t_sim,species_pos[:,gfp_sim_idx[idx]]*1e6,species_neg[:,gfp_sim_idx[idx]]*1e6,color=shade,alpha=0.6,linewidth=0)


    else
        if plotline == 1
            plot(t_sim,species_mean[:,gfp_sim_idx[idx]]*1e3,color=lcolor)
            axis([0,16.0,0,30])
            xticks([0,4,8,12,16])
        end
        fill_between(t_sim,species_pos[:,gfp_sim_idx[idx]]*1e3,species_neg[:,gfp_sim_idx[idx]]*1e3,color=shade,alpha=0.6,linewidth=0)
    end
    # axis([0,16.0,0,30])
    # # xticks([0,4,8,12,16])
    # yticks([0,10,20,30])
    # if names(gfp)[idx] == Symbol("mRNA (nM)")
    #     axis([0,16.0,0,1000])
    #     yticks([0,300,600,900])
    # end
    xlabel("Time (h)")
end
tight_layout()
if plotsave == 1
    savefig("$(filedir)/$(savecase)_protein.pdf")
    close()
end
