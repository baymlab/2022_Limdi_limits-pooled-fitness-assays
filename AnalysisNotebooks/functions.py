#!/usr/bin/env python3

### Script with a lot of commonly used functions for simulating and reanalyzing fitness assay data

import numpy as np
from scipy import stats
import pandas as pd

def search_gene(locations,ta_sites,i):
    """
    Input: 
    locations: gene coordinates, start, end and orientation
    ta_sites: list of TA sites
    i: query gene
    
    Returns:
    mask corresponding to the interior of the 'i'th gene in the E. coli genome
    """
    start = locations[i, 0]
    end = locations[i, 1]
    length = end - start
    #if the gene is on the forward strand
    if locations[i,2]==1:
        search_area = (ta_sites > start+length*frac5p)&(ta_sites < end - length*frac3p)
    #if the gene is on the reverse strand
    elif locations[i,2]==-1:
        search_area = (ta_sites < start+length*frac5p)&(ta_sites > end - length*frac3p)
    return search_area

#calculating fitnesses
def fitness_estimate(counts_red, counts_green, locs, genes_lost, initial_depth, min_sites, min_frac, max_reads, t_start, t_end):
    """
    data:
    - counts_red, counts_green: trajectories for both replicates
    - exclude: filtering out genes that are lost/have too few sites from analysis as they're indicated as -1
    params:
    - locs: start and end points of genes
    - initial_depth: exclude sites below a certain threshold in fitness estimation
    - min_sites: minimum number of sites required for estimating fitness for each replicate
    - min_frac: fraction of TA sites that need to considered in estimating fitness
    - site_thresh: minimum number of reads for a trajectory to be included (strictly greater than)
    - max_reads: maximum reads that can be used for inverse variance estimate
    
    """
    green = np.sum(counts_green,axis=1)
    red = np.sum(counts_red, axis=1)
    #fitness effect (best estimate)
    fitness_inverse_var = np.zeros([len(genes_lost),2])-1
    #three different metrics of measurement error
    stderr_fitness_inverse_var = np.zeros([len(genes_lost)])-1
    #defining time interval based on gene essentiality status:
    time_gens = np.linspace(t_start,t_end,t_end-t_start+1)*6.64
    
    #iterating over all genes that are not lost over evolution
    for gene in np.where(genes_lost==0)[0] :
        #TA sites within this gene
        search_area = search_gene(locations=locs, ta_sites=ta_sites, i=gene)
        #weights
        weights_green = []
        weights_red = []
        #fitnesses
        s_green = []
        s_red = []
        
        ta = np.sum(search_area==True)
        #now, iterating over sites with at least tm1_depth reads at the start of the assay:
        for j in np.where((search_area==True))[0]:
            if (counts_green[t_start,j]/green[t_start]*10**7>initial_depth) & (counts_red[t_start,j]/red[t_start]*10**7>initial_depth):
                #extracting the trajectories
                traj_green = counts_green[t_start:t_end+1,j]/green[t_start:t_end+1]*10**7
                traj_red = counts_red[t_start:t_end+1,j]/red[t_start:t_end+1]*10**7
                n_before_g = traj_green[t_start]
                n_before_r = traj_red[t_start]
                n_after_g = traj_green[t_start+1]
                n_after_r = traj_red[t_start+1]
                
                #green replicate first
                if np.min(traj_green):
                    sg = np.polyfit(time_gens, np.log(traj_green), 1)[0]
                    wg = min(inverse_var_weight(n_before_g, n_after_g), inverse_var_weight(max_reads, max_reads))
                    #sotring the data
                    weights_green.append(wg)
                    s_green.append(sg)
                
                #now the red replicate
                if np.min(traj_red):
                    sr = np.polyfit(time_gens, np.log(traj_red), 1)[0]
                    wr = min(inverse_var_weight(n_before_r, n_after_r), inverse_var_weight(max_reads, max_reads))
                    #sotring the data
                    weights_red.append(wr)
                    s_red.append(sr)
        
        #converting to numpy arrays
        s_green = np.array(s_green)
        s_red = np.array(s_red)
        weights_green = np.array(weights_green)
        weights_red = np.array(weights_red)
        
        #pooling over all sites, require minimum number of sites, and fraction of sites represented
        if min(len(s_green),len(s_red))>=min_sites and min(len(s_green),len(s_red))/ta>=min_frac:

            #weighing based on inverse variacne:
            fitness_inverse_var[gene,0] = np.average(s_green, weights=weights_green)
            fitness_inverse_var[gene,1] = np.average(s_red, weights=weights_red)
            
            #inverse variance weighted standard error measurement
            stacked_s = np.hstack([s_green,s_red])
            stacked_w = np.hstack([weights_green,weights_red])
            average = np.average(stacked_s, weights=stacked_w)
            #inverse variance weighted sem
            sem_inv = np.sqrt(np.average((stacked_s-average)**2, weights=stacked_w)/(len(stacked_s)-1))
            #saving the inverse variance wieghted measurements.
            stderr_fitness_inverse_var[gene] = sem_inv
            
    return fitness_inverse_var, stderr_fitness_inverse_var


#calculating fitnesses
def fitness_site(counts_red, counts_green, locs, genes_lost, initial_depth, t_start, t_end):
    #select insertion mutations with at least a minimum number of reads at the start:
    green = np.sum(counts_green,axis=1)
    red = np.sum(counts_red, axis=1)
    #initial depth is above threshold
    condition1 = counts_green[t_start,:]/green[t_start]*10**7>initial_depth
    condition2 = counts_red[t_start,:]/red[t_start]*10**7>initial_depth
    #trajectories do not die out
    condition3 = np.min(counts_green[t_start:t_end+1, :], axis=0) > 0
    condition4 = np.min(counts_red[t_start:t_end+1, :], axis=0) > 0
    nonzero_sites = np.where(condition1 & condition2 & condition3 & condition4)[0]
#     print(nonzero_sites, len(nonzero_sites))
    
    fitnesses = []
    
    #defining time interval based on gene essentiality status:
    time_gens = np.linspace(t_start,t_end,t_end-t_start+1)*6.64
    #iterate over each site and calculate fitness:
    for site in nonzero_sites:
        traj_green = counts_green[:t_end+1, site]
        traj_red = counts_red[:t_end+1, site]
#         print(traj_green)
        fitnesses.append((np.polyfit(time_gens[:t_end+1], np.log(traj_green), 1)[0] + np.polyfit(time_gens[:t_end+1], np.log(traj_red), 1)[0])*0.5)
    
    return fitnesses
    

def downsample(data, scale):
    """
    Inputs: 
    - data: counts matrix for bulk fitness assay
    - scale: scaling factor for downsampling, must be greater than 1.
    
    Process:
    - downsample number of reads mapping to an insertion site as follows (for each time point, here: by row)
    - use np.repeat to get an list with every insertion site repeated N times, where N is the number of mapped reads
    - use np.shuffle to rearrange this list
    - pick the first 1/scale fraction of this list
    - use np.unique to which sites are represented, and how frequently after downsampling.
    
    Output:
    - data_scaled: same shape as data but each row of the matrix downsampled by the scaling factor
    """
    assert scale >= 1, f"downsampling scale factor must be greater than/equal 1"
    
    if scale == 1: #do not downsample the data at all:
        return data
    
    else:
        data_scaled = np.zeros_like(data)

        for t in range(data.shape[0]): #there are 5 timepoints in the data
            #this is the key step in the process, every TA site is repeated as many times as number of reads mapping to it
            explicit_data = np.repeat(np.arange(0,data.shape[1]), data[t,:].astype('int'))
            #this list is then shuffled
            np.random.shuffle(explicit_data)
            #and a subset of this list becomes the new data
            N_ds = int(data.sum(axis=1)[t]/scale)
            #as we shuffled the data, taking the first N_ds reads is equivalent to taking a 1/scale random subset of the data
            downsampled = explicit_data[:N_ds]
            #getting the counts and unique TA sites represented after downsampling
            unique, counts = np.unique(downsampled, return_counts=True)
            data_scaled[t,unique] = counts
        
        return data_scaled
    
def simulate_assay(initial_counts, s_mut, t, gens):
    """
    Input:
    initial_counts: read counts before fitness assay
    s_mut: fitness effect of disrupting gene (or more generally of the mutation)
    t: number of timepoints for assay
    gens: number of generations per timepoint
    
    Output:
    simulates a fitness trajectory assuming that expected number of counts at any given time point is Poisson
    distributed with mean specified by mutant fitness
    """
    sim_traj = np.zeros([initial_counts.shape[0], t])
    #initialize:
    sim_traj[:,0] = initial_counts
    #running the simulation
    for t in range(1, t):
        expected_counts = sim_traj[:,t-1]*np.exp(s_mut*gens)
        #actual counts are drawn from a Poisson distribution
        sim_traj[:, t] = np.random.poisson(expected_counts)
    
    return sim_traj

def fitness_calculator(trajectories, n_gens):
    """
    Calculates fitness from simulated mutant trajectory data
    """
    sites = trajectories.shape[0]
    time = trajectories.shape[1]
    time_range = np.linspace(0, time-1, time)*n_gens
    
    s_site = [np.polyfit(time_range, np.log(trajectories[i, :]), 1)[0] for i in range(sites) if np.min(trajectories[i,:]>0)]
    return np.mean(s_site), np.std(s_site)

def fitness_calculator_replicates(trajectories, n_gens, reps, uncalculated_count=None):
    """
    trajectories: simulated data
    n_gens: number of generations of selection per day
    reps: number of replicates to average over for estimating fitness
    """
    time = trajectories.shape[1]
    time_range = np.linspace(0, time-1, time)*n_gens
    mutations = int(trajectories.shape[0]/reps)
    s_mutation = np.zeros(mutations)-1
    err_mutation = np.zeros(mutations)-1
    
    reps_used = []
    
    for i in range(mutations):
        splice = trajectories[i:i+reps, :]
        s_rep = [np.polyfit(time_range, np.log(splice[k, :]), 1)[0] for k in range(reps) if np.min(splice[k,:]>0)]
        if np.size(s_rep)>1:
            reps_used.append(np.size(s_rep))
            err_mutation[i] = np.std(s_rep)/np.sqrt(len(s_rep)-1)
            s_mutation[i] = np.mean(s_rep)
    
    num_uncalculated=np.sum(s_mutation==-1)
    
    if uncalculated_count==None:
        return np.mean(s_mutation[s_mutation>-1]), np.mean(err_mutation[err_mutation>-1])
    else:
        return np.mean(s_mutation[s_mutation>-1]), np.mean(err_mutation[err_mutation>-1]), num_uncalculated, reps_used

def bounds_estimator(coverage, t_library, t_recovery, t_assay, fitness_range):
    """
    Inputs:
    coverage: what is the average number of counts per gene (in experiments)
    t_library: number of generations of selection while constructing library
    t_recovery: number of generations of recovery from freezer stock
    t_assay: number of generations in the fitness assay
    fitness_range: the range for which we want error bounds in fitness measurements
    
    Output:
    upper and lower bound in fitness for each measurement
    
    """
    
    #t0 is the start of the fitness assay. 
    #Assumption that all mutants exhibit exponential growth  
    coverage_site_t0 = coverage*np.exp(fitness_range*(t_library+t_recovery))
    #After the fitness assay
    coverage_site_t1 = coverage_site_t0*np.exp(fitness_range*t_assay)
    #bounds
    upper = np.log((coverage_site_t1 + coverage_site_t1**0.5)/((coverage_site_t0 - coverage_site_t0**0.5)))/t_assay
    lower = np.log((coverage_site_t1 - coverage_site_t1**0.5)/((coverage_site_t0 + coverage_site_t0**0.5)))/t_assay
    
    #returning bounds
    
    return upper, lower