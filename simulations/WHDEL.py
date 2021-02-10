#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np 
import pandas as pd

import multiprocessing as mp 
import matplotlib as mpl 

import matplotlib.pyplot as plt
import sys 
import os

from scipy.integrate import odeint
from scipy.sparse import csr_matrix, lil_matrix, dok_matrix, vstack, issparse
from scipy.stats import binom, poisson, multinomial


# In[ ]:


#parameters 
# 1 generation = entry to shedding  
hours_per_gen = 6 # McCrone, bioRxiv 2020
R0 = 5 # Within-host R0 (Hadjichrysanthou et al., 2016)
C_max = 4e8 # max number of cells that can be targeted (Hadjichrysanthouet al., 2016;)
d = 2*(hours_per_gen/24) # per-capita virus decay rate per generation (clearance)

# variables 
bottleneck = 1 # starting number of (V)irions 
r = 100 # branching factor (average number of virions produced by each cell)
days_to_analyze = 14

# computed variables 
days_to_analyze = 14
hours_per_gen = 6
T = int(days_to_analyze*24/hours_per_gen) # number of replication cycles  
t_vec = np.arange(T+1)
beta = (R0*d)/(C_max*r)


# In[ ]:


# Find number of viruses within host using within-host ODE
def within_host_ODE(x, t): 
    c, v = x 
    nc = beta*c*v
    dc = -nc
    dv = nc*r - d*v
    return dc, dv

# initialize 
x_0 = C_max, bottleneck 
def solve_path(t_vec, x_init=x_0):
    F = lambda x, t: within_host_ODE(x, t)
    c_path, v_path = odeint(F, x_init, t_vec).transpose() 
    return c_path, v_path

c_path, v_path = solve_path(t_vec)
v_path = np.rint(v_path).astype('int64')


# In[ ]:


def binomial_robust(NN, p, size=None):
    NNp = NN*p
    if NN > 20 and NNp < 5:
        return np.random.poisson(NNp, size)
    else:
        return np.random.binomial(NN, p, size)

def multinomial_robust(NN, p, size=None):
    if NN < 1000: 
        return multinomial.rvs(NN, p)
    else: 
        results = np.array([binomial_robust(NN, pi, size) for pi in p])
        last_entry = int(NN) - results[:-1].sum(0)
        while last_entry < 0: 
            results = np.array([binomial_robust(NN, pi, size) for pi in p])
            last_entry = int(NN) - results[:-1].sum(0)
        return np.rollaxis(results, 0, results.ndim)

def generate_randint(N, high_lim_arr, low_lim=0): 
    #print (high_lim_arr)
    if (high_lim_arr <= 0).any(): 
        print (high_lim_arr)
        raise Exception
        
    result = np.random.randint(low=0, high=high_lim_arr)
    result_sum = result.sum()
    exit_counter = 0
    while result_sum == 0: 
        result = np.random.randint(low=0, high=high_lim_arr)
        result_sum = result.sum()
        exit_counter += 1
        if exit_counter > 10: 
            break 
    
    if result_sum > 0: 
        result = np.rint(N*(result/result_sum)).astype(np.int)
    
    xs = result.sum()-N
    #print (xs, '*')
    while (abs(xs) != 0): 
        if abs(xs) < 5: 
            rand_add = xs
        else: 
            if xs < 0: 
                rand_add = np.random.randint(low=xs, high=0)
            else: 
                rand_add = np.random.randint(low=0, high=xs)
        
        result_idx = np.random.choice(np.arange(len(result)), 100 if len(result) >= 100 else np.int(0.9*len(result)), replace=False)
        np.random.shuffle(result_idx)
        for idx in result_idx: 
            value = result[idx]
            if value - rand_add < high_lim_arr[idx] and value - rand_add >= 0: 
                #print (value, rand_add)
                result[idx] = value - rand_add
                break 
                
        xs = result.sum()-N
        #print (xs, "**t")
    return result

def remove_extinct_genotypes(g_dist, g_boolean):
    # remove extinct genotype 
    extinct_g = np.where(g_dist<1)[0]
    if len(extinct_g) > 0: 
        g_dist = g_dist[g_dist>0]
        mask = np.ones(g_boolean.shape[0], dtype=bool)
        mask[extinct_g] = False
        g_boolean = g_boolean[mask]
        return g_dist, g_boolean
    else: 
        return g_dist, g_boolean
    
def combine_identical_genotypes(g_dist, g_boolean):

    indices = g_boolean.indices
    indptr = g_boolean.indptr
    
    mpos_to_g = [tuple(indices[indptr[g]:indptr[g+1]]) for g in np.arange(g_dist.size)]
    #print ("step 1 done, {}".format(len(mpos_to_g)))
    if len(mpos_to_g) == len(set(mpos_to_g)): 
        return g_dist, g_boolean
    else: 
        mpos_to_glist = {}
        for g, mpos in enumerate(mpos_to_g): 
            #mpos_count = mpos_to_g.count(mpos)
            #if mpos_count > 1: 
            try: 
                mpos_to_glist[mpos].append(g)
            except:
                mpos_to_glist[mpos] = [g]
        
        discard_g = []
        for mpos, glist in mpos_to_glist.items(): 
            if len(glist) > 1:
                discard_g += glist[1:]
                g_dist[glist[0]] += g_dist[glist[1:]].sum() 
        
        mask = np.ones(g_dist.size, dtype=bool)
        mask[discard_g] = 0
        
        g_dist, g_boolean = g_dist[mask], g_boolean[mask]
        return g_dist, g_boolean

def add_mutations(g_dist, g_boolean, m_v, loci):
    g_mutate = generate_randint(m_v, g_dist) # randomly select individuals from genotype classes (index = genotype)
    g_dist -= g_mutate # remove mutated individuals from original genotype 

    # convert csr g_boolean to lil 
    g_boolean = g_boolean.tolil()
    # get genotype rows to copy 
    abv0_g_boolean = g_boolean[np.where(g_mutate>0)[0]]
    abv0_g_mutate = g_mutate[g_mutate>0]

    #print ("..reconstruct new genotype boolean array and add mutations..")
    new_g_boolean_rows, new_g_boolean_cols = g_boolean.nonzero()
    new_g_boolean_rows = list(new_g_boolean_rows)
    new_g_boolean_cols = list(new_g_boolean_cols)

    new_g_rowcounter = g_boolean.shape[0]
    for g, nind_to_mutate in enumerate(abv0_g_mutate):
        existing_g_cols = list(abv0_g_boolean[g].nonzero()[-1])
        # random select locus to mutate for each individual with genotype
        indpos_to_mutate = np.random.choice(loci, size=nind_to_mutate, replace=True) 
        for i in range(nind_to_mutate):
            # append to rows and cols
            new_pos_to_mutate = indpos_to_mutate[i]
            if new_pos_to_mutate in existing_g_cols: # back mutation 
                new_g_boolean_cols += list(set(existing_g_cols)-set([new_pos_to_mutate]))
                new_g_boolean_rows += [new_g_rowcounter+i]*(len(existing_g_cols)-1)
            else: # forward mutation
                new_g_boolean_cols += existing_g_cols + [new_pos_to_mutate]
                new_g_boolean_rows += [new_g_rowcounter+i]*(len(existing_g_cols)+1)
        # update new_g_rowcounter
        new_g_rowcounter += nind_to_mutate

    # print (len(new_g_boolean_rows) == len(new_g_boolean_cols))
    # make new boolean array 
    m, n = g_boolean.shape[0]+abv0_g_mutate.sum(), g_boolean.shape[-1]
    g_boolean = csr_matrix(([1]*len(new_g_boolean_rows), (new_g_boolean_rows, new_g_boolean_cols)), 
                           shape=(m, n), dtype=np.int8)

    #print ("..reconstruct genotype distribution array..")
    # concatenate old and new g dist and boolean 
    g_dist = np.hstack([g_dist, np.ones(abv0_g_mutate.sum(), dtype=np.int)]) # concatenate 

    #print ("..remove extinct genotypes..")
    # remove any extinct genotypes 
    g_dist, g_boolean = remove_extinct_genotypes(g_dist, g_boolean)
        
    return g_dist, g_boolean


# In[ ]:


def simulate(sim, abs_del_s, frac_neu_nonsyn, seed, syn_l=390, nonsyn_l=420, mu=3e-6, exp_s_dist=0, s_mode=1, N_arr=v_path, max_N=1e7, t_vec=t_vec, N0=bottleneck, sample_Nmax=5e3, verbose=0):
    
    np.random.seed(seed)
    
    L = syn_l + nonsyn_l # total length
    
    # make loci map 
    if verbose > 0: 
        print ("initialise monoclonal population..")
    all_loci = np.arange(L)
    syn_loci = np.random.choice(all_loci, syn_l, replace=False)
    nonsyn_loci = np.setdiff1d(all_loci, syn_loci)
    
    neu_nonsyn_loci = np.random.choice(nonsyn_loci, np.round(frac_neu_nonsyn*nonsyn_l).astype(np.int), replace=False)
    del_nonsyn_loci = np.setdiff1d(nonsyn_loci, neu_nonsyn_loci)
    print (len(neu_nonsyn_loci), len(del_nonsyn_loci))
    
    # make fitness array of mutant 
    mut_fitness_arr = np.zeros(L) # first index = 0 
    if exp_s_dist > 0:
        mut_fitness_arr[del_nonsyn_loci] = s_mode*np.random.exponential(abs_del_s, size=len(del_nonsyn_loci)) # del_s
    else: 
        mut_fitness_arr[del_nonsyn_loci] = s_mode*abs_del_s
    
    # initialise monoclonal population 
    rate_arr = np.zeros((2, len(t_vec)-2))
    curr_g_dist = np.zeros(1) + N0 # distribution of genotypes
    curr_g_boolean = csr_matrix((N0, L), dtype=np.int8) # genotype boolean array (if < 0 = all WT)
    
    # expand monoclonal population 
    N = max_N if N_arr[1] > max_N else N_arr[1]  # hard upper limit on virion pop size
    # get fitness distribution of previous gen 
    curr_f_dist = np.multiply(np.exp(1+curr_g_boolean.dot(mut_fitness_arr)), curr_g_dist)
    pop_f = curr_f_dist.sum() # fitness of population
    curr_f_dist = curr_f_dist/pop_f # normalize 
    curr_g_dist = multinomial_robust(N, curr_f_dist) # resample based on fitness distribution to get new curr_g_dist
    
    for t in t_vec[1:-1]: 
            
        # add synonymous mutations 
        m_v = np.random.poisson(N*mu*syn_l) 
        if m_v > 0:
            if verbose > 0:
                print ("add synonymous mutations..%i"%(m_v))
            curr_g_dist, curr_g_boolean = add_mutations(curr_g_dist, curr_g_boolean, m_v, syn_loci)
            if verbose > 0:
                print ("combine identical genotypes..")
            curr_g_dist, curr_g_boolean = combine_identical_genotypes(curr_g_dist, curr_g_boolean)
        
        # add nonsynonymous mutations 
        if len(neu_nonsyn_loci) > 0:
            m_v = np.random.poisson(N*mu*len(neu_nonsyn_loci)) # randomly select m_v number of individuals to mutate 
            if m_v > 0:
                if verbose > 0:
                    print ("add neutral nonsynonymous mutations..%i"%(m_v))
                curr_g_dist, curr_g_boolean = add_mutations(curr_g_dist, curr_g_boolean, m_v, neu_nonsyn_loci)
                
                if verbose > 0:
                    print ("combine identical genotypes..")
                curr_g_dist, curr_g_boolean = combine_identical_genotypes(curr_g_dist, curr_g_boolean)
        
        if len(del_nonsyn_loci) > 0:
            m_v = np.random.poisson(N*mu*len(del_nonsyn_loci)) 
            if m_v > 0:
                if verbose > 0:
                    print ("add deleterious nonsynonymous mutations..%i"%(m_v)) 
                curr_g_dist, curr_g_boolean = add_mutations(curr_g_dist, curr_g_boolean, m_v, del_nonsyn_loci)
            
                if verbose > 0:
                    print ("combine identical genotypes..")
                curr_g_dist, curr_g_boolean = combine_identical_genotypes(curr_g_dist, curr_g_boolean)
        
        # resampling 
        N = max_N if N_arr[t+1] > max_N else N_arr[t+1]  # hard upper limit on virion pop size
        if verbose > 0:
            print ("{} ({:.2f} days), {:.2e} virions".format(t, (t*hours_per_gen)/24, N), "resampling..")
        # get fitness distribution of previous gen 
        curr_f_dist = np.multiply(np.exp(1+curr_g_boolean.dot(mut_fitness_arr)), curr_g_dist)
        pop_f = curr_f_dist.sum() # fitness of population
        curr_f_dist = curr_f_dist/pop_f # normalize 
        # resample based on fitness distribution to get new curr_g_dist
        curr_g_dist = multinomial_robust(N, curr_f_dist)
        if verbose > 0:
            print ("remove extinct genotypes..")
        # remove extinct genotypes 
        curr_g_dist, curr_g_boolean = remove_extinct_genotypes(curr_g_dist, curr_g_boolean)
        
        # compute evo rates 
        sample_g_dist = multinomial_robust(sample_Nmax, curr_g_dist/curr_g_dist.sum())
        syn_r = curr_g_boolean[:,syn_loci].T.dot(sample_g_dist)/sample_g_dist.sum()
        nonsyn_r = curr_g_boolean[:,nonsyn_loci].T.dot(sample_g_dist)/sample_g_dist.sum()

        syn_r = syn_r.sum()/L/(1+((t+1)*hours_per_gen)/24)
        nonsyn_r = nonsyn_r.sum()/L/(1+((t+1)*hours_per_gen)/24)
        
        if ((t+1)*hours_per_gen)/24%1 == 0:
            print (sim, "{:.2f} days".format(((t+1)*hours_per_gen)/24), "{:.4e}".format(syn_r), "{:.4e}".format(nonsyn_r), "{:.2f}".format(nonsyn_r/syn_r))
        rate_arr[0,t-1] = syn_r
        rate_arr[1,t-1] = nonsyn_r
    return rate_arr             
#simulate(sim=0, abs_del_s=0.01, frac_neu_nonsyn=0., s_mode=-1, exp_s_dist=0, verbose=0)


# In[ ]:


abs_del_s = float(sys.argv[1])
frac_neu_nonsyn = float(sys.argv[2])
sim_N = int(sys.argv[3])
threadnum = int(sys.argv[4])


# In[ ]:

if os.path.isdir("./sim_results"): 
    os.mkdir("./sim_results")

all_rate_arr = np.zeros((sim_N, 2, len(t_vec[1:])-1))

numpy_seeds = np.random.choice(np.arange(sim_N*100), sim_N, replace=False)

pool = mp.Pool(processes=threadnum)
results = [pool.apply_async(simulate, args=(sim, abs_del_s, frac_neu_nonsyn, numpy_seeds[sim],)) for sim in range(sim_N)]
for sim, p in enumerate(results): 
    all_rate_arr[sim] = p.get()
pool.close()

np.savez("./sim_results/mean_pos_s{}_frac_neu_nonsyn{}.npz".format(abs_del_s, frac_neu_nonsyn), all_rate_arr=all_rate_arr)

