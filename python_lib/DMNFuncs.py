
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
import re
import os
from scipy.special import loggamma, digamma, polygamma

def get_var_haplotype_df(vcf_df, mapped_reads_df):
    """
    # Haplotype reconstruction inference

    We can reconstruct the distribution of the most parsimonious alleles that make up the virus population of a sample
    using a maximum-likelihood approach ([Ghafari et al., 2020](https://jvi.asm.org/content/early/2020/04/09/JVI.00014-20)).

    We can think of the short reads covering a subset of loci of the genome $l$ as observed partial haplotype counts
    ($X_{l}^{P}$). The true (unobserved) proportions underlying these data is denoted by $q_{l}^{P}$. The proportions
    ($q_{l}^{P}$) of these partial haplotypes ($h_l^P$) can be mapped to the actual haplotypes ($h_i$ where $i$ indexes each
    distinct haplotype) proportions ($q_{i}$) by $q_l^P = T_l\times q_i$. For instance, if the partial haplotypes consists of
    {A-, G-} and the hypothetical true haplotypes are {AT, AC, GT, GC}:

    $
    \begin{bmatrix}
    q_{A-}^{P}\\
    q_{G-}^{P}
    \end{bmatrix}=
    \begin{bmatrix}
    1 & 1 & 0 & 0\\
    0 & 0 & 1 & 1
    \end{bmatrix}
    \begin{bmatrix}
    q_{AT}\\
    q_{AC}\\
    q_{GT}\\
    q_{GC}
    \end{bmatrix}
    $

    Assuming that $X_{l}^{P} \sim DirMultinomial(N_l^P, \pi_l^P, \varphi)$ where $N_l^P$ is the total number of reads
    covering loci subset $l$, we reconstruct the most parsimonious distribution of haplotypes by making guesses of $q_i$
    (and assume some level of overdisperion $\varphi$) and maximise liklihood  of the partial haplotype Diriclet multinomial
    model.
    """

    try:
        vcf_df = vcf_df.reset_index()
    except:
        pass
    vcf_df = vcf_df.drop_duplicates(['nucpos', 'nuc_var'])
    vcf_df = vcf_df.set_index("nucpos")

    try:
        mapped_reads_df = mapped_reads_df.reset_index()
    except:
        pass

    # filter reads with polymorphic position
    readidx_to_nucpos = [{"idx":read_idx, "nucpos":nucpos} for nucpos in vcf_df.index for read_idx in mapped_reads_df[(mapped_reads_df['start_nucaln_pos']<=nucpos)&(mapped_reads_df['end_nucaln_pos_plus1']>nucpos)].index]
    readidx_to_nucpos = pd.DataFrame.from_dict(readidx_to_nucpos).set_index("idx").sort_index()

    # variants must be called in vcf
    nucpos_to_variants  = {}
    for nucpos in set(vcf_df.index):
        if isinstance(vcf_df.loc[nucpos], pd.Series):
            nucpos_to_variants[nucpos] = [vcf_df.loc[nucpos,'nuc_consensus']]
            nucpos_to_variants[nucpos].append(vcf_df.loc[nucpos,'nuc_var'])
        else:
            nucpos_to_variants[nucpos] = [vcf_df.loc[nucpos, 'nuc_consensus'].iloc[0]]
            nucpos_to_variants[nucpos] += list(vcf_df.loc[nucpos,'nuc_var'])

    # generate observed partial haplotype DataFrame
    partial_haplotype_tally = {}
    for read_idx in set(readidx_to_nucpos.index):
        nucpos_covered = readidx_to_nucpos.loc[read_idx, "nucpos"]

        if isinstance(nucpos_covered, pd.Series):
            nucpos_covered = sorted(set(nucpos_covered))
        else:
            nucpos_covered = [nucpos_covered]

        # get alleles of nucpos covered
        read_seq = list(mapped_reads_df.loc[read_idx, "seq"])
        start_nucaln_pos = mapped_reads_df.loc[read_idx, "start_nucaln_pos"]
        alleles = [read_seq[nucpos-start_nucaln_pos] for nucpos in nucpos_covered]

        # filter out gaps and variants that failed quality filters
        remove_nucpos = []
        alleles_covered = []
        for _, nucpos in enumerate(nucpos_covered) :
            if alleles[_] in nucpos_to_variants[nucpos]:
                alleles_covered.append(alleles[_])
            else:
                remove_nucpos.append(nucpos)
        nucpos_covered = sorted(set(nucpos_covered)-set(remove_nucpos))

        if len(nucpos_covered) == 0:
            continue

        try:
            partial_haplotype_tally[tuple(nucpos_covered)][tuple(alleles_covered)] += 1
        except:
            try:
                partial_haplotype_tally[tuple(nucpos_covered)][tuple(alleles_covered)] = 1
            except:
                partial_haplotype_tally[tuple(nucpos_covered)] = {tuple(alleles_covered):1}

    # generate count dataframe
    partial_haplotype_df = []
    for nucpos_covered in sorted(partial_haplotype_tally.keys()):
        for alleles, alleles_count in partial_haplotype_tally[nucpos_covered].items():
            # compute relative frequency for each position set
            partial_haplotype_df.append({"nucpos":nucpos_covered, "var":alleles, "count":alleles_count})

    partial_haplotype_df = pd.DataFrame.from_dict(partial_haplotype_df)

    return partial_haplotype_df

def generate_all_possible_haplotypes(var_haplotype_df):
    # compute variant coverage
    nucpos_to_alleles_to_count = {}
    for r, row in var_haplotype_df.iterrows():
        nucpos_covered = list(row["nucpos"])
        alleles = list(row["var"])

        for _, nucpos in enumerate(nucpos_covered):
            try:
                nucpos_to_alleles_to_count[nucpos][alleles[_]] += row['count']
            except:
                try:
                    nucpos_to_alleles_to_count[nucpos][alleles[_]] = row['count']
                except:
                    nucpos_to_alleles_to_count[nucpos] = {alleles[_]:row['count']}

    # genrate all possible haplotypes
    sorted_nucpos = sorted(set(nucpos_to_alleles_to_count.keys()))
    all_possible_haplotypes = list(itertools.product(*[list(nucpos_to_alleles_to_count[nucpos].keys()) for nucpos in sorted_nucpos]))

    # rank possible haplotypes by expected joint prob computed from single-site frequencies
    haplotype_to_ss_score = {}
    for haplotype in all_possible_haplotypes:
        for nidx, nucpos in enumerate(sorted_nucpos):
            allele = haplotype[nidx]
            try:
                haplotype_to_ss_score[haplotype].append(nucpos_to_alleles_to_count[nucpos][allele]/sum(nucpos_to_alleles_to_count[nucpos].values()))
            except:
                haplotype_to_ss_score[haplotype] = [nucpos_to_alleles_to_count[nucpos][allele]/sum(nucpos_to_alleles_to_count[nucpos].values())]
        haplotype_to_ss_score[haplotype] = np.product(haplotype_to_ss_score[haplotype])

    ss_ranked_haplotypes = sorted(haplotype_to_ss_score.keys(), key=haplotype_to_ss_score.get)[::-1]

    # rank haplotypes by each unique read_type
    # skip if <10 reads for read_type read_type covered
    read_type_list = []
    multi_mutant_boolean = 0
    for nucpos_list in set(var_haplotype_df['nucpos']):
        if  var_haplotype_df[var_haplotype_df['nucpos']==nucpos_list]['count'].sum() < 10:
            continue

        if len(nucpos_list) > 1:
            # check if we have multi-mutants in read
            for var in set(var_haplotype_df[var_haplotype_df['nucpos']==nucpos_list]['var']):
                if sum([1 if nucpos_to_alleles_to_count[nucpos_list[vidx]][v]/sum(nucpos_to_alleles_to_count[nucpos_list[vidx]].values()) < 0.5 else 0 for vidx, v in enumerate(var)]) > 1:
                    multi_mutant_boolean = 1

        read_type_list.append(nucpos_list)

    if multi_mutant_boolean == 0:
        # if only single mutants found, ranked by single site joint prob
        sorted_haplotypes = ss_ranked_haplotypes[:]
    else:
        read_type_hierarchy = {}
        top_read_type = []

        # sort read type by number of sites covered (higer first)
        for nucpos_list in sorted(read_type_list, key=lambda nucpos_list: len(nucpos_list))[::-1]:

            add_boolean = 1
            for prev_nucpos_list in top_read_type:
                # absorb any read_type with positions that were subset to larger-sites (top) read_type
                if set(nucpos_list) <= set(prev_nucpos_list):
                    read_type_hierarchy[prev_nucpos_list].append(nucpos_list)
                    add_boolean = 0
                    break

            # if shorted read type is unique -> becomes a top read_type
            if add_boolean == 1:
                read_type_hierarchy[nucpos_list] = []
                top_read_type.append(nucpos_list)
                continue

        # identifying read type with nonoverlapping sites to all other top read types
        non_overlapping_top_site = []
        for nucpos_list in top_read_type:
            overlap = 0
            for j in top_read_type:
                if j == nucpos_list:
                    continue
                elif len(set(nucpos_list)&(set(j))) > 0:
                    overlap = 1

            if overlap == 0:
                non_overlapping_top_site.append(nucpos_list)

        haplotype_to_ms_score = {}
        haplotype_to_disregard = []
        for top_nucpos_list in top_read_type:
            top_var_haplotype_df = var_haplotype_df[var_haplotype_df['nucpos']==top_nucpos_list]
            topvar_to_count = pd.Series(top_var_haplotype_df['count'].values, index=top_var_haplotype_df['var']).to_dict()
            topvar_to_count = {"".join(k):v for k, v in topvar_to_count.items()}

            # absorb any lower hierarchy alleles
            if len(read_type_hierarchy[top_nucpos_list]) > 0:
                for btm_nucpos_list in sorted(read_type_hierarchy[top_nucpos_list], key=lambda nucpos_list: var_haplotype_df[var_haplotype_df['nucpos']==nucpos_list]['count'].sum(), reverse=True):
                    btm_var_haplotype_df = var_haplotype_df[var_haplotype_df['nucpos']==btm_nucpos_list]

                    for btm_var in btm_var_haplotype_df['var']:
                        btm_var_regex = ''.join([btm_var[btm_nucpos_list.index(top_nucpos)] if top_nucpos in btm_nucpos_list else '[A-Z]' for top_nucpos in top_nucpos_list])

                        matched = []
                        for top_var in topvar_to_count:
                            if re.search(btm_var_regex, top_var):
                                matched.append(top_var)

                        btm_var_count = btm_var_haplotype_df[btm_var_haplotype_df['var']==btm_var]['count'].iloc[0]
                        if len(matched) == 0:
                            topvar_to_count[btm_var_regex] = btm_var_count
                        elif len(matched) == 1:
                            topvar_to_count[matched[0]] += btm_var_count
                        else:
                            matched_total = sum([topvar_to_count[top_var] for top_var in matched])
                            for top_var in matched:
                                topvar_to_count[top_var] += btm_var_count*(topvar_to_count[top_var]/matched_total)

            # compute joint probability
            for haplotype in all_possible_haplotypes:
                hap_var = "".join([haplotype[sorted_nucpos.index(top_nucpos)] for top_nucpos in top_nucpos_list])

                score = 0
                for top_var_regex, count in topvar_to_count.items():
                    if re.search(top_var_regex, hap_var):
                        score = count/sum(topvar_to_count.values())
                        break

                if score == 0 and top_nucpos_list in non_overlapping_top_site:
                    haplotype_to_disregard.append(haplotype)
                    continue

                try:
                    haplotype_to_ms_score[haplotype].append(score)
                except:
                    haplotype_to_ms_score[haplotype] = [score]

        haplotype_to_ms_score = {haplotype:np.prod(score) for haplotype, score in haplotype_to_ms_score.items()}
        ms_ranked_haplotypes = sorted(haplotype_to_ms_score, key=haplotype_to_ms_score.get)[::-1]
        ms_ranked_haplotypes += haplotype_to_disregard

        sorted_haplotypes = {haplotype:np.prod([ms_ranked_haplotypes.index(haplotype), ss_ranked_haplotypes.index(haplotype)]) for haplotype in all_possible_haplotypes if haplotype not in haplotype_to_disregard}
        sorted_haplotypes = sorted(sorted_haplotypes, key=sorted_haplotypes.get)

    # generate transpose matrix
    T_matrix = []

    ss_var_haplotype_df = []
    for nucpos, allele_to_count in nucpos_to_alleles_to_count.items():
        for allele, count in allele_to_count.items():
            ss_var_haplotype_df.append({'nucpos':tuple([nucpos]), 'var':tuple([allele]), 'count':count})
    ss_var_haplotype_df = pd.DataFrame.from_dict(ss_var_haplotype_df)

    for haplotype in sorted_haplotypes:
        T_row = []

        try:
            print (haplotype, haplotype_to_ss_score[haplotype], haplotype_to_ms_score[haplotype])
        except:
            print (haplotype, haplotype_to_ss_score[haplotype])

        """for r, row in var_haplotype_df.iterrows():"""
        for r, row in ss_var_haplotype_df.iterrows():
            if [haplotype[_] for _, nucpos in enumerate(sorted_nucpos) if nucpos in list(row["nucpos"])] != list(row['var']):
                T_row.append(0)
            else:
                T_row.append(1)

        T_matrix.append(list(T_row))

    return sorted_haplotypes, np.array(T_matrix), sorted_nucpos, ss_var_haplotype_df

def initialise_possible_haplotype_frequency(stopidx, X_haplotype_boolean):
    # initialise haplotype frequencies
    haplo_freq = np.random.uniform(size=stopidx)

    if X_haplotype_boolean == 1:
        # unknown haplotype freq (<1%)
        unknown_freq = np.random.uniform()/100
        haplo_freq = (haplo_freq/np.sum(haplo_freq))*(1-unknown_freq) # normalise
        haplo_freq = np.append(haplo_freq, 1.-np.sum(haplo_freq))
    else:
        haplo_freq = haplo_freq/np.sum(haplo_freq) # normalise

    return haplo_freq

def convert_possible_haplo_freq_to_partial_haplo_freq(T, haplo_freq):
    partial_PI = np.dot(T.T, haplo_freq)
    partial_PI = partial_PI/partial_PI.sum() # normalise
    return partial_PI

def ComputeLL(X_df, PI, C):
    ALPHA = PI*C
    ll = 0
    for nucpos in set(X_df['nucpos']):
        nucpos_X_df = X_df[X_df['nucpos']==nucpos]
        X = np.array(nucpos_X_df["count"])
        alpha_arr = ALPHA[nucpos_X_df.index]
        ll += loggamma(alpha_arr.sum()) - loggamma(X.sum()+alpha_arr.sum()) + (loggamma(X+alpha_arr)-loggamma(alpha_arr)).sum()
    return ll

def sa_proposal(curr_freq, sigma, X_haplotype_boolean):

    if X_haplotype_boolean == 1:
        # propose new X freq
        new_X_freq = curr_freq[-1] + np.random.normal(0, 0.005)
        new_X_freq = 1e-10 if new_X_freq < 0. else new_X_freq
        new_X_freq = np.random.uniform()/100 if new_X_freq > 0.01 else new_X_freq

        # gaussian kernel on previous haplo frequencies
        new_freq = curr_freq[:-1] + np.random.normal(0, sigma, len(curr_freq[:-1]))
        new_freq[new_freq<0.] = 1e-10 # always positive
        new_freq = (new_freq/np.sum(new_freq))*(1-new_X_freq) # normalise

        # combine with new X freq
        new_freq = np.append(new_freq, new_X_freq)
    else:
        # gaussian kernel on previous haplo frequencies
        new_freq = curr_freq + np.random.normal(0, sigma, len(curr_freq))
        new_freq[new_freq<0.] = 1e-10 # always positive
        new_freq = new_freq/np.sum(new_freq) # normalise

    return new_freq

def sa_acceptance(i, new_freq, curr_freq, curr_LL, obsX, curr_T, C):

    # convert full haplotype to partial haplotype frequencies (partial pi)
    partial_PI = convert_possible_haplo_freq_to_partial_haplo_freq(curr_T, new_freq)

    # compute initial LL
    new_LL = ComputeLL(obsX, partial_PI, C)

    # compute acceptance likelihood (probability)
    if new_LL-curr_LL > np.log(np.random.uniform()):
        return new_freq, new_LL
    else:
        # no change
        return curr_freq, curr_LL

def compute_BIC(param_num, data_num, LL):
    return param_num*np.log(data_num) - 2*LL

def optimise_haploset_frequencies(obsX, curr_haploset, X_haplotype_boolean, curr_T, C, maxiter, save_state, searches, sigma, verbose):

    state_df = []
    for _sidx, s in enumerate(np.random.choice(range(10000), searches)):

        np.random.seed(s) # set seed

        # initialise full haplotype frequencies
        curr_freq = initialise_possible_haplotype_frequency(len(curr_haploset), X_haplotype_boolean)

        # convert full haplotype to partial haplotype frequencies (partial pi)
        partial_PI = convert_possible_haplo_freq_to_partial_haplo_freq(curr_T, curr_freq)

        # compute initial LL
        curr_LL = ComputeLL(obsX, partial_PI, C)

        state_df.append({"search":_sidx+1, "state":0, "freq":curr_freq, "LL":curr_LL})

        for i in range(maxiter):
            proposed_freq = sa_proposal(curr_freq, sigma, X_haplotype_boolean)

            curr_freq, curr_LL = sa_acceptance(i, proposed_freq, curr_freq, curr_LL, obsX, curr_T, C)

            if i == 0 or (i+1)%save_state == 0:
                state_df.append({"search":_sidx+1, "state":i+1, "freq":curr_freq, "LL":curr_LL})
                if verbose == 1:
                    if i == 0 or np.random.uniform() > 0.9:
                        np.set_printoptions(precision=5)
                        print ("{:<2d}{:>6d}\t{}\t{:>.5f}".format(_sidx+1, i if i == 0 else i+1, curr_freq, curr_LL))

    state_df = pd.DataFrame.from_dict(state_df).set_index(["search", "state"])

    return state_df

def get_expected_read_frequencies(var_haplotype_df):
    """
    compute empirical read coverage frequencies to correct for parital PI
    """
    nucpos_to_coverage = {}
    nucposlist_to_coverage = {}
    for r, row in var_haplotype_df.iterrows():
        # read coverage
        try:
            nucposlist_to_coverage[row['nucpos']] += row['count']
        except:
            nucposlist_to_coverage[row['nucpos']] = row['count']

        # individual site coverage
        for nucpos in row['nucpos']:
            try:
                nucpos_to_coverage[nucpos] += row['count']
            except:
                nucpos_to_coverage[nucpos] = row['count']

    for r, row in var_haplotype_df.iterrows():
        var_haplotype_df.at[r, 'prob'] = np.prod([nucposlist_to_coverage[row['nucpos']]/nucpos_to_coverage[nucpos] for nucpos in row['nucpos']])

    return var_haplotype_df

def reconstruct_haplotype(var_haplotype_df, log_prefix, maxopt_iter=5000, save_state=25, searches=2, overdispersion=1e-3, sigma=0.05, verbose=1, show_state_plt=1):

    # compute C factor from overdispersion
    C = (1-overdispersion)/overdispersion

    # generate all possible hypothetical haplotypes and transpose matrix
    all_possible_haplotypes, T_matrix, sorted_nucpos, ss_var_haplotype_df = generate_all_possible_haplotypes(var_haplotype_df)

    """
    # get partial prob of read
    var_haplotype_df = get_expected_read_frequencies(var_haplotype_df)
    """

    # for ploting state_df
    fig, ax = plt.subplots(1, 1, figsize=(8.3, 5.8))
    all_states_df = []
    deltabic_abv_zero = 0

    for chs_idx, curr_haploset_stopidx in enumerate(range(1, len(all_possible_haplotypes)+1)):

        # set up curr haplotype set (all consensus)
        curr_haploset = all_possible_haplotypes[:curr_haploset_stopidx]
        curr_T = T_matrix[:curr_haploset_stopidx]
        if (curr_T.sum(axis=1)==0).any():
            break

        # add a further haplotype "X" describing cloud of potential haplotypes not already described by current haploset
        if ((1-curr_T.sum(axis=0))<=0).all() == False:
            X_haplotype_boolean = 1
        else:
            X_haplotype_boolean = 0

        if X_haplotype_boolean == 1:
            # add X haplotype to transpose matrix
            curr_T = np.append(curr_T, [1-curr_T.sum(axis=0)], axis=0)
            curr_T[curr_T<0] = 0

        if verbose == 1:
            print (curr_haploset, "+ X"+str(X_haplotype_boolean))
            print (curr_T)

        """
        # correct transpose matrix by expected read frequencies
        curr_T = curr_T*np.array(var_haplotype_df["prob"])
        curr_T = (curr_T.T/curr_T.sum(axis=1)).T
        """

        # search for optimised haploset frequencies with the maxmimum likelihood of explaining count data
        """
        state_df = optimise_haploset_frequencies(var_haplotype_df[["nucpos", "count"]], curr_haploset, X_haplotype_boolean, curr_T, C, int(np.ceil(maxopt_iter/2)) if chs_idx == 0 else maxopt_iter, save_state, 1 if chs_idx == 0 else searches, sigma, verbose)
        """

        state_df = optimise_haploset_frequencies(ss_var_haplotype_df[["nucpos", "count"]], curr_haploset, X_haplotype_boolean, curr_T, C, int(np.ceil(maxopt_iter/2)) if chs_idx == 0 else maxopt_iter, save_state, 1 if chs_idx == 0 else searches, sigma, verbose)

        # plot ML search trace
        for search in sorted(set(state_df.index.get_level_values(0))):
            search_state_df = state_df.loc[search]
            ax.plot(search_state_df.index, search_state_df['LL'], label="S%i-H%i(X%i)"%(search, curr_haploset_stopidx, X_haplotype_boolean))

        state_df = state_df.reset_index()
        state_df["haplo_len"] = curr_haploset_stopidx
        state_df["X"]= X_haplotype_boolean
        state_df["haplo_set"] = str(curr_haploset)

        # compute BIC based on best LL
        best_state = state_df.loc[state_df['LL'].idxmax()]
        curr_best_LL = best_state["LL"]
        curr_bic = compute_BIC(curr_haploset_stopidx+X_haplotype_boolean, var_haplotype_df['count'].sum(), curr_best_LL)
        state_df['best_BIC'] = curr_bic
        if chs_idx > 0:
            delta_bic = curr_bic-prev_bic
            print ("BIC: %.3f (delta=%.3f)\n"%(curr_bic, delta_bic))

            # break if delta BIC is positive for two consecutive searches
            if delta_bic > 0:
                deltabic_abv_zero += 1
                if deltabic_abv_zero == 2:
                    break
                continue
            else:
                # reset if deltabic_abv_zero was > 0 previously
                if deltabic_abv_zero > 0:
                    deltabic_abv_zero = 0

        else:
            print ("BIC: %.3f\n"%(curr_bic))
        prev_bic = curr_bic

        # append to all_states_df
        all_states_df.append(state_df)

    # concat all state_df
    all_states_df = pd.concat(all_states_df, ignore_index=True)
    all_states_df.to_csv(log_prefix+".csv", index=False)

    # plot ML-search trace
    ax.set_xlabel("States")
    ax.set_ylabel("LL")
    ax.set_title(log_prefix)
    plt.legend(loc='center left', bbox_to_anchor=(1., 0.5))
    plt.tight_layout()
    plt.savefig(log_prefix+".pdf", bbox_inches='tight', pad_inches=0.5)
    if show_state_plt == 1:
        plt.show()
    # close figure
    plt.close(fig)

    # choose best freq and haploset based on lowest BIC
    best_BIC = min(all_states_df['best_BIC'])
    best_all_states_df = all_states_df[all_states_df['best_BIC']==best_BIC]
    best_state = best_all_states_df.loc[best_all_states_df['LL'].idxmax()]
    best_freq = best_state["freq"]
    best_haploset = eval(best_state['haplo_set'])
    best_X_boolean = best_state['X']
    best_LL = best_state['LL']

    return best_freq, best_haploset, best_X_boolean, best_LL, best_BIC

def generate_ML_haplotype_freq_dist(ML_freq, ML_haploset, X_boolean, partial_haplotype_df, outdir, overdispersion=1e-3, N=100):

    # compute alpha from computed ML_freq
    alpha_0 = 1/overdispersion - 1
    alpha_k = ML_freq*alpha_0

    # generate N number of random full haplotype frequencies
    if X_boolean == 1:
        pi_N = np.random.dirichlet(alpha_k[:-1], N)
    else:
        pi_N = np.random.dirichlet(alpha_k, N)

    all_nucpos_list = list(set(partial_haplotype_df['nucpos']))
    sorted_nucpos = sorted(set([k for v in all_nucpos_list for k in v]))

    ML_haplotype_freq_dist = np.zeros((N, len(ML_haploset)))
    Nidx = 0
    for pi in pi_N:
        # recreate partial haplotype df
        sample_partial_haplotype_df = []
        for nucpos_list in all_nucpos_list:

            # generate partial haploset
            nucpos_idx = [sorted_nucpos.index(nucpos) for nucpos in nucpos_list]
            partial_haploset = [tuple([haplotype[idx] for idx in nucpos_idx]) for haplotype in ML_haploset]

            # multinomial sampling based on pi
            n = partial_haplotype_df[partial_haplotype_df['nucpos']==nucpos_list]['count'].sum()
            X = np.random.multinomial(n, pi)

            # count partial haplotype
            phap_to_count = {}
            for cidx, count in enumerate(X):
                if count == 0:
                    continue
                try:
                    phap_to_count[partial_haploset[cidx]] += count
                except:
                    phap_to_count[partial_haploset[cidx]] = count

            for phap, count in phap_to_count.items():
                sample_partial_haplotype_df.append({'nucpos':nucpos_list, 'var': phap, 'count':count})

        # generate sample partial haplotype df
        sample_partial_haplotype_df = pd.DataFrame.from_dict(sample_partial_haplotype_df)
        sample_freq, sample_haploset, sample_X_boolean, sample_LL, sample_BIC = reconstruct_haplotype(partial_haplotype_df,
                                                                                                      "%s/temp_sample"%(outdir),
                                                                                                      verbose=0, searches=1,
                                                                                                      show_state_plt=0)

        # append sample_freq to distribution
        for hidx, haplotype in enumerate(ML_haploset):
            if haplotype in sample_haploset:
                ML_haplotype_freq_dist[(Nidx, hidx)] = sample_freq[sample_haploset.index(haplotype)]
        Nidx += 1

    return ML_haplotype_freq_dist

class DirichletMultiNomModel():
    """
    ## Aim: To get a maximum-likelihood estimate of noise - i.e. overdispersion $\varphi$ and variant
    proportions of each locus (Steven H Wu et al., Bioinformatics 2017)

    For each locus, we assume that the observed base counts in sample (replicate) $i$
    ($\mathbf{X_i}=\{x_a\}$ where $a$ denotes the residue type) follows a Dirchlet multinomial distribution:


    $\mathbf{X_i}\sim Multinomial(N_i,\mathbf{p})$ and $\mathbf{p} \sim Dirchlet(\mathbf{\alpha})$

    where $\mathbf{p}$ is a vector of proportions $\{p_a\}$ and $N_i$ is the sample size. The resulting
    likelihood function would then be:

    $
    \mathcal{L_i} = \log\left(\frac{\Gamma\left(\sum_a\alpha_{a}\right)}{\Gamma\left(\sum_a\alpha_{a}+N_i\right)}\right) + \log\left(\prod_a\frac{\Gamma\left(x_{a,i} + \alpha_{a}\right)}{\Gamma\left(\alpha_{a}\right)}\right)
    $

    Given that the expectation of the Dirchlet distribution $\pi_{a} = \frac{\alpha_{a}}{\sum_a\alpha_{a}}$,
    the expectation and variance of $x_{a,i}$ would be:

    $E(X_{a,i})=N_i\frac{\alpha_{a}}{\sum_a\alpha_{a}}=N_i\pi_{a}$
    $Var\left(X_{a,i}\right)=N_i\pi_{a}\left(1-\pi_{a}\right)\left(\frac{N_i+\sum_a\alpha_{a}}{1+\sum_a\alpha_{a}}\right)$

    Defining $\varphi=\frac{1}{1+\sum_a\alpha_{a}}$, we can reparameterise $\alpha_{a} = \frac{1-\varphi}{\varphi}\pi_{a}$ and as a result:

    $
    Var(X_{a,i}) = N\pi_{a}\left(1-\pi_{a}\right)\left(1+\left(N_i-1\right)\right)\varphi
    $


    If $\varphi\to0$, the Dirchlet multinomial model converges to a multinomial. In other words, $\varphi$ is a
    measure how overdispersed a multinomial distribution is (i.e. how much greater variability) and this can be
    a measure of noise when we fit the model to the dataset.

    Here, we estimate the maximum-likelihood model for each locus by finding the parameters that maximises
    $\mathcal{L}\left(X|\varphi,\pi\right)$ using the fixed-point Newton-Rhapson iteration method proposed by
    [Minka, 2012](https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf).
    """
    def __init__(self, x, max_iter=1000000, tol=1e-6, verbose=1):
        self.x = x
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose

    def d_ll(self, alpha_arr, x_arr):
        d_ll_row = digamma(alpha_arr.sum(axis=1, keepdims=True)) - digamma(x_arr.sum(axis=1, keepdims=True) + alpha_arr.sum(axis=1, keepdims=True)) + digamma(x_arr+alpha_arr) - digamma(alpha_arr)
        return d_ll_row.sum(axis=0, keepdims=True)

    def ll(self, alpha_arr, x_arr):
        ll_row = loggamma(alpha_arr.sum(axis=1, keepdims=True)) - loggamma(x_arr.sum(axis=1, keepdims=True) + alpha_arr.sum(axis=1, keepdims=True)) + (loggamma(x_arr+alpha_arr) - loggamma(alpha_arr)).sum(axis=1, keepdims=True)
        return ll_row.sum()

    def initial_guess(self, x_arr):
        return np.array([np.mean(x_arr/x_arr.sum(axis=1, keepdims=True), axis=0)])

    # NR optimisation to find pi and overdispersion
    def optimise(self):
        # add nominal tol to zero counts
        x = self.x+self.tol

        # initial guess of alpha
        alpha = self.initial_guess(x)
        #print ('initial', alpha)

        i = 1 # iteration counter
        step = int(np.ceil(self.max_iter/20)) # step for printing progress
        while i <= self.max_iter:
            #  solution found
            #print (self.d_ll(alpha, x)[0][0])
            if (self.d_ll(alpha, x)<self.tol).all() == True:
                if self.verbose == 1:
                    print ('Solution found at iteration {:<4d} (loglikelihood: {:.5f}).'.format(i, self.ll(alpha, x)))
                pi = alpha[0]/alpha[0].sum()
                dispersion = 1/(1+alpha[0].sum())
                # results rounded to 3 decimal places
                return (pi, dispersion)

            # update alpha
            num = (digamma(x + alpha) - digamma(alpha)).sum(axis=0, keepdims=True)
            den = (digamma(x.sum(axis=1, keepdims=True)+alpha.sum(axis=1, keepdims=True)) - digamma(alpha.sum(axis=1, keepdims=True))).sum()
            alpha = alpha*(num/den)
            i += 1

            # print progress
            if self.verbose == 1 and i%step==0:
                print ("Current iteration: {:<4d} (loglikelihood={:.5f})".format(i, self.ll(alpha, x)))

        # no solution found
        if self.verbose == 1:
            print ('Maximum iterations ({}) reached, solution not found.'.format(self.max_iter))
        return None

class DirichletMultiNomModel_ODEst():
    """
    The above estimates true variant proportions alongside the amount of overdispersion given
    repeated observations of a sample.

    Here, we compute the overdispersion parameter if we assume the mean variant proportions of the
    repeated observations adequately estimates true proportions.
    """
    def __init__(self, x, max_iter = 1000000, tol=1e-10, verbose=1):
        self.x = x
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose

    def d_ll(self, C, pi_arr, x_arr):
        d_ll_row = digamma(C) - digamma(x_arr.sum(axis=1, keepdims=True)+C) + (pi_arr*digamma(x_arr+(C*pi_arr)) - pi_arr*digamma(C*pi_arr)).sum(axis=1, keepdims=True)
        return d_ll_row.sum()

    def ll(self, C, pi_arr, x_arr):
        alpha_arr = C*pi_arr
        ll_row = loggamma(alpha_arr.sum(axis=1, keepdims=True)) - loggamma(x_arr.sum(axis=1, keepdims=True) + alpha_arr.sum(axis=1, keepdims=True)) + (loggamma(x_arr + alpha_arr) - loggamma(alpha_arr)).sum(axis=1, keepdims=True)
        return ll_row.sum()

    def optimise(self):
        # initial guess of overdispersion
        C = 200

        # data
        X = self.x+self.tol
        pi = X.sum(axis=0, keepdims=True)/X.sum()

        i = 1
        step = int(np.ceil(self.max_iter/20)) # step for printing progress

        while i <= self.max_iter:
            # satisfy tol
            if np.isclose(self.d_ll(C, pi, X), 0):
                if self.verbose == 1:
                    print ('Solution found at iteration {:<4d} (loglikelihood: {:.5f}).'.format(i, self.ll(C, pi, X)))
                return (1/(C+1))

            num = (pi*digamma(X+(C*pi)) - pi*digamma(C*pi)).sum()
            den = (digamma(X.sum(axis=1, keepdims=True)+C) - digamma(C)).sum()
            C = C*(num/den)
            i += 1

            # print progress
            if self.verbose == 1 and i%step==0:
                print ("Current iteration: {:<4d} (loglikelihood: {:.5f})".format(i, self.ll(C, pi, X)))

        # no solution found
        if self.verbose == 1:
            print ('Maximum iterations ({}) reached, solution not found.'.format(max_iter))
        return None
