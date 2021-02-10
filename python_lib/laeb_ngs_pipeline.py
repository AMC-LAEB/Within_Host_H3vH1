from Bio import SeqIO
import pandas as pd
import numpy as np
import subprocess
import os
import re
import time
import random
import itertools
import gzip
import json
import platform
import ast
import multiprocessing as mp
from multiprocessing import Manager
from os.path import expanduser
from importlib.machinery import SourceFileLoader
from scipy.stats import binom

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

def barcode_errors(sorted_sampid):

    base_map = {"A":1, "T":2, "G":3, "C":4, "-":0}

    for sampid in sorted_sampid:
        # get mapped reads with true aln pos
        mapped_reads_df = pd.read_csv("./results/barcode_sort/mapped_%s_true_nucpos.csv"%(sampid)).set_index(["amplicon_nr", "barcode"])
        for amplicon_nr in mapped_reads_df.index.get_level_values(0):
            amp_mapped_reads_df = mapped_reads_df.loc[amplicon_nr]

            # for each unique barcode in amplicon
            for barcode in set(amp_mapped_reads_df.index):
                bc_amp_mapped_reads_df = amp_mapped_reads_df.loc[barcode]
                # make mini map_nuc_results
                map_nuc_results = {}
                #if isinstance(bc_amp_mapped_reads_df, pd.Series):

                print (bc_amp_mapped_reads_df)
                aln_block = np.zeros((len(bc_amp_mapped_reads_df), bc_amp_mapped_reads_df["end_nucaln_pos_plus1"].max()))
                for r, (barcode, row) in enumerate(bc_amp_mapped_reads_df.iterrows()):
                    print (aln_block[(r,)])
                    print (list(map(lambda res:base_map[res], list(row["seq"]))))
                raise Exception

    return

def assess_perr(var_freq, depth, q, err_tol):
    i = 0
    pe = 10**(-q/10)
    pErr = 1.
    while i <= var_freq:
        pErr -= binom.pmf(i, depth, pe)
        i += 1

    if (np.isnan(pErr)):
        raise Exception('pErr = nan; check freq/coverage/q = {}, {}, {}'.format(var_freq, coverage, q))

    return pErr < err_tol

def variant_calling_worker(sampid, codon_table, all_bases, variant_call_df, sorted_refnames, base_qual_threshold, min_cov, min_var_prop, gene_to_proteinorf, err_tol, alnrefnun_to_truerefnum, ha_numbering_conversion, HAnum_subtype, min_var_freq):

    # parse map_nuc_results and map_codon_results
    try:
        map_nuc_results = pd.read_csv('./results/map_nuc_results_{}.csv'.format(sampid), keep_default_na=False).set_index(['Gene', 'Position']).sort_index()
    except:
        #print ('WARNING: No mapped reads for %s.'%sampid)
        return

    map_codon_results = pd.read_csv('./results/map_codon_results_{}.csv'.format(sampid), keep_default_na=False).set_index(['Gene', 'Position']).sort_index()

    for gene in sorted_refnames:
        # section map_nuc_results for each gene segment
        try:
            gene_map_nuc_results = map_nuc_results.loc[gene].copy().sort_index()
        except:
            #print ('WARNING: No mapped reads for gene segment %s for %s.'%(gene, sampid))
            continue

        ## -------------------------
        #   Variant calling
        ## -------------------------

        # only consider nucleotide above minimum coverage
        gene_map_nuc_results = gene_map_nuc_results[gene_map_nuc_results['Coverage']>=min_cov]
        if len(gene_map_nuc_results) == 0:
            # lack of coverage for gene segment
            #print ('WARNING: Lack of coverage for gene segment %s for %s.'%(gene, sampid))
            continue

        # filter for nucleotide positions which are not conserved
        polymorphic_nuc_positions = list(gene_map_nuc_results[gene_map_nuc_results[all_bases].max(axis=1)!=gene_map_nuc_results['Coverage']].index)
        if len(polymorphic_nuc_positions) == 0: # no polymorphic nucleotide positions
            continue

        # first consolidate all nucleotide variants passing quality filters
        nucpos_to_varinfo = {}
        for nucpos in polymorphic_nuc_positions:

            # nucpos should be within cdr
            try:
                gene_coord_nucrow = gene_to_proteinorf.loc[gene, nucpos]
            except:
                # non-coding region ignored
                continue

            # get map_nuc_results for nucpos in gene
            nucpos_gene_map_nuc_results = gene_map_nuc_results.loc[nucpos]

            # consensus base
            nuc_consensus = nucpos_gene_map_nuc_results['Consensus']
            nuc_consensus_freq = nucpos_gene_map_nuc_results[nuc_consensus]

            # nuc_depth for position
            nuc_depth = nucpos_gene_map_nuc_results['Coverage']

            # for any base > 0 counts
            for base in (nucpos_gene_map_nuc_results[all_bases][nucpos_gene_map_nuc_results[all_bases]>0].index):
                if base == nuc_consensus:
                    continue

                # get nucleotide frequency for each variant base
                var_nuc_freq = nucpos_gene_map_nuc_results[base]

                # check that it is above min frequency
                if var_nuc_freq < min_var_freq:
                    continue

                # check that it is above min proportion
                var_nuc_prop = var_nuc_freq/nuc_depth
                if var_nuc_prop < min_var_prop:
                    continue

                # check if var_nuc_freq could be attributed to error
                if assess_perr(var_nuc_freq, nuc_depth, base_qual_threshold, err_tol) == False:
                    continue

                # save variants
                try:
                    nucpos_to_varinfo[nucpos]['variants'][base] = {'var_nuc_freq':var_nuc_freq, 'var_nuc_prop':var_nuc_prop}
                except:
                    # ...alongside other information in order of consensus, nuc_depth
                    nucpos_to_varinfo[nucpos] = {'nuc_consensus':nuc_consensus,
                                                 'nuc_consensus_freq':nuc_consensus_freq,
                                                 'nuc_depth':nuc_depth,
                                                 'variants':{base:{'var_nuc_freq':var_nuc_freq, 'var_nuc_prop':var_nuc_prop}}}

        for nucpos in nucpos_to_varinfo.keys():

            variant_info = nucpos_to_varinfo[nucpos]

            # nucleotide consensus and depth
            nuc_consensus = variant_info['nuc_consensus']
            nuc_consensus_freq = variant_info['nuc_consensus_freq']
            nuc_depth = variant_info['nuc_depth']

            # get protein coordinates for nuc_pos
            try:
                gene_coord_nucrow = gene_to_proteinorf.loc[gene, nucpos]
            except:
                continue

            # for each protein that the nucpos in the gene encode
            for protein in gene_coord_nucrow.index:
                # position of protein that nucleotide translates to
                aapos = gene_coord_nucrow.loc[protein]['aa']

                # convert HA numbering to desired subtype numbering
                HA_num_type = None
                if protein == "HA" and isinstance(ha_numbering_conversion, pd.DataFrame):
                    if aapos > ha_numbering_conversion.index.max():
                        continue
                    else:
                        ha_aapos = ha_numbering_conversion.loc[aapos, HAnum_subtype]

                        if pd.isna(ha_numbering_conversion.loc[aapos, HAnum_subtype]):
                            HA_num_type = ha_numbering_conversion.index.name
                        else:
                            aapos = int(ha_aapos)
                            HA_num_type = HAnum_subtype

                # codon position of nucleotide
                frame = gene_coord_nucrow.loc[protein]['frame']

                # get starting position of codon
                codon_start_pos = nucpos-(frame-1)

                # expected aa consensus based on nucleotide consensus
                try:
                    expected_aa_consensus = codon_table[''.join(list(gene_map_nuc_results.loc[codon_start_pos:codon_start_pos+2]['Consensus']))]
                except:
                    # not enough reads on nucleotide sites constituting the codon
                    for base, curr_var_dict in variant_info['variants'].items():
                        var_nuc_freq = curr_var_dict['var_nuc_freq']
                        var_nuc_prop = curr_var_dict['var_nuc_prop']

                        print_nuc_pos = alnrefnun_to_truerefnum[gene][nucpos] if alnrefnun_to_truerefnum != None else nucpos
                        variant_call_df.append({'sampid':sampid,
                                                'gene':gene,
                                                'nucpos':print_nuc_pos, 'nuc_var':base, 'nuc_prop':var_nuc_prop, 'nuc_freq':var_nuc_freq,
                                                'nuc_consensus':nuc_consensus, 'nuc_consensus_freq':nuc_consensus_freq, 'nuc_coverage':nuc_depth,
                                                'protein':protein,
                                                'aapos':aapos, "HA_num_type":HA_num_type, 'aa_var':None, 'aa_prop':None, 'aa_freq':None,
                                                'expected_aa_consensus':None, 'aa_consensus':None, 'aa_consensus_freq':None, 'aa_coverage':None,
                                                'codon_pos':frame, 'codon':None, 'nonsyn':None})
                    continue

                # get the actual codon reads
                codon_map_codon_results = map_codon_results.loc[gene, codon_start_pos]
                codon_map_codon_results = codon_map_codon_results[codon_map_codon_results.gt(0)]

                # get depth of aa based on sum of codons
                aa_depth = codon_map_codon_results.sum()

                # ignore premature stop codon
                stop_codon_aapos = max(gene_to_proteinorf.xs(protein, level=2)['aa'])
                # only allowed nucpos to have stop codons
                stop_codon_nucpos = list(gene_to_proteinorf.xs(protein, level=2)[gene_to_proteinorf.xs(protein, level=2)['aa'] == stop_codon_aapos].index.get_level_values(1))

                # only consider variant codon where bases for site in qn is a significant variant or consensus
                variant_codon_list = [codon for codon in codon_map_codon_results.index if codon[frame-1] in [nuc_consensus]+list(variant_info['variants'].keys())]
                if nucpos not in stop_codon_nucpos:
                   # remove premature stop codons
                    variant_codon_list = [codon for codon in variant_codon_list if codon_table[codon] != '*']

                # continue if no more variants left
                if len(variant_codon_list) == 0:
                    continue

                # filter out any potentially errorneous nucleotide variants in other codon sites
                for other_frame in range(1, 4):
                    if other_frame == frame:
                        continue
                    other_nucpos = nucpos+(other_frame-frame)

                    try:
                        allowed_other_nuc_variants = [nucpos_to_varinfo[other_nucpos]['nuc_consensus']] + list(nucpos_to_varinfo[other_nucpos]['variants'].keys())
                    except:
                        # other_nucpos has no significant variants
                        allowed_other_nuc_variants = [gene_map_nuc_results.loc[other_nucpos, 'Consensus']]

                    variant_codon_list = [codon for codon in variant_codon_list if codon[other_frame-1] in allowed_other_nuc_variants]

                # go through each variant base
                for base, curr_var_dict in variant_info['variants'].items():
                    var_nuc_freq = curr_var_dict['var_nuc_freq']
                    var_nuc_prop = curr_var_dict['var_nuc_prop']

                    # get aa variants and consensus
                    aa_consensus_freq = -1
                    varcodon_to_count = {}
                    varbase_accounted_boolean = 0

                    for codon in variant_codon_list:
                        codon_count = codon_map_codon_results[codon]
                        aa = codon_table[codon]

                        # observed codon accounted for nuc consensus
                        if (codon[frame-1]) == nuc_consensus:
                            if codon_count > aa_consensus_freq:
                                aa_consensus_freq = codon_count
                                aa_consensus = aa
                                aa_consensus_codon = codon
                                continue

                        # observed codon accounted for var base
                        elif (codon[frame-1]) == base:
                            varcodon_to_count[codon] = (aa, codon_count)
                            varbase_accounted_boolean = 1
                            continue

                    if aa_consensus_freq < 0:
                        # don't know what is codon consensus by nuc consensus position
                        if varbase_accounted_boolean == 0:
                            # cant find varbase codon either
                            print_nuc_pos = alnrefnun_to_truerefnum[gene][nucpos] if alnrefnun_to_truerefnum != None else nucpos
                            variant_call_df.append({'sampid':sampid,
                                                    'gene':gene,
                                                    'nucpos':print_nuc_pos, 'nuc_var':base, 'nuc_prop':var_nuc_prop, 'nuc_freq':var_nuc_freq,
                                                    'nuc_consensus':nuc_consensus, 'nuc_consensus_freq':nuc_consensus_freq, 'nuc_coverage':nuc_depth,
                                                    'protein':protein,
                                                    'aapos':aapos, "HA_num_type":HA_num_type, 'aa_var':None, 'aa_prop':None, 'aa_freq':None,
                                                    'expected_aa_consensus':None, 'aa_consensus':None, 'aa_consensus_freq':None, 'aa_coverage':None,
                                                    'codon_pos':frame, 'codon':None, 'nonsyn':None})
                        else:
                            # ... but variant codon accounted for
                            # - use expected nucleotide-inferred codon
                            for codon, (aa, codon_count) in varcodon_to_count.items():
                                # boolean to determine if aa is syn, nonsyn or stop codon
                                if aa == expected_aa_consensus:
                                    nonsyn_bool = 0
                                elif aa == "*":
                                    nonsyn_bool = -1 # stop codon
                                else:
                                    nonsyn_bool = 1
                                print_nuc_pos = alnrefnun_to_truerefnum[gene][nucpos] if alnrefnun_to_truerefnum != None else nucpos
                                variant_call_df.append({'sampid':sampid,
                                                        'gene':gene,
                                                        'nucpos':print_nuc_pos, 'nuc_var':base, 'nuc_prop':var_nuc_prop, 'nuc_freq':var_nuc_freq,
                                                        'nuc_consensus':nuc_consensus, 'nuc_consensus_freq':nuc_consensus_freq, 'nuc_coverage':nuc_depth,
                                                        'protein':protein,
                                                        'aapos':aapos, "HA_num_type":HA_num_type, 'aa_var':aa, 'aa_prop':codon_count/aa_depth, 'aa_freq':codon_count,
                                                        'expected_aa_consensus':expected_aa_consensus, 'aa_consensus':None, 'aa_consensus_freq':None, 'aa_coverage':aa_depth,
                                                        'codon_pos':frame, 'codon':codon, 'nonsyn':nonsyn_bool})
                        continue
                    else:
                        # consensus codon accounted for
                        if varbase_accounted_boolean == 0:
                            # but cant find varbase codon - use consensus codon
                            hypothetical_varbase_codon = list(aa_consensus_codon)
                            hypothetical_varbase_codon[frame-1] = base
                            hypothetical_varbase_codon = "".join(hypothetical_varbase_codon)

                            codon_count = None
                            aa = codon_table[hypothetical_varbase_codon]

                            varcodon_to_count[hypothetical_varbase_codon] = (aa, codon_count)

                        for codon, (aa, codon_count) in varcodon_to_count.items():
                            # boolean to determine if aa is syn, nonsyn or stop codon
                            if aa == aa_consensus:
                                nonsyn_bool = 0
                            elif aa == "*":
                                nonsyn_bool = -1 # stop codon
                            else:
                                nonsyn_bool = 1
                            # note that nucpos is corrected to start from ATG
                            print_nuc_pos = alnrefnun_to_truerefnum[gene][nucpos] if alnrefnun_to_truerefnum != None else nucpos
                            variant_call_df.append({'sampid':sampid,
                                                    'gene':gene,
                                                    'nucpos':print_nuc_pos, 'nuc_var':base, 'nuc_prop':var_nuc_prop, 'nuc_freq':var_nuc_freq,
                                                    'nuc_consensus':nuc_consensus, 'nuc_consensus_freq':nuc_consensus_freq, 'nuc_coverage':nuc_depth,
                                                    'protein':protein,
                                                    'aapos':aapos, "HA_num_type":HA_num_type, 'aa_var':aa, 'aa_prop':None if codon_count == None else codon_count/aa_depth, 'aa_freq':codon_count,
                                                    'expected_aa_consensus':expected_aa_consensus, 'aa_consensus':aa_consensus, 'aa_consensus_freq':aa_consensus_freq, 'aa_coverage':aa_depth,
                                                    'codon_pos':frame, 'codon':codon, 'nonsyn':nonsyn_bool})

    return

def variant_calling(sorted_sampid, sorted_refnames, base_qual_threshold, min_cov, min_var_prop, gene_to_proteinorf, err_tol, alnrefnun_to_truerefnum=None, ha_numbering_conversion=None, HAnum_subtype=None, min_var_freq=0, threadnum=4, reanalyze_bool=0):
    """
    Call variants
    """
    print ("\nTallying minority variants with minimum coverage of %i at >%.1f%% with base calling error tolerance at %.1f%%..."%(min_cov, 100*min_var_prop, 100*err_tol))
    # parameters
    codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
    all_bases = ['A', 'T', 'G', 'C']

    varcall_fname = './results/variant_call_MinCoV{}_MinProp{}_MinFreq{}_ErrTol{}.csv'.format(min_cov, min_var_prop, min_var_freq, err_tol)

    if os.path.isfile(varcall_fname) and reanalyze_bool == 0:
        variant_call_df = pd.read_csv(varcall_fname)
        variant_call_df['protein'] = variant_call_df['protein'].fillna('NA')
        variant_call_df = variant_call_df.set_index(['sampid', 'gene', 'nucpos']).sort_index()
    else:
        variant_call_df = Manager().list()

        pool = mp.Pool(processes=threadnum)
        results = [pool.apply_async(variant_calling_worker, args=(sampid, codon_table, all_bases, variant_call_df, sorted_refnames, base_qual_threshold, min_cov, min_var_prop, gene_to_proteinorf, err_tol, alnrefnun_to_truerefnum, ha_numbering_conversion, HAnum_subtype, min_var_freq,)) for sampid in sorted_sampid]
        output = [p.get() for p in results]

        # change it back to python list type
        variant_call_df = list(variant_call_df)
        variant_call_df = pd.DataFrame.from_dict(variant_call_df)
        variant_call_df = variant_call_df.set_index(['sampid', 'gene', 'nucpos']).sort_index()
        variant_call_df.to_csv(varcall_fname)

    return variant_call_df

def sort_reads_worker(sampid, barcode_amp_stats, gene_to_proteinorf, gene_amplicon_to_unique_actual_pos, alnrefnun_to_truerefnum, barcode_coordinates, min_barcode_filter):

    more_than_one_possible_amp_count = 0

    # get all mapped reads for sampid
    mapped_reads_fname = './results/mapped/mapped_%s.csv'%(sampid)
    try:
        mapped_reads_df = pd.read_csv(mapped_reads_fname)
    except:
        print ('No mapped reads found for {}'.format(sampid))
        return

    # skip if already checked for amplicon_nr
    headers_to_check = ['amplicon_nr']
    if barcode_coordinates != None:
        headers_to_check.append("barcode")
    if len(set(list(mapped_reads_df))&set(headers_to_check)) == len(headers_to_check):
        #print ("HAHA")
        return

    # write new mapped reads with amplicon_nr and barcode info (if any)
    new_mapped_reads = []
    for r, row in mapped_reads_df.iterrows():
        new_row = row.to_dict().copy()

        gene = row["gene"]
        start_pos = row["start_nucaln_pos"]
        end_pos = row["end_nucaln_pos_plus1"]
        read_nuc_seq = row["seq"]
        read_pos_range = range(start_pos, end_pos)

        # convert refnum pos range to true nuc pos if needed
        if alnrefnun_to_truerefnum == None:
            true_pos_range = read_pos_range
        else:
            true_pos_range = [alnrefnun_to_truerefnum[gene][p] for p in read_pos_range if p in alnrefnun_to_truerefnum[gene]]

        # get amplicon_nr of read
        found_amp_nr = None
        cand_amp_nr = {amplicon_nr:len(set(true_pos_range)&set(gene_amplicon_to_unique_actual_pos[gene][amplicon_nr])) for amplicon_nr in gene_amplicon_to_unique_actual_pos[gene].keys() if len(set(true_pos_range)&set(gene_amplicon_to_unique_actual_pos[gene][amplicon_nr])) > 0}

        if len(cand_amp_nr) == 1:
            found_amp_nr = max(cand_amp_nr, key=cand_amp_nr.get)
        # more than one possible candidate amplicon when searched by unique sites to ampliccon - does not make sense so skip
        elif len(cand_amp_nr) > 1:
            more_than_one_possible_amp_count += 1
            continue

        # find barcodes of mapped reads
        if barcode_coordinates != None:
            barcode = None
            gene_barcode_coordinates = barcode_coordinates[gene]
            for idx, bc_range in enumerate(gene_barcode_coordinates):
                bc_range = range(bc_range[0], bc_range[1])
                if set(bc_range) <= set(read_pos_range):
                    barcode = "".join([read_nuc_seq[idx] for idx, read_pos in enumerate(read_pos_range) if read_pos in bc_range])

                    if re.search("^\-+$", barcode): # overlapping amplicon
                        continue
                    else:
                        # ambiguous amp nr can also solved by barcode order if each amplicon has only one barcodes
                        if found_amp_nr == None:
                            if len(gene_barcode_coordinates) == len(gene_amplicon_to_unique_actual_pos[gene].keys()):
                                found_amp_nr = idx+1
                        else:
                            # barcode order does not ties with amplicon_nr based on unique sites - again, skip.
                            if idx+1 != found_amp_nr:
                                continue
                        break
            # barcode with poor quality bases (gaps) = no barcode
            if barcode != None and re.search("-", barcode):
                barcode = None
            # discard reads with no barcode
            if barcode == None:
                continue
            # save barcode
            new_row['barcode'] = barcode

        # save amplicon_nr
        new_row['amplicon_nr'] = found_amp_nr
        # save new read
        new_mapped_reads.append(new_row)

    new_mapped_reads = pd.DataFrame.from_dict(new_mapped_reads)

    # filter reads where barcode < min_barcode_filter
    if min_barcode_filter != None:

        if not os.path.isdir("./results/barcode_sort"):
            os.mkdir("./results/barcode_sort")

        # find all barcodes to discard

        new_mapped_reads_with_truenucpos = [] # new mapped_reads_df for sampid with cdr sequences and true nuc pos positions
        index_to_remove = []
        new_mapped_reads = new_mapped_reads.set_index(["amplicon_nr", "barcode"]).sort_index()

        for amplicon_nr in sorted(set(new_mapped_reads.index.get_level_values(0))):
            barcode_to_discard = []

            amp_new_mapped_reads = new_mapped_reads.loc[amplicon_nr]

            for barcode in set(amp_new_mapped_reads.index):
                barcode_new_mapped_reads = amp_new_mapped_reads.loc[barcode]
                # only one read per barcode (and if min_barcode_filter > 1)
                if isinstance(barcode_new_mapped_reads, pd.Series) and min_barcode_filter > 1:
                    barcode_to_discard.append(barcode)
                # number of reads with barcode < min_barcode_filter
                elif len(barcode_new_mapped_reads) < min_barcode_filter:
                    barcode_to_discard.append(barcode)

            # write result to files
            # barcode statistics
            barcode_amp_stats.append({"sampid":sampid, "gene":gene, "amplicon_nr":amplicon_nr,
                                      "total_reads":len(amp_new_mapped_reads), # total number of reads for each amplicon
                                      "num_of_reads_with_barcode_abv_thres":len(amp_new_mapped_reads[~amp_new_mapped_reads.index.isin(barcode_to_discard)]), # number of reads with unique barcode satisfying min_barcode_filter
                                      "num_of_unique_barcodes":len(set(amp_new_mapped_reads.index)), # total number of uniquie barcodes
                                      "num_of_unique_barcoes_abv_thres":len(set(amp_new_mapped_reads.index)-set(barcode_to_discard))}) # number of unique barcodes above min_barcode_filter

            # write to file all barcodes
            with open("./results/barcode_sort/%s_%i_allbar.txt"%(sampid, amplicon_nr), "w") as output:
                output.write("\n".join(list(amp_new_mapped_reads.index)))
            # all barcodes satisfying thres
            with open("./results/barcode_sort/%s_%i_allminbar.txt"%(sampid, amplicon_nr), "w") as output:
                output.write("\n".join(list([barcode for barcode in amp_new_mapped_reads.index if barcode not in barcode_to_discard])))
            # all reads with barcodes abv thres
            with open("./results/barcode_sort/%s_%i_allminread.txt"%(sampid, amplicon_nr), "w") as output:
                # discard reads with < min bcc filter
                amp_new_mapped_reads = amp_new_mapped_reads[~amp_new_mapped_reads.index.isin(barcode_to_discard)]
                amp_new_mapped_reads_with_truenucpos = []

                for barcode, row in amp_new_mapped_reads.iterrows():
                    new_row = row.to_dict().copy()
                    new_row["amplicon_nr"] = amplicon_nr
                    new_row["barcode"] = barcode

                    gene = row["gene"]
                    start_pos = row["start_nucaln_pos"]
                    end_pos = row["end_nucaln_pos_plus1"]
                    read_nuc_seq = row["seq"]
                    read_pos_range = range(start_pos, end_pos)

                    # write sequence within coding region only
                    cdr_range = list(gene_to_proteinorf.loc[gene].index.get_level_values(0))
                    cdr_seq = []
                    true_nuc_pos = []
                    for idx, ref_pos in enumerate(read_pos_range):
                        if ref_pos in cdr_range:
                            cdr_seq.append(read_nuc_seq[idx])
                            true_nuc_pos.append(alnrefnun_to_truerefnum[gene][ref_pos])
                    if len(cdr_seq) == 0: # read does not cover any cdr region
                        continue
                    cdr_seq = "".join(cdr_seq)

                    output.write("%i %s\n"%(amplicon_nr, cdr_seq))



                    new_row["seq"] = cdr_seq
                    new_row["start_nucaln_pos"] = min(true_nuc_pos)
                    new_row["end_nucaln_pos_plus1"] = max(true_nuc_pos)+1
                    amp_new_mapped_reads_with_truenucpos.append(new_row)

                amp_new_mapped_reads_with_truenucpos = pd.DataFrame.from_dict(amp_new_mapped_reads_with_truenucpos)
                new_mapped_reads_with_truenucpos.append(amp_new_mapped_reads_with_truenucpos)

            # add amplicon, barcode sets to discard
            index_to_remove += [(amplicon_nr, barcode) for barcode in barcode_to_discard]

        new_mapped_reads_with_truenucpos = pd.concat(new_mapped_reads_with_truenucpos, ignore_index=True)
        new_mapped_reads_with_truenucpos = new_mapped_reads_with_truenucpos.set_index(['gene', 'amplicon_nr', 'barcode']).sort_index()
        new_mapped_reads_with_truenucpos.to_csv("./results/barcode_sort/mapped_%s_true_nucpos.csv"%(sampid))

        # discard reads with < min bcc filter
        new_mapped_reads = new_mapped_reads[~new_mapped_reads.index.isin(index_to_remove)]
        new_mapped_reads = new_mapped_reads.reset_index()

    # save new_mapped_reads
    new_idx_names = ["gene", "amplicon_nr"]
    if min_barcode_filter != None:
        new_idx_names.append("barcode")
    new_mapped_reads = new_mapped_reads.set_index(new_idx_names).sort_index()
    new_mapped_reads.to_csv(mapped_reads_fname)
    return

def sort_reads(sorted_sampid, primer_coords, gene_to_proteinorf, threadnum=4, alnrefnun_to_truerefnum=None, barcode_coordinates=None, min_barcode_filter=None):

    print ("Sort reads by amplicons %s..."%("and barcodes" if barcode_coordinates!=None else ""))

    primer_coords_df = pd.read_csv(primer_coords).set_index(["gene", "amplicon_nr"])

    # assign amplicon by matching unique residues that can only be found each amplicon
    gene_amplicon_to_unique_actual_pos = {}
    # for each gene

    for gene in set(primer_coords_df.index.get_level_values(0)):
        gene_primer_coords_df = primer_coords_df.loc[gene]
        # get list of amplicons
        amp_list = sorted(set(gene_primer_coords_df.index))
        # get all positions covered by amplicon
        amp_to_full_range = {amplicon_nr:range(gene_primer_coords_df.loc[amplicon_nr, "coord_start"].min(), gene_primer_coords_df.loc[amplicon_nr, "coord_end"].max()) for amplicon_nr in amp_list}

        # get unique positions covered by amplicon
        for amplicon_nr in amp_list:
            try:
                gene_amplicon_to_unique_actual_pos[gene][amplicon_nr] = list(set(amp_to_full_range[amplicon_nr]) - set([k for v in [amp_to_full_range[amp2] for amp2 in amp_list if amp2 != amplicon_nr] for k in v]))
            except:
                gene_amplicon_to_unique_actual_pos[gene] = {amplicon_nr:list(set(amp_to_full_range[amplicon_nr]) - set([k for v in [amp_to_full_range[amp2] for amp2 in amp_list if amp2 != amplicon_nr] for k in v]))}

    # to save barcode statistics
    # define multiprocessing list
    barcode_amp_stats = Manager().list()
    # set up worker threads
    pool = mp.Pool(processes=threadnum)
    results = [pool.apply_async(sort_reads_worker, args=(sampid, barcode_amp_stats, gene_to_proteinorf, gene_amplicon_to_unique_actual_pos, alnrefnun_to_truerefnum, barcode_coordinates, min_barcode_filter,)) for sampid in sorted_sampid]
    output = [p.get() for p in results]
    # change it back to python list type
    barcode_amp_stats = list(barcode_amp_stats)
    pool.close()

    print ("...done.")

    if len(barcode_amp_stats) > 0:
        barcode_amp_stats = pd.DataFrame.from_dict(barcode_amp_stats).set_index(["sampid", "gene", "amplicon_nr"]).sort_index()
        barcode_amp_stats.to_csv("./results/barcode_amplicon_stats.csv")
        print (barcode_amp_stats)
    elif os.path.isfile("./results/barcode_amplicon_stats.csv"):
        barcode_amp_stats = pd.read_csv("./results/barcode_amplicon_stats.csv").set_index(["sampid", "gene", "amplicon_nr"]).sort_index()
        print (barcode_amp_stats)

    return

def parse_sam_worker(sampid, mapping_stats, sorted_refnames, data_folder, base_qual_threshold, max_indel_count, max_indel_freq, nucpos_shift, sampid_prefix, sampid_suffix, standardised_gene_name, save_unmapped_reads_bool):

    mapped_reads_fname = './results/mapped/mapped_%s.csv'%(sampid)
    unmapped_reads = {}
    fail_qual_reads = {}
    mapped_read_count = 0
    map_nuc_results = {}
    map_codon_results = {}
    mapped_reads = []

    sam_fname = '%s/align/%s%s%s.sam.gz'%(data_folder, sampid, sampid_prefix, sampid_suffix)
    fhandle = gzip.open(sam_fname, 'rb').readlines()

    reference_names = []
    for l in range(len(fhandle)):

        line = fhandle[l].decode('utf-8').strip()
        # skip file header lines
        if re.search('^\@', line):
            continue

        # sort fields
        fields = line.split('\t')
        read_id = fields[0]
        flag = int(fields[1])
        # convert gene name to standardize name
        try:
            refname = standardised_gene_name[fields[2]]
        except:
            refname = fields[2]
        start_pos = int(fields[3])
        map_qual = int(fields[4])
        cigar = fields[5]
        #template_len = int(fields[8])
        read_seq = fields[9]
        read_qual = fields[10]

        # exclude unmapped and non-primary read alignments (flag >= 256)
        if (flag >= 256):
            continue

        # unmapped flags for single-end reads
        if (flag == 4) or (flag == 20):
            unmapped_reads[read_id] = read_seq
            continue

        # exclude reads that did not map based on CIGAR string
        if cigar == '*':
            unmapped_reads[read_id] = read_seq
            continue

        # exclude reads with too many indels
        if re.search("[ID]", cigar):
            indel_count = {'I':0, 'D':0}
            for expr in re.findall("(\d+)([ID])", cigar):
                indel_count[expr[-1]] += int(expr[0])

            match_count = sum(map(int, re.findall("(\d+)M", cigar)))
            indel_count = sum(indel_count.values())

            if max_indel_count == None:
                if indel_count > match_count: # discard reads where indels > match
                    fail_qual_reads[read_id] = read_seq
                    continue
            else:
                if indel_count > max_indel_count:
                    fail_qual_reads[read_id] = read_seq
                    continue

            if max_indel_freq == None:
                if indel_count/match_count > 1.: # discard reads where there are more indels than matched
                    fail_qual_reads[read_id] = read_seq
                    continue
            else:
                if indel_count/match_count > max_indel_freq:
                    fail_qual_reads[read_id] = read_seq
                    continue

        # exclude read that did not map to reference
        if refname not in sorted_refnames:
            raise Exception("%s was mapped to %s which was not found in reference FASTA."%(read_id, refname))
            #unmapped_reads[read_id] = read_seq
            #continue

        # shift start_pos to cdr-numbering
        try:
            start_pos -= nucpos_shift[refname] # based on cdr-numbering
        except:
            pass

        # expand cigar
        ExpandedCIGAR = []
        for _ in re.finditer('(\d+)([MIDNSHP\=X])', cigar):
            ExpandedCIGAR += int(_.group(1))*[_.group(2)]

        curr_aln_num = -1 # reset reference base # addition
        rbase_idx = -1 # reset base index in read_seq
        read_nuc_seq = []
        insertion_ridx = []

        # add to mapped_read_count
        mapped_read_count += 1

        for char in ExpandedCIGAR:
            # H/P consumes neither query or reference
            if re.search('[HP]', char):
                continue

            # I/S consumes only query
            elif re.search('[IS]', char): # insertion
                rbase_idx += 1
                insertion_ridx.append(rbase_idx)
            else:
                # add to reference_position for M=X and D/N
                curr_aln_num += 1
                ref_pos = start_pos + curr_aln_num

                # D/N consumes only reference
                if re.search('[DN]', char):
                    if char == "D":
                        read_nuc_seq.append('-') # deletion
                    else:
                        read_nuc_seq.append('N') # skipped region from reference
                    continue
                else:
                    rbase_idx += 1

            # only tally ATGC
            read_base = read_seq[rbase_idx]
            if re.search("[^ATGC]", read_base):
                read_nuc_seq.append('N')
                continue

            # only tally bases > required quality threshold and isn't ambiguous
            base_qual = ord(read_qual[rbase_idx])-33

            if base_qual < base_qual_threshold:
                read_nuc_seq.append('N')
                continue
            else:
                # add read_base to read_nuc_seq
                read_nuc_seq.append(read_base)

        mapped_reads.append({'gene':refname, 'start_nucaln_pos':start_pos, 'end_nucaln_pos_plus1':ref_pos+1, "insertion":insertion_ridx, 'seq':''.join(read_nuc_seq)})

    mapping_stats.append({'sampid':sampid, 'mapped':mapped_read_count,
                          'unmapped':len(unmapped_reads), 'fail_qual':len(fail_qual_reads)})

    # check if there were any mapped reads
    if len(mapped_reads) == 0:
        print ('No mapped reads found for {}'.format(sampid))
        return

    # save mapped reads
    mapped_reads = pd.DataFrame.from_dict(mapped_reads).set_index(['gene', 'start_nucaln_pos']).sort_index()
    mapped_reads.to_csv(mapped_reads_fname)

    # save unampped reads
    if save_unmapped_reads_bool == 1 and len(unmapped_reads) > 0:
        with open('./results/unmapped/{}_unmapped.fasta'.format(sampid), 'w') as output:
            if len(unmapped_reads) > 100:
                sample_unmapped = random.sample(unmapped_reads.keys(), 100)
            else:
                sample_unmapped = unmapped_reads.keys()
            for read_id in sample_unmapped:
                output.write('>{}\n{}\n'.format(read_id, unmapped_reads[read_id]))

    return

def parse_sam(sorted_sampid, sorted_refnames, data_folder, base_qual_threshold, max_indel_count=None, max_indel_freq=None, nucpos_shift=None, sampid_prefix="", sampid_suffix="", standardised_gene_name=None, save_unmapped_reads_bool=0, threadnum=4, reanalyze_bool=0, plt_show=0):
    """
    Parse output SAM files from bowtie2
    """
    print ("\nParsing SAM files...")

    # create results directories
    if not os.path.isdir('./results/'):
        os.mkdir('./results')

    if save_unmapped_reads_bool == 1 and not os.path.isdir('./results/unmapped'):
        os.mkdir('./results/unmapped')

    if not os.path.isdir('./results/mapped'):
        os.mkdir('./results/mapped')

    if os.path.isfile('./results/mapping_stats.csv') and reanalyze_bool == 0:
        mapping_stats = pd.read_csv('./results/mapping_stats.csv')
    else:
        # define mp list type
        mapping_stats = Manager().list()

        # set up worker threads
        pool = mp.Pool(processes=threadnum)
        results = [pool.apply_async(parse_sam_worker, args=(sampid, mapping_stats, sorted_refnames, data_folder, base_qual_threshold, max_indel_count, max_indel_freq, nucpos_shift, sampid_prefix, sampid_suffix, standardised_gene_name, save_unmapped_reads_bool,)) for sampid in sorted_sampid]
        output = [p.get() for p in results]

        # change it back to python list type
        mapping_stats = list(mapping_stats)
        pool.close()

        mapping_stats = pd.DataFrame.from_dict(mapping_stats)
        mapping_stats.to_csv('./results/mapping_stats.csv', index=False)

    # plot mapping stats
    # create figures folder
    if not os.path.isdir('./results/figures'):
        os.mkdir('./results/figures')

    # plot mapping proportions
    mapping_stats = mapping_stats.set_index('sampid')
    mapping_stats = mapping_stats.sort_index()
    mapping_stats['total'] = mapping_stats[['fail_qual', 'mapped', 'unmapped']].sum(axis=1)
    mapping_stats['mapped_prop'] = mapping_stats['mapped']/mapping_stats['total']
    mapping_stats['unmapped_prop'] = mapping_stats['unmapped']/mapping_stats['total']
    mapping_stats['fail_qual_prop'] = mapping_stats['fail_qual']/mapping_stats['total']

    with plt.style.context("default"):
        ax = mapping_stats[['mapped_prop', 'fail_qual_prop', 'unmapped_prop']].plot.bar(figsize=(18, 10), stacked=True, color=['#8da0cb', '#66c2a5', '#fc8d62'])
        fig = ax.get_figure()
        ax.tick_params(axis='y', which='major', labelsize=16)
        ax.tick_params(axis='x', which='major', labelsize=12)
        ax.set_ylabel('Frequency')
        ax.set_ylim((0., 1.))
        ax.yaxis.label.set_fontsize(16)
        ax.set_xlabel('Samples')
        ax.xaxis.label.set_fontsize(16)
        ax.tick_params(axis='both', which='major', labelsize=8)

        plt.legend(bbox_to_anchor=(0.5, 1.03),
                   loc='center',
                   ncol=3, facecolor='white',
                   labels=['Mapped', 'Low quality', 'Unmapped'], fontsize=14)

        plt.tight_layout()
        plt.savefig('./results/figures/mapping_results.pdf', bbox_inches='tight', pad_inches=0)
        if plt_show == 1:
            plt.show()
        # close figure
        plt.close(fig)

    # plot number of reads
    with plt.style.context("default"):
        ax = mapping_stats['total'].plot.bar(figsize=(18, 10), color='#8da0cb')
        fig = ax.get_figure()
        ax.tick_params(axis='y', which='major', labelsize=16)
        ax.tick_params(axis='x', which='major', labelsize=12)
        ax.set_ylabel('Number of reads')
        ax.set_ylim((1, 10**np.ceil(np.log10(mapping_stats['total'].max()))))
        ax.set_yscale('log')
        ax.yaxis.label.set_fontsize(16)
        ax.set_xlabel('Samples')
        ax.xaxis.label.set_fontsize(16)
        ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.savefig('./results/figures/read_counts.pdf', bbox_inches='tight', pad_inches=0)
        if plt_show == 1:
            plt.show()
        plt.close(fig)

    print ("...done.")

    return mapping_stats

def tally_bases_worker(sampid, codon_table, all_bases, no_map_codon_results, ignore_indels_bool, reanalyze_bool):

    if ignore_indels_bool > 0:
        all_bases = ['A', 'T', 'G', 'C']
    else:
        all_bases = ['A', 'T', 'G', 'C', '-', 'N']

    map_nuc_results_fname = './results/map_nuc_results_%s.csv'%(sampid)
    map_codon_results_fname = './results/map_codon_results_%s.csv'%(sampid)
    mapped_reads_fname = './results/mapped/mapped_%s.csv'%(sampid)

    if reanalyze_bool == 0 and os.path.isfile(map_nuc_results_fname):
        return

    map_nuc_results = {}
    map_codon_results = {}

    # get all mapped reads for sampid
    try:
        mapped_reads_df = pd.read_csv(mapped_reads_fname)
    except:
        print ('No mapped reads found for {}'.format(sampid))
        return

    for r, row in mapped_reads_df.iterrows():
        refname = row["gene"]
        start_pos = row["start_nucaln_pos"]
        end_pos = row["end_nucaln_pos_plus1"]
        read_nuc_seq = row["seq"]

        try:
            insertion_ridx = ast.literal_eval(row["insertion"])
        except:
            insertion_ridx = []
        insertion_memory = {}

        refpos_nuc_seq = [] # (don't care about insertions for now)

        ref_pos = start_pos-1
        for idx, read_base in enumerate(read_nuc_seq):
            # skip if ignore indels
            if ignore_indels_bool > 0:
                if idx in insertion_ridx or re.search("[^ATGC]", read_base, re.I):
                    continue

            # read base in an insertion
            if idx in insertion_ridx:
                try:
                    print_pos = insertion_memory[ref_pos][-1] + 0.00001
                    insertion_memory[ref_pos].append(print_pos)
                except:
                    print_pos = ref_pos + 0.00001
                    insertion_memory[ref_pos] = [print_pos]
            else:
                ref_pos += 1
                print_pos = ref_pos
                refpos_nuc_seq.append(read_base)

            # add nuc count
            try:
                map_nuc_results[refname][print_pos][read_base] += 1
            except:
                try:
                    map_nuc_results[refname][print_pos] = {base:1 if base == read_base else 0 for base in all_bases}
                except:
                    map_nuc_results[refname] = {print_pos:{base:1 if base == read_base else 0 for base in all_bases}}

        if no_map_codon_results > 0: # skip codon counts
            continue

        # count codon (don't care about insertions for now)
        refpos_read_len = len(refpos_nuc_seq)
        for _b in range(refpos_read_len):
            if _b+3 >= refpos_read_len:
                continue

            ref_pos = start_pos+_b
            read_codon = "".join(refpos_nuc_seq[_b:_b+3])

            if re.search('[^ATGC]', read_codon):
                continue

            # add nuc count (except for gaps)
            try:
                map_codon_results[refname][ref_pos][read_codon] += 1
            except:
                try:
                    map_codon_results[refname][ref_pos] = {codon:1 if codon == read_codon else 0 for codon in codon_table.keys()}
                except:
                    map_codon_results[refname] = {ref_pos:{codon:1 if codon == read_codon else 0 for codon in codon_table.keys()}}

    map_nuc_results = pd.DataFrame.from_dict({(i,j): map_nuc_results[i][j]
                                          for i in map_nuc_results.keys()
                                          for j in map_nuc_results[i].keys()},
                                         orient='index')
    map_nuc_results.index.names = ['Gene', 'Position']

    # calculate coverage
    map_nuc_results['Coverage'] = map_nuc_results[all_bases].sum(axis=1)

    # get consensus
    nuc_consensus = list(map_nuc_results[all_bases].idxmax(axis=1))
    map_nuc_results['Consensus'] = nuc_consensus

    # print to csv
    map_nuc_results = map_nuc_results.sort_index()
    map_nuc_results.to_csv(map_nuc_results_fname)

    # reorder map_codon_results and print csv
    if len(map_codon_results) > 0:
        map_codon_results = pd.DataFrame.from_dict({(i,j): map_codon_results[i][j]
                                              for i in map_codon_results.keys()
                                              for j in map_codon_results[i].keys()},
                                             orient='index')
        map_codon_results.index.names = ['Gene', 'Position']

        # print to csv
        map_codon_results = map_codon_results.sort_index()
        map_codon_results.to_csv(map_codon_results_fname)

    return

def tally_bases(sorted_sampid, threadnum=4, no_map_codon_results=0, ignore_indels_bool=1, reanalyze_bool=0):
    """
    Tally base and codon counts based on saved mapped reads
    """
    print ("\nTallying bases and codons...")
    # parameters
    codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
            }

    pool = mp.Pool(processes=threadnum)
    results = [pool.apply_async(tally_bases_worker, args=(sampid, codon_table, no_map_codon_results, ignore_indels_bool, reanalyze_bool,)) for sampid in sorted_sampid]
    output = [p.get() for p in results]
    pool.close()

    print ("..done.")
    return

def read_mapping_bt2(sorted_sampid, data_folder, gene_len, reffasta, threadnum=4, bowtie2_fpath="bowtie2"):
    """
    Read mapping using bowtie2
    """
    # index reference sequence
    print ("\nIndex reference sequences...")
    ref_key = re.sub('(^.+/|\.[^\.+]$)', '', reffasta)
    cmd = ['%s-build'%(bowtie2_fpath), reffasta, ref_key] # map to barcode ref fasta
    subprocess.call(cmd)

    # align sequences
    if not os.path.isdir("%s/align"%(data_folder)):
        os.mkdir("%s/align"%(data_folder))

    timestamp = time.strftime("%Y%m%d-%H%M%S") # timestamp for log file
    print ("\nMapping reads with bowtie2...")
    for sampid in sorted_sampid:
        input_fpath = '%s/merge/%s.extendedFrags.fastq.gz'%(data_folder, sampid)
        if not os.path.isfile(input_fpath):
            raise Exception("Unable to find input FASTQ file %s for mapping"%(input_fpath))

        if not os.path.isfile('%s/align/%s.sam.gz'%(data_folder, sampid)):
            # mapping with bt2
            with open('./data/bt_aln_%s.log'%(timestamp), 'a') as output:
                output.write('%s\n'%(sampid))
                cmd = [bowtie2_fpath,
                       '-x', ref_key,
                       '-X', str(max(gene_len.values())),
                       '-k', '2',
                       '--very-sensitive-local',
                       '-p', str(threadnum),
                       '-U', input_fpath,
                       '-S', '%s/align/%s.sam'%(data_folder, sampid)]
                subprocess.call(cmd, stderr=subprocess.STDOUT, stdout=output)
                output.write('\n')

            # gzip output sam file
            cmd = ['gzip', '%s/align/%s.sam'%(data_folder, sampid)]
            subprocess.call(cmd)

    print ("...done.")
    return

def merge_reads(dat_df, data_folder, sorted_sampid, flash_fpath="flash"):
    if not os.path.isdir('{}/merge'.format(data_folder)):
         os.mkdir('{}/merge'.format(data_folder))

    if os.path.isfile("{}/merge_df.csv".format(data_folder)):
        merge_df = pd.read_csv("{}/merge_df.csv".format(data_folder))
    else:
        merge_df = []

        for sampid in sorted_sampid:

            r1_fpath = '{}/trimmed/{}_R1.paired.fastq.gz'.format(data_folder, sampid)
            r2_fpath = '{}/trimmed/{}_R2.paired.fastq.gz'.format(data_folder, sampid)

            # merge trimmed, paired reads
            cmd = [flash_fpath, '--output-directory={}/merge'.format(data_folder),
                   '--output-prefix={}'.format(sampid),
                   '-m', '10',
                   '-M {}'.format(int(dat_df.loc[sampid]['end_pos'].max())),
                   '--compress',
                   r1_fpath, r2_fpath]

            p = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
            output, err = p.communicate()
            if p.returncode == 1:
                raise Exception('Error merging {}'.format(sampid))

            total = int(re.search('Total pairs:\s+(\d+)', output.decode('utf-8')).group(1))
            combined = int(re.search('Combined pairs:\s+(\d+)', output.decode('utf-8')).group(1))
            merge_df.append({'sampid':sampid, 'total':total, 'combined':combined})

        merge_df = pd.DataFrame.from_dict(merge_df)
        merge_df['percent_combined'] = merge_df['combined']/merge_df['total']
        merge_df = merge_df[['sampid', 'total', 'combined', 'percent_combined']]
        merge_df.to_csv("{}/merge_df.csv".format(data_folder), index=False)

    print ("Overall mean prop of pairs merged: {:.2%} (SD: {:.2%})".format(np.mean(merge_df['percent_combined']), np.std(merge_df['percent_combined'], ddof=1)))
    return merge_df.set_index('sampid')

def trim_raw_fastq(dat_df, sorted_sampid, data_folder, trimmomatic_fpath, threadnum=4, single_end=0, adapter_to_clip=None):
    """
    Trimming raw FASTQ sequences with trimmomatic
    """

    print ("\nTrimming raw FASTQ sequences...")
    if not os.path.isdir('{}/trimmed'.format(data_folder)):
        os.mkdir('{}/trimmed'.format(data_folder))

    if os.path.isfile('{}/trimmomatic_stats.csv'.format(data_folder)):
        trimmed_df = pd.read_csv('{}/trimmomatic_stats.csv'.format(data_folder))
    else:
        trimmed_df = []
        for sampid in sorted_sampid:
            sampid_fdf = dat_df.loc[sampid]

            if single_end > 0: # single end modde
                if os.path.isfile('{}/trimmed/{}_trimmed.fastq.gz'.format(data_folder, sampid)):
                    continue
                else:
                    cmd = ["trimmomatic", "SE", "-phred33",
                           "-threads", str(threadnum),
                           sampid_fdf.loc['R1']['fpath'],
                           '{}/trimmed/{}_trimmed.fastq.gz'.format(data_folder, sampid)]
                    if adapter_to_clip != None:
                        cmd += ["ILLUMINACLIP:{}/adapters/{}:2:3:10".format(trimmomatic_fpath, adapter_to_clip)]
                    cmd += ["LEADING:3", "TRAILING:3",
                    "MAXINFO:40:0.4", "CROP:{}".format(int(sampid_fdf['end_pos'].max())),
                    "MINLEN:30"]

                    p = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
                    output, err = p.communicate()
                    if p.returncode == 1:
                        raise Exception('Error trimming {}'.format(sampid))

                    total = int(re.search('Input Reads: (\d+)', output.decode('utf-8')).group(1))
                    survived = int(re.search('Surviving: (\d+)', output.decode('utf-8')).group(1))
                    dropped = int(re.search('Dropped: (\d+)', output.decode('utf-8')).group(1))

                    trimmed_df.append({'sampid':sampid, 'total':total, 'survived':survived,'dropped':dropped})
                    #print (trimmed_df)
                    #raise Exception

            else:
                # trim low quality ends with trimmomatic
                ## settings: LEADING/TRAILING

                ## MAXINFO:<targetLength>:<strictness>
                ## CROP:<length> - crop AFTER <length> nt from 5' end
                ## MINLEN:<length> - drop reads that are < <length> nt

                ## ILLUMINACLIP:
                ## <fastaWithAdaptersEtc>: NexteraPE-PE.fa
                ## <seed mismatches>: 2 mismatch counts allowed
                ## <palindrome clip threshold>: 30 (accuracy of match between paired adapter ligated reads)
                ## <simple clip threshold>: 10 accuracy of match between any adapter
                ## <minAdapterLength>: 2 min adapter length to detect (can be low for palindrome mode; low false positive rate) TruSeq3-PE-2.fa:2:30:10:2
                ## <keepBothReads>: false - drop reverse read

                if os.path.isfile('{}/trimmed/{}_R1.paired.fastq.gz'.format(data_folder, sampid)):
                    continue
                else:
                    cmd = ["trimmomatic".format(trimmomatic_fpath), "PE", "-phred33",
                           "-threads", str(threadnum),
                           sampid_fdf.loc['R1']['fpath'],
                           sampid_fdf.loc['R2']['fpath'],
                           '{}/trimmed/{}_R1.paired.fastq.gz'.format(data_folder, sampid),
                           '{}/trimmed/{}_R1.unpaired.fastq.gz'.format(data_folder, sampid),
                           '{}/trimmed/{}_R2.paired.fastq.gz'.format(data_folder, sampid),
                           '{}/trimmed/{}_R2.unpaired.fastq.gz'.format(data_folder, sampid)]
                    if adapter_to_clip != None:
                        cmd += ["ILLUMINACLIP:{}/adapters/{}:2:3:10".format(trimmomatic_fpath, adapter_to_clip)]
                    cmd += ["LEADING:3", "TRAILING:3",
                    "MAXINFO:40:0.4", "CROP:{}".format(int(sampid_fdf['end_pos'].max())),
                    "MINLEN:30"]
                    #,
                    p = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
                    output, err = p.communicate()
                    if p.returncode == 1:
                        raise Exception('Error trimming {}'.format(sampid))

                    #print  (output.decode('utf-8'))
                    total_pairs = int(re.search('Read Pairs: (\d+)', output.decode('utf-8')).group(1))
                    both_survive = int(re.search('Both Surviving: (\d+)', output.decode('utf-8')).group(1))
                    forward_survive = int(re.search('Forward Only Surviving: (\d+)', output.decode('utf-8')).group(1))
                    reverse_survive = int(re.search('Reverse Only Surviving: (\d+)', output.decode('utf-8')).group(1))
                    dropped = int(re.search('Dropped: (\d+)', output.decode('utf-8')).group(1))

                    trimmed_df.append({'sampid':sampid, 'total_pairs':total_pairs, 'both':both_survive,
                                       'forward':forward_survive, 'reverse':reverse_survive, 'dropped':dropped})

        trimmed_df = pd.DataFrame.from_dict(trimmed_df)
        if single_end > 0:
            trimmed_df['survived_prop'] = trimmed_df['survived']/trimmed_df['total']
            trimmed_df = trimmed_df[['sampid', 'total', 'survived', 'survived_prop', 'dropped']]
            #display (trimmed_df)
            print ("Overall reads surived: {:.3%} (SD: {:.3%})".format(np.mean(trimmed_df['survived_prop']), np.std(trimmed_df["survived_prop"])))
        else:
            trimmed_df['both_prop'] = trimmed_df['both']/trimmed_df['total_pairs']
            trimmed_df = trimmed_df[['sampid', 'total_pairs', 'both', 'both_prop', 'forward', 'reverse', 'dropped']]
            print ("Overall mean prop of pairs retained: {:.3%} (SD: {:.3%})".format(np.mean(trimmed_df['both_prop']), np.std(trimmed_df['both_prop'], ddof=1)))

        trimmed_df.to_csv("{}/trimmomatic_stats.csv".format(data_folder), index=False)

    return trimmed_df.set_index('sampid')

def pretrim_fastqc(dat_df, data_folder, base_qual_threshold, fastqc_fpath="fastqc"):
    print ("\nPerform pre-trim FASTQC with minimum base quality %i (change with --base_qual_threshold if needed)..."%(base_qual_threshold))

    if not os.path.isdir('{}/fastqc'.format(data_folder)):
        os.mkdir('{}/fastqc'.format(data_folder))

    def analyse_fastqc(fname):
        fhandle = open(fname, 'r').readlines()
        module = ''
        df_dict = {}

        for line in fhandle:
            if re.search('>>Per base sequence quality', line):
                module = 'pbsq'
                continue
            if re.search('>>Adapter Content', line):
                # check whether there is adapter content
                module = 'ac'
                continue
            if re.search('>>Per sequence quality scores', line):
                module = 'psqs'
                continue
            if module != '' and re.search('>>END_MODULE', line):
                df_dict[module] = pd.DataFrame(df)
                module = ''

            if module != '':
                if re.search('^#', line):
                    headers = re.sub('#', '', line.strip()).split('\t')
                    df = {h:[] for h in headers}
                else:
                    values = line.strip().split('\t')
                    for _h, h in enumerate(headers):
                        df[h].append(values[_h])

        # % of reads with acceptable phred score > base_qual_threshold
        df = df_dict['psqs']
        df['Quality'] = pd.to_numeric(df['Quality'], errors='coerce')
        df['Count'] = pd.to_numeric(df['Count'], errors='coerce')
        percent_abv_qualthres = sum(df[df['Quality']>=base_qual_threshold]['Count'])/sum(df['Count'])

        # max position where median phred score > base_qual_threshold
        df = df_dict['pbsq']
        df['Median'] = pd.to_numeric(df['Median'], errors='coerce')

        try:
            start_pos = int(re.search('^(\d+)-*', df[df['Median']>=base_qual_threshold]['Base'].iloc[0]).group(1))
        except:
            start_pos = None
        try:
            end_pos = int(re.search('-*(\d+)$', df[df['Median']>=base_qual_threshold]['Base'].iloc[-1]).group(1))
        except:
            end_pos = None

        df = df_dict['ac'].drop('Position', axis=1)
        ac = {}
        for adapter in list(df):
            val = np.sum(pd.to_numeric(df[adapter], errors='coerce'))/len(df[adapter])
            if adapter in ac:
                if val > ac[adapter]:
                    ac[adapter] = val
            else:
                if val > 0.:
                    ac[adapter] = val

        return percent_abv_qualthres, ac, start_pos, end_pos

    adapters_to_rm = {}
    for (sampid, read) in dat_df.index:

        # run fastqc
        if not os.path.isfile('{}/fastqc/{}_{}_fastqc.html'.format(data_folder, sampid, read)):
            fpath = dat_df.loc[sampid, read]['fpath']
            if re.search('Darwin', platform.platform()):
                cmd = ['zcat', '<', fpath, '|', fastqc_fpath, 'stdin:{}_{}'.format(sampid, read), '--outdir={}/fastqc'.format(data_folder)]
            else:
                cmd = ['zcat', fpath, '|', fastqc_fpath, 'stdin:{}_{}'.format(sampid, read), '--outdir={}/fastqc'.format(data_folder)]
            subprocess.call(' '.join(cmd), shell=True)

            ## unzip files
            cmd = ['unzip', '{}/fastqc/{}_{}_fastqc.zip'.format(data_folder, sampid, read),
                   '-d', '{}/fastqc/'.format(data_folder)]
            subprocess.call(cmd)

        # analyse results
        percent_abv_qualthres, ac, start_pos, end_pos = analyse_fastqc('{}/fastqc/{}_{}_fastqc/fastqc_data.txt'.format(data_folder, sampid, read))
        if len(ac) > 0:
            for k, v in ac.items():
                if k in adapters_to_rm:
                    if v > adapters_to_rm[k]:
                        adapters_to_rm[k] = v
                else:
                    if v > 0.:
                        adapters_to_rm[k] = v

        dat_df.at[(sampid, read), 'percent_abv_qualthres'] = percent_abv_qualthres
        dat_df.at[(sampid, read), 'start_pos'] = start_pos
        dat_df.at[(sampid, read), 'end_pos'] = end_pos

    if len(adapters_to_rm) > 0:
        print ('\n#-- Presence of adapter sequence (max. proportion of reads) --#')
        for k, v in adapters_to_rm.items():
            print ("{}: {:.2f}%".format(k, v))

    return dat_df.sort_index()

def generate_raw_fastq_df(sorted_sampid, data_folder,):
    print ('\nGenerate dataframe of input FASTQ files...')
    # check for raw subdirectory in data_folder
    if not os.path.isdir("%s/raw"%(data_folder)):
        raise Exception("Raw FASTQ files must be stored in %s/raw subdirectory."%(data_folder))

    # get path to raw FASTQ files sorted by read direction
    dat_df = []
    for sampid in sorted_sampid:
        both_reads_found = 0
        for fname in os.listdir('{}/raw'.format(data_folder)):
            if re.search('{}.*\.gz'.format(sampid), fname):

                try:
                    read = re.search('R(1|2)', fname, re.I).group().upper()
                except:
                    raise Exception("Paired read files must be denoted as \'R1\' or \'R2\'  in their file names.")
                dat_df.append({'sampid':sampid, 'read':read, 'fpath':'{}/raw/{}'.format(data_folder, fname)})
                both_reads_found += int(read[-1])

        # fastq raw datafiles not found
        if (both_reads_found) < 3:
            raise Exception("Missing %s file(s) for sample %s."%("both R1 and R2" if both_reads_found == 0 else "R%i"%(both_reads_found), sampid))

    dat_df = pd.DataFrame.from_dict(dat_df)
    dat_df = dat_df.set_index(['sampid', 'read'])
    dat_df = dat_df.sort_index()
    return dat_df

def initialisation(cds_coords, reffasta, laeb_lib_fpath, nucpos_shift=None):
    """
    Generate protein open reading frame dataframe, sorted gene names, and any shifts of nucletotide positions to CDS.
    """
    print ('\nInitialising CDS coordinates...')

    # load flu common libraries
    fc = SourceFileLoader('fc', "%s/flu_common.py"%(laeb_lib_fpath)).load_module()

    # read cds_coords
    if re.search("^(H3N2_Bris07|H1N1pdm09_Cali09)$", cds_coords):
        cds_coords_df = pd.read_csv("%s/reference/CDS_%s.csv"%(laeb_lib_fpath, cds_coords))
        cds_coords_df['protein'] = cds_coords_df['protein'].fillna("NA")
        cds_coords_df = cds_coords_df.set_index(['gene', 'protein']).sort_index()
        if cds_coords == "H1N1pdm09_Cali09":
            nucpos_shift = pd.read_csv("%s/reference/CDS_shift_%s.csv"%(laeb_lib_fpath, cds_coords))
            nucpos_shift = nucpos_shift.set_index("gene").T.to_dict('records')[0]
        else:
            nucpos_shift = None
    else:
        # non-preset sequences
        cds_coords_df = pd.read_csv(cds_coords)
        cds_coords_df['protein'] = cds_coords_df['protein'].fillna("NA")
        cds_coords_df = cds_coords_df.set_index(['gene', 'protein']).sort_index()
        try:
            nucpos_shift = pd.read_csv(nucpos_shift)
            nucpos_shift = nucpos_shift.set_index("gene").T.to_dict('records')[0]
        except:
            nucpos_shift = None

    # parameters
    codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }

    # read reference fasta
    ref_fdat = fc.parsefasta(reffasta, rtype='nuc', aln=0)
    sorted_refnames = sorted(ref_fdat.keys())

    # get dataframe of protein ORF in gene
    gene_to_proteinorf = []
    print ('\nCheck translated protein sequences...')
    for (gene, protein), row in cds_coords_df.iterrows():

        coordinates = re.findall("\d+,\s*[\d\'\*]+", str(row.coordinates))

        try:
            refseq = list(ref_fdat[gene])
        except:
            continue

        refnum_dict = {"nuc":[], "aa":[], "frame":[]}

        cdr_nucrange = []
        for coor in coordinates:
            coor = coor.replace(" ", "").split(",")
            start, end = coor
            start = int(start)
            if re.search('\*', end):
                end = "*"
            else:
                end = int(end)

            start -= 1
            if end != '*':
                cdr_nucrange += list(range(start, end))
            else:
                cdr_nucrange += list(range(start, len(refseq)))

        refaaseq = []
        refaa_num = 0
        cdr_refseq = [nuc for i, nuc in enumerate(refseq) if i in cdr_nucrange]

        for _i in range(0, len(cdr_refseq), 3):
            codon = cdr_refseq[_i:_i+3]
            aa = codon_table["".join(codon).upper()]

            refaaseq.append(aa)
            refaa_num += 1

            for _j in range(3):
                refnum_dict["aa"].append(refaa_num)
                # remove nucpos shift for cdf numbering if necessary
                try:
                    num_shift = nucpos_shift[gene]
                except:
                    num_shift = 0
                refnum_dict["nuc"].append(cdr_nucrange[_i+_j]+1-num_shift)
                refnum_dict['frame'].append(_j+1)

            if aa == '*':
                break

        refnum_dict["gene"] = [gene]*len(refnum_dict["nuc"])
        refnum_dict["protein"] = [protein]*len(refnum_dict["nuc"])
        refnum_dict = pd.DataFrame(refnum_dict)
        gene_to_proteinorf.append(refnum_dict)

        # for checking protein sequence
        print (protein, "".join(refaaseq))

    gene_to_proteinorf = pd.concat(gene_to_proteinorf)
    gene_to_proteinorf = gene_to_proteinorf.set_index(['gene', 'nuc', 'protein']).sort_index()

    # determine CDS gene length
    gene_len = {}
    for gene in sorted_refnames:
        cds_nuc_positions = list(gene_to_proteinorf.loc[gene].index.get_level_values(0))
        gene_len[gene] = max(cds_nuc_positions)-min(cds_nuc_positions)+1

    return gene_to_proteinorf, gene_len, sorted_refnames, nucpos_shift
