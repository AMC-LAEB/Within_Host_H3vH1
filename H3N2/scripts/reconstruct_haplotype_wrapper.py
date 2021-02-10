from importlib.machinery import SourceFileLoader
import multiprocessing as mp
import pandas as pd
import numpy as np
import argparse
import itertools
import os
import re

def multiproc_haplorecon_wrapper(sampid, gene, variant_call_df, DMNFuncs_fpath, mapped_reads_fpath, outdir):
    # load Dirichlet Multinomial-related functions
    dm = SourceFileLoader('DMNFuncs', DMNFuncs_fpath).load_module()

    # get sample variant calls
    try:
        sample_vcf = variant_call_df.loc[sampid].copy()
    except:
        # no variant calls for sample
        return {"sampid":sampid, "gene":None,
                "nuc_pos":None, "haplotype":None,
                "X":None, "ML_freq":None,
                "LL":None, "BIC":None}

    # get mapped reads
    sample_mapped_reads = pd.read_csv("%s/mapped_%s.csv"%(mapped_reads_fpath, sampid)).set_index("gene")

    # reconstruct haplotypes for each gene
    try:
        gene_sample_vcf = sample_vcf.loc[gene]
    except:
        # no variant found in gene for sample
        return {"sampid":sampid, "gene":gene,
                "nuc_pos":None, "haplotype":None,
                "X":None, "ML_freq":None,
                "LL":None, "BIC":None}

    if isinstance(gene_sample_vcf, pd.Series):
        gene_sample_vcf = gene_sample_vcf.to_frame().T

    gene_sample_vcf = gene_sample_vcf.drop_duplicates(["nucpos", "nuc_var"])
    nucpos_list = list(gene_sample_vcf['nucpos'])

    if len(gene_sample_vcf) == 1:
        # only one variant found - no need to apologise
        best_haploset = [tuple(gene_sample_vcf["nuc_var"].iloc[0]), tuple(gene_sample_vcf['nuc_consensus'].iloc[0])]
        best_freq = np.array([gene_sample_vcf['nuc_prop'].iloc[0],
                              gene_sample_vcf['nuc_consensus_freq'].iloc[0]/gene_sample_vcf['nuc_coverage'].iloc[0]])
        return {"sampid":sampid, "gene":gene,
                "nuc_pos":nucpos_list, "haplotype":best_haploset,
                "X":0, "ML_freq":best_freq,
                "LL":None, "BIC":None}

    # ML search already completed
    if os.path.isfile("%s/HaploRecon_%s_%s.csv"%(outdir, sampid, gene)):
        all_states_df = pd.read_csv("%s/HaploRecon_%s_%s.csv"%(outdir, sampid, gene))

        best_BIC = min(all_states_df['best_BIC'])
        best_all_states_df = all_states_df[all_states_df['best_BIC']==best_BIC]
        best_state = best_all_states_df.loc[best_all_states_df['LL'].idxmax()]
        best_freq = np.fromstring(re.sub("[\[\]]", "",  best_state["freq"]), sep=' ')
        best_haploset = eval(best_state['haplo_set'])
        best_X_boolean = best_state['X']
        best_LL = best_state['LL']

        if not os.path.isfile("%s/HaploFreqDist_%s_%s.npy"%(outdir, sampid, gene)):
            # get partial_haplotype_df from mapped reads
            gene_sample_mapped_reads = sample_mapped_reads.loc[gene]
            partial_haplotype_df = dm.get_var_haplotype_df(gene_sample_vcf, gene_sample_mapped_reads)

            # generate  ML haplotype frequency distribution by resampling based on best_freq
            print ("generate ML freq dist...")
            best_freq_dist = dm.generate_ML_haplotype_freq_dist(best_freq, best_haploset, best_X_boolean, partial_haplotype_df, outdir)
            # save distribution to file
            with open("%s/HaploFreqDist_%s_%s.npy"%(outdir, sampid, gene), "wb") as fhandle:
                np.save(fhandle, best_freq_dist)

        return {"sampid":sampid, "gene":gene,
                "nuc_pos":nucpos_list, "haplotype":best_haploset,
                "X":best_X_boolean, "ML_freq":best_freq,
                "LL":best_LL, "BIC":best_BIC}

    gene_sample_mapped_reads = sample_mapped_reads.loc[gene]

    # get partial_haplotype_df from mapped reads
    partial_haplotype_df = dm.get_var_haplotype_df(gene_sample_vcf, gene_sample_mapped_reads)

    # reconstruct most parsimonious haplotypes
    best_freq, best_haploset, best_X_boolean, best_LL, best_BIC = dm.reconstruct_haplotype(partial_haplotype_df, "%s/HaploRecon_%s_%s"%(outdir, sampid, gene), verbose=0, show_state_plt=0)

    # generate  ML haplotype frequency distribution by resampling based on best_freq
    print ("generate ML freq dist...")
    best_freq_dist = dm.generate_ML_haplotype_freq_dist(best_freq, best_haploset, best_X_boolean, partial_haplotype_df, outdir)
    # save distribution to file
    with open("%s/HaploFreqDist_%s_%s.npy"%(outdir, sampid, gene), "wb") as fhandle:
        np.save(fhandle, best_freq_dist)

    # save to dataframe
    return {"sampid":sampid, "gene":gene,
            "nuc_pos":nucpos_list, "haplotype":best_haploset,
            "X":best_X_boolean, "ML_freq":best_freq,
            "LL":best_LL, "BIC":best_BIC}

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Haplotype Reconstruction wrapper script', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-v", '--vcf', type=str, required=True, help='Variant calling CSV generated by LAEB NGS pipeline.')
    required_args.add_argument("-f", '--meta', type=str, required=True, help='Sample metadata file used as input for LAEB NGS pipeline.')
    required_args.add_argument("-m", '--mapped_reads_fpath', type=str, required=True, help='Path to directory where mapped reads are stored.')

    analyses_options = parser.add_argument_group('Analyses options')
    analyses_options.add_argument('--ncpu', type=int, default=8, help='nCPUs')
    analyses_options.add_argument('--outdir', type=str, default="./haplo_reconstruct", help='Output directory path')
    analyses_options.add_argument("--DMNFuncs_fpath", type=str, default="/data/home/xhan/python_lib/DMNFuncs.py", help="Path to DMNFuncs.py")

    params = parser.parse_args()

    # create output directory if not present
    if not os.path.isdir(params.outdir):
        os.mkdir(params.outdir)

    # read variant call df
    variant_call_df = pd.read_csv(params.vcf)
    sorted_refnames = sorted(set(variant_call_df["gene"]))
    variant_call_df = variant_call_df.set_index(["sampid", "gene"])

    # read meta df to get sampid
    meta_df = pd.read_csv(params.meta)
    sorted_sampid = sorted(set(meta_df['sampid']))

    # reconstruct haplotypes
    pool = mp.Pool(processes=params.ncpu)
    results = [pool.apply_async(multiproc_haplorecon_wrapper, args=(sampid, gene, variant_call_df, params.DMNFuncs_fpath, params.mapped_reads_fpath, params.outdir,)) for (sampid, gene) in itertools.product(sorted_sampid, sorted_refnames)]
    haplotype_reconstruction_df = [p.get() for p in results]
    haplotype_reconstruction_df = pd.DataFrame.from_dict(haplotype_reconstruction_df)
    haplotype_reconstruction_df.to_csv("%s/haplotype_reconstruct_df.csv"%(params.outdir), index=False)

    return

if __name__ == "__main__":
    main()
