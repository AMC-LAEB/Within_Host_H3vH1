from Bio import SeqIO
from os.path import expanduser
import pandas as pd
import numpy as np
import re
import os
import ete3
import subprocess

# create protein-guided nucleotide alignment
codon_table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',}
all_aas = list(sorted(set(codon_table.values()))) + ["X"]

def translate_codon_nuc_to_aa(sequence):
    pro_seq = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        try:
            aa = codon_table[codon.upper()]
        except:
            aa = "X"
        pro_seq.append(aa)
    return "".join(pro_seq)

def sort_passage_annotation(annotation):
    """
    Sort passage annotations
    """
    if annotation == "":
        return "u"

    # check if annotation has been previously sorted
    previously_sorted_passage = pd.read_csv(expanduser("~/Dropbox/LAEB-NGS/reference/passage_annotations.csv"))
    previously_sorted_passage = previously_sorted_passage.drop_duplicates('passage')
    previously_sorted_passage = previously_sorted_passage.set_index("passage")

    try:
        return previously_sorted_passage.loc[annotation, "ptype"]
    except:
        #print ('previously unsorted! - %s'%(annotation))
        annotation = annotation.lower()
        # unknown
        if re.search("(^(nvd\d+|llc|nvd\d+/e\d+|llc-mk\d+_\d+pass|rna|nc2|llc-mk\d+-mdck\d+/siat\d+|mdck33016pf_2|e7_is2_cln25|a/alborz/153084/2019|r-m1x1\,mdck1|i-cs|rnacs|passage_details:_na_extract|\.|passaged_in_mouse_lung|mek|mus_musculus|different_isolation_sources|\d+_mouse_passages|mek\d+|passage_details:_x\,_*mdck\d+|dee_\d+_challenge_virus_stock|llc-mk\d+_\d+_\+\d+|y|na|\d+_day_\d+|passage_details\:_(x|p\d+))$|rmix|vero|tmk\d+|\?|unknown|^x|^rm|qmc2|[pr]h*mk|^(p(\d+|x)|pi)$|rii|caco|^r(\d+|x)|[^a-z]r(\d+|x))", annotation):
            print (annotation, "u")
            return "u"

        passage_type = []

        if re.search("^(am-*\d+|s31n_amantadine_resistance_passage_details:_e1)$", annotation):
            print (annotation, "e")
            return "e"

        if re.search("^(cs|cell|c|c_\d+|c_\d+_\+\d+)$", annotation):
            print (annotation, "c")
            return "c"

        if re.search("^(not_passaged|passage_details:_no_passage|clincal_specimen|lung-\d+|nasal_swab|nasopharyngeal_epithelium)$", annotation):
            print (annotation, "o")
            return "o"

        # original
        if re.search("(^(org$|or|p0_\(h3n2\))$|ori|o[rt]*i*gini*al|initial|clinical|direct|primary)", annotation):
            print (annotation, "o")
            passage_type.append("o")

        # cell-based
        if re.search("(qmc|beas-2b|ax-4|hck|mdck|siat|^(c|s|m)(\d+|x)|[^a-z](c|s|m)(\d+|x))", annotation):
            print (annotation, "c")
            passage_type.append("c")

        # egg-based
        if re.search("(^e$|spfck|egg|amniotic|am\d+al\d+|^(e)(\d+|x)|[^a-z](e)(\d+|x))", annotation):
            print (annotation, "e")
            passage_type.append("e")

        if len(passage_type) == 0:
            print (annotation)
            raise Exception
        """else:
            print (annotation, passage_type)"""

        if "e" in passage_type:
            return "e"
        elif "c" in passage_type:
            return "c"
        else:
            return "o"

def safe_ln(x):
    """
    safe natural logarithm
    """
    if isinstance(x, np.ndarray):
        if (x > 0).all() == True:
            return np.log(x)
        else:
            x[x>0] = np.log(x[x>0])
            x[x==0] = -np.Infinity
            return x
    else:
        if x > 0:
            return np.log(x)
        else:
            return -np.Infinity

def toYearFraction(date, date_format='%Y-%m-%d', discard_year_only=0):
    """
    Convert string date to decimal date.
    """
    from datetime import datetime as dt
    import time

    # partial date strings on GISAID syntax
    if re.search("\d+_\(Month_and_day_unknown\)", date) or re.search("^\d\d\d\d$", date):
        if discard_year_only > 0:
            return False
        else:
            return float(re.search("(\d+)_\(Month_and_day_unknown\)", date).group(1)) + 0.5

    if re.search("\d+-\d+_\(Day_unknown\)", date) or re.search("^\d\d\d\d-\d\d$", date):
        year, month = re.search("^(\d+)-(\d+)", date).group(1,2)
        date = "{}-{}-15".format(year, month)

    date = dt.strptime(date, date_format)
    def sinceEpoch(date): # returns seconds since epoch
        epoch = dt(1970, 1, 1)
        t = dt(date.year, date.month, date.day)
        diff = t-epoch
        return diff.days * 24 * 3600 + diff.seconds
        #return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def filter_same_sname_by_passage_type(fdat_df, ptype_preference=["o", "c", "e", "u"]):
    try:
        fdat_df = fdat_df.reset_index()
    except:
        pass

    iid_to_remove = []
    for sname in set(fdat_df['sname']):
        sname_fdat_df = fdat_df[fdat_df['sname']==sname]
        if len(set(sname_fdat_df['iid'])) > 1:
            # prefer "o" > "c" > "e" > "u"
            for preferred_ptype in ptype_preference:
                p_sname_fdat_df = sname_fdat_df[sname_fdat_df['ptype']==preferred_ptype]
                if len(p_sname_fdat_df) > 0:
                    iid_to_remove += (list(set(sname_fdat_df['iid'])-set([np.random.choice(list(p_sname_fdat_df['iid']))])))
                    break

    print ("Number of input sequences: %i"%(len(fdat_df)))
    fdat_df = fdat_df[~fdat_df['iid'].isin(iid_to_remove)]
    fdat_df = fdat_df.set_index("iid")
    print ("Number of output sequences: %i"%(len(fdat_df)))
    return fdat_df

def generate_gisaid_fasta_df(fname, rtype="nuc", ambiguous_tol=0.01, len_tol=0.9):
    """
    Generate pandas dataframe for sequences downloaded from GISAID
    """
    fdat_df  = []
    standardise_gene_name = {"PB2":1, "PB1":2, "PA":3, "HA":4, "NP":5, "NA":6, "MP":7, "NS":8}
    subtype_to_influenza_gene_len = {"A":{'1-PB2': 2280, '2-PB1': 2274, '3-PA': 2151, '4-HA': 1701, '5-NP': 1497, '6-NA': 1410, '7-M': 982, '8-NS': 838}, "B":{'1-PB2': 2259, '2-PB1': 2313, '3-PA': 2178, '4-HA': 1749, '5-NP': 1683, '6-NA': 1398, '7-M': 1076, '8-NS': 1024}}

    fasta_dat = parsefasta(fname, rtype=rtype)
    print ("Number of input sequences: %i"%(len(fasta_dat)))
    amb_count = 0
    len_count = 0
    for header, sequence in fasta_dat.items():
        sname, gene, iid, date, passage, subtype = header.split("|")
        gene = "%i-%s"%(standardise_gene_name[gene], gene)
        flu_type = re.search("^(A|B)", subtype).group()

        # uncount sequences with > amb_tol of amb res
        amb_res = "n" if rtype == "nuc" else "X"
        if sequence.count(amb_res)/len(sequence) > ambiguous_tol:
            amb_count += 1
            continue
        # min sequence length
        if len(sequence) < len_tol*subtype_to_influenza_gene_len[flu_type][gene]:
            len_count += 1
            continue

        date = toYearFraction(date)
        sname = re.sub("(\(h\dn\d\)|[^a-z0-9\-\.\/_])", "", sname.lower())
        fdat_df.append({"sname":sname, "gene":gene, "iid":iid, "subtype":subtype, "date":date, "passage":passage, "seq":sequence})

    fdat_df = pd.DataFrame.from_dict(fdat_df).set_index("iid")
    print ("Number of output sequences: %i"%(len(fdat_df)))
    print ("Removed because AMM(<%.2f)/LEN(>%.2f) = %i/%i"%(ambiguous_tol, len_tol, amb_count, len_count))
    print ("Number of unique iid: %i"%(len(set(fdat_df.index))))

    return fdat_df

def parsefasta(fname, rtype="nuc", aln=0, DasDel=0):
    """
    Parse FASTA files
    """
    if rtype == "nuc":
        # nucleotide sequences
        if aln == 0:
            if DasDel > 0:
                # not alignment
                return {rec.description:re.sub("[^atgcnd]", "n", str(rec.seq).lower().replace("-", "")) for rec in SeqIO.parse(fname, "fasta")}
            else:
                # not alignment
                return {rec.description:re.sub("[^atgcn]", "n", str(rec.seq).lower().replace("-", "")) for rec in SeqIO.parse(fname, "fasta")}
        else:
            if DasDel > 0:
                return {rec.description:re.sub("[^\-atgcnd]", "n", str(rec.seq).lower()) for rec in SeqIO.parse(fname, "fasta")}
            else:
                return {rec.description:re.sub("[^\-atgcn]", "n", str(rec.seq).lower()) for rec in SeqIO.parse(fname, "fasta")}
    elif rtype == "pro":
        # protein sequences
        if aln == 0:
            return {rec.description:re.sub("[^%s]"%("".join(all_aas)), "X", str(rec.seq).upper().replace("-", "")) for rec in SeqIO.parse(fname, "fasta")}
        else:
            return {rec.description:re.sub("[^\-%s]"%("".join(all_aas)), "X", str(rec.seq).upper()) for rec in SeqIO.parse(fname, "fasta")}
    else:
        raise Exception("Unable to recognise rtype (nuc or pro).")

def proaln_to_nucaln(pro_aln, nuc_fasta):
    """
    Create nucleotide alignment dictionary based on protein alignment
    """
    """# create protein-guided nucleotide alignment
    codon_table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',}
    all_aas = list(sorted(set(codon_table.values()))) + ["X"]
    codon_table = pd.DataFrame.from_dict(codon_table.items())
    codon_table.columns = ['codon', 'aa']"""

    # must have matching headers between pro_aln and nuc_fasta
    if set(pro_aln.keys()) != set(nuc_fasta.keys()):
        raise Exception("Protein alignment and nucleotide fasta must have identical matching headers.")

    header_to_nucaln = {}
    codon_table = codon_table.set_index("aa")
    for header, proaln in pro_aln.items():
        nucseq = nuc_fasta[header]

        regex_nuc = []
        proidx_to_fill = []
        proidx_to_codon = {}
        for i, prores in enumerate(proaln):
            if prores == "-":
                proidx_to_codon[i] = '---'
                continue
            elif prores == "X":
                encoding_codons = "[atgcn]"*3
            else:
                encoding_codons = codon_table.loc[prores, "codon"]
                if isinstance(encoding_codons, pd.Series):
                    encoding_codons = "(%s)"%("|".join(map(lambda _: _.lower(), list(encoding_codons))))
                else:
                    encoding_codons = encoding_codons.lower()

            regex_nuc.append(encoding_codons)
            proidx_to_fill.append(i)

        # add stop codon at the end (if there is one)
        encoding_codons = codon_table.loc["*", "codon"]
        encoding_codons = "(%s)*"%("|".join(map(lambda _: _.lower(), list(encoding_codons))))
        regex_nuc.append(encoding_codons)
        proidx_to_fill.append(i+1)

        # get codon_nuc_seq
        try:
            codon_nuc_seq = re.search("".join(regex_nuc), nucseq).group()
        except:
            print ("Unable to search for codon nuc seq for %s"%(header))
            continue

        # codon align nucleotide sequence
        nucaln_seq = [codon_nuc_seq[idx:idx+3] for idx in range(0, len(codon_nuc_seq), 3)]
        for idx, codon in enumerate(nucaln_seq):
            proidx_to_codon[proidx_to_fill[idx]] = codon

        header_to_nucaln[header] = "".join([proidx_to_codon[proidx] for proidx in sorted(proidx_to_codon)])

    return header_to_nucaln

def parse_pdb(pdb_fname):
    fhandle = open(pdb_fname, "r").readlines()
    pdb_coordinates = []
    three_letter_aa_conversion = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Glu":"E","Gln":"Q","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
    for k, v in three_letter_aa_conversion.items():
        del three_letter_aa_conversion[k]
        three_letter_aa_conversion[k.upper()] = v

    for line in fhandle:
        try:
            atom_num, res, chain, pos = re.search("^ATOM\s+(\d+)\s+CA\s+([A-Z][A-Z][A-Z])\s+([A-Z])\s+(\d+)", line).group(1, 2, 3, 4)
            pdb_coordinates.append({"chain":chain, "atom_num":atom_num, "pos":int(pos), "aa":three_letter_aa_conversion[res]})
        except:
            continue
    return pd.DataFrame.from_dict(pdb_coordinates).set_index(["chain", "pos"]).sort_index()

class foldx():
    """
    class of foldx analyses
    """
    def __init__(self, pdb_structure):
        self.pdb = pdb_structure
        # parse pdb for coordinates
        self.pdb_coordinates = parse_pdb(pdb_structure)
        return

    def repair_pdb(self):
        """
        Repair structures before you do any modelling with FoldX.
        RepairPDB identify those residues which have bad torsion angles,
        or VanderWaals' clashes, or total energy, and repairs them.

        Outputs for input file name PDB.pdb - PDB_Repair.pdb and
        PDB_Repair.fxout (energies of repaired residues)
        """
        # make directory for outputs
        if os.path.isdir("./foldx") == False:
            os.mkdir("./foldx")

        cmd = ["foldx", "--command=RepairPDB", "--pdb=%s"%(self.pdb), "--output-dir=./foldx"]
        subprocess.call(cmd)
        return

    def generate_mutant_inputs(self, mutant_df, aln_seq):
        """
        Align WT PDB sequence to structure and generate mutant-file input
        for foldx. Note that the current script assumes all chains in PDB
        structure are monomers of the aln_seq!
        """

        # align sequence to pdb structure
        with open("temp.fasta", "w") as output:
            for chain in set(self.pdb_coordinates.index.get_level_values(0)):
                chain_pdb_coords = self.pdb_coordinates.loc[chain].sort_index() # sort sequence by position
                chain_seq = "".join(list(chain_pdb_coords['aa']))
                output.write(">%s\n%s\n"%(chain, chain_seq))
            output.write(">ref\n%s\n"%(aln_seq))

        # align using ginsi
        cmd = ["ginsi", "temp.fasta", ">", "temp_mafft.fasta"]
        subprocess.call(" ".join(cmd), shell=True)

        # read alignment
        position_conversion = []
        for chain, chain_seq in parsefasta("temp_mafft.fasta", rtype="pro", aln=1).items():
            if chain == "ref":
                continue

            chain_pdb_coords = self.pdb_coordinates.loc[chain].copy().sort_index().reset_index()

            chain_seq = list(chain_seq)
            for idx, res in enumerate(chain_seq):
                if res == "-":
                    continue
                ref_pos = idx+1 # reference position
                chain_seq_sofar = list(filter(("-").__ne__, chain_seq[:idx+1]))
                chain_res_idx = len(chain_seq_sofar)-1
                chain_res = chain_seq_sofar[-1]
                position_conversion.append({"ref_pos":ref_pos, "chain":chain, "chain_pos":chain_pdb_coords.iloc[chain_res_idx]["pos"], "chain_res":chain_res})
        position_conversion = pd.DataFrame.from_dict(position_conversion).set_index("ref_pos").sort_index()

        # remove temporary files
        cmd = ["rm", "temp*"]
        subprocess.call(cmd)

        # write mutant-file for foldx
        mutant_input_dict = {}
        premut_list = [] # if chain_res != wt
        mutant_df = mutant_df.set_index(["comut_idx", "pos"]).sort_index()
        for comut_idx in mutant_df.index.get_level_values(0):
            idx_mutant_df = mutant_df.loc[comut_idx]
            if isinstance(idx_mutant_df, pd.Series):
                raise Exception("idx_mutant_df is a Series! XUETETY!")

            poswise_comut_list = []
            if len(idx_mutant_df) > 1:
                ref_mut_key_memory = []

            for ref_pos in sorted(set(idx_mutant_df.index)):
                wt = idx_mutant_df.loc[ref_pos, 'wt']
                mt = idx_mutant_df.loc[ref_pos, 'mt']

                try:
                    ref_pos_position_conversion = position_conversion.loc[ref_pos].copy()
                except:
                    print ("WARNING: Position %i is not present in given structure!"%(ref_pos))
                    continue

                if isinstance(ref_pos_position_conversion, pd.Series):
                    ref_pos_position_conversion = ref_pos_position_conversion.to_frame().T
                    #raise Exception("position_conversion.loc[ref_pos] is a Series! DOUBLE XUETETY!")

                ref_pos_position_conversion = ref_pos_position_conversion.set_index("chain")

                chainwise_comut_list = []
                for chain in sorted(ref_pos_position_conversion.index):
                    # do single mutation first
                    chain_pos = ref_pos_position_conversion.loc[chain, "chain_pos"]
                    chain_res = ref_pos_position_conversion.loc[chain, "chain_res"] # if chain_res != wt ???

                    if chain_res != wt:
                        premut_list.append("%s%s%i%s"%(chain_res, chain, chain_pos, wt))

                    chainwise_comut_list.append("%s%s%i%s"%(wt, chain, chain_pos, mt))

                mutant_input_dict["".join([wt, str(ref_pos), mt])] = "%s;"%(",".join(chainwise_comut_list))
                if len(idx_mutant_df) > 1:
                    ref_mut_key_memory.append("".join([wt, str(ref_pos), mt]))
                # add to poswise
                poswise_comut_list += chainwise_comut_list

            if len(set([re.search("\d+", mut).group() for mut in poswise_comut_list])) > 1:
                # now do multisite mutant if any
                mutant_input_dict[",".join(ref_mut_key_memory)] = "%s;"%(",".join(poswise_comut_list))

        return mutant_input_dict, sorted(set(premut_list))

    def available_presets(self):
        """
        Return list of available preset AA sequences for alignment to
        PDB structure
        """
        return list(parsefasta(expanduser("~/Dropbox/LAEB-NGS/reference/foldx_preset_sequences.fasta")).keys())

    def parse_foldx_output(self, pdf_fname, premut_boolean, numberOfRuns):
        foldx_output_df = pd.read_csv("./foldx/Dif_" + pdf_fname + ".fxout", skiprows=8, sep="\t")

        if premut_boolean == 0:
            mutidx = 0

        output_df = []

        for idx, row in foldx_output_df.iterrows():
            if premut_boolean == 1:
                output_df.append({"mutation":self.sorted_mutant_list[idx], "ddG":row["total energy"]})
            else:
                if idx > 0 and idx%numberOfRuns == 0:
                    mutidx += 1
                output_df.append({"mutation":self.sorted_mutant_list[mutidx], "ddG":row["total energy"]})

        return pd.DataFrame.from_dict(output_df)

    def mutate(self, mutant_df_fname, aln_seq, repair_pdb_boolean=1, numberOfRuns=5):
        """
        Perform mutation and compute ddG
        - mutant_df must be a csv file with heading: "wt,pos,mt,comut_idx"
          "comut_idx" is a unique index for mutations to be applied together.
          Here, we will mutate each mutation individually and then do another
          run combining them together.
        - aln_seq must either be a sequence or AA fasta file of the WT
        """
        aln_presets = parsefasta(expanduser("~/Dropbox/LAEB-NGS/reference/foldx_preset_sequences.fasta"), rtype="pro")
        try:
            aln_seq = aln_presets[aln_seq]
            print ("Using preset sequence as aln_seq...")
        except:
            if re.search("\.fa(sta)*$", aln_seq, re.I):
                print ("aln_seq given as a FASTA...")
                aln_seq = parsefasta(aln_seq, rtype="pro")
                aln_seq = list(aln_seq.values())[0]
            else:
                print ("raw aln_seq given...")

        # read mutant_df
        mutant_df = pd.read_csv(mutant_df_fname)
        if set(["wt", 'pos', 'mt', 'comut_idx']) != set(list(mutant_df)):
            raise Exception("mutant_df must be a csv file with heading: \"wt,pos,mt,comut_idx\"")
        # generate required mutant inputs
        mutant_input_dict, premut_list = self.generate_mutant_inputs(mutant_df, aln_seq)

        # make directory for outputs
        if os.path.isdir("./foldx") == False:
            os.mkdir("./foldx")

        # repair pdb first
        if repair_pdb_boolean == 1:
            print ("Repair PDB structure...")
            pdb_repaired_fname = re.sub("\.pdb$", "_Repair.pdb", self.pdb)
            pdb_repaired_fname = re.sub("^[^/]*/", "", pdb_repaired_fname)
            if os.path.isfile("./foldx/%s"%(pdb_repaired_fname)) == False:
                self.repair_pdb()
        else:
            print ("Warning: Input PDB structure not repaired - make sure it already was!")

        if len(premut_list) > 0:
            # write mutant input file
            with open("./foldx/individual_list.txt", "w") as output:
                output.write("%s;\n"%(",".join(premut_list)))

            # perform stability analyses
            print ("Perform pre-mutations...")
            cmd = ["foldx", "--command=BuildModel", "--pdb=%s"%(pdb_repaired_fname), "--pdb-dir=./foldx/",  "--mutant-file=./foldx/individual_list.txt", "--numberOfRuns=%i"%(numberOfRuns), "--output-dir=./foldx/"]
            subprocess.call(cmd)

            # write mutant input file
            self.sorted_mutant_list = list(mutant_input_dict.keys())
            with open("./foldx/individual_list.txt", "w") as output:
                for mut_set in self.sorted_mutant_list:
                    output.write("%s\n"%(mutant_input_dict[mut_set]))

            print ("Perform stability analyses...")
            all_output_df = []
            for run in range(numberOfRuns):
                cmd = ["foldx", "--command=BuildModel", "--pdb=%s_1_%i.pdb"%(re.sub(".pdb$", "", pdb_repaired_fname), run), "--pdb-dir=./foldx/",  "--mutant-file=./foldx/individual_list.txt", "--numberOfRuns=1", "--output-dir=./foldx/"]
                subprocess.call(cmd)

                all_output_df.append(self.parse_foldx_output("%s_1_%i"%(re.sub(".pdb$", "", pdb_repaired_fname), run), 1, 1))
            all_output_df = pd.concat(all_output_df, ignore_index=True)

        else:
            # write mutant input file
            self.sorted_mutant_list = list(mutant_input_dict.keys())
            with open("./foldx/individual_list.txt", "w") as output:
                for mut_set in self.sorted_mutant_list:
                    output.write("%s\n"%(mutant_input_dict[mut_set]))

            # perform stability analyses
            print ("Perform stability analyses...")
            cmd = ["foldx", "--command=BuildModel", "--pdb=%s"%(pdb_repaired_fname), "--pdb-dir=./foldx/",  "--mutant-file=./foldx/individual_list.txt", "--numberOfRuns=%i"%(numberOfRuns), "--output-dir=./foldx/"]
            subprocess.call(cmd)

            all_output_df = self.parse_foldx_output(re.sub(".pdb$", "", pdb_repaired_fname), 0, numberOfRuns)

        # save output to file
        all_output_df.to_csv("./foldx/foldx_output.csv", index=False)
        return all_output_df
