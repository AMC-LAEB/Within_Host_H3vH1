# Custom scripts and files for Han et al. (2021)  

## Citation  
Within-host evolutionary dynamics of seasonal and pandemic human influenza A viruses in young children  
Alvin X. Han, Zandra C. Felix Garza, Matthijs R. A. Welkers, René M. Vigeveno, Tran Nhu Duong, Le Thi Quynh Mai, Pham Quang Thai, Dang Dinh Thoang, Tran Thi Ngoc Anh, Ha Manh Tuan, Nguyen Thanh Hung, Le Quoc Thinh, Le Thanh Hai, Hoang Thi Bich Ngoc, Kulkanya Chokephaibulkit, Pilaipan Puthavathana, Nguyen Van Vinh Chau, Nghiem My Ngoc, Nguyen Van Kinh, Dao Tuyet Trinh, Tran Tinh Hien, Heiman F. L. Wertheim, Peter Horby, Annette Fox, H. Rogier van Doorn, Dirk Eggink, Menno D. de Jong, Colin A. Russell  
_bioRxiv_ 2021.03.20.436248; doi: https://doi.org/10.1101/2021.03.20.436248  

For any questions/queries with regards to the scripts/notebooks, please email [Alvin Han](x.han@amsterdamumc.nl).  

## Sequencing data
All raw sequence data have been deposited at NCBI sequence read archive under BioProject Accession number PRJNA722099.  

## Installation  
You will need Python 3 and several dependencies to run the Jupyter notebooks/scripts and reproduce our analyses. We have only ran our scripts/notebooks on MacOS. Linux should work just as well. We recommend using ```conda``` to download and install all dependencies:  

1. Clone this repository.  
2. ```cd Within_Host_H3vH1```  
3. ```conda env create -f environment.yml```  
4. Prior to running the notebooks/scripts: ```conda activate h3vh1```  

## Running the analyses  
1. Separate folders were created for the different subtypes (```H3N2``` and ```H1N1pdm09```). Move the raw sequencing data files thtat you have downloaded into the respective ```<subtype>/data/raw``` subfolder.  
3. *Read mapping* and *variant calling* were performed by running the ```SEA_<subtype>_01_NGS_Mapping.ipynb``` jupyter notebook.  
4. For *H3N2*, we were able to analyse the precision of our variant calls (```SEA_H3N2_02a_Variant_Precision.ipynb```).  
5. Phylogenetic analyses of the concatenated sequences were performed to identify samples with potential mixed infections and/or were contaminated (```SEA_<subtype>_02b*_Phylogenetic_Analyses.ipynb```). You may use ```R``` script, available in ```H3N2/scripts/ggtree001-2.R``` to create the Supplementary figures S9 and S10. You will need to change the working directory and output filename accordingly as stipulated in the script.    
6. *Within-host genetic diversity* and *evolutionary dynamics* analyses were performed using ```SEA_<subtype>_03_Genetic_Diversity_Analyses.ipynb```.
7. *Haplotype reconstruction* was performed by running ```H3N2/scripts/reconstruct_haplotype_wrapper.py``` using the following commands (for *H3N2* as an example):  
```
cd H3N2  
python ./scripts/reconstruct_haplotype_wrapper.py -v ./results/variant_call_MinCoV100_MinProp0.02_MinFreq0_ErrTol0.01.csv -f ./results/metadata_wo_mixed_infections.csv -m ./results/mapped --DMNFuncs_fpath ../python_lib/DMNFuncs.py  
```
8. For *H1N1pdm09*, *transmmission bottleneck size* is estimated by running ```SEA_H1N1pdm09_04_Transmission_Analyses.ipynb```.
9. *Within-host simulations* was performed by running ```simulations/scripts/WHDEL.py```. A bash script is enclosed to run multiple simulations in parallel:  
```
cd simulations
bash WHDEL_wrapper.sh
```
