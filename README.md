# The identification of tracrRNAs in genomic assemblies


### Installation
To install all dependencies required by these notebooks and scripts, create a named environment using conda and the environment.yml in this repository. If you do not have access to conda, you are doing bioinformatics, it is time to download it [here](https://docs.anaconda.com/anaconda/install/)! Otherwise, here are the steps:

1. Clone the github repository
	
	``git clone https://github.com/skDooley/TRACR_RNA.git``

1. Navigate into the directory the TRACR_RNA directory

   ``cd TRACR_RNA``

1. Create an environment from the environment.yml file

	``conda env create -f environment.yml``

1. activate the environment

	``conda activate cas``

Once these steps are complete, you are ready to launch Jupyter Notebook server and begin finding CRISPR-Cas systems and tracrRNAs.


### Usage
The python notebooks here were created to work together (See the <font style="color:red">note</font>font> below if you want to pluggin results from a different CRISPR-Cas finding tool). Within the notebooks themselves, all user input is identified by a coding cell descripton with the red background of jupyter notebook markdown warnings. All cells with a green background above the coding cells can be run without user input.

**CRISPR_CasFinding.ipynb** takes a directory (or directories) and recursively walks through the directory to identify any files with the extensions: fasta, fna, or fa. The notebook then provides code to create bash scripts to search for CRISPR arrays and proteins that are similary to those found in data/proteins/DiverseCas9s.faa. Additionally, CRISPR_CasFinding.ipynb also uses data/hmm/phi_domains.hmm to identify the 3 RuvC subdomains and the HNH domain of all the assemblies found using the hmm.

**IdentifingTracrRNAs.ipynb** uses the fasta and pickled files created with CRISPR_CasFinding.ipynb to find potential tracrRNA sequence/structural homologs. No additional input is required, but the only input cell can be changed to changed to point to a different fasta file. 


<font style="color:red">
 It is possible to use notebooks separately, but data needs to be processed into objects that are recognized by the code in the notebook(s) you want to use. If you choose to not use the CRISPR-Cas finding options, determine how to plug your results into the objects in scripts/CRISPRtools.py. Once your results are parsed into the CRISPRtools' objects, pickle the CasOperons objects. See "Process the CRISPR results" and "Read the HMM results and create protein and nucleotide sequence files" to see how I parsed  my data into the objects.
</font>


## If you use my code to find tracrRNAs, please cite the paper [Identification and Evolution of Cas9 tracrRNAs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8336202/). If you use parts of my code in research that isn't on identifying tracrRNAs, please cite this Github repository. If you have any questions or run into any issues, please feel free to [log them](https://github.com/skDooley/TRACR_RNA/issues) in the repository.