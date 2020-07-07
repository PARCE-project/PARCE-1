#!/usr/bin/python3

"""
PARCE: Protocol for Amino acid Refinement through Computational Evolution

From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
Computer Physics Communications 
Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio
Year: 2020

Third-party tools required:

- Scwrl4: http://dunbrack.fccc.edu/scwrl4/license/index.html
- Gromacs 5.1.4 (tested version): http://manual.gromacs.org/documentation/5.1.4/download.html
- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python3-rdkit

Data required:
- A PDB file containing the starting system (including protein, peptide/small protein and solvent)
- A gro file of the PDB template structure
- The topology files of the structure chains
- (Optional) Files containing atoms restrained during the simulations
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Miguel A. Soler","Alessandro Laio", "Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################
from src import scoring
from src import mutation
from src import general
import argparse
import yaml
import os

########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':
    
    # Script arguments
    parser = argparse.ArgumentParser(description='PeptideDesigner: an open source protocol to design peptide bound to protein targets')
    parser.add_argument('-c', dest='config_file', type=argparse.FileType(mode='r'), required=True,
                        help='File containing all the necessary parameters to run the protocol') 
    
    #####################################################################################
    # Assignment of parameters
    #####################################################################################
    args = parser.parse_args()
    if args.config_file:
        data = yaml.load(args.config_file)
        delattr(args, 'config_file')
        arg_dict = args.__dict__
        for key, value in data.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].extend(v)
            else:
                arg_dict[key] = value
    else:
        print("A config file is necessary to run the protocol. Exiting ...")
        exit()
    
    # Check the arguments
    if args.folder:
        folder=args.folder
    else:
        print("The parameter 'folder' is required for the analysis. Exiting ...")
        exit()
    if args.src_route:
        src_route=args.src_route
    else:
        print("The parameter 'src_route' is required for the analysis. Exiting ...")
        exit()
    if args.mode in ("start","restart","nothing"):
        mode=args.mode
    else:
        print("The parameter 'mode' is required for the analysis, or an option should be selected from: start, restart and nothing. Exiting ...")
        exit()
    if args.peptide_reference:
        peptide_reference=args.peptide_reference
    else:
        print("The parameter 'peptide_reference' is required for the analysis. Exiting ...")
        exit()
    if args.pdbID:
        pdbID=args.pdbID
    else:
        print("The parameter 'pdbID' is required for the analysis. Exiting ...")
        exit()
    if args.chain:
        chain=args.chain
    else:
        print("The parameter 'chain' is required for the analysis. Exiting ...")
        exit()
    if args.sim_time:
        sim_time=args.sim_time
    else:
        print("The parameter 'sim_time' in nanoseconds is required for the analysis. Exiting ...")
        exit()
    if args.num_mutations:
        num_mutations=args.num_mutations
    else:
        print("The parameter 'num_mutations' is required for the analysis. Exiting ...")
        exit()
    if args.try_mutations:
        try_mutations=args.try_mutations
    else:
        print("The parameter 'try_mutations' is required for the analysis. Exiting ...")
        exit()
    if args.half_flag:
        half_flag=args.half_flag
    else:
        half_flag=0
    if args.residues_mod:
        residues_mod=args.residues_mod.split(",")
        residues_mod = list(map(int, residues_mod))
    else:
        print("The parameter 'residues_mod' is required for the analysis. Exiting ...")
        exit()
    if args.md_route:
        md_route=args.md_route
    else:
        print("The parameter 'md_route' is required for the analysis. Exiting ...")
        exit()
    if args.md_original:
        md_original=args.md_original
    else:
        print("The parameter 'md_original' is required for the analysis. Exiting ...")
        exit()
    if args.score_list:
        score_list=args.score_list.split(",")
    else:
        print("The parameter 'score_list' is required for the analysis. Exiting ...")
        exit()
    if args.threshold:
        threshold=args.threshold
    else:
        print("The parameter 'threshold' is required for the analysis. Exiting ...")
        exit()
    if args.mutation_method in ("faspr","scwrl4"):
        mutation_method=args.mutation_method
    else:
        print("The parameter 'mutation_method' is required for the analysis, or an option should be selected from: faspr and scwrl4. Exiting ...")
        exit()
    try:    
        if args.scwrl_path:
            scwrl_path=args.scwrl_path
    except:
        print("The Scwrl4 path has not been provided. In that case the software will be called from the system path")
        scwrl_path="Scwrl4"
    try:    
        if args.gmxrc_path:
            os.system(". {}".format(args.gmxrc_path))
            os.environ['GMX_MAXBACKUP'] = "-1"
    except:
        print("The Gromacs path has not been provided. In that case the software will be called using the system path")
        os.environ['GMX_MAXBACKUP'] = "-1"
    
    ####################################################################################
    # Starting the design
    ####################################################################################
    # Start some variables
    iteration=0
    score_dictionary={}
    score_mode="all"
    
    # Create complex object with the basic information to start
    protein_complex=general.complex(chain,pdbID,iteration,score_list,threshold,num_mutations,score_mode,sim_time,mode,try_mutations)
    protein_complex.configure_folder(folder,src_route,md_route,md_original)
    protein_complex.setup(folder,src_route)
    # Run the first NPT simulation based on the provided data
    print("Starting first npt simulation ...")
    protein_complex.run_npt(initial=True,run_md=True)
    print("Getting molecules from simulation ...")
    protein_complex.get_molecules_after_md()
    
    # Score the firs run
    print("Scoring the system ...")
    protein_complex.score_complex(score_list,half_flag)
    print(protein_complex.sc_dict_complete)
    
    # Write the scores in the report document
    mutation_document=open("design_output/"+folder+"/mutation_report.txt","w")
    mutation_document.write("Structure: {}\n".format(protein_complex.pdbID))
    score_sentence=""
    score_dictionary[0]={}
    for key in protein_complex.sc_dict_complete:
        score_sentence=score_sentence+key+":"+str(protein_complex.sc_dict_complete[key])+" "
        score_dictionary[0][key]=float(protein_complex.sc_dict_complete[key])
    mutation_document.write("Iteration_{}: Original - Accepted Score: {} Sequence:{}\n".format(iteration,score_sentence,peptide_reference))
    
    # Start the mutation of random amino acids
    #protein_complex.mutation_random(residues_mod,mutation_document,score_dictionary,half_flag,mutation_method,scwrl_path)
