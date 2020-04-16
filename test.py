#!/usr/bin/python

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
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Miguel A. Soler", "Alessandro Laio", "Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################
from src import scoring
from src import mutation
from src import general
# BioPython
from Bio.PDB import *

########################################################################################
# Test
########################################################################################

if __name__ == '__main__':
    
    # List of variables that should be defined to create the class
    chain="B"
    pdbID="1ppg_AAPAAAPP"
    score_list=["bach","pisa","zrank","irad","bmf-bluues","firedock"]
    threshold=3
    half_flag=0
    sim_time="5"
    num_mutations=5
    temperature="310"
    
    # Additional variables for the main script
    residues_mod=[1,2,3,4,5,6,7,8]
    
    # Variables sensible to the mode: start or restart
    folder="test_result"
    mode="start"
    peptide_reference="AAPAAAPP"
    md_route="./design_output/1ppg_init"
    md_original="1ppg_AAPAAAPP"
    
    iteration=0
    score_dictionary={}
    score_mode="all"
    
    # Variables for mutation test
    scwrl_path="/usr/local/bin/scwrl4/Scwrl4"
    peptide_mutated="AAPATAPP"
    
    # File containing the results
    report=open("report_test.txt","w")
    flag_step1=0
    flag_step2=0
    flag_step3=0
    
    # 1. Test Gromacs functionality
    try:
        print("1. Starting the test of Gromacs functionalities ...")
        protein_complex=general.complex(chain,pdbID,iteration,score_list,threshold,num_mutations,score_mode,sim_time,mode,temperature)
        protein_complex.configure_folder(folder,md_route,md_original)
        protein_complex.setup(folder)
        print("####################################")
        print("The Gromacs test passed successfully")
        report.write("The Gromacs test passed successfully\n")
        print("####################################")
        flag_step1=1
    except:
        print("####################################")
        print("The Gromacs test failed. Please verify that Gromacs is installed, as well as the GromacsWrapper package in python")
        report.write("The Gromacs test failed. Please verify that Gromacs is installed, as well as the GromacsWrapper package in python\n")
        print("####################################")
    
    # 2. Test Scwrl4 functionality
    try:
        print("2. Starting the test of Scwrl4 functionalities ...")
        parser = PDBParser()
        reference = parser.get_structure('REF',"design_output/{}/{}.pdb".format(folder,pdbID))
        local_mutation=mutation.mutate_peptide("design_output/{}".format(folder),peptide_reference,peptide_reference,5,chain,reference,"A","T",0,["A"],0,scwrl_path)
        local_mutation.replace_amino_acid()
        local_mutation.mutate_scwrl()
        print("####################################")
        print("The Scwrl4 test passed successfully")
        report.write("The Scwrl4 test passed successfully\n")
        print("####################################")
        flag_step2=1
    except:
        print("####################################")
        print("The Scwrl4 test failed. Please verify that Scwrl4 is correctly installed, and the permissions to execute are activated")
        report.write("The Scwrl4 test failed. Please verify that Scwrl4 is correctly installed, and the permissions to execute are activated\n")
        print("####################################")
    
    # 3. Test the scoring functions    
    try:
        print("3. Starting the test of scoring functions ...")
        # Call the scoring functions
        sc=scoring.score_protein_protein("complex","design_output/{}".format(folder),["A"],chain)
        # Calculate the designated scores
        for s in score_list:
            if s=="pisa":
                sc.computePisa()    
                print(float(sc.pisa_score))
            if s=="bach":
                sc.computeBach()    
                print(float(sc.bach_score))
            if s=="firedock":
                sc.computeFiredock()    
                print(float(sc.firedock_score))
            if s=="zrank":
                sc.computeZrank()    
                print(float(sc.zrank_score))
            if s=="irad":
                sc.computeIrad()    
                print(float(sc.irad_score))
            if s=="bmf-bluues":
                sc.computeBMF()    
                print(float(sc.totalbmfbluues_score))
        print("####################################")
        print("The scoring functions test passed successfully")
        report.write("The scoring functions test passed successfully\n")
        report.write("Pisa score: {}\n".format(sc.pisa_score))
        report.write("BACH score: {}\n".format(sc.bach_score))
        report.write("Firedock score: {}\n".format(sc.firedock_score))
        report.write("ZRANK score: {}\n".format(sc.zrank_score))
        report.write("IRAD score: {}\n".format(sc.irad_score))
        report.write("BMF-BLUUES score: {}\n".format(sc.totalbmfbluues_score))
        print("####################################")
        flag_step3=1
    except:
        print("####################################")
        print("The scoring functions test failed. Please verify that the required scoring files are present in the src/scores folder")
        report.write("The scoring functions test failed. Please verify that the required scoring files are present in the src/scores folder\n")
        print("####################################")
    
    if flag_step1==1 and flag_step2==1 and flag_step3==1:
        print("####################################")
        print("Everything is ready. You can start the protocol :)")
        print("####################################")
        report.write("####################################\n")
        report.write("Everything is ready. You can start the protocol :)\n")
        report.write("####################################\n")
        
    # Close the report
    report.close()
