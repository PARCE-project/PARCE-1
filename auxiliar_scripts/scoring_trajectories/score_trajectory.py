#!/usr/bin/python3

"""
PARCE: Protocol for Amino acid Refinement through Computational Evolution
Script to score an MD trajectory in PDB format

From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
Computer Physics Communications 
Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio
Year: 2020

NOTE: TO RUN THE SCRIPT A SET OF VARAIBLES SHOULD BE MODIFIED AT THE END BASED ON YOUR LOCAL PATHS

Explanation:

This script calculates the scores for a given MD trajectory or single file in PDB format, printing the score per frame
and the final average for each score. To run the script first update it with the path of the PARCE source folder and
the path of the trajectory, together with the name of the file. NOTE: The file should be in PDB format

The variables to modify are:
- src_route: Path of the PARCE src folder
- traj_path: Path where the PDB trajectory or single file is located
- pdb_name: Name of the PDB file containing the single complex or the MD trajectory frames

To run the command use the instruction 'python3 score_trajectory.py'. The outputs are single files per scoring function
with the calculated scores for all the frames included in the PDB file. With that output, the user can easily implement
a consensus by vote based on the trajectories they require to compare. 
"""

# Import local modules
import os
import sys
from statistics import mean
from statistics import stdev

########################################################################################
# Functions
########################################################################################
def read_score_trajectory(pdb_file,score_list,src_route,path,chain_join,binder):
    
    # Read the pdb information from the trajectory
    PDB_traj=[x.strip() for x in open(pdb_file)]
    model_number = 1
    new_file_text = ""
    
    score_dictionary={}
    for sc in score_list: score_dictionary[sc]=0.0
    
    total_score={}
    for sc in score_list: total_score[sc]=[]
    
    # Iterate over all the MD frames
    for i,line in enumerate(PDB_traj):
        if line == "ENDMDL" or i==len(PDB_traj)-1:
            
            # Save file with file number in name
            output_file = open("model"+str(model_number)+".pdb", "w")
            output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
            output_file.close()
            
            # Change residue names
            os.system("sed -i 's/HIP/HIS/g' model"+str(model_number)+".pdb")
            os.system("sed -i 's/HID/HIS/g' model"+str(model_number)+".pdb")
            os.system("sed -i 's/GLH/GLU/g' model"+str(model_number)+".pdb")
            os.system("sed -i 's/ASH/ASP/g' model"+str(model_number)+".pdb")
            os.system("sed -i 's/CD  ILE/CD1 ILE/g' model"+str(model_number)+".pdb")
            
            # Model ID
            p="model"+str(model_number)
            
            # Call the function to score
            sc=scoring.score_protein_protein(p,path,src_route,chain_join,binder)
            
            # Calculate the designated scores
            for s in score_list:
                
                if s=="pisa":
                    sc.computePisa()    
                    total_score[s].append(float(sc.pisa_score))
                if s=="bach":
                    sc.computeBach()    
                    total_score[s].append(float(sc.bach_score))
                if s=="firedock":
                    sc.computeFiredock()    
                    total_score[s].append(float(sc.firedock_score))
                if s=="zrank":
                    sc.computeZrank()    
                    total_score[s].append(float(sc.zrank_score))
                if s=="irad":
                    sc.computeIrad()    
                    total_score[s].append(float(sc.irad_score))
                if s=="bmf-bluues":
                    sc.computeBMF()    
                    total_score["bmf-bluues"].append(float(sc.totalbmfbluues_score))
                    #total_score["bmf"].append(float(sc.bmf_score))
                    #total_score["bluues"].append(float(sc.bluues_score))

            # Reset everything for next model
            model_number += 1
            new_file_text = ""
            os.system("rm {}.pdb".format(p))
        elif not line.startswith("MODEL"):
            new_file_text += line + '\n'
            
    # Store the score of the target
    for s in score_list:
        if s=="bmf-bluues":
            score_dictionary["bmf-bluues"]=mean(total_score["bmf-bluues"])
            #score_dictionary["bmf"]=mean(total_score["bmf"])
            #score_dictionary["bluues"]=mean(total_score["bluues"])
        else:
            score_dictionary[s]=mean(total_score[s])
        
    # Return the dictionary with the calculated scores
    return score_dictionary,total_score

########################################################################################
def score_complex(score_list,pdb_name,src_route,local_path,traj_path,chain_join,binder):
        
    # Calculate the scores based on the trajectories
    sc_dict_complex,totalAll=read_score_trajectory("{}/{}.pdb".format(traj_path,pdb_name),score_list,src_route,local_path,chain_join,binder)
    
    # Write the scores per frame in a file for each scoring function
    for score in totalAll:
        report=open("score_{}_{}.txt".format(score,pdb_name),"w")
        print("Processing {} score ...".format(score))
        for i,valAll in enumerate(totalAll[score]):
            valComplete=float(valAll)
            report.write("{:.3f}\n".format(float(valAll)))
        report.close()
    
    # Create a summary file with all the scoring function averages
    score_document=open("score_summary_{}.txt".format(pdb_name),"w")
    score_sentence=""
    for key in sc_dict_complex:
        score_sentence=score_sentence+key+":"+str(sc_dict_complex[key])+" "
    score_document.write("Score: {}\n".format(score_sentence))

########################################################################################        
# Main execution
########################################################################################
if __name__ == '__main__':
    
    # MAIN PARAMETERS TO MODIFY
    src_route="/home/PARCE-1" # Local path of your PARCE installation folder
    sys.path.append(src_route)
    from src import scoring
    
    # Local path defined automatically
    local_path=os.getcwd()
    # List of scoring functions
    score_list=["bach","pisa","zrank","irad","firedock","bmf-bluues"] 
    # Path of the trajectory or single PDB file
    traj_path="/home/PARCE-1/auxiliar_scripts/scoring_trajectories"
    # Name of the PDB trajectory file or single file
    pdb_name="example_trajectory"
    # List with the ids of the target chains
    target_chains=["A"]
    # Chain id of the peptide
    binder_chain="B" 
    
    # Main function to score the trajectory or PDB file
    score_complex(score_list,pdb_name,src_route,local_path,traj_path,target_chains,binder_chain)
