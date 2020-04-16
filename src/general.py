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

# Gromacs Wrapper
import gromacs
import gromacs.environment
from gromacs.fileformats.ndx import NDX
from gromacs.fileformats.mdp import MDP

# BioPython
from Bio.PDB import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import BiopythonWarning

# General modules
from statistics import mean
from statistics import stdev
import random
from random import shuffle
from random import randint
import warnings
import logging
import os
import glob
import subprocess
import math
import re
import numpy as np

# Local modules
from src import scoring
from src import mutation

# To activate depending on the logging requirements
warnings.simplefilter('always', gromacs.BadParameterWarning)
gromacs.environment.flags['capture_output'] = True
#gromacs.environment.flags['capture_output'] = "file"
#gromacs.environment.flags['capture_output_filename'] = 'gromacs_general_output.txt'
warnings.simplefilter('ignore', BiopythonWarning)
gromacs.start_logging()

# Classes and functions
class complex:
    ########################################################################################
    def __init__(self,binder,pdbID,iteration,score_list,consensus_threshold,num_mutations,score_mode="all",sim_time="10",mode="start",try_mutations=10):
        """
        Initialization of the class complex with a set of arguments:
        
        Arguments:
        binder -- chain identifier of the peptide that will be modified
        pdbID -- name of the PDB file containing the complex after a MD simulation
        iteration -- number of the iteration the run will start. Useful to control posterior restart of the process
        score_list -- list of scoring functions that will be calculated from the trajectory
        consensus_threshold -- threshold of number of scoring functions in agreement to accept the consensus
        num_mutations -- number of mutations defined for the protocol
        score_mode -- mode to calculate the scoring functions (Default: all)
        sim_time -- time that wil be configured to run the steps in nanoseconds (Default: 10)
        mode -- general mode of the running that have two values: start and restart (Default: start)
        try_mutations -- number of mutations tried after having minimization or equilibration problems (Default: 10)
        """
        self.binder=binder
        self.iteration=iteration
        self.pdbID=pdbID
        self.score_mode=score_mode
        self.score_list=score_list
        self.consensus_threshold=consensus_threshold
        self.mode=mode
        self.sim_time=sim_time
        self.num_mutations=num_mutations
        self.try_mutations=try_mutations
    
    ########################################################################################
    def configure_folder(self,folder_name,md_route,md_original_pdb_name):
        """
        Function to configure the folder with the required information to start a running
        
        Arguments:
        folder_name -- name of the folder that will contain all the information
        md_route -- route of an MD simulation where the specific data will be acquired to start/restart the protocol
        md_original_pdb_name -- name of the pdb file that will be used as starting point for the design
        
        Required:
        - PDB file that will be used as starting point
        - GRO file from the PDB file used as starting point
        - Topol files, including the general and the topol files per protein chain
        
        Output:
        - Folder with all the required information for running
        """
        
        # Check the mode the protocol will run
        if self.mode=="start":
            # Create folder and copy the initial standard files
            os.system("rm -r design_output/{folder}; mkdir design_output/{folder}".format(folder=folder_name))
            os.system("cp -r src/start/mdp src/start/tip3p_mod.itp design_output/{}".format(folder_name))
            
            # Modify the npt file time
            mdp_npt=MDP()
            mdp_npt.read("design_output/{}/mdp/md-NPT.mdp".format(folder_name))
            mdp_npt['nsteps']=int(float(self.sim_time)*500000)
            mdp_npt.write()
            
            # Create folders to store the information
            os.system("mkdir design_output/{folder}/binder design_output/{folder}/complexP\
                      design_output/{folder}/log_npt design_output/{folder}/log_nvt design_output/{folder}/score_trajectory\
                      design_output/{folder}/solvent design_output/{folder}/system design_output/{folder}/target\
                      design_output/{folder}/trajectory design_output/{folder}/iterations".format(folder=folder_name))
            
            # Copy the files from the MD
            os.system("cp {}/{}.pdb design_output/{}/{}.pdb".format(md_route,md_original_pdb_name,folder_name,self.pdbID))
            os.system("cp {}/{}.gro design_output/{}/{}.gro".format(md_route,md_original_pdb_name,folder_name,self.pdbID))
            os.system("cp {}/topol_Protein_chain_* {}/topol.top design_output/{}".format(md_route,md_route,folder_name))
            os.system("cp {}/posre_Selection_Protein_chain_* design_output/{}".format(md_route,folder_name))
        
        # Restart from a given iteration
        if self.mode=="restart":
            # Create folder and copy the initial standard files
            os.system("rm -r design_output/{folder}; mkdir design_output/{folder}".format(folder=folder_name))
            os.system("cp -r src/start/mdp src/start/tip3p_mod.itp design_output/{}".format(folder_name))
            
            # Modify the npt file time
            mdp_npt=MDP()
            mdp_npt.read("design_output/{}/mdp/md-NPT.mdp".format(folder_name))
            mdp_npt['nsteps']=int(float(self.sim_time)*500000)
            mdp_npt.write()
            
            # Copy the files from the MD
            os.system("cp {}/system/system_{}.pdb design_output/{}/{}.pdb".format(md_route,str(self.iteration),folder_name,self.pdbID))
            # Generate the gro file
            rc,sout,serr=gromacs.editconf(f="design_output/{}/{}.pdb".format(folder_name,self.pdbID), o="design_output/{}/{}.gro".format(folder_name,self.pdbID), stdout=False)
            
            # Copy the topol files and restrain atom files
            os.system("cp {}/topol_Protein_chain_* {}/topol.top design_output/{}".format(md_route,md_route,folder_name))
            os.system("cp {}/posre_Selection_Protein_chain_* design_output/{}".format(md_route,folder_name))
            os.system("rm design_output/{}/topol_Protein_chain_{}.itp".format(folder_name,self.binder))
            
            # Generate a novel itp file for the chain
            os.system("cp {}/binder/binder_{}.pdb design_output/{}".format(md_route,str(self.iteration),folder_name))
            rc,sout,serr=gromacs.pdb2gmx(f="design_output/{}/binder_{}.pdb".format(folder_name,str(self.iteration)), p="design_output/{}/binder.top".format(folder_name), o="design_output/{}/binder_{}.gro".format(folder_name,str(self.iteration)), stdout=False, input=('6','6'))
            os.system("sed -i '/forcefield/d' design_output/{}/binder.top".format(folder_name))
            os.system("sed -i '/\[ system \]/,$d' design_output/{}/binder.top".format(folder_name))
            os.system("head -n -18 design_output/{}/binder.top > design_output/{}/topol_Protein_chain_{}.itp".format(folder_name,folder_name,self.binder))
            os.system("rm design_output/{}/binder.top design_output/{}/binder_{}.pdb design_output/{}/binder_{}.gro".format(folder_name,folder_name,self.iteration,folder_name,self.iteration))
            os.system("rm posre.itp")
            
            # Create folders to store the information
            os.system("mkdir design_output/{folder}/binder design_output/{folder}/complexP\
                      design_output/{folder}/log_npt design_output/{folder}/log_nvt design_output/{folder}/score_trajectory\
                      design_output/{folder}/solvent design_output/{folder}/system design_output/{folder}/target\
                      design_output/{folder}/trajectory design_output/{folder}/iterations".format(folder=folder_name))
        
        if self.mode=="nothing":
            pass

    ########################################################################################    
    def setup(self,folder_path):
        """
        Configure the protein for starting the simulation. Based on the input files it generates the index file with the respective groups
        
        Arguments:
        folder_path -- folder with the structure information to start
        
        Output:
        - molecules.ndx -- index containing the required groups of the system
        """
                
        # Store the path of the folder with the information
        self.folder_path=folder_path
        self.path="design_output/{}".format(folder_path)
        
        # Read the structure
        parser = PDBParser()
        
        # Get the chains of the initial structure        
        warnings.simplefilter("ignore")
        self.structure = parser.get_structure('PEP', self.path+"/{}.pdb".format(self.pdbID))
        self.models = self.structure[0]
        self.chains=[]
        self.chain_pdbs={}
        self.chain_join=[]
        for chain in self.models:
            chLetter=chain.get_id()
            if chLetter != ' ':
                self.chains.append(chLetter)
                self.chain_pdbs[chLetter]=chain
                if chLetter != self.binder: self.chain_join.append(chLetter)
            else:
                self.chain_pdbs["SOL"]=chain
            
        # The last index reference
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/{}.pdb".format(self.pdbID), o=self.path+"/reference.ndx", stdout=False, input=('q'))
        ref_ndx = NDX()
        ref_ndx.read(self.path+"/reference.ndx")
        self.index_ref=len(ref_ndx)-1
        gromacs.utilities.unlink_gmx(self.path+"/reference.ndx")

        # Create the index from the gro file
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/{}.gro".format(self.pdbID), o=self.path+"/molecules_1.ndx", stdout=False, input=('q'))
        ndx1 = NDX()
        ndx1.read(self.path+"/molecules_1.ndx")
            
        # Create the index used for the rest of the simulation
        input_for_ndx=()
        counter=self.index_ref
        for i,ch in enumerate(self.chains):
            input_for_ndx+=('chain {}'.format(ch),); counter+=1
            input_for_ndx+=('{} & 1'.format(counter),); counter+=1
            if ch != self.binder: input_for_ndx+=('name {} chain{}'.format(counter,ch),) # Changed target by chainA
            else: input_for_ndx+=('name {} binder'.format(counter),) # Create the binder group in the index file
        
        # Generate the sentence that join the chains belonging to the target
        sentence=""
        for i,ch in enumerate(self.chain_join):
            if i==0: sentence=sentence+"chain {}".format(ch)
            else: sentence=sentence+" | chain {}".format(ch)
        
        # Create the target group in the index file
        target_chains="-".join(self.chain_join)
        input_for_ndx+=(sentence,); counter+=1
        input_for_ndx+=('name {} target'.format(counter),)
        input_for_ndx+=('q',)
            
        # Create the index file
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/{}.pdb".format(self.pdbID), o=self.path+"/molecules_2.ndx", stdout=False, input=input_for_ndx)
        
        # This part was made to avoid the reposition of the SOL group in the index file after the first mutation, avoiding errors after that
        ndx2 = NDX()
        ndx2.read(self.path+"/molecules_2.ndx")
        for i,group in enumerate(ndx2):
            if i>self.index_ref:
                ndx1[group]=ndx2[group]
        ndx1.write(self.path+'/molecules.ndx')
        
        # Delete the temp files
        os.system("rm {path}/molecules_1.ndx {path}/molecules_2.ndx".format(path=self.path))
        
    ########################################################################################
    def run_minim(self,pdbID,run_minim=False):
        """
        Function to run the final minimization of the complex
        
        Arguments:
        pdbID -- name of the the pdb file that will be minimized
        run_minim -- boolean flag that control if the minimization will be run (Default: False)
        
        Output
        - Run the minimization and save the output file with the name provided by pdbID
        """
        
        # Prepare the files for running the minimization
        rc,sout,serr=gromacs.grompp(f=self.path+"/mdp/minim.mdp", o=self.path+"/min.tpr", p=self.path+"/{}.top".format(pdbID),
                                    c=self.path+"/{}.gro".format(pdbID), stdout=False)
        gromacs.utilities.unlink_gmx("mdout.mdp")
        
        if run_minim:
            # Run the minimization
            rc,sout,serr=gromacs.mdrun(deffnm=self.path+"/min", stdout=False)
            
            # Rename the files with the final names
            os.system("mv {path}/min.gro {path}/{pdb}.gro".format(path=self.path,pdb=pdbID))
            rc,sout,serr=gromacs.trjconv(f=self.path+"/{}.gro".format(pdbID),s=self.path+"/min.tpr", o=self.path+"/{}.pdb".format(pdbID),stdout=False,input=("0"))
            os.system("sed -i 's/Cl/CL/g' {}/{}.pdb".format(self.path,pdbID)) 
            
            # Delete temporary files
            os.system("rm {path}/min.edr {path}/min.trr {path}/min.log {path}/min.tpr".format(path=self.path))
            
    ########################################################################################
    def run_nvt(self,run_md=False):
        """
        Run NVT simulations based on the mdp file configured for that purpose
        
        Arguments:
        run_md -- boolean flag to indicate if we want to run the MD, or just create the files necessary to run (Default:False)
        
        Output:
        - Files derived from the MD, as well as the respective tpr files necessary to run the simulations
        """
        
        rc,sout,serr=gromacs.grompp(f=self.path+"/mdp/md-NVT.mdp", o=self.path+"/md-nvt.tpr", p=self.path+"/system.top",
                                    c=self.path+"/system.gro", n=self.path+"/molecules.ndx", stdout=False)
        gromacs.utilities.unlink_gmx("mdout.mdp")
        ref_mdp_nvt=MDP()
        ref_mdp_nvt.read(self.path+"/mdp/md-NVT.mdp")
        
        if run_md: rc,sout,serr=gromacs.mdrun(deffnm=self.path+"/md-nvt", stdout=False)
        
    ########################################################################################
    def run_npt(self,initial=False,run_md=False):
        """
        Run NPT simulations based on the mdp file configured for that purpose
        
        Arguments:
        initial -- boolean flag to indicate if this NPT is made to the initial PDB file. Otherwise it will use the files obtained after each mutation (Default:False)
        run_md -- boolean flag to indicate if we want to run the MD, or just create the files necessary to run (Default:False)
        
        Output:
        - Files derived from the MD, as well as the respective tpr files necessary to run the simulations
        """
        if initial:
            # Prepare the NPT
            rc,sout,serr=gromacs.grompp(f=self.path+"/mdp/md-NPT.mdp", o=self.path+"/md-npt.tpr", p=self.path+"/topol.top",
                                        c=self.path+"/{}.gro".format(self.pdbID), n=self.path+"/molecules.ndx", stdout=False)
        else: 
            rc,sout,serr=gromacs.grompp(f=self.path+"/mdp/md-NPT.mdp", o=self.path+"/md-npt.tpr", p=self.path+"/system.top",
                                        c=self.path+"/md-nvt.gro", n=self.path+"/molecules.ndx", stdout=False)
        
        gromacs.utilities.unlink_gmx("mdout.mdp")
        ref_mdp_npt=MDP()
        ref_mdp_npt.read(self.path+"/mdp/md-NPT.mdp")
        self.time_ref=int(ref_mdp_npt['nsteps'])/500
        
        if run_md: rc,sout,serr=gromacs.mdrun(deffnm=self.path+"/md-npt", stdout=False)
    
    ########################################################################################    
    def get_molecules_after_md(self,md_name="md-npt"):
        """
        Function to obtain the molecules after the MD and avoid Periodic Boundary Conditions from the trajectory
        
        Arguments:
        md_name -- file code name where the simulations are stored
        
        Output:
        - Physical files with the molecules obtained from the latest NPT simulation
        """
        
        # Center the xtc file using the target and the protein as reference
        rc,sout,serr=gromacs.trjconv(f=self.path+"/{}.xtc".format(md_name),s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",pbc="nojump",
                                     o=self.path+"/NOJUMP.xtc",stdout=False,input=("0"))
        rc,sout,serr=gromacs.trjconv(f=self.path+"/NOJUMP.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",
                                     o=self.path+"/CENTER.xtc",center=True,stdout=False,input=("Protein \n 0"))
        rc,sout,serr=gromacs.trjconv(f=self.path+"/CENTER.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",pbc="mol",
                                     ur="compact",o=self.path+"/npt-pbc.xtc",stdout=False,input=("0"))
        rc,sout,serr=gromacs.trjconv(f=self.path+"/CENTER.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",pbc="mol",
                                     b="{}".format(str(self.time_ref)),ur="compact",o=self.path+"/npt-pbc.gro",stdout=False,input=("0"))
 
        # After having the final xtc file centered without jumps, we obtain the corresponding files    
        # Get the target
        rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",b="{}".format(str(self.time_ref)),
                                     o=self.path+"/target/target_"+str(self.iteration)+".pdb",stdout=False,input=("target"))
        
        # Get the binder
        rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",b="{}".format(str(self.time_ref)),
                                     o=self.path+"/binder/binder_"+str(self.iteration)+".pdb",stdout=False,input=("binder"))
        
        # Get the complex
        rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",b="{}".format(str(self.time_ref)),
                                     o=self.path+"/complexP/complexP_"+str(self.iteration)+".pdb",stdout=False,input=("Protein"))
        
        # Get the solvent
        rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",b="{}".format(str(self.time_ref)),
                                     o=self.path+"/solvent/solvent_"+str(self.iteration)+".pdb",stdout=False,input=("Water_and_ions"))
        
        # Get the system
        rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",b="{}".format(str(self.time_ref)),
                                     o=self.path+"/system/system_"+str(self.iteration)+".pdb",stdout=False,input=("System"))
 
        # Get the trajectory
        rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",
                                     o=self.path+"/all_npt_traj.pdb",stdout=False,input=("1"))
        
        # Save the system_trajectory file and the log files of the NVT and npt if that is the case
        os.system("mv {} {}/trajectory/system_trajectory_{}.pdb".format(self.path+"/all_npt_traj.pdb",self.path,self.iteration))
        os.system("cp {} {}/log_nvt/md-nvt_{}.log".format(self.path+"/md-nvt.log",self.path,self.iteration))
        os.system("cp {} {}/log_npt/md-npt_{}.log".format(self.path+"/md-npt.log",self.path,self.iteration))
        
        # Delete unnecessary xtc files
        os.system("rm {path}/NOJUMP.xtc {path}/CENTER.xtc".format(path=self.path))
        
    ########################################################################################
    def read_score_trajectory(self,pdb_file,score_list,half_flag,only_complex=False):
        """
        Function to score each snapshot of the MD trajectory for the scoring function selected
        
        Arguments:
        pdb_file -- file containing the trajectory with pdb files
        score_list -- list of scoring functions that will be used to select the mutations
        half_flag -- flag  to use half the trajectory or not during the scoring
        only_complex -- boolean variable that check if the calculation is with the complex only or taking into account the parts as well (Default: False)
        
        Output:
        Score dictionaries and list containing the calculated scores
        """
        
        # Read the pdb information from the trajectory
        PDB_traj=[x.strip() for x in open(pdb_file)]
        len_trajectory=self.time_ref/100
        model_number = 1
        new_file_text = ""
        
        # Dictionaries where the scores will be stored
        score_dictionary={'pisa':0.0,'bach':0.0,'firedock':0.0,'irad':0.0,'totalbmfbluues':0.0,'bmf':0.0,'bluues':0.0,'zrank':0.0}
        total_score={'pisa':[],'bach':[],'firedock':[],'irad':[],'totalbmfbluues':[],'bmf':[],'bluues':[],'zrank':[]}
        
        # Iterate over all the MD frames
        for line in PDB_traj:
            if line == "ENDMDL":
                # Save file with file number in name
                p="model"+str(model_number)
                ref_value=0
                if half_flag==1: ref_value=len_trajectory/2
                
                # Check the flag to select the number of snapshots
                if model_number > ref_value:
                    output_file = open(self.path+"/model" + str(model_number) + ".pdb", "w")
                    output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
                    output_file.close()
                 
                    # Change residue names
                    os.system("sed -i 's/HIP/HIS/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/HID/HIS/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/GLH/GLU/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/ASH/ASP/g' "+self.path+"/model" + str(model_number) + ".pdb")    
                 
                    bash = "grep OC1 {}/model{}.pdb | head -n 1 | tail -n 1 | awk '{{print $4}}'".format(self.path,str(model_number))
                    resA = subprocess.check_output(['bash','-c', bash])
                    bash = "grep OC1 {}/model{}.pdb | head -n 2 | tail -n 1 | awk '{{print $4}}'".format(self.path,str(model_number))
                    resB = subprocess.check_output(['bash','-c', bash])
                    bash = "grep OC1 {}/model{}.pdb | head -n 3 | tail -n 1 | awk '{{print $4}}'".format(self.path,str(model_number))
                    resC = subprocess.check_output(['bash','-c', bash])
            
                    os.system("sed -i 's/CD  ILE/CD1 ILE/g' {}/model{}.pdb".format(self.path,str(model_number)))
                    os.system("sed -i 's/OC1 {} A/O   {} A/g' {}/model{}.pdb".format(resA.strip(),resA.strip(),self.path,str(model_number)))
                    os.system("sed -i 's/OC2 {} A/OXT {} A/g' {}/model{}.pdb".format(resA.strip(),resA.strip(),self.path,str(model_number)))
                    os.system("sed -i 's/OC1 {} B/O   {} B/g' {}/model{}.pdb".format(resB.strip(),resB.strip(),self.path,str(model_number)))
                    os.system("sed -i 's/OC2 {} B/OXT {} B/g' {}/model{}.pdb".format(resB.strip(),resB.strip(),self.path,str(model_number)))
                    os.system("sed -i 's/OC1 {} C/O   {} C/g' {}/model{}.pdb".format(resC.strip(),resC.strip(),self.path,str(model_number)))
                    os.system("sed -i 's/OC2 {} C/OXT {} C/g' {}/model{}.pdb".format(resC.strip(),resC.strip(),self.path,str(model_number)))
                    
                    # Call the scoring functions
                    sc=scoring.score_protein_protein(p,self.path,self.chain_join,self.binder)
                    
                    # Calculate the designated scores. If only complex the scores are assigned 0 values
                    for s in score_list:
                        if s=="pisa":
                            if only_complex:
                                total_score[s].append(0.0)
                            else:
                                sc.computePisa()    
                                total_score[s].append(float(sc.pisa_score))
                        if s=="bach":
                            sc.computeBach()    
                            total_score[s].append(float(sc.bach_score))
                        if s=="firedock":
                            if only_complex:
                                total_score[s].append(0.0)
                            else:
                                sc.computeFiredock()    
                                total_score[s].append(float(sc.firedock_score))
                        if s=="zrank":
                            if only_complex:
                                total_score[s].append(0.0)
                            else:
                                sc.computeZrank()    
                                total_score[s].append(float(sc.zrank_score))
                        if s=="irad":
                            if only_complex:
                                total_score[s].append(0.0)
                            else:
                                sc.computeIrad()    
                                total_score[s].append(float(sc.irad_score))

                        if s=="goap":
                            sc.computeGoap()    
                            total_score["totalgoap"].append(float(sc.totalgoap_score))
                            total_score["dfire"].append(float(sc.dfire_score))
                            total_score["goap"].append(float(sc.goap_score))
                        if s=="bmf-bluues":
                            if only_complex:
                                total_score["totalbmfbluues"].append(0.0)
                                total_score["bmf"].append(0.0)
                                total_score["bluues"].append(0.0)
                            else:
                                sc.computeBMF()    
                                total_score["totalbmfbluues"].append(float(sc.totalbmfbluues_score))
                                total_score["bmf"].append(float(sc.bmf_score))
                                total_score["bluues"].append(float(sc.bluues_score))
                
                # Reset everything for next model
                model_number += 1
                new_file_text = ""
                os.system("rm {}/{}.pdb".format(self.path,p))
            elif not line.startswith("MODEL"):
                new_file_text += line + '\n'
                
        # Store the scores in the dictionary
        for s in score_list:
            if s=="goap":
                score_dictionary["totalgoap"]=mean(total_score["totalgoap"])
                score_dictionary["dfire"]=mean(total_score["dfire"])
                score_dictionary["goap"]=mean(total_score["goap"])
            elif s=="bmf-bluues":
                score_dictionary["totalbmfbluues"]=mean(total_score["totalbmfbluues"])
                score_dictionary["bmf"]=mean(total_score["bmf"])
                score_dictionary["bluues"]=mean(total_score["bluues"])
            else:
                score_dictionary[s]=mean(total_score[s])
        
        # Delete the trajectory file
        os.system("rm {}".format(pdb_file))
        
        # Return the dictionary with the calculated scores
        return score_dictionary,total_score
            
    ########################################################################################
    def score_complex(self,score_list,half_flag,md_name="md-npt"):
        """
        Function that will score the system based on the scoring list and the mode selected
        
        Arguments:
        score_list -- list containing the scoring functions that will be used to run the scoring
        md_name -- name of the trajectory file that will be used to calculate the score averages
        
        Output:
        Files with the scores per trajectory and the score dictionary
        """
        
        # If is selected to score the full complex
        if self.score_mode=="all":
            # Convert the last xtc file to a pdb trajectory
            rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",
                                         o=self.path+"/all_npt_traj.pdb",stdout=False,input=("1"))
            # Target trajectory
            rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",
                                         o=self.path+"/target_npt_traj.pdb",stdout=False,input=("target"))
            # Binder trajectory
            rc,sout,serr=gromacs.trjconv(f=self.path+"/npt-pbc.xtc",s=self.path+"/{}.tpr".format(md_name),n=self.path+"/molecules.ndx",
                                         o=self.path+"/binder_npt_traj.pdb",stdout=False,input=("binder"))
            
            
            # Calculate the scores based on the trajectories
            self.sc_dict_complex,totalAll=self.read_score_trajectory(self.path+"/all_npt_traj.pdb",score_list,half_flag)
            self.sc_dict_target,totalTarget=self.read_score_trajectory(self.path+"/target_npt_traj.pdb",score_list,half_flag,True)
            self.sc_dict_binder,totalBinder=self.read_score_trajectory(self.path+"/binder_npt_traj.pdb",score_list,half_flag,True)
            
            # Calculate the substracted score from complex and individual parts
            self.sc_dict_complete={'pisa':0.0,'bach':0.0,'firedock':0.0,'irad':0.0,'totalbmfbluues':0.0,'bmf':0.0,'bluues':0.0,'zrank':0.0}
            for s in score_list:
                if s=="goap":
                    self.sc_dict_complete["totalgoap"]=self.sc_dict_complex["totalgoap"]-self.sc_dict_target["totalgoap"]-self.sc_dict_binder["totalgoap"]
                    self.sc_dict_complete["dfire"]=self.sc_dict_complex["dfire"]-self.sc_dict_target["dfire"]-self.sc_dict_binder["dfire"]
                    self.sc_dict_complete["goap"]=self.sc_dict_complex["goap"]-self.sc_dict_target["goap"]-self.sc_dict_binder["goap"]
                if s=="bmf-bluues":
                    self.sc_dict_complete["totalbmfbluues"]=self.sc_dict_complex["totalbmfbluues"]-self.sc_dict_target["totalbmfbluues"]-self.sc_dict_binder["totalbmfbluues"]
                    self.sc_dict_complete["bmf"]=self.sc_dict_complex["bmf"]-self.sc_dict_target["bmf"]-self.sc_dict_binder["bmf"]
                    self.sc_dict_complete["bluues"]=self.sc_dict_complex["bluues"]-self.sc_dict_target["bluues"]-self.sc_dict_binder["bluues"]
                else:
                    self.sc_dict_complete[s]=self.sc_dict_complex[s]-self.sc_dict_target[s]-self.sc_dict_binder[s]
            
            # Create the score trajectory files
            for score in totalAll:
                report=open(self.path+"/score_trajectory/score_trajectory_{}_{}.txt".format(score,self.iteration),"w")
                print("Processing {} scores ...".format(score))
                for i,valAll in enumerate(totalAll[score]):
                    valComplete=float(valAll)-float(totalTarget[score][i])-float(totalBinder[score][i])
                    report.write("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(float(valAll),float(totalTarget[score][i]),float(totalBinder[score][i]),float(valComplete)))
                report.close()
    
    ########################################################################################
    def mutation_random(self,residues_mod,mutation_document,score_dictionary,half_flag,scwrl_path):
        """
        Function to mutate any amino acid of a peptide sequence randomly
        
        Arguments:
        residues_mod -- list of positions of the peptide sequence that can be mutated
        mutation_document -- text file with the report of the scores and the mutations
        score_dictionary -- general dictionary where the scores are stored per iteration
        
        Output:
        - Loop that will run the mutations until the number determined is achieved
        """
        
        # Get the sequence of the binder chain
        aminoacids={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                    "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
        self.sequence_binder=""
        self.initial_res_pos_binder=0
        parser = PDBParser()
        reference = parser.get_structure('REF',"{}/binder/binder_{}.pdb".format(self.path,self.iteration))
        for chain in reference[0]:
            for i,residue in enumerate(chain):
                seq=aminoacids[residue.get_resname()]
                self.sequence_binder=self.sequence_binder+str(seq)
                if i==0: self.initial_res_pos_binder=residue.get_full_id()[3][1]
        
        # Check which amino acids should mutate
        self.number_aa_used=len(residues_mod)
        self.aa_list=["A","D","E","F","H","I","K","L","M","G","N","P","Q","R","S","T","V","W","Y","C"]
        self.total_aa_used=len(self.aa_list)
        
        # Define the initial settings of the reference peptide
        mutation_list=[]
        last_good_iteration=self.iteration
        pep_single_mutation_1=self.sequence_binder
        flag_restart_mutation=0
        
        # Save the good iterations in a list
        good_iterations=[last_good_iteration]
        peptide_reference=[pep_single_mutation_1]
        
        # Start iteration over the number of mutations
        for i in range(0,self.num_mutations):
            
            # Loop to check if the mutation was minimized without issues
            counter_try=0
            while flag_restart_mutation==0:
                # Sum the number of trys
                counter_try+=1

                # Check to change the last good iteration and finish the process if the protocol is stuck
                if counter_try==self.try_mutations:
                    size_iter=len(good_iterations)
                    if size_iter==1:
                        mutation_document.write("Terminated because of iteration at first attempt\n")
                        quit()
                    if size_iter>1:
                        last_good_iteration=good_iterations[size_iter-2]
                        pep_single_mutation_1=peptide_reference[size_iter-2]
                        mutation_document.write("Second try with iteration {} based on last {}\n".format(last_good_iteration,good_iterations[size_iter-1]))

                if counter_try==self.try_mutations+10: 
                     mutation_document.write("Terminated because of iteration at second attempt\n")
                     quit()
                

                # Obtain the random positions on the peptide chain and in the set of AA that will be used to do the mutation
                position_pep_prev=randint(0,self.number_aa_used-1)
                position_peptide=int(residues_mod[position_pep_prev])-1
                position_mutation=randint(0,self.total_aa_used-1)
                
                # Check if the old AA is equal to the new one, in such cases it will modify it to keep them different
                old_aa=pep_single_mutation_1[position_peptide]
                new_aa=self.aa_list[position_mutation]
                if new_aa==old_aa:
                    if new_aa=="C":
                        new_aa=self.aa_list[position_mutation-1]
                    else:
                        new_aa=self.aa_list[position_mutation+1]
                position=position_peptide+1
                
                # Create the string with the novel sequence
                temporal_pep = list(pep_single_mutation_1)
                temporal_pep[position-1] = new_aa
                pep_single_mutation_2=''.join(temporal_pep)
                 
                # Read the new binder
                parser = PDBParser()
                reference = parser.get_structure('REF',"{}/complexP/complexP_{}.pdb".format(self.path,last_good_iteration))
                
                # Call the mutation module and run the required steps before running the simulation
                local_mutation=mutation.mutate_peptide(self.path,pep_single_mutation_1,pep_single_mutation_2,self.initial_res_pos_binder+position-1,
                                                       self.binder,reference,old_aa,new_aa,self.iteration,self.chain_join,last_good_iteration,scwrl_path)
                
                # Mutate and run the local minimizations
                local_mutation.replace_amino_acid()
                local_mutation.mutate_scwrl()
                try:
                    local_mutation.run_minim_complex(True)
                    local_mutation.add_box_mutation()
                    # Run minimization
                    print("Running full minimization ...")
                    self.run_minim("system",True)
                    steps=glob.glob("step*")
                    mutation_document.write("Attempt mutation {}{}{}{} Steps {}\n".format(old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa,len(steps)))
                    flag_restart_mutation=1
                    os.system("rm step*")

                except:
                    print("Minimization problem")
                    mutation_document.write("Minimization problem with {}{}{}{}\n".format(old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa))
                    flag_restart_mutation=0
                    os.system("rm step*")
                    
                # Try for the NVT
                if flag_restart_mutation==1:
                    try:
                        # Increase the iteration
                        self.iteration+=1
            
                        # Configure the new system
                        mutated_system=complex(self.binder,"system",self.iteration,self.score_list,self.consensus_threshold,self.num_mutations,self.score_mode,self.sim_time,self.mode,self.try_mutations)
                        mutated_system.setup(self.folder_path)
                        print("Running NVT equilibration ...")
                        mutated_system.run_nvt(True)
                        flag_restart_mutation=1
                        os.system("rm step*")
                    except:
                        print("NVT problem")
                        self.iteration=self.iteration-1
                        mutation_document.write("NVT problem with {}{}{}{}\n".format(old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa))
                        flag_restart_mutation=0
                        os.system("rm step*")
            
            print("Running NPT simulation with new complex ...")
            mutated_system.run_npt(initial=False,run_md=True)
            mutated_system.get_molecules_after_md()

            # Calculate the score
            print("Scoring the new complex ...")
            mutated_system.score_complex(mutated_system.score_list,half_flag)
            print(mutated_system.sc_dict_complex)
            print(mutated_system.sc_dict_target)
            print(mutated_system.sc_dict_binder)
            print(mutated_system.sc_dict_complete)
            
            # After running the protocol, write the mutation in a document
            score_sentence=""
            score_dictionary[self.iteration]={}
            for key in mutated_system.sc_dict_complete:
                score_sentence=score_sentence+key+":"+str(mutated_system.sc_dict_complete[key])+" "
                score_dictionary[self.iteration][key]=float(mutated_system.sc_dict_complete[key])
            
            # Execute the consensus
            accept=mutated_system.consensus_criteria(score_dictionary,last_good_iteration,mutated_system.score_list)
            if accept==1:
                mutation_document.write("Iteration_{}: {}{}{}{} - Accepted Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa,score_sentence,pep_single_mutation_2))
            else:
                mutation_document.write("Iteration_{}: {}{}{}{} - Rejected Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa,score_sentence,pep_single_mutation_2))
            
            # Save additional file with the scores per iteration
            iteration_file=open("{}/iterations/score_iterations_{}.txt".format(self.path,self.iteration),"w")
            if accept==1:
                iteration_file.write("Iteration_{}: {}{}{}{} - Accepted Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa,score_sentence,pep_single_mutation_2))
            else:
                iteration_file.write("Iteration_{}: {}{}{}{} - Rejected Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,self.initial_res_pos_binder+position-1,new_aa,score_sentence,pep_single_mutation_2))
            iteration_file.close()

            # Print the mutation and update the peptide sequence if the consensus is accepted            
            print(pep_single_mutation_1,pep_single_mutation_2)
            if accept==1:
                pep_single_mutation_1=pep_single_mutation_2
                last_good_iteration=self.iteration
            print(accept,last_good_iteration)

            # Store the iterations that were good
            good_iterations.append(last_good_iteration)
            peptide_reference.append(pep_single_mutation_1)
            # Configure again the restart
            flag_restart_mutation=0

    ########################################################################################   
    def consensus_criteria(self,score_dictionary,last_iteration,score_list):
        """
        Function to apply the rank-by-vote strategy using consensus threshold
        
        Arguments:
        score_dictionary -- general dictionary containing the scores calculated by previous iterations
        last_iteration -- iteraction with the latest scores accepted that will be used to runthe differences
        
        Output:
        acceptance or not of the mutation
        """
        
        # First test using complex score from mutation report
        sc1_old=0.0;sc1_new=0.0
        sc2_old=0.0;sc2_new=0.0
        if len(score_list)>2:
            sc3_old=0.0;sc3_new=0.0
        if len(score_list)>3:
            sc4_old=0.0;sc4_new=0.0
        if len(score_list)>4:
            sc5_old=0.0;sc5_new=0.0
        if len(score_list)>5:
            sc6_old=0.0;sc6_new=0.0
        
        # Check in the dictionary the calculated scores
        for key in score_dictionary:
            if key==last_iteration:
                for scores in score_dictionary[key]:
                    if scores==score_list[0]: sc1_old=score_dictionary[key][scores]
                    if scores==score_list[1]: sc2_old=score_dictionary[key][scores]
                    if len(score_list)>2:
                        if scores==score_list[2]: sc3_old=score_dictionary[key][scores]
                    if len(score_list)>3:
                        if scores==score_list[3]: sc4_old=score_dictionary[key][scores]
                    if len(score_list)>4:
                        if scores==score_list[4]: sc5_old=score_dictionary[key][scores]
                    if len(score_list)>5:
                        if scores==score_list[5]: sc6_old=score_dictionary[key][scores]
            if key==self.iteration:
                for scores in score_dictionary[key]:
                    if scores==score_list[0]: sc1_new=score_dictionary[key][scores]
                    if scores==score_list[1]: sc2_new=score_dictionary[key][scores]
                    if len(score_list)>2:
                        if scores==score_list[2]: sc3_new=score_dictionary[key][scores]
                    if len(score_list)>3:
                        if scores==score_list[3]: sc4_new=score_dictionary[key][scores]
                    if len(score_list)>4:
                        if scores==score_list[4]: sc5_new=score_dictionary[key][scores]
                    if len(score_list)>5:
                        if scores==score_list[5]: sc6_new=score_dictionary[key][scores]
        
        # Apply the the rank-by-vote consensus approach  
        counter_threshold=0
        
        if sc1_new-sc1_old < 0.0: counter_threshold+=1
        if sc2_new-sc2_old < 0.0: counter_threshold+=1
        if len(score_list)>2:
            if sc3_new-sc3_old < 0.0: counter_threshold+=1
        if len(score_list)>3:
            if sc4_new-sc4_old < 0.0: counter_threshold+=1
        if len(score_list)>4:
            if sc5_new-sc5_old < 0.0: counter_threshold+=1
        if len(score_list)>5:
            if sc6_new-sc6_old < 0.0: counter_threshold+=1
        
        # Acceptance flag
        acceptance=0
        if counter_threshold>=self.consensus_threshold: acceptance=1
        
        # Return the acceptance or not of the mutation
        return acceptance
