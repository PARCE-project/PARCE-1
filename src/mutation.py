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
import os
import subprocess
import numpy as np
import warnings

# BioPython
from Bio.PDB import *
from Bio import BiopythonWarning

# Gromacs Wrapper
import gromacs
from gromacs.fileformats.ndx import NDX

# To activate depending on the logging requirements
warnings.simplefilter('always', gromacs.BadParameterWarning)
gromacs.environment.flags['capture_output'] = True
#gromacs.environment.flags['capture_output'] = "file"
#gromacs.environment.flags['capture_output_filename'] = 'gromacs_mutation_output.txt'
warnings.simplefilter('ignore', BiopythonWarning)
gromacs.start_logging()

# Class and functions
class mutate_peptide:
    ########################################################################################
    def __init__(self,path,pep_original,pep_reference,pep_position,pep_chain,pep_pdb,old_aa,new_aa,iteration,chain_join,last_iteration,scwrl_path):
        """
        Initializer of the mutation protocol
        
        Arguments:
        path -- path where the files required for the mutation are stored
        pep_original -- peptide sequence before the single mutation
        pep_reference -- peptide sequence with the generated mutation
        pep_position -- position of the peptide that was mutated
        pep_chain -- chain of the peptide that is modified
        pep_pdb -- pdb file containing the complex that will be used to run the mutation
        old_aa -- the amino acid that wil lbe changed
        new_aa -- the new amino acid mutated
        iteration -- the current iteration
        chain_join -- the chains of the target
        last_iteration -- the last iteration wth an accepted mutation
        scwrl_path -- path to the Scwrl4 program
        """
        self.path=path
        self.pep_original=pep_original
        self.pep_reference=pep_reference
        self.pep_position=pep_position
        self.pep_chain=pep_chain
        self.pep_pdb=pep_pdb # Protein structure
        self.old_aa=old_aa
        self.new_aa=new_aa
        self.iteration=iteration
        self.chain_join=chain_join
        self.last_iteration=last_iteration
        self.scwrl_path=scwrl_path
    
    ########################################################################################   
    def replace_amino_acid(self):
        """
        Function that remove the amino acid that will be mutated
        
        Output
        pre-mutated.pdb -- file without the side chain of the old amino acid
        """
        aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                    "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}

        # Read the PDB file
        residues=self.pep_pdb.get_residues()
        chain=self.pep_pdb[0][self.pep_chain]
        
        # Report the mutation made
        message="The residue {} in chain {} and position {} will be changed by {}".format(self.old_aa,self.pep_chain,self.pep_position,self.new_aa)
        
        # Rename the residue
        chain[self.pep_position].resname=aminoacids[self.new_aa]
        
        # Delete the other atoms leaving only the atoms of the backbone
        ids=[]
        for a in chain[self.pep_position]:
            atomId=a.id
            if atomId not in ("N","CA","O","C","OC1"): ids.append(atomId)
        for i in ids: chain[self.pep_position].detach_child(i)
            
        # Saving the new structure
        io = PDBIO()
        io.set_structure(self.pep_pdb)
        io.save(self.path+"/pre-mutated.pdb")

########################################################################################
    def mutate_scwrl(self):
        """
        Function to generate the new amino acid side chain using scwrl
        
        Output:
        complex.pdb -- file containing the complex with the respective mutation
        """
        aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                    "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
        
        # Get the last amino acid of the peptide
        sequence=self.pep_reference
        s = list(sequence)
        seq_size=len(s)
        lastAA=aminoacids[s[-1]]
        
        # Modify atom name of the last AA to avoid errors
        if seq_size>=10:
            if seq_size>=100:
                os.system("sed -i 's/OC1 {} {} {}/O   {} {} {}/g' {}/pre-mutated.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/OC2 {} {} {}/OXT {} {} {}/g' {}/pre-mutated.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/CD  ILE/CD1 ILE/g' {}/pre-mutated.pdb".format(self.path))
            else:
                os.system("sed -i 's/OC1 {} {}  {}/O   {} {}  {}/g' {}/pre-mutated.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/OC2 {} {}  {}/OXT {} {}  {}/g' {}/pre-mutated.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/CD  ILE/CD1 ILE/g' {}/pre-mutated.pdb".format(self.path))
        else:
            os.system("sed -i 's/OC1 {} {}   {}/O   {} {}   {}/g' {}/pre-mutated.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
            os.system("sed -i 's/OC2 {} {}   {}/OXT {} {}   {}/g' {}/pre-mutated.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
            os.system("sed -i 's/CD  ILE/CD1 ILE/g' {}/pre-mutated.pdb".format(self.path))

        # Separate the chains
        os.system("python3 src/scores/get_chains.py {}/pre-mutated.pdb {}".format(self.path,self.path))
        
        # Obtain the region of the peptide before the position that will be mutated
        if self.pep_position>=10:
            bash = "grep -n \"N   {} {}  {}\" {}/pre-mutated_{}.pdb | grep -Eo '^[^:]+'".format(aminoacids[self.new_aa],self.pep_chain,self.pep_position,self.path,self.pep_chain)
        else:
            bash = "grep -n \"N   {} {}   {}\" {}/pre-mutated_{}.pdb | grep -Eo '^[^:]+'".format(aminoacids[self.new_aa],self.pep_chain,self.pep_position,self.path,self.pep_chain)

        nini = subprocess.check_output(['bash','-c', bash])
        os.system("head -n{} {}/pre-mutated_{}.pdb > {}/start".format(int(nini)-1,self.path,self.pep_chain,self.path))
        
        # Obtain the region of the peptide after the position that will be mutated
        if self.pep_position>=10:
            bash = "grep -n \"O   {} {}  {}\" {}/pre-mutated_{}.pdb | grep -Eo '^[^:]+'".format(aminoacids[self.new_aa],self.pep_chain,self.pep_position,self.path,self.pep_chain)
        else:
            bash = "grep -n \"O   {} {}   {}\" {}/pre-mutated_{}.pdb | grep -Eo '^[^:]+'".format(aminoacids[self.new_aa],self.pep_chain,self.pep_position,self.path,self.pep_chain)

        nfin = subprocess.check_output(['bash','-c', bash])
        bash = "wc -l {}/pre-mutated_{}.pdb | awk '{{print $1}}'".format(self.path,self.pep_chain)
        tot = subprocess.check_output(['bash','-c', bash])
        os.system("tail -n{} {}/pre-mutated_{}.pdb > {}/end".format(int(tot)-int(nfin),self.path,self.pep_chain,self.path))
        
        # Obtain the sequence of the peptide with the new position in upper case
        sequence=self.pep_reference.lower() 
        s = list(sequence)
        s[self.pep_position-1]=self.new_aa
        sequence="".join(s)
        print(sequence)
        os.system("echo {} > {}/mutant.seq".format(sequence,self.path))
        
        # Join the target chains in one file
        for i,ch in enumerate(self.chain_join):
            if i==0:            
                os.system("grep ATOM {}/pre-mutated_{}.pdb > {}/final_target.pdb".format(self.path,ch,self.path))
                os.system("echo 'TER' >> {}/final_target.pdb".format(self.path))
            else:
                os.system("grep ATOM {}/pre-mutated_{}.pdb >> {}/final_target.pdb".format(self.path,ch,self.path))
                os.system("echo 'TER' >> {}/final_target.pdb".format(self.path))
        
        # Run the Scwrl4 program
        os.system("{} -i {}/pre-mutated_{}.pdb -f {}/final_target.pdb -o {}/pre-mutated_{}_mod.pdb -s {}/mutant.seq -h".format(self.scwrl_path,self.path,self.pep_chain,self.path,self.path,self.pep_chain,self.path))
        
        # Join the parts again, and the peptide with the target, creating the complex.pdb file
        if self.pep_position>=10:
            os.system("grep \"{} {}  {}\" {}/pre-mutated_{}_mod.pdb > {}/mut".format(aminoacids[self.new_aa],self.pep_chain,self.pep_position,self.path,self.pep_chain,self.path))
        else:
            os.system("grep \"{} {}   {}\" {}/pre-mutated_{}_mod.pdb > {}/mut".format(aminoacids[self.new_aa],self.pep_chain,self.pep_position,self.path,self.pep_chain,self.path))

        os.system("cat {}/start {}/mut {}/end | grep -v END > {}/mutated_binder.pdb".format(self.path,self.path,self.path,self.path))
        os.system("cat {}/final_target.pdb {}/mutated_binder.pdb > {}/complex.pdb".format(self.path,self.path,self.path))
        
        if seq_size>=10:
            if seq_size>=100:
                os.system("sed -i 's/O   {} {} {}/OC1 {} {} {}/g' {}/complex.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/OXT {} {} {}/OC2 {} {} {}/g' {}/complex.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/CD1 ILE/CD  ILE/g' {}/complex.pdb".format(self.path))
            else:
                os.system("sed -i 's/O   {} {}  {}/OC1 {} {}  {}/g' {}/complex.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/OXT {} {}  {}/OC2 {} {}  {}/g' {}/complex.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
                os.system("sed -i 's/CD1 ILE/CD  ILE/g' {}/complex.pdb".format(self.path))
        else:
            os.system("sed -i 's/O   {} {}   {}/OC1 {} {}   {}/g' {}/complex.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
            os.system("sed -i 's/OXT {} {}   {}/OC2 {} {}   {}/g' {}/complex.pdb".format(lastAA,self.pep_chain,str(seq_size),lastAA,self.pep_chain,str(seq_size),self.path))
            os.system("sed -i 's/CD1 ILE/CD  ILE/g' {}/complex.pdb".format(self.path))

       
        # Delete temporary files
        os.system("rm {path}/final_target.pdb {path}/mutated_binder.pdb {path}/pre-mutated.pdb {path}/pre-mutated_* {path}/start {path}/mut {path}/end {path}/chains.seq {path}/mutant.seq".format(path=self.path))
        
    ########################################################################################
    def run_minim_complex(self,run_minim=False):
        """
        Function to run a local minimization on the side chain that was mutated
        
        Arguments:
        run_minim -- boolean flag that will control if the minimization is run or not
        
        Output:
        complex.pdb -- new complex pdb with the minimization and the new itp files
        """
        
        # Get the chain with the peptide to generate a novel itp file
        os.system("python3 src/scores/get_chains.py {}/complex.pdb {}".format(self.path,self.path))
        rc,sout,serr=gromacs.pdb2gmx(f=self.path+"/complex_"+self.pep_chain+".pdb", p=self.path+"/binder.top", o=self.path+"/complex_"+self.pep_chain+".gro", stdout=False, input=('6','6'))
        os.system("sed -i '/forcefield/d' {}/binder.top".format(self.path))
        os.system("sed -i '/\[ system \]/,$d' {}/binder.top".format(self.path))
        os.system("mv {}/binder.top {}/complex_Protein_chain_{}.itp".format(self.path,self.path,self.pep_chain))
        rc,sout,serr=gromacs.editconf(f=self.path+"/complex_"+self.pep_chain+".gro", o=self.path+"/complex_"+self.pep_chain+".pdb", stdout=False)
        
        # Fix the amino acid nomenclature
        os.system("for i in ASP ARG HIS HIE HID HIP LYS GLU SER THR ASN GLN CYS CYX GLY PRO ALA VAL ILE LEU MET PHE TYR TRP; do sed -i s/\"$i  \"/\"$i {}\"/g {}/complex_{}.pdb; done".format(self.pep_chain,self.path,self.pep_chain))
        for i,ch in enumerate(self.chain_join):
            if i==0:            
                os.system("grep ATOM {}/complex_{}.pdb > {}/complex.pdb".format(self.path,ch,self.path))
                os.system("echo 'TER' >> {}/complex.pdb".format(self.path))
            else:
                os.system("grep ATOM {}/complex_{}.pdb >> {}/complex.pdb".format(self.path,ch,self.path))
                os.system("echo 'TER' >> {}/complex.pdb".format(self.path))
        
        # Get the new complex.pdb and the peptide chain itp file and delete temporal files
        os.system("grep ATOM {}/complex_{}.pdb >> {}/complex.pdb".format(self.path,self.pep_chain,self.path))
        os.system("echo 'TER' >> {}/complex.pdb".format(self.path))
        os.system("rm {}/complex_*.pdb".format(self.path))
        os.system("rm {}/complex_*.gro".format(self.path))
        os.system("rm {}/chains.seq".format(self.path))
        os.system("head -n -18 {}/complex_Protein_chain_{}.itp > {}/temp; mv {}/temp {}/complex_Protein_chain_{}.itp".format(self.path,self.pep_chain,self.path,self.path,self.path,self.pep_chain))
        
        # Copy the topol files of the target chains, which are the same always
        for ch in self.chain_join:
            os.system("cp {}/topol_Protein_chain_{}.itp {}/complex_Protein_chain_{}.itp".format(self.path,ch,self.path,ch))
        
        # Copy the topol.top to complex.top and delete all the additional atoms
        os.system("cp {}/topol.top {}/complex.top".format(self.path,self.path))
        os.system("sed -i '/Ion/d' {}/complex.top".format(self.path))
        os.system("sed -i '/SOL/d' {}/complex.top".format(self.path))
        os.system("sed -i '/NA/d' {}/complex.top".format(self.path))
        os.system("sed -i '/CL/d' {}/complex.top".format(self.path))
        os.system("sed -i '/solvent/d' {}/complex.top".format(self.path))
        os.system("sed -i 's/topol_/complex_/g' {}/complex.top".format(self.path))
        
        # Get a pdb of the complex where an index will be created
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/complex.pdb", o=self.path+"/reference.ndx", stdout=False, input=('q'))
        ref_ndx = NDX()
        ref_ndx.read(self.path+"/reference.ndx")
        index_ref=len(ref_ndx)-1
        gromacs.utilities.unlink_gmx(self.path+"/reference.ndx")
        
        # Create the side chain index
        input_for_ndx=()
        counter=index_ref
        input_for_ndx+=('chain {}'.format(self.pep_chain),); counter+=1
        input_for_ndx+=('name {} binder'.format(counter),);
        input_for_ndx+=('"SideChain" & "binder"'+' &  r {}'.format(self.pep_position),); counter+=1
        input_for_ndx+=('"System" &! {}'.format(counter),); counter+=1
        input_for_ndx+=('name {} scmut'.format(counter),);
        sentence=""
        for i,ch in enumerate(self.chain_join):
            if i==0: sentence=sentence+"chain {}".format(ch)
            else: sentence=sentence+" | chain {}".format(ch)
        input_for_ndx+=(sentence,); counter+=1
        input_for_ndx+=('name {} target'.format(counter),)
        input_for_ndx+=('q',)
        
        # Generate the index file
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/complex.pdb", o=self.path+"/scmut.ndx", stdout=False, input=input_for_ndx)
        
        # Generate the gro file
        rc,sout,serr=gromacs.editconf(f=self.path+"/complex.pdb", o=self.path+"/complex.gro", stdout=False)
        
        # Add a small box for the residues
        os.system("sed -i '$ d' {}/complex.gro".format(self.path))
        os.system('echo "   20.0   20.0   20.0" >> {path}/complex.gro'.format(path=self.path))
        
        # Prepare the files for the minimization
        rc,sout,serr=gromacs.grompp(f=self.path+"/mdp/minim_scmut.mdp", o=self.path+"/complex.tpr", p=self.path+"/complex.top", n=self.path+"/scmut.ndx", c=self.path+"/complex.gro", stdout=False)
        gromacs.utilities.unlink_gmx("mdout.mdp")
        
        # Run the minimization of the side chain alone and the residues around it
        if run_minim:
            # Run the minimization
            print("Running first minimization ...")
            rc,sout,serr=gromacs.mdrun(deffnm=self.path+"/complex", stdout=False)
            
            # Get the complex pdb file
            rc,sout,serr=gromacs.trjconv(f=self.path+"/complex.gro",s=self.path+"/complex.tpr", o=self.path+"/min_complex.pdb",stdout=False,input=("1"))
            os.system("rm posre.itp {path}/complex.tpr {path}/complex.top; grep -v ENDMDL {path}/min_complex.pdb | grep -v MODEL > {path}/complex.pdb; rm {path}/min_complex.pdb {path}/complex.log {path}/complex.trr {path}/complex.edr {path}/scmut.ndx".format(path=self.path))
            
    ########################################################################################  
    def add_box_mutation(self):
        """
        Function to add the solvent and run a minimization of the local amino acid with the waters included
        
        Output:
        system.gro -- file containing the complete system minimized after the mutation
        """
        
        # Read the solvent
        os.system("grep ATOM {path}/solvent/solvent_{itera}.pdb | grep -v ENDMDL > {path}/solvent.pdb".format(path=self.path,itera=self.last_iteration))
        
        # Concatenate complex and system
        os.system("cat {path}/complex.pdb {path}/solvent.pdb > {path}/system.pdb".format(path=self.path))
        rc,sout,serr=gromacs.editconf(f=self.path+"/system.pdb", o=self.path+"/system_mod.pdb", stdout=False)
        os.system("mv {}/system_mod.pdb {}/system.pdb".format(self.path,self.path))

        # Make an index of the system
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/system.pdb", o=self.path+"/index.ndx", stdout=False, input=('chain {} \n q'.format(self.pep_chain)))
        
        # Copy the topol files of the target which are always the same
        for ch in self.chain_join:
            os.system("cp {}/topol_Protein_chain_{}.itp {}/system_Protein_chain_{}.itp".format(self.path,ch,self.path,ch))
        
        # Copy the topol.top to system.top
        os.system("cp {}/topol.top {}/system.top".format(self.path,self.path))
        os.system("cp {}/complex_Protein_chain_{}.itp {}/system_Protein_chain_{}.itp".format(self.path,self.pep_chain,self.path,self.pep_chain))
        os.system("sed -i 's/topol_/system_/g' {}/system.top".format(self.path))
       
        # Select water and ions within 0.2 distance of the residue
        rc,sout,serr=gromacs.make_ndx(f=self.path+"/system.pdb", o=self.path+"/index.ndx", stdout=False, input=('chain {}'.format(self.pep_chain),'q'))
        rc,sout,serr=gromacs.select(f=self.path+"/system.pdb", n=self.path+"/index.ndx", s=self.path+"/system.pdb", on=self.path+"/index_sol.ndx", stdout=False, select="group Water_and_ions and same residue as within 0.2 of (group ch{} and resnr {})".format(self.pep_chain,self.pep_position))
        
        # Solve the issue with atoms overlapped with the selected residue
        values=[x.strip() for x in open(self.path+"/index_sol.ndx")]
        atomsSOL=[]
        for i,v in enumerate(values):
            if i!=0:
                info=v.split()
                atomsSOL=atomsSOL+info
        # List of the atoms that will be deleted
        atomsDelete=[]
         
        # Check the overlapped atoms
        for a in atomsSOL:
            # Obtain the list of atoms from the index file
            bash = "awk '$2 == '{}' {{print $6','$7','$8}}' {}/system.pdb".format(a,self.path)
            coordinates = subprocess.check_output(['bash','-c', bash])
            comp=coordinates.strip().split()
            comparison=[]
            for c in comp: comparison.append(float(c))
            ndComp=np.array(comparison)
            
            distancesSOL=[]
            
            # Read the structure in biopython
            parser = PDBParser()
            structure = parser.get_structure('PEP', self.path+"/system.pdb")
            model = structure[0]
            
            # Check the distances with all the atoms from the selected residue
            for residue in model[self.pep_chain]:
                resC=residue.get_resname()
                resNumber=residue.get_full_id()[3][1]
                if resNumber==self.pep_position:
                    for atom in residue:
                        idAtom = atom.get_id()
                        if idAtom[0].isdigit() == False:
                            if resC=="ILE" and idAtom=="CD": idAtom="CD1"
                            
                            diff = atom.coord - ndComp
                            diffValue=np.sqrt(np.sum(diff * diff))
                            distancesSOL.append(float(diffValue))
                            
            # Threshold to determine which atoms can be overlapped
            if min(distancesSOL)<1.0:
                if a not in atomsDelete: atomsDelete.append(a)
         
        # Selection of the final atoms that will be included in the index
        final_index=[]
        for element in atomsSOL:
            flag=0
            for delete in atomsDelete:
                if abs(int(element)-int(delete))<=2: flag=1
            if flag==0:
                final_index.append(element)
                
        # Update of the index sol file
        new_index=open(self.path+"/index_sol2.ndx","w")
        new_index.write("{}\n".format(values[0]))
        group=[]
        counter=1
        for ele in final_index:
            if counter <15:
                group.append(ele)
                counter+=1
            else:
                group.append(ele)
                new_index.write(" ".join(group)+" \n")
                counter=1
                group=[]
        new_index.write(" ".join(group)+" ")
        new_index.close()
        
        # Update the file
        os.system("mv {}/index_sol2.ndx {}/index_sol.ndx".format(self.path,self.path))
                
        ref_ndx = NDX()
        ref_ndx.read(self.path+"/index.ndx")
        index_ref=len(ref_ndx)-1
                                          
        # Create the side chain index in a template file
        os.system("echo 'name 0 overlap' > %s/template" %self.path)
        os.system("echo '\"SideChain\" & \"ch{}\" & r {}' >> {}/template".format(self.pep_chain,str(self.pep_position),self.path))
        os.system("echo '\"overlap\" | \"SideChain_&_ch{}_&_r_{}\"' >> {}/template".format(self.pep_chain,str(self.pep_position),self.path))
        os.system("echo '\"System\" &! \"overlap_SideChain_&_ch{}_&_r_{}\"' >> {}/template".format(self.pep_chain,str(self.pep_position),self.path))
        os.system("echo 'q' >> {}/template".format(self.path))
        
        # Create an index joining both created before
        os.system("gmx -quiet make_ndx -f {path}/system.pdb -n {path}/index_sol.ndx {path}/index.ndx -o {path}/total_index.ndx < {path}/template".format(path=self.path))
        os.system("sed -i 's/System_&_\!overlap_SideChain_&_ch{}_&_r_{}/to_block/g' {}/total_index.ndx".format(self.pep_chain,str(self.pep_position),self.path))
        
        # Generate the gro file
        rc,sout,serr=gromacs.editconf(f=self.path+"/system.pdb", o=self.path+"/system.gro", stdout=False)
        
        # Prepare the files for the minimization and run
        rc,sout,serr=gromacs.grompp(f=self.path+"/mdp/minim_overlap.mdp", o=self.path+"/systemNEW.tpr", p=self.path+"/system.top", n=self.path+"/total_index.ndx", c=self.path+"/system.gro", stdout=False)
        gromacs.utilities.unlink_gmx("mdout.mdp")
        print("Running second minimization ...")
        rc,sout,serr=gromacs.mdrun(deffnm=self.path+"/systemNEW", stdout=False)

        # Copy the system.gro file that will be used to run the last minimization
        os.system("cp {}/systemNEW.gro {}/system.gro".format(self.path,self.path))
        os.system("sed -i '$ d' {}/system.gro".format(self.path))
        os.system("tail -n1 {path}/npt-pbc.gro >> {path}/system.gro".format(path=self.path)) 
        
        # Delete temporal files
        os.system("rm {path}/complex.pdb {path}/solvent.pdb {path}/systemNEW* {path}/template {path}/index.ndx {path}/index_sol.ndx {path}/total_index.ndx *.itp".format(path=self.path))
        
