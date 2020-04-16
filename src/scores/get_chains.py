#!/usr/bin/python

from Bio.PDB import *
import sys

protein=sys.argv[1]
path=sys.argv[2]

report=open("%s/chains.seq" %path,"w")

parser = PDBParser()
structure = parser.get_structure('PEP',protein)
for chain in structure[0]:
    ppb = CaPPBuilder()
    chLetter=chain.get_id()
    name=protein.split(".")[0]
    for pp in ppb.build_peptides(chain):
        seq=pp.get_sequence()
    io = PDBIO()
    io.set_structure(chain)
    io.save('{}_{}.pdb'.format(name,chLetter))

    report.write(chain.get_id()+"\n")
    
report.close()
