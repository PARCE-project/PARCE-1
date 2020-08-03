#!/bin/bash


#PARCE: Protocol for Amino acid Refinement through Computational Evolution

#From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
#Computer Physics Communications
#Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio
#Year: 2020


########################################################################################
# Ubuntu 16.04 tested commands
########################################################################################

sudo apt-get update
sudo apt-get install pdb2pqr
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-yaml
sudo apt-get install python3-tk
python3 -m pip install GromacsWrapper==0.7 numpy==1.18 scipy==1.4 matplotlib==3.0
