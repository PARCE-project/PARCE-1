# PARCE

## Protocol for Amino acid Refinement through Computational Evolution

* From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
* Computer Physics Communications, 2020
* Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio


## Purpose

Here we present PARCE, an open source Protocol for Amino acid Refinement through Computational Evolution that implements an advanced and promising method for the design of peptides and proteins. The protocol performs a random mutation in the binder sequence, then samples the bound conformations using molecular dynamics simulations, and evaluates the protein-protein interactions from multiple scoring. Finally, it accepts or rejects the mutation by applying a consensus criterion based on binding scores. The procedure is iterated with the aim to explore efficiently novel sequences with potential better affinities toward their targets. We also provide a tutorial for running and reproducing the methodology.

## Third-party tools required:

- Scwrl4: http://dunbrack.fccc.edu/scwrl4/license/index.html
- Gromacs 5.1.4 (tested version): http://manual.gromacs.org/documentation/5.1.4/download.html

**NOTE: Path to both executables can be provided in the configuration file**

Scwrl4 can be installed freely after filling a form available in the website to obtain an academic license. **Please verify the permissions to run the program**. Gromacs 5.1.4 **(version tested in the protocol)** can be compiled and installed using the source code. The scoring functions are provided in the **src** folder and configured to run the analysis.

The BioPython and additional python modules can be installed directly from the OS repositories. An example in Ubuntu 16.04 is:

```
sudo apt-get install pdb2pqr
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-tk
sudo apt-get install python3-yaml
pip3 install GromacsWrapper
```

**NOTE: A `install_dependencies.sh` file is provided to automatize the installation of dependencies in the Linux (Ubuntu) operating system**

## Input files required

To run the protocol, **it is required to previously run a simulation in Gromacs with the system of interest**. After that, the following input files are required:

- A PDB file containing the starting system, including the protein, peptide and solvent. **Ideally renumber the chains to start from position 1 at each chain**
- A gro file of the PDB template structure
- The topology files of the structure chains
- (Optional) Files with itp extensions that selects the atoms restrained during the simulations

## MDP files for Gromacs MD simulations

The protocol has included a set of mdp files (with fixed names) to run the multiple minimization steps, the NVT equilibrations and the NPT production stages. The parameters have been optimized for the protocol efficiency and accuracy. **However, if these parameters wants to be modified by the user, the source files can be found in the folder src/start/mdp. The default temperature of the system is 310K.**

## Graphical summary

![summary](src/pictures/summary_PARCE.png)

## How to run the protocol script

**The protocol has been created and tested using Python3.5**. The basic command line to run the script is:

`python3 run_protocol.py [-h] -c CONFIG_FILE`
                                       
where the arguments are:

```
arguments:
  -h, --help      show this help message and exit
  -c CONFIG_FILE  File containing all the necessary parameters to run the
                  protocol
 ```
 
The configuration file describe all the parameters required to run the protocol. A full detailed explanation is provided in the next section.

## Configuration file

A configuration file is required with the following parameters:

- **folder**: Name of the folder that has all the input and output files of the protocol.
- **mode**: The design mode, which has three possible options: *start* (start the protocol from zero), *restart* (start from a particular iteration of a previous run) and *nothing* (just run without modifying existing files).
- **peptide_reference**: The sequence of the peptide, or protein fragment that will be modified.
- **pdbID**: Name of the structure that is used as input, which contains the protein, the peptide/sprotein and the solvent molecules.
- **chain**: Chain id of the peptide/protein in the structural complex
- **sim_time**: Time in nanoseconds that will be used to sample the complex after each mutation. *Recommended 5 or 10 ns*
- **num_mutations**: Number of mutations that will be attempted
- **try_mutations**: Number of mutations tried after having minimization or equilibration problems. *These issues can be found depending on the system and the previous sampling of the complex before starting the protocol*
- **residues_mod**: These are the specific positions of the residues that want to be modified. This depends on the peptide/protein length and the numbering in the PDB file. *It is recommended to renumber the input structure to associate the first position to residue number 1, looking to avoid errors in posterior stages related with the renumbering of the coordinate files*
- **md_route**: Path to the folder containing the input files, which are the files used during the previous MD sampling of the system.
- **md_original**: Name of the system file located in the folder containing the previous MD sampling.
- **score_list**: List of the scoring functions that will be used to calculate the consensus. Currently the package only has available the code for six of them: BACH, Pisa, ZRANK, IRAD, BMF-BLUUES and FireDock. *At least two should be selected.*
- **half_flag**: Flag that controls which part of the trajectory is used to obtain the average score. If *0*, the full trajectory is used, if *1*, only the last half.
- **threshold**: Threshold used for the consensus. If the number of scoring functions in agreement are equal or greater than the threshold, then the mutation is accepted.

*Optional arguments:*
-  **scwrl_path**: Provide the path to Scwrl4 in case it is not installed in a PATH folder. By default the system will use the system path to call the program
-  **gmxrc_path**: Provide the path to GMXRC in case Gromacs was not included previously in the system path

The following is an example of the configuration file *(config_peptide.txt)* for a protease-peptide complex provided in the code:

```
folder: 1ppg_design
mode: start
peptide_reference: AAPAAAPP
pdbID: 1ppg_AAPAAAPP
chain: B
sim_time: 5
num_mutations: 30
try_mutations: 10
residues_mod: 1,2,3,4,5,6,7,8
md_route: ./design_output/peptide_protein
md_original: 1ppg_AAPAAAPP
score_list: bach,pisa,zrank,irad,bmf-bluues,firedock
half_flag: 0
threshold: 3
scwrl_path: /usr/local/bin/scwrl4/Scwrl4
gmxrc_path: /usr/local/gromacs/bin/GMXRC
```
If any of these parameters are missing, the protocol stops and prints a warning messsage to the user.

In addition, another configuration file called *(config_protein.txt)* can be used to run a protein-protein example based on a nanobody protein interaction. The configuration file contains all the required information to run the analysis based on the starting structural data provided in the folder `design_output/protein_protein`

## Folder content
When a design run start, an initial folder is created with the required input files, and the folders that will store the outputs step-by-step. The following is a list of the folders created, and the specific content stored:

- **binder**: Store the peptide/protein structure after each mutation attempt
- **target**: Store the target structure after each mutation attempt
- **complexP**: Store the target-peptide/protein structure after each mutation attempt
- **solvent**: Store the solvent box after each mutation attempt
- **system**: Store the complete target-peptide/protein-solvent complex after each mutation attempt
- **trajectory**: Store the MD trajectory of the previous mutations
- **score_trajectory**: Store the average scores for each snapshot from the trajectories. *The file is split into four columns. The first column is the score of the complex. The second and third are the scores for the receptor and peptide alone. The fourth column is the total score after doing the difference between the complex and each component*
- **log_npt**: Store the log file from each npt run to verify possible errors
- **log_nvt**: Store the log file from each nvt run to verify possible errors

The design protocol results are summarized in the output file called `mutation_report.txt`, which contains details per mutation step like the type of mutation, the average scores, the peptide/small protein sequence and if the mutation was accepted or not. In addition, the report includes failed attempts based on minimization or equilibration problems. The latest can happen depending on the MD result. To overcome these issues, the protocol automatically attempt a number of mutations using the immediately accepted structure. If after that number the system keeps failing, the new mutations will use the accepted structure before the last structure available. If the problem persist during a number of mutations, the system will stop.
In addition, a file named `gromacs.log` stores the logging of all the Gromacs commands, and the file `gromacs_general_output.txt` store the latest Gromacs command used, just to track the progress during the run. To guarantee that the additional tools and dependencies are functioning, a set of tests are provided.

## Tests
A number of tests are provided to check the PARCE functionalities of the third-party tools. These are:

- Call to Gromacs functions to configure an example input file
- Attempt a mutation using the Scwrl4 program
- Calculate all the scoring functions using the initial system provided

The test can be run using the following command: `python3 test.py`. A report with the results per test is generated in the main folder with the name `report_test.txt`.

**NOTE: Please change in the test.py script the paths to Scwrl4 and Gromacs based on your local installation. The variables are: scwrl_path and gmxrc_path**

## Post-analysis of the results

With the final `mutation_report.txt`, it is possible to check and select the accepted sequences, and plot the scores to verify that the energies are getting minimized after the mutation steps. The results are numbered per iteration step, and the folders content facilitates locating the information. Examples of the analysis are provided in the original manuscript.

## Docker details

To run the docker image, first you require to install Docker in the operating system of interest. A guide of the distributions for Docker Desktop (Windows and Mac) and Docker server (Linux) can be found here: https://docs.docker.com/engine/install/

After verifying the correct installation, and checking if `sudo` is required or not, the image can be downloaded as follows:

```docker pull rochoa85/parce-1:latest```

To create the container and start playing with the protocol use the following command:

```docker run -it rochoa85/parce-1 /bin/bash```

After that, you can find the code in the folder: `/home/PARCE` and start after following the instructions provided in the README

To exit, just enter the command `exit`. Then you can check the container created with the command `sudo docker ps -a`. The **container-id** is in the first column, which will be used to access later the docker container. To achieve that, first activate the **container-id** with:

```docker start container-id```

and then open the bash environment with the following command

```docker exec -it container-id /bin/bash```

**Note: To have instructions about how to use docker in different OS, follow these tutorials:**

- Windows: https://docs.docker.com/docker-for-windows/
- Mac: https://docs.docker.com/docker-for-mac/
- Linux distributions: https://docs.docker.com/engine/install/ubuntu/

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
