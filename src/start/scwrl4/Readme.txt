SCWRL4.0
(c) 2009 Georgii G. Krivov, Maxim V. Shapovalov, Roland L. Dunbrack
Fox Chase Cancer Center, 333 Cottman Avenue, Philadelphia PA 19111
Roland.Dunbrack@fccc.edu
http://dunbrack.fccc.edu/scwrl4


SCWRL4 is a program for predicting side-chain conformations for a given protein backbone.
It runs as a console or command-line application and is parameterized with a configuration file.
Data is transferred within files which are specified via command-line arguments.

You can always get a QUICK REFERENCE by launching the program without any arguments.

Two command-line arguments are mandatory:

-i <backbone_IN>

		This specifies the backbone upon which side-chains should be packed.

		This must be a well-formed PDB flat text file (extension ".ent" or ".pdb")
		Please see http://www.wwpdb.org/documentation/format23/v2.3.html
		for elaborate description of the format.
		Only ATOM records from this file will be considered for structure construction.

		Atoms are grouped into residues based on their IDs (chainID+resNumber+insertionCode).

		Backbone sites will be resolved solely from main conformation of the atoms labeled
		as N,C,CA and O. If any of these atoms is absent then the corresponding residue will
		be completely ignored. Supplementary conformations (e.g., alt_loc of "B") of the
		backbone are ignored.

		The peptide bond connectivity between backbone sites will be determined based on
		backbone C-N atom-atom distances. More than one peptide bond for any N or C atom
		will be treated as invalid input resulting in execution termination.

		OXT atom at any C-terminus will be modeled only if the corresponding residue in the
		input contains an OXT atom (atom is rebuilt; input OXT coordinates are ignored).

		Three hydrogen atoms (instead of one) at any N atom will be modeled in either of two
		cases: this nitrogen does not participate in any peptide bond and its resNumber = 1,
		or the corresponding residue in the input contains "H2" atom (its coordinates are ignored).

-o <model_OUT>

		Specifies the local path and filename where the predicted model should be written.

		The format is the same as for input (flat text, PDB).

		Output will contain ATOM records which form one or several chains. All residues will 
		have only one conformation and their chain IDs and residue numbers will match the input.

		If the input PDB contained residues with incomplete backbone conformation of the backbone,
		then these residues will not be present in the output.

		The order in which residues are given is guided by the actual peptide bond
		connectivity and thus may be different from the input.

		If disulfide bonds were detected then corresponding SSBOND records will be placed
		at the beginning of the output file. If the input file contained CRYST1 and SCALE* records
		they will also be written in the beginning of the output.

		TER record will be placed in the end of each chain.


Also there are several optional arguments:

-s <sequence>

		Specifies the location of the sequence file.

		This enables side-chain prediction against an arbitrary subset of residues
		and also provides a convenient way to switch amino-acid type of side-chains.

		Content of the file should represent a sequence of single-letter amino acid codes
		specifying the amino-acid type for the corresponding residue in the input file.
		Letters in the sequence are associated with the residues according to the order
		in which the former are observed in the input file.

		An UPPER-CASE letter specifies that the corresponding residue's side-chain conformation
		should be predicted for the amino acid type given in the sequence file.

		Lower case letter means that the input conformation of a side-chain should be preserved.
		In this case a side-chain in the input must be complete and its actual amino-acid type
		must match the one specified by the letter in the sequence. If this holds then
		the dihedral chi-angles of the input side-chain will be computed and used for rotamer
		construction and the variances for flexible rotamer model will be extracted from the closest
		rotamer in the rotamer library. Side-chains may lack hydrogen atoms. Hydroxyl groups
		will be modeled in all possible conformations and the optima will be placed in output.
		If the residue is multiconformational then the described procedure is applied to each
		of the conformation whose amino-acid type is consistent with the amino acid in the
		sequence.

		Letters that do not stand for any amino-acid type are replaced with 'G' and the
		corresponding warning is printed.

		All non-letter characters (numbers, whitespaces, etc) are ignored.

-p <params>

		Forces the program to use values of the parameters from the specified file.

		A widespread ".ini" file format (http://en.wikipedia.org/wiki/INI_file) is used
		for configuration files. Values are named and arranged in sections.
		Single-line comments should start with '#' symbol.

		[RotLib] section contains the path to the binary rotamer library. This must be
		consistent with the actual location of the rotamer library.

		[Params_ALL] section may contain a numeric value for "frmSigmaBoost".
		This parameter adds extra multiplier to the weights that are used in front of
		variances in the flexible rotamer model. Setting it to zero will disable
		flexible rotamer model.

-#

		Switches on the prediction in crystal mode.

		CRYST1 record from the input file is used to resolve the symmetry.
		If SCALEn records are present then the scale matrix in the input file is
		checked for consistency with the one derived from the CRYST1 record.
		All MTRIXn records are ignored which means that only the actual content
		of input file is considered.

-% <sym_ops>

		Forces the program to account for symmetry that is defined by the spatial movement
		operators from the specified file. This enables accounting for arbitrary symmetry.

		Each operator defines a movement which is performed on the original input structure
		to obtain virtual symmetric copies. Interaction of the original structure with these
		copies is accounted for during the energy calculations.

		A special format is used. Each spatial movement operator is defined by one line.
		The first character is a name of the operator. It is followed by six
		floating point numbers separated by whitespaces. The first number specifies
		the angle of rotation. The second and the third specify a direction around
		which this rotation occurs. This direction is given in the spherical coordinates
		system by angle theta and angle phi. All angles are in degrees.
		The last three numbers represent the translation in the original coordinates.

		Covalent bonds between original structure and symmetry copies are not considered.

-f <frame.pdb>

		Load atoms from the specified file and treat them as static frame atoms.

		Both ATOM and HETATM records are considered as source. Chemical elements are
		resolved from its short name which should be present in element field (columns 77-78).
		All charges are ignored and atoms are treated as neutral.
		Hydrogen bonds between the model and atoms within static frame are not accounted for.

		If symmetric closure is being involved, static frame is not subject to virtual clone
		generation induced by symmetry operators. If a symmetric static frame should be modeled
		a user should prepare corresponding frame file manually. This enables a general case: when
		inner symmetry of the static frame differs from one used for model.

-v
		Ignore chi variances in the rotamer library.
		
		This option forces to use Rigid Rotamer Model. Flexible subrotamers will not be involved
		in the energy calculations. Typically this yields faster but less accurate prediction.

-h
		This will omit hydrogen atoms in the output.

-t
		This disables capping of the terminal residues.



For the sake of better usability the default priority of SCWRL4 process is 'BelowNormal'.
The following two options will override this value.

-2
		Run the programm with 'Idle' priority
-0
		Run the program with 'Normal' priority


The last two options are of special interest and are not likely to be commonly used

-g <graph>

		This will save the interaction graph into the specified text file.
		In the begining of this file all residues are listed and self energies of all active
		rotamers are provided. In the rest part of this file pairwise energies are listed.

-w <workspace>

		Saves precomputed workspace into specified binary file.
		This may be helpful for reporting bugs during the optimization.


