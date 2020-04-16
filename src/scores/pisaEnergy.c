	/* 
	Last updated Sept 27, 2011
	*/ 
	
	#include<stdio.h>
	#include<string.h>
	#include<math.h>
	#include<stdlib.h>
	#include<time.h>
	
	#define PDB_FILENAME_LEN 100
	#define MAX_NUM_ATOMS 60000
	#define MAX_LINE_LEN 200
	#define NUM_ATOM_TYPES 32
	#define FIELD_LEN 10
	#define MAX_DISTANCE 8.0
	#define MIN_DISTANCE 2.0
	#define NUM_BINS 3 
	#define NUM_BINS_SPLINE 5
	
	const double r1Spline[NUM_BINS_SPLINE]={2.0,3.0,4.0,4.5,6.0};
	const double r2Spline[NUM_BINS_SPLINE]={3.0,4.0,4.5,6.0,8.0};
	
	const double r1[NUM_BINS]={2.0,3.5,5.0};
	const double r2[NUM_BINS]={3.5,5.0,8.0};
	
	typedef struct atom_info {
	//long atom_id;
	char atom_name[10];
	char residue_name[4];
	char chain;
	int residue_number;
	double x;
	double y;
	double z;
	short atom_type;
	} atom;
	
	short getAtomType32(char *atom_name,char *res_name)
	{  
	
 	if(strcmp(res_name,"LYS")==0 && strcmp(atom_name,"NZ")==0)
		return 1;
			
	else if(strcmp(atom_name,"N")==0)
		return 2;
		
	else if((strcmp(atom_name,"C")==0) || (strcmp(res_name,"ASN")==0 && strcmp(atom_name,"CG")==0) || (strcmp(res_name,"GLN")==0&& strcmp(atom_name,"CD")==0))
		return 3;
		
	else if((strcmp(atom_name,"O")==0)|| (strcmp(res_name,"ASN")==0 && strcmp(atom_name,"OD1")==0) || (strcmp(res_name,"GLN")==0 && strcmp(atom_name,"OE1")==0) || (strcmp(res_name,"GLN")==0 && strcmp(atom_name,"OE")==0))
		return 4;
	
	else if(strcmp(atom_name,"CA")==0 && strcmp(res_name,"PRO")!=0)
		return 5;
	
	else if((strcmp(res_name,"ALA")==0 && strcmp(atom_name,"CB")==0) || (strcmp(res_name,"ILE")==0 && (strcmp(atom_name,"CG2")==0|| strcmp(atom_name,"CD1")==0)) || (strcmp(res_name,"LEU")==0&& (strcmp(atom_name,"CD1")==0|| strcmp(atom_name,"CD2")==0)) || (strcmp(res_name,"LEU")==0&& (strcmp(atom_name,"CD")==0|| strcmp(atom_name,"CE")==0)) || (strcmp(res_name,"THR")==0&& strcmp(atom_name,"CG2")==0) || (strcmp(res_name,"VAL")==0&& (strcmp(atom_name,"CG1")==0|| strcmp(atom_name,"CG2")==0)))
		return 6;
		
	else if((strcmp(res_name,"ARG")==0&& (strcmp(atom_name,"CB")==0|| strcmp(atom_name,"CG")==0)) || (strcmp(res_name,"ASN")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"GLN")==0&& (strcmp(atom_name,"CB")==0|| strcmp(atom_name,"CG")==0)) ||  (strcmp(res_name,"GLU")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"HIS")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"ILE")==0&& (strcmp(atom_name,"CB")==0|| strcmp(atom_name,"CG1")==0)) || (strcmp(res_name,"LEU")==0&& (strcmp(atom_name,"CB")==0|| strcmp(atom_name,"CG")==0)) || (strcmp(res_name,"LYS")==0&& (strcmp(atom_name,"CB")==0|| strcmp(atom_name,"CG")==0|| strcmp(atom_name,"CD")==0)) || (strcmp(res_name,"MET")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"PHE")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"PRO")==0&& (strcmp(atom_name,"CB")==0|| strcmp(atom_name,"CG")==0)) || (strcmp(res_name,"TRP")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"TYR")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"VAL")==0&& strcmp(atom_name,"CB")==0))
		return 7;
		
	else if((strcmp(res_name,"PHE")==0&& (strcmp(atom_name,"CG")==0|| strcmp(atom_name,"CD1")==0|| strcmp(atom_name,"CD2")==0|| strcmp(atom_name,"CE1")==0|| strcmp(atom_name,"CE2")==0|| strcmp(atom_name,"CZ")==0)) || (strcmp(res_name,"TRP")==0&& (strcmp(atom_name,"CE3")==0|| strcmp(atom_name,"CZ2")==0|| strcmp(atom_name,"CZ3")==0|| strcmp(atom_name,"CH2")==0)) || (strcmp(res_name,"TYR")==0&& (strcmp(atom_name,"CG")==0|| strcmp(atom_name,"CD1")==0|| strcmp(atom_name,"CD2")==0|| strcmp(atom_name,"CE1")==0|| strcmp(atom_name,"CE2")==0)) || (strcmp(res_name,"TYR")==0&& (strcmp(atom_name,"CC")==0|| strcmp(atom_name,"CD")==0|| strcmp(atom_name,"CE")==0|| strcmp(atom_name,"CG")==0|| strcmp(atom_name,"CH")==0)))
		return 8;
	
	else if(strcmp(res_name,"TYR")==0&& (strcmp(atom_name,"CZ")==0|| strcmp(atom_name,"CF")==0))
		return 9;
		
	else if((strcmp(res_name,"SER")==0&& strcmp(atom_name,"OG")==0) || (strcmp(res_name,"THR")==0&& strcmp(atom_name,"OG1")==0) || (strcmp(res_name,"TYR")==0&& strcmp(atom_name,"OH")==0))
		return 10;
		
	else if(strcmp(res_name,"TRP")==0&& (strcmp(atom_name,"CG")==0|| strcmp(atom_name,"CD2")==0))
		return 11;
		
	else if(strcmp(res_name,"TRP")==0&& (strcmp(atom_name,"CD1")==0|| strcmp(atom_name,"CE2")==0))
		return 12;
	
	else if(strcmp(res_name,"TRP")==0&& strcmp(atom_name,"NE1")==0)	
		return 13;
	
	else if(strcmp(res_name,"MET")==0&& strcmp(atom_name,"CG")==0)
		return 14;
		
	else if(strcmp(res_name,"MET")==0&& (strcmp(atom_name,"SD")==0|| strcmp(atom_name,"S")==0|| strcmp(atom_name,"SE")==0))
		return 15;
		
	else if (strcmp(res_name,"MET")==0&& strcmp(atom_name,"CE")==0)
		return 16;
	
	else if ((strcmp(res_name,"LYS")==0&& strcmp(atom_name,"CE")==0) || (strcmp(res_name,"LYS")==0&& strcmp(atom_name,"CZ")==0))
		return 17;
				
	else if ((strcmp(res_name,"SER")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"THR")==0&& strcmp(atom_name,"CB")==0))
		return 18;
		
	else if (strcmp(res_name,"PRO")==0&& (strcmp(atom_name,"CD")==0|| strcmp(atom_name,"CA")==0))
		return 19;
		
	else if (strcmp(res_name,"CYS")==0&& strcmp(atom_name,"CB")==0)
		return 20;
	
	else if (strcmp(res_name,"CYS")==0&& (strcmp(atom_name,"SG")==0|| strcmp(atom_name,"S")==0|| strcmp(atom_name,"SE")==0))
		return 21;
	
	else if (strcmp(res_name,"HIS")==0&& (strcmp(atom_name,"CG")==0|| strcmp(atom_name,"CD2")==0))
		return 22;
		
	else if (strcmp(res_name,"HIS")==0&& (strcmp(atom_name,"ND1")==0|| strcmp(atom_name,"NE2")==0))
		return 23;
		
	else if(strcmp(res_name,"HIS")==0&& strcmp(atom_name,"CE1")==0)
		return 24;
	
	else if(strcmp(res_name,"ARG")==0&& strcmp(atom_name,"CD")==0)
		return 25;
	
	else if(strcmp(res_name,"ARG")==0&& strcmp(atom_name,"NE")==0)
		return 26;
		
	else if(strcmp(res_name,"ARG")==0&& strcmp(atom_name,"CZ")==0)
		return 27;
		
	else if(strcmp(res_name,"ARG")==0&& (strcmp(atom_name,"NH1")==0|| strcmp(atom_name,"NH2")==0))
		return 28;
	
	else if((strcmp(res_name,"ASN")==0&& strcmp(atom_name,"ND2")==0) || (strcmp(res_name,"GLN")==0&& strcmp(atom_name,"NE2")==0))
		return 29;
		
	else if((strcmp(res_name,"ASP")==0&& strcmp(atom_name,"CB")==0) || (strcmp(res_name,"GLU")==0&& strcmp(atom_name,"CG")==0))
		return 30;
		
	else if((strcmp(res_name,"ASP")==0&& strcmp(atom_name,"CG")==0) || (strcmp(res_name,"GLU")==0&& strcmp(atom_name,"CD")==0) || (strcmp(res_name,"GLU")==0&& strcmp(atom_name,"CD1")==0))
		return 31;
		
	else if((strcmp(atom_name,"OXT")==0)|| (strcmp(res_name,"ASP")==0&& (strcmp(atom_name,"OD1")==0|| strcmp(atom_name,"OD2")==0)) || (strcmp(res_name,"GLU")==0&& (strcmp(atom_name,"OE1")==0|| strcmp(atom_name,"OE2")==0)) || (strcmp(res_name,"GLU")==0&& (strcmp(atom_name,"OE11")==0|| strcmp(atom_name,"OE21")==0)) || (strcmp(res_name,"GLU")==0&& strcmp(atom_name,"OE")==0) || strcmp(atom_name,"OX2")==0 )
		return 32;
	
	else
		return -1;		
	
	}
	
	
	int isStdAtomType32(int type)
	{       if(type>=1&&type<=32)
			return 1;
		else
			return 0;
	}
	
	int charBelongstoString(char str[],char c)
	{
	  int i;
	  for(i=0;i<strlen(str);i++)
	  {
	       if(str[i]==c)
		    return 1;
	  }
	  return 0;  
	  
	} 
	
	int getDistanceBin(double dist)
	{       int i;
		if(dist>=MAX_DISTANCE||dist<MIN_DISTANCE)  //erroneous values
			return -1;
		i=NUM_BINS_SPLINE-1; 
		while(!(dist>=r1Spline[i]&&dist<r2Spline[i])) // going backwards because the lower bins are less populated
			i--;
		return(i);
	   
	}
	
	int main(int argc, char *argv[])
	{

	char pdbFile[PDB_FILENAME_LEN],pdbReceptorChains[FIELD_LEN],pdbLigandChains[FIELD_LEN],paramFile[PDB_FILENAME_LEN]; 
		
	char str_buf[MAX_LINE_LEN],line[MAX_LINE_LEN],flname[PDB_FILENAME_LEN],ch;
	char atom_check[FIELD_LEN],atom_id[FIELD_LEN],atm_name[FIELD_LEN],res_name[FIELD_LEN],xcoord[FIELD_LEN],ycoord[FIELD_LEN],zcoord[FIELD_LEN],res_num_string[FIELD_LEN];
	long atomID;
	
	double x_pos,y_pos,z_pos,dist,xdif,ydif,zdif; 
	
	short typ, itype,jtype,temp;
	    
	FILE *fp;
	
	int resNum,i,j,r;
	
	int curr_atom_count,currProtein=0; /* currProtein is 0 or 1 according to whether current ATOM line is that of receptor or ligand. curr_atom_count is the number of atoms currently in the current protein list: receptor / ligand */	 
	
	int rec_atom_count=-1,lig_atom_count=-1; 
	/* rec_atom_count refers to the receptor atom count */
	
	double paramsAtomPot[NUM_ATOM_TYPES][NUM_ATOM_TYPES][NUM_BINS];  
	/* parameters of the atomic potential */
	
	double numContacts[NUM_ATOM_TYPES][NUM_ATOM_TYPES][NUM_BINS]; 
	/* number of 2-body contacts of each type */
	
	atom models[2][MAX_NUM_ATOMS]; 
	/* models[0] contains model receptor and models[1] contains ligand.  */
	
	double energy=0.0;
	/* Final energy value of the structure */
				
	if(argc!=5)
	{ printf("Error. usage: ./a.out <location of model> <receptor chains of model> <ligand chains of model> <location of parameter file>\n"); 
	  return(0);
        }
        
        /* Parse command line parameters */	
        strcpy(pdbFile,argv[1]);
	strcpy(pdbReceptorChains,argv[2]);
	strcpy(pdbLigandChains,argv[3]);
	strcpy(paramFile,argv[4]);
	
	/* Initialize param values and bin counts */
	for(i=0;i<NUM_ATOM_TYPES;i++)
	{	for(j=0;j<NUM_ATOM_TYPES;j++)
		{	
		     for(r=0;r<NUM_BINS;r++)
		     {
		        paramsAtomPot[i][j][r]=0.0;
			numContacts[i][j][r]=0.0;
		     }		
				
		}
	}
	
	/* Get the training parameters */
	if((fp = fopen(paramFile, "r"))== NULL) 
	{
		printf("Cannot open training data file.\n");
		return 0;
	}
	
	// Initialize the parameters matrix	
	for(i=0;i<NUM_ATOM_TYPES;i++)
	{	for(j=i;j<NUM_ATOM_TYPES;j++)
		{
		     for(r=0;r<NUM_BINS;r++)
		     {
		       fscanf(fp,"%lf",&paramsAtomPot[i][j][r]);
		     }  
		}
	}	
	
	fclose(fp);
	
	 /* parse the model file with their corresponding receptor and ligand chains */
	if((fp = fopen(pdbFile, "r"))== NULL) 
	{
	       printf("Cannot open file %s\n",pdbFile);
	       return 1;
	}
		    
	 /* Parse the PDB file of the given model */	       
	 while(fgets(str_buf, MAX_LINE_LEN, fp)!=NULL) 
	 {   
		    strcpy(line,str_buf); // all splicing operations are done on the string, line
		    
		    strncpy(atom_check,line,4);
		    atom_check[4]='\0';
		    
		    if(strcmp(atom_check,"ATOM")==0)
		    {         
			      /* get atom ID */
// 			      for(i=0,j=6;j<11;i++,j++)  
// 				   atom_id[i]=line[j];
// 			      atom_id[i]='\0';
// 			      
// 			      atomID=atol(atom_id);
// 			      
			      /* get atom name */
			      for(i=0,j=12;j<16;j++)
			      {	       if(line[j]!=' ')
					{	atm_name[i]=line[j];
						  i++;
					}
			      } 
			      atm_name[i]='\0';
			      
			      /* get residue name */
			      for(i=0,j=17;j<20;i++,j++)
				   res_name[i]=line[j];
			      res_name[i]='\0';
									   
			      /* get residue number */
			      for(i=0,j=22;j<26;i++,j++)
				   res_num_string[i]=line[j];
			      res_num_string[i]='\0';
			      
			      resNum=atoi(res_num_string);
			      
			      /* get atom type */
			      typ=getAtomType32(atm_name,res_name);	
			      				   
			      if(isStdAtomType32(typ)==0) /* check if it is one of the 32 heavy atom types */
				   continue;

						  
			      /* check which chain the atom belongs to: receptor or ligand, insert it to the appropriate set of atoms */
			       
			      ch=line[21];
			      
			      if (charBelongstoString(pdbReceptorChains,ch))
			      {    rec_atom_count+=1;
				   curr_atom_count=rec_atom_count; 
				   currProtein=0;			      
			      }
				   
				   
			      else if(charBelongstoString(pdbLigandChains,ch))
			      {    lig_atom_count+=1;
				   curr_atom_count=lig_atom_count; 
				   currProtein=1;
			      }   
			      
			      else  /* unnecessary atom */
				   continue;
									   
			      /* get the X,Y,Z coordinates */ 
			      for(i=0,j=30;j<38;i++,j++)
				   xcoord[i]=line[j];
			      xcoord[i]='\0';
			      
			      x_pos=atof(xcoord);
			      
			      for(i=0,j=38;j<46;i++,j++)
				   ycoord[i]=line[j];
			      ycoord[i]='\0';
			      
			      y_pos=atof(ycoord);
			      
			      for(i=0,j=46;j<54;i++,j++)
				   zcoord[i]=line[j];
			      zcoord[i]='\0';
			      
			      z_pos=atof(zcoord);
			      
			      //models[currProtein][curr_atom_count].atom_id=atomID;
			      strcpy(models[currProtein][curr_atom_count].atom_name,atm_name);
			      strcpy(models[currProtein][curr_atom_count].residue_name,res_name);
			      models[currProtein][curr_atom_count].chain=ch;
			      models[currProtein][curr_atom_count].residue_number=resNum;
			      models[currProtein][curr_atom_count].x=x_pos;
			      models[currProtein][curr_atom_count].y=y_pos;
			      models[currProtein][curr_atom_count].z=z_pos;
			      models[currProtein][curr_atom_count].atom_type=typ;
			      
			     		    
		    } //end ATOM line
	       
	 }  //end parsing in while loop
	       
	 fclose(fp);
	       
	 rec_atom_count+=1;
	 lig_atom_count+=1;
	 
	 /* Now get the contacts for the current model */	 
	       
	 for(i=0;i<rec_atom_count;i++)
	 {
		for(j=0;j<lig_atom_count;j++)
		{       			     
				xdif=models[0][i].x-models[1][j].x;
				ydif=models[0][i].y-models[1][j].y;
				zdif=models[0][i].z-models[1][j].z;
				
				if(abs(xdif)>MAX_DISTANCE || abs(ydif)>MAX_DISTANCE || abs(zdif)>MAX_DISTANCE  ) /* No need to calculate square root etc expensively if this condition is met */
				     continue;
				
				dist=sqrt(pow(xdif,2)+pow(ydif,2)+pow(zdif,2));	
				
				if(dist<MIN_DISTANCE || dist>MAX_DISTANCE)
				     continue;
				
				r=getDistanceBin(dist);  // distance bin numbering starts from 0. bin is  {0,1}  
						
					
				itype=models[0][i].atom_type;
				jtype=models[1][j].atom_type;
				
							
				if(jtype<itype) // swap itype and jtype, so that itype is always <= jtype
				{	temp=jtype;
					jtype=itype;
					itype=temp;
				}
								

				switch(r)
				{
				     case 0: { numContacts[itype-1][jtype-1][0]+=1.0; break; }
				     case 1: { numContacts[itype-1][jtype-1][0]+=4.0-dist; numContacts[itype-1][jtype-1][1]+=dist-3.0; break; }
				     case 2: { numContacts[itype-1][jtype-1][1]+=1.0; break;}
 				     case 3: { numContacts[itype-1][jtype-1][1]+=4.0-(dist/3.0); numContacts[itype-1][jtype-1][2]+=(dist/3.0)-3.0; break; }
				     case 4: { numContacts[itype-1][jtype-1][2]+=1.0; break; }
				     default: break;
				}
		}	
	} 


	for(i=0;i<NUM_ATOM_TYPES;i++)
	{	for(j=i;j<NUM_ATOM_TYPES;j++)
		{	
		     for(r=0;r<NUM_BINS;r++) 
		     {
			 energy+=numContacts[i][j][r]*paramsAtomPot[i][j][r];
		      }
		}
	}
	
	printf("%.4lf",energy);
	
	return 0;

	}
