#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


struct global_data
{
    int NpartX;
    double FactorZ;
    double BoxSize;
    double Density;
        
    char *output_fname;
    char *input_fname;
    
    int NX;
    int NY;
    int NZ;

    double dx;
    double dy;
    double dz;
        
} All;



struct io_header_1
{
    int npart[6];			/*!< number of particles of each type in this file */
    double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
                         stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;			/*!< time of snapshot file */
    double redshift;		/*!< redshift of snapshot file */
    int flag_sfr;			/*!< flags whether the simulation was including star formation */
    int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
    unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
                                 different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;		/*!< flags whether cooling was included  */
    int num_files;		/*!< number of files in multi-file snapshot */
    double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
    double Omega0;		/*!< matter density in units of critical density */
    double OmegaLambda;		/*!< cosmological constant parameter */
    double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
    int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
    int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
                             particles */
    unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
    int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
    int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */
    
    int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                   or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                   For snapshots files, the value informs whether the simulation was evolved from
                                   Zeldoch or 2lpt ICs. Encoding is as follows:
                                   FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                   FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                   FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                   FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                   FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                   All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                   */
    float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */
    char fill[18];		/*!< fills to 256 Bytes */    
    char names[15][2];
} header;



int     NumPart, Ngas, NumParttot, Ngastot, Nhalo, Nstar;
double  Mgas, Mhalo;


struct particle_data 
{
	float  Pos[3];
	float  Vel[3];
	float  Mass;
	int    Type;
	
	float  Rho, U, Temp, Ne, Nh, Hsml, Sfr, SAge, Z, Pot;
} *P;

int *Id;

void read_parameterfile(char *fname);


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
    
  All.input_fname = malloc(200*sizeof(char));
  All.output_fname = malloc(200*sizeof(char));
 


  if (argc < 2) {
    printf("Enter input filename: \n");
    scanf("%s", All.input_fname);
	//printf("Enter filename for output: \n");
    //scanf("%s", All.output_fname);
  }
  else {
    //printf("#reading from %s\n",argv[1]);
    sprintf(All.input_fname, "%s", argv[1]);
	//sprintf(All.output_fname, "%s", argv[1]);
  }  
  //printf("Enter filename for output: \n");
  //scanf("%s", All.output_fname);

  //load_file(input_fname, files);
  //printf("Particles: %d\n",NumPart);
  
    load_particles(All.input_fname);

    int i, j;
    j=0;

    // for my enzo, need M x y z vx vy vz
    //   in units of   Msun, pc, pc, pc, km/s, km/s, km/s
    //      convert kpc to pc below
//    for(i=Ngas+1;i<=Ngas+Nhalo;i++){ // get the DM
//    for (i=1; i <= Ngas; i ++){

    // write three files
    //   1) gas
    //   2) stars
    //   3) DM

    FILE *output;
    const double mu = 100.0 / 1.0122738558968995e-08; // mass units in solar masses - could not figure out how to get this natively from dataset - will need to verify that this works

    output = fopen("gas.dat", "w");
    for (i = 1; i <= Ngas; i++){
      fprintf(output, "%f %f %f %f %f %f %.8e\n", P[i].Pos[0]*1000.0, P[i].Pos[1]*1000.0, P[i].Pos[2]*1000.0,
                                                        P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], header.mass[0]*mu);
    }
    fclose(output);

    output = fopen("halo.dat", "w");
    for (i=Ngas+1; i <= Ngas+Nhalo; i++){
      fprintf(output, "%f %f %f %f %f %f %.8e \n", P[i].Pos[0]*1000.0, P[i].Pos[1]*1000.0, P[i].Pos[2]*1000.0, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], header.mass[1]*mu);
    }
    fclose(output);

    output = fopen("disk.dat", "w");
    for (i=Ngas+Nhalo+1; i <= Ngas+Nhalo+Nstar;i++){
      fprintf(output, "%f %f %f %f %f %f %.8e \n", P[i].Pos[0]*1000.0, P[i].Pos[1]*1000.0, P[i].Pos[2]*1000.0, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], header.mass[2]*mu);
    }
    fclose(output);






    //read_parameterfile(All.input_fname);
    
    //get_values();
    
    //set_particles();
    
    //define_header();

    //modify();
    //header.flag_entropy_instead_u = 0;


    //save_particles(All.output_fname);

}





/* define the values in the header
 */
int define_header(void)
{
    
    int i;

    header.npart[0] = NumPart;
    header.npart[1] = 0;
    header.npart[2] = 0;
    header.npart[3] = 0;
    header.npart[4] = 0;
    header.npart[5] = 0;
    
    header.mass[0] = 0.0;
    //header.mass[0] = Mgas;
    header.mass[1] = 0.0;
    header.mass[2] = 0.0;
    header.mass[3] = 0.0;
    header.mass[4] = 0.0;
    header.mass[5] = 0.0;
    
    header.time = 0.0;
    header.redshift = 0.0;
    header.flag_sfr = 0;
    header.flag_feedback = 0;

    header.npartTotal[0] = NumPart;
    header.npartTotal[1] = 0;
    header.npartTotal[2] = 0;
    header.npartTotal[3] = 0;
    header.npartTotal[4] = 0;
    header.npartTotal[5] = 0;
    
    header.flag_cooling = 0;
    header.num_files = 1;
    header.BoxSize = All.BoxSize;
    header.Omega0 = 0.0;
    header.OmegaLambda = 0.0;
    header.HubbleParam = 1.0;
    header.flag_stellarage = 0;
    header.flag_metals = 0;

    header.npartTotalHighWord[0] = 0;//(unsigned int) (NumPart >> 32);
    header.npartTotalHighWord[1] = 0;
    header.npartTotalHighWord[2] = 0;
    header.npartTotalHighWord[3] = 0;
    header.npartTotalHighWord[4] = 0;
    header.npartTotalHighWord[5] = 0;
    
    
    header.flag_entropy_instead_u = 0;
    header.flag_doubleprecision = 0;
    header.flag_ic_info = 0;
    header.lpt_scalingfactor = 1.0;

    return 0;
}



/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
	
  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");

    return 0;
}




int load_particles(char *fname)
{
	FILE *fd;
	int i,d,k;
	int nwithmass;
	float xyz[3];
	double t;
	int blklen;
#define SKIP fread(&blklen, sizeof(blklen), 1, fd);
	
	
	if(!(fd=fopen(fname,"r")))
    {
		printf("error opening file %s\n",fname);
		exit(0);
    }
	
	// printf("loading initial conditions from file `%s'\n\n",fname);
	
	blklen=sizeof(header);
	SKIP;
	fread(&header, sizeof(header), 1, fd);
	SKIP;
	
	NumPart = header.npart[0]+header.npart[1]+header.npart[2]+header.npart[3]+header.npart[4]+header.npart[5];
        printf("%.16e %.16e %.16e %.16e %.5e %.5e\n", header.mass[0], header.mass[1], header.mass[2], header.mass[3], header.mass[4], header.mass[5]);
	printf("NumPart = %d\n", NumPart);
	Ngas = header.npart[0];
	printf("Ngas = %d\n", Ngas);
	Nhalo = header.npart[1];
	printf("Nhalo = %d\n", Nhalo);
        Nstar = header.npart[2];
        printf("Nstar = %d\n", Nstar);

	P=malloc(NumPart*sizeof(struct particle_data));
	P--;   /* start with offset 1 */

	blklen=3*(NumPart)*sizeof(float);
	SKIP;
	for(i=1;i<=NumPart;i++)
    {
		//xyz[0]=P[i].Pos[0];
		//xyz[1]=P[i].Pos[1];
		//xyz[2]=P[i].Pos[2];
		//fwrite(xyz,sizeof(float),3,fd);
		fread(P[i].Pos, sizeof(float), 3, fd);
		//printf("%g\n",P[i].Pos[0]);
    }      
	SKIP;
	
	blklen=3*(NumPart)*sizeof(float);
	SKIP;
	for(i=1;i<=NumPart;i++)
    {
		//xyz[0]=P[i].Vel[0];
		//xyz[1]=P[i].Vel[1];
		//xyz[2]=P[i].Vel[2];
		//fwrite(xyz,sizeof(float),3,fd);
		fread(P[i].Vel, sizeof(float), 3, fd);
		//printf("%g\n",P[i].Vel[0]);
    }
	SKIP;

        blklen=(NumPart)*sizeof(float);
        SKIP;
        for(i=1;i<=NumPart;i++)
    {
                //xyz[0]=P[i].Vel[0];
                //xyz[1]=P[i].Vel[1];
                //xyz[2]=P[i].Vel[2];
                //fwrite(xyz,sizeof(float),3,fd);
                fread(&P[i].Mass, sizeof(float), 1, fd);
                //printf("%g\n",P[i].Vel[0]);
    }
        SKIP;

	
	
	blklen=(NumPart)*sizeof(int);
	SKIP;
	for(i=1;i<=NumPart;i++)
    {
		//fwrite(&i,sizeof(int),1,fd);  // ID
		fread(xyz, sizeof(int), 1, fd);
		
    }
	SKIP;
	
	
	nwithmass=0;
	for(k=0; k<6; k++)
    {
		if (header.mass[k]==0) nwithmass+= header.npart[k];
    }
	
	blklen=nwithmass*sizeof(float);
	if(nwithmass>0) SKIP;
	
	if(header.mass[0]==0)
    {
		for(i=1;i<=Ngas;i++)
        {
			//xyz[0]=P[i].Mass;
			fread(&P[i].Mass,sizeof(float),1,fd);
			//printf("Mass=%g\n !!!",P[i].Mass);
        }
    }
        if(header.mass[1]==0)
    {
                for(i=Ngas+1;i<=Ngas+Nhalo;i++)
        {
                        //xyz[0]=P[i].Mass;
                        fread(&P[i].Mass,sizeof(float),1,fd);
                        //printf("Mass=%g\n !!!",P[i].Mass);
        }
    }
        if(header.mass[2]==0)
    {
                for(i=Ngas+Nhalo+1;i<=Ngas+Nhalo+Nstar;i++)
        {
                        //xyz[0]=P[i].Mass;
                        fread(&P[i].Mass,sizeof(float),1,fd);
                        //printf("Mass=%g\n !!!",P[i].Mass);
        }
    }

	
	
	if(nwithmass>0) SKIP;
	
	if(Ngas)
    {
		blklen=(Ngas)*sizeof(float);
		SKIP;
		for(i=1;i<=Ngas;i++)
		{
			//xyz[0]= P[i].U;
			//fwrite(xyz, sizeof(float), 1, fd);
			fread(&(P[i].U), sizeof(float), 1, fd);
			//printf("%g\n",P[i].U);
		}
		SKIP;
		
		SKIP;
		for(i=1;i<=Ngas;i++)
		{
			//xyz[0]= P[i].Rho;
			//fwrite(xyz, sizeof(float), 1, fd);
			fread(&(P[i].Rho), sizeof(float), 1, fd);
			//printf("%g\n",P[i].Rho);
		}
		SKIP;
		
		if(header.flag_cooling){
			
			SKIP;
			for(i=1;i<=Ngas;i++)
			{
				//xyz[0]= P[i].Ne;
				//fwrite(xyz, sizeof(float), 1, fd);
				fread(&(P[i].Ne), sizeof(float), 1, fd);
			}
			SKIP;
			
			SKIP;
			for(i=1;i<=Ngas;i++)
			{
				//xyz[0]= P[i].Nh;
				//fwrite(xyz, sizeof(float), 1, fd);
				fread(&(P[i].Nh), sizeof(float), 1, fd);
			}
			SKIP;
			
		}
		
		SKIP;
		for(i=1;i<=Ngas;i++)
		{
			//xyz[0]= P[i].Hsml;
			//fwrite(xyz, sizeof(float), 1, fd);
			fread(&(P[i].Hsml), sizeof(float), 1, fd);
			//printf("%g\n",P[i].Hsml);
		}
		SKIP;
		
		if(header.flag_sfr){
			
			SKIP;
			for(i=1;i<=Ngas;i++)
			{
				//xyz[0]= P[i].Sfr;
				//fwrite(xyz, sizeof(float), 1, fd);
				fread(&(P[i].Sfr), sizeof(float), 1, fd);
			}
			SKIP;
			
		}
		
    }
	

     
	fclose(fd);
    
    return 0;
}






int modify(){

  int i, j;
  float gamma = 5./3.;
  float Mach = 0.1;
  float P0 = 1.0 / (gamma * Mach*Mach);

  j=0;
  for(i=1;i<=Ngas;i++){    
    float du = (P0 - 5.0) / (gamma - 1.0) / P[i].Rho;
    j++;
    printf("u = %g  du = %g  \n", P[i].U, du);
    P[i].U += du;
    printf("unew = %g \n", P[i].U);

  }
  printf("j = %d\n",j);

  return 0;
}


int save_particles(char *fname)
{
	FILE *fd;
	int i,d,k;
	int nwithmass;
	float xyz[3];
	double t;
	int blklen;
#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);
	
	
	if(!(fd=fopen(fname,"w")))
    {
		printf("error opening file %s\n",fname);
		exit(0);
    }
	
	printf("saveing initial conditions to file `%s'\n\n",fname);
	
	blklen=sizeof(header);
	BLKLEN;
	fwrite(&header, sizeof(header), 1, fd);
	BLKLEN;
	
	
	blklen=3*(NumPart)*sizeof(float);
	BLKLEN;
	for(i=1;i<=NumPart;i++)
    {
		xyz[0]=P[i].Pos[0];
		xyz[1]=P[i].Pos[1];
		xyz[2]=P[i].Pos[2];
		fwrite(xyz,sizeof(float),3,fd);
    }      
	BLKLEN;
	
	blklen=3*(NumPart)*sizeof(float);
	BLKLEN;
	for(i=1;i<=NumPart;i++)
    {
		xyz[0]=P[i].Vel[0];
		xyz[1]=P[i].Vel[1];
		xyz[2]=P[i].Vel[2];
		fwrite(xyz,sizeof(float),3,fd);
    }
	BLKLEN;
	
	
	blklen=(NumPart)*sizeof(int);
	BLKLEN;
	for(i=1;i<=NumPart;i++)
    {
		fwrite(&i,sizeof(int),1,fd);  // ID
    }
	BLKLEN;
	
	
	nwithmass=0;
	for(k=0; k<6; k++)
    {
		if (header.mass[k]==0) nwithmass+= header.npart[k];
    }
	
	blklen=nwithmass*sizeof(float);
	if(nwithmass>0) BLKLEN;
	
	if(header.mass[0]==0)
    {
		for(i=1;i<=Ngas;i++)
        {
	  if(P[i].Pos[1] < 64. && P[i].Pos[1] > -64.)  
	    xyz[0]=P[i].Mass * 4.;
	  else
	    xyz[0]=P[i].Mass;
	  //printf("mass=%g\n",xyz[0]);
			fwrite(xyz,sizeof(float),1,fd);
        }
    }

	/*
	if(header.mass[1]==0)
    {
		for(i=Ngas+1;i<=Ngas+Nhalo;i++)
        {
			xyz[0]=P[i].Mass;
			fwrite(xyz,sizeof(float),1,fd);
        }
    }
	
	if(header.mass[2]==0)
    {
		for(i=Ngas+Nhalo+1;i<=Ngas+Nhalo+Ndisk;i++)
        {
			xyz[0]=P[i].Mass;
			fwrite(xyz,sizeof(float),1,fd);
        }
    }
	
	if(header.mass[3]==0)
    {
		for(i=Ngas+Nhalo+Ndisk+1;i<=Ngas+Nhalo+Ndisk+Nbulge;i++)
        {
			xyz[0]=P[i].Mass;
			fwrite(xyz,sizeof(float),1,fd);
        }
    }
	
	if(header.mass[4]==0)
    {
		for(i=Ngas+Nhalo+Ndisk+Nbulge+1;i<=Ngas+Nhalo+Ndisk+Nbulge+Nstars;i++)
        {
			xyz[0]=P[i].Mass;
			fwrite(xyz,sizeof(float),1,fd);
        }
    }
     */
	
	if(nwithmass>0) BLKLEN;
	
	if(Ngas)
    {
		blklen=(Ngas)*sizeof(float);
		BLKLEN;
		for(i=1;i<=Ngas;i++)
		{
		  xyz[0]= P[i].U;
		  fwrite(xyz, sizeof(float), 1, fd);
		}
		BLKLEN;
		
		BLKLEN;
		for(i=1;i<=Ngas;i++)
		{
			xyz[0]= P[i].Rho;
			fwrite(xyz, sizeof(float), 1, fd);
		}
		BLKLEN;
		
		if(header.flag_cooling){
			
			BLKLEN;
			for(i=1;i<=Ngas;i++)
			{
				xyz[0]= P[i].Ne;
				fwrite(xyz, sizeof(float), 1, fd);
			}
			BLKLEN;
			
			BLKLEN;
			for(i=1;i<=Ngas;i++)
			{
				xyz[0]= P[i].Nh;
				fwrite(xyz, sizeof(float), 1, fd);
			}
			BLKLEN;
			
		}
		
		BLKLEN;
		for(i=1;i<=Ngas;i++)
		{
			xyz[0]= P[i].Hsml;
			fwrite(xyz, sizeof(float), 1, fd);
		}
		BLKLEN;
		
		if(header.flag_sfr){
			
			BLKLEN;
			for(i=1;i<=Ngas;i++)
			{
				xyz[0]= P[i].Sfr;
				fwrite(xyz, sizeof(float), 1, fd);
			}
			BLKLEN;
			
		}
		
    }
	
    /*
	if(Nstars)
    {
		blklen=(Nstars)*sizeof(float);
		BLKLEN;
		for(i=1+Ngas+Nhalo+Ndisk+Nbulge;i<=Ngas+Nhalo+Ndisk+Nbulge+Nstars;i++)
		{
			xyz[0]= P[i].SAge;
			fwrite(xyz, sizeof(float), 1, fd);
		}
		BLKLEN;
	}	
	
	if((Ngas)||(Nstars))
    {
		blklen=(Ngas+Nstars)*sizeof(float);
		BLKLEN;
		for(i=1;i<=Ngas;i++)
		{
			xyz[0]= P[i].Z;
			fwrite(xyz, sizeof(float), 1, fd);
		}
		for(i=1+Ngas+Nhalo+Ndisk+Nbulge;i<=Ngas+Nhalo+Ndisk+Nbulge+Nstars;i++)
		{
			xyz[0]= P[i].Z;
			fwrite(xyz, sizeof(float), 1, fd);
		}		
		BLKLEN;
	}
     */
     
	fclose(fd);
    
    return 0;
}

void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300
    
    FILE *fd;
    char buf[200], buf1[200], buf2[200], buf3[200];
    int i, j, nt;
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];
    int errorFlag = 0;
    
    // read parameter file on all processes for simplicty
    
    nt = 0;
    
    printf("\nReading parameter file... ");
    
    strcpy(tag[nt], "NpartX");
    addr[nt] = &All.NpartX;
    id[nt++] = INT;
    
    strcpy(tag[nt], "FactorZ");
    addr[nt] = &All.FactorZ;
    id[nt++] = FLOAT;
        
    strcpy(tag[nt], "BoxSize");
    addr[nt] = &All.BoxSize;
    id[nt++] = FLOAT;
        
    strcpy(tag[nt], "Density");
    addr[nt] = &All.Density;
    id[nt++] = FLOAT;
    
    strcpy(tag[nt], "OutputFile");
    addr[nt] = All.output_fname;
    id[nt++] = STRING;
    
    if((fd = fopen(fname, "r")))
    {
        while(!feof(fd))
        {
            buf[0] = 0;
            fgets(buf, 200, fd);
            
            if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                continue;
            
            if(buf1[0] == '%')
                continue;
            
            for(i = 0, j = -1; i < nt; i++)
                if(strcmp(buf1, tag[i]) == 0)
                {
                    j = i;
                    tag[i][0] = 0;
                    break;
                }
            
            if(j >= 0)
            {
                switch (id[j])
                {
                    case FLOAT:
                        *((double *) addr[j]) = atof(buf2);
                        break;
                    case STRING:
                        strcpy(addr[j], buf2);
                        break;
                    case INT:
                        *((int *) addr[j]) = atoi(buf2);
                        break;
                }
            }
            else
            {
                fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
                errorFlag = 1;
            }
        }
        fclose(fd);
        
    }
    else
    {
        fprintf(stdout, "Parameter file %s not found.\n", fname);
        errorFlag = 1;
    }
    
    
    for(i = 0; i < nt; i++)
    {
        if(*tag[i])
        {
            fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
            errorFlag = 1;
        }
    }
    
    if(errorFlag)
    {
        exit(1);
    }
    
    printf("done\n");
    
    
#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}



