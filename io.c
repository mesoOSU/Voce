#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "evp.h"

//typedef enum{N_I, N_R} VType;
typedef enum{N_I, N_R, N_S} VType;

typedef struct{
	char *vName;
	void *vPtr;
	VType vType;
	int vLen;
	int vStatus;
}NameList;

#define NameI(x) {#x, &x, N_I, sizeof(x)/sizeof(int)}
#define NameR(x) {#x, &x, N_R, sizeof(x)/sizeof(real)}
#define NameS(x) {#x, &x, N_S, 1}
#define NP_I ((int *)(nameList[k].vPtr) + j)
#define NP_R ((real *)(nameList[k].vPtr) + j)
#define NP_S ((char *)(nameList[k].vPtr) + j)

NameList nameList[] = {
	NameI(CellDim),
	NameI(N_phases),
	NameI(Type_phases),
	NameR(TimeStep),
	NameI(N_steps),
	NameI(Update_Flag),
	NameI(Hard_Flag),
	NameI(Tex_Flag),
	NameI(IterMax),
	NameI(PrintControl),
	NameR(Err),
	NameI(VelGrad_BC_Flag),
	NameR(VelGrad_BC),
	NameI(Stress_BC_Flag),
	NameR(Stress_BC),
	NameR(ElastConst_PhaseI),
	NameR(ElastConst_PhaseII),
	NameS(Slip_PhaseI),
	NameS(Slip_PhaseII),
	NameS(initial_ms),
	NameR(self),
	NameR(coplane),
	NameR(cross),
	NameR(glissile),
	NameR(hirth),
	NameR(lomer),
};

int GetNameList(char **argv){
	int j, k, match, ok;
	char buff[80], *token;
	FILE *fp;
	strcpy(buff, argv[1]);
	if((fp = fopen(buff, "r")) == 0){
		printf("PE#%d: Input file not found!\n", rank);
		return(1);
	}
	for(k=0; k<sizeof(nameList)/sizeof(NameList); k++)
		nameList[k].vStatus = 0;
	ok = 1;
	while(1){
		fgets(buff, 80, fp);
		if(feof(fp))
			break;
		if(buff[0]=='#'||buff[0]=='\n')
			continue;
		token = strtok(buff, " \t\n");
		if(!token)
			break;
		match = 0;
		for(k = 0; k<sizeof(nameList)/sizeof(NameList); k++){
			if(strcmp(token, nameList[k].vName) == 0){
					match = 1;
					if(nameList[k].vStatus == 0){
						nameList[k].vStatus = 1;
						for(j=0; j<nameList[k].vLen; j++){
							token = strtok(NULL, ", \t\n");
							if(token){
								switch(nameList[k].vType){
									case N_I:
										*NP_I = atol(token);
										break;
									case N_R:
										*NP_R = atof(token);
										break;
									case N_S:
										strcpy(NP_S,token);
										break;
								}
							}
							else{
								nameList[k].vStatus = 2;
								ok = 0;
							}
						}
						token = strtok(NULL, ", \t\n");
						if(token){
							nameList[k].vStatus = 3;
							ok = 0;
						}
						break;
					}
					else{
						nameList[k].vStatus = 4;
						ok = 0;
					}
			}
		}
		if(!match) ok = 0;
	}
	fclose(fp);
	for(k=0; k<sizeof(nameList)/sizeof(NameList); k++){
		if(nameList[k].vStatus != 1)
			ok = 1;
	}
	return(ok);
}/* end GetNameList() */

void PrintNameList(){
	int j, k;

	printf("===================================================\n");
	printf("-----------NameList for input data-----------------\n\n");
	for(k=0; k<sizeof(nameList)/sizeof(NameList); k++){
		printf("%s\t", nameList[k].vName);
		if(strlen(nameList[k].vName) < 8) printf("\t");
		if(nameList[k].vStatus > 0){
			for(j=0; j<nameList[k].vLen; j++){
				switch(nameList[k].vType){
					case N_I:
						printf("%d ", *NP_I);
						break;
					case N_R:
						printf("%#e ", *NP_R);
						break;
					case N_S:
						printf("%s ", NP_S);
						break;
				}
			}
		}
		switch(nameList[k].vStatus){
			case 0:
				printf("** no data");
				break;
			case 1:
				break;
			case 2:
				printf("** missing data");
				break;
			case 3:
				printf("** extra data");
				break;
			case 4:
				printf("** multiply defined");
				break;
		}
		printf("\n");
	}
	printf("\n-------------------End of NameList-----------------\n");
	printf("===================================================\n\n");

	return;
}/* end PrintNameList() */

void OpenFiles(void)
{
	if(rank!=0){
		PError("ERROR: Only rank#0 is allowed to modify global files",28);
	}
	fp_vm = fopen("vm.out","w+");
	fp_err = fopen("err.out","w+");

	return;
}/*end OpenFiles()*/

void CloseFiles(void)
{
	if(rank!=0){
		PError("ERROR: Only rank#0 is allowed to modify global files",28);
	}
	fclose(fp_vm);
	fclose(fp_err);

	return;
}/*end CloseFiles()*/

void WriteEpsMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Eps[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Eps[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Eps[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Eps[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%04d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteEpsMPI()*/
void WriteWorkEMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;
	voigt els;
	ten2nd els33;
	voigt66 aux_66; 
	ten4th aux_3333;

	AllocMem(tmpVector, lsize, real);
	

		local_loop{
		    tmpVector[pIDX]=0.0;
		    	for(i=0;i<6;i++){
			if(i<3){
				els[i] = DisGrad[pIDX][i][i] - Eps[pIDX][i][i];
			}
			else if(i==3){
				els[i] = (DisGrad[pIDX][1][2]+DisGrad[pIDX][2][1])/2.0 - Eps[pIDX][1][2];
			}
			else if(i==4){
				els[i] = (DisGrad[pIDX][0][2]+DisGrad[pIDX][2][0])/2.0 - Eps[pIDX][0][2];
			}
			else if(i==5){
				els[i] = (DisGrad[pIDX][0][1]+DisGrad[pIDX][1][0])/2.0 - Eps[pIDX][0][1];
			}
		}
	    	chg_basis(els,els33,aux_66,aux_3333,1);
	    	T2_loop{
	    	tmpVector[pIDX]+=0.5*Sig[pIDX][mi][mj]*els33[mi][mj];
	    	}
		}
		sprintf(fname, "%s_S%06d.iout", s,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	

	free(tmpVector);

  return;
}/* WriteWorkEMPI()*/

void WriteElsMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = DisGrad[pIDX][i][i] - Eps[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = (DisGrad[pIDX][1][2]+DisGrad[pIDX][2][1])/2.0 - Eps[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = (DisGrad[pIDX][0][2]+DisGrad[pIDX][2][0])/2.0 - Eps[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = (DisGrad[pIDX][0][1]+DisGrad[pIDX][1][0])/2.0 - Eps[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%06d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

  return;
}/* WriteElsMPI()*/



void WriteSigMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Sig[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Sig[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Sig[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Sig[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%04d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteSigMPI()*/

void WriteEdotMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Edot[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Edot[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Edot[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Edot[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%04d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteEdotMPI()*/
void WriteRhoCRMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	

		local_loop{
				tmpVector[pIDX] = rho_crss[pIDX];
}
		
	
		sprintf(fname, "%s_S%04d.iout", s ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
//	}

	free(tmpVector);

	return;
}/*end WriteRhoMPI()*/
void WriteRhoRMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	

		local_loop{
				tmpVector[pIDX] = rho_rss[pIDX];
}
		
	
		sprintf(fname, "%s_S%04d.iout", s ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
//	}

	free(tmpVector);

	return;
}/*end WriteRhoMPI()*/

/*void WriteSLIPMPI(char *s, int step)
{
	typedef struct{
		//real angle[3];	// Euler angles
		//int coord[3];		// coordinates
		int jgr;			// grain type
		//int jph;			//phase type
		//real shear;
                real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;
                // buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	//real t1,ph,t2;
	//ten2nd sa2xt;
	char fname[100] = {0};
	MPI_File fp;
	int i;

	int count = 25;	
	int block_lens[25] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[25];
	MPI_Datatype old_types[25] = {MPI_INT,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	
	//MPI_Address(&value.coord, &indices[0]);
	//MPI_Address(value.coord, &indices[1]);
	MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[1]);
        MPI_Address(&value.SS2, &indices[2]);
        MPI_Address(&value.SS3, &indices[3]);
        MPI_Address(&value.SS4, &indices[4]);
        MPI_Address(&value.SS5, &indices[5]);
        MPI_Address(&value.SS6, &indices[6]);
        MPI_Address(&value.SS7, &indices[7]);
        MPI_Address(&value.SS8, &indices[8]);
        MPI_Address(&value.SS9, &indices[9]);
        MPI_Address(&value.SS10, &indices[10]);
        MPI_Address(&value.SS11, &indices[11]);
        MPI_Address(&value.SS12, &indices[12]);
        MPI_Address(&value.SS13, &indices[13]);
        MPI_Address(&value.SS14, &indices[14]);
        MPI_Address(&value.SS15, &indices[15]);
        MPI_Address(&value.SS16, &indices[16]);
        MPI_Address(&value.SS17, &indices[17]);
        MPI_Address(&value.SS18, &indices[18]);
        MPI_Address(&value.SS19, &indices[19]);
        MPI_Address(&value.SS20, &indices[20]);
        MPI_Address(&value.SS21, &indices[21]);
        MPI_Address(&value.SS22, &indices[22]);
        MPI_Address(&value.SS23, &indices[23]);
        MPI_Address(&value.SS24, &indices[24]);
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		//T2_loop{
		//	sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		//}
		//EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		//tmpVector[pIDX].angle[0] = t1;
		//tmpVector[pIDX].angle[1] = ph;
		//tmpVector[pIDX].angle[2] = t2;
		//tmpVector[pIDX].coord[0] = px+lxs+1;
		//tmpVector[pIDX].coord[1] = py+1;
		//tmpVector[pIDX].coord[2] = pz+1;
		tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = gammatotal[pIDX][0];
                tmpVector[pIDX].SS2 = gammatotal[pIDX][1];
                tmpVector[pIDX].SS3 = gammatotal[pIDX][2];
                tmpVector[pIDX].SS4 = gammatotal[pIDX][3];
                tmpVector[pIDX].SS5 = gammatotal[pIDX][4];
                tmpVector[pIDX].SS6 = gammatotal[pIDX][5];
                tmpVector[pIDX].SS7 = gammatotal[pIDX][6];
                tmpVector[pIDX].SS8 = gammatotal[pIDX][7];
                tmpVector[pIDX].SS9 = gammatotal[pIDX][8];
                tmpVector[pIDX].SS10 = gammatotal[pIDX][9];
                tmpVector[pIDX].SS11 = gammatotal[pIDX][10];
                tmpVector[pIDX].SS12 = gammatotal[pIDX][11];
                tmpVector[pIDX].SS13 = gammatotal[pIDX][12];
                tmpVector[pIDX].SS14 = gammatotal[pIDX][13];
                tmpVector[pIDX].SS15 = gammatotal[pIDX][14];
                tmpVector[pIDX].SS16 = gammatotal[pIDX][15];
                tmpVector[pIDX].SS17 = gammatotal[pIDX][16];
                tmpVector[pIDX].SS18 = gammatotal[pIDX][17];
                tmpVector[pIDX].SS19 = gammatotal[pIDX][18];
                tmpVector[pIDX].SS20 = gammatotal[pIDX][19];
                tmpVector[pIDX].SS21 = gammatotal[pIDX][20];
                tmpVector[pIDX].SS22 = gammatotal[pIDX][21];
                tmpVector[pIDX].SS23 = gammatotal[pIDX][22];
                tmpVector[pIDX].SS24 = gammatotal[pIDX][23];
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}*//*end WriteTextureMPI()*/
void WriteSVMMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	

		local_loop{
				tmpVector[pIDX] = SVM[pIDX];
}
		
	
		sprintf(fname, "%s_S%04d.iout", s ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
//	}

	free(tmpVector);

	return;
}/*end WriteRhoMPI()*/
void WriteEVMMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	

		local_loop{
				tmpVector[pIDX] = EVM[pIDX];
}
		
	
		sprintf(fname, "%s_S%04d.iout", s ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
//	}

	free(tmpVector);

	return;
}/*end WriteEVMMPI()*/
void WriteWorkPMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;
	

	AllocMem(tmpVector, lsize, real);
	
	
		local_loop{
		tmpVector[pIDX]=work_p[pIDX];
		}
		sprintf(fname, "%s_S%06d.iout", s ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, rank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	

	free(tmpVector);

	return;
}/*end WriteWorkPMPI()*/
//cp
void WriteSLIPMPI(char *s, int step)
{
	typedef struct{
			// Euler angles
			// coordinates
		int jgr;
		 real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SS12;//SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;// grain type
		// real buff;//SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;		// buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	char fname[100] = {0};
	MPI_File fp;
	int i;

	int count = 13;	
	int block_lens[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[13];
	MPI_Datatype old_types[13] = {MPI_INT,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	

	MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[1]);
        MPI_Address(&value.SS2, &indices[2]);
        MPI_Address(&value.SS3, &indices[3]);
        MPI_Address(&value.SS4, &indices[4]);
        MPI_Address(&value.SS5, &indices[5]);
        MPI_Address(&value.SS6, &indices[6]);
        MPI_Address(&value.SS7, &indices[7]);
        MPI_Address(&value.SS8, &indices[8]);
        MPI_Address(&value.SS9, &indices[9]);
        MPI_Address(&value.SS10, &indices[10]);
        MPI_Address(&value.SS11, &indices[11]);
        MPI_Address(&value.SS12, &indices[12]);
	
	
	
	
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = gammatotal[pIDX][0];
                tmpVector[pIDX].SS2 = gammatotal[pIDX][1];
                tmpVector[pIDX].SS3 = gammatotal[pIDX][2];
                tmpVector[pIDX].SS4 = gammatotal[pIDX][3];
                tmpVector[pIDX].SS5 = gammatotal[pIDX][4];
                tmpVector[pIDX].SS6 = gammatotal[pIDX][5];
                tmpVector[pIDX].SS7 = gammatotal[pIDX][6];
                tmpVector[pIDX].SS8 = gammatotal[pIDX][7];
                tmpVector[pIDX].SS9 = gammatotal[pIDX][8];
                tmpVector[pIDX].SS10 = gammatotal[pIDX][9];
                tmpVector[pIDX].SS11 = gammatotal[pIDX][10];
                tmpVector[pIDX].SS12 = gammatotal[pIDX][11];
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteTextureMPI()*/
void WritetauxMPI(char *s, int step)
{
	typedef struct{
			// Euler angles
			// coordinates
		int jgr;
		 real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SS12;//SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;// grain type
		// real buff;//SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;		// buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	char fname[100] = {0};
	MPI_File fp;
	int i;

	int count = 13;	
	int block_lens[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[13];
	MPI_Datatype old_types[13] = {MPI_INT,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	

	MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[1]);
        MPI_Address(&value.SS2, &indices[2]);
        MPI_Address(&value.SS3, &indices[3]);
        MPI_Address(&value.SS4, &indices[4]);
        MPI_Address(&value.SS5, &indices[5]);
        MPI_Address(&value.SS6, &indices[6]);
        MPI_Address(&value.SS7, &indices[7]);
        MPI_Address(&value.SS8, &indices[8]);
        MPI_Address(&value.SS9, &indices[9]);
        MPI_Address(&value.SS10, &indices[10]);
        MPI_Address(&value.SS11, &indices[11]);
        MPI_Address(&value.SS12, &indices[12]);
	
	
	
	
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = tau_x[pIDX][0];
                tmpVector[pIDX].SS2 = tau_x[pIDX][1];
                tmpVector[pIDX].SS3 = tau_x[pIDX][2];
                tmpVector[pIDX].SS4 = tau_x[pIDX][3];
                tmpVector[pIDX].SS5 = tau_x[pIDX][4];
                tmpVector[pIDX].SS6 = tau_x[pIDX][5];
                tmpVector[pIDX].SS7 = tau_x[pIDX][6];
                tmpVector[pIDX].SS8 = tau_x[pIDX][7];
                tmpVector[pIDX].SS9 = tau_x[pIDX][8];
                tmpVector[pIDX].SS10 = tau_x[pIDX][9];
                tmpVector[pIDX].SS11 = tau_x[pIDX][10];
                tmpVector[pIDX].SS12 = tau_x[pIDX][11];
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteTextureMPI()*/

void WriteNewPositionMPI(char *s, int step)
{
	typedef struct{
        real SS1,SS2,SS3;
                // buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	char fname[200] = {0};
	MPI_File fp;
	int i;

	int count = 3;	
	int block_lens[3] = {1,1,1};
	MPI_Aint indices[3];
	MPI_Datatype old_types[3] = {MPI_real,MPI_real,MPI_real};
	
        MPI_Address(&value.SS1, &indices[0]);
        MPI_Address(&value.SS2, &indices[1]);
        MPI_Address(&value.SS3, &indices[2]);
       
       
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{

                tmpVector[pIDX].SS1 = displacement_fluct[pIDX][0];
                tmpVector[pIDX].SS2 = displacement_fluct[pIDX][1];
                tmpVector[pIDX].SS3 = displacement_fluct[pIDX][2];
               
               
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end NewPosition()*/


void WriteCRSSMPI(char *s, int step)
{
	typedef struct{
			// Euler angles
			// coordinates
		int jgr;
		 real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SS12;//SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;// grain type
		// real buff;//SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,SS21,SS22,SS23,SS24;		// buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	char fname[100] = {0};
	MPI_File fp;
	int i;

	int count = 13;	
	int block_lens[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[13];
	MPI_Datatype old_types[13] = {MPI_INT,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	

	MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[1]);
        MPI_Address(&value.SS2, &indices[2]);
        MPI_Address(&value.SS3, &indices[3]);
        MPI_Address(&value.SS4, &indices[4]);
        MPI_Address(&value.SS5, &indices[5]);
        MPI_Address(&value.SS6, &indices[6]);
        MPI_Address(&value.SS7, &indices[7]);
        MPI_Address(&value.SS8, &indices[8]);
        MPI_Address(&value.SS9, &indices[9]);
        MPI_Address(&value.SS10, &indices[10]);
        MPI_Address(&value.SS11, &indices[11]);
        MPI_Address(&value.SS12, &indices[12]);
	
	
	
	
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = crss[pIDX][0][0];
                tmpVector[pIDX].SS2 = crss[pIDX][1][0];
                tmpVector[pIDX].SS3 = crss[pIDX][2][0];
                tmpVector[pIDX].SS4 = crss[pIDX][3][0];
                tmpVector[pIDX].SS5 = crss[pIDX][4][0];
                tmpVector[pIDX].SS6 = crss[pIDX][5][0];
                tmpVector[pIDX].SS7 = crss[pIDX][6][0];
                tmpVector[pIDX].SS8 = crss[pIDX][7][0];
                tmpVector[pIDX].SS9 = crss[pIDX][8][0];
                tmpVector[pIDX].SS10 = crss[pIDX][9][0];
                tmpVector[pIDX].SS11 = crss[pIDX][10][0];
                tmpVector[pIDX].SS12 = crss[pIDX][11][0];
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteTextureMPI()*/
//

void WriteTextureMPI(char *s, int step)
{
	typedef struct{
		real angle[3];	// Euler angles
		int coord[3];		// coordinates
		int jgr;			// grain type
		int jph;			//phase type
		int buff;			// buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	real t1,ph,t2;
	ten2nd sa2xt;
	char fname[100] = {0};
	MPI_File fp;
	int i;

	int count = 5;	
	int block_lens[5] = {3,3,1,1,1};
	MPI_Aint indices[5];
	MPI_Datatype old_types[5] = {MPI_real,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
	
	MPI_Address(&value, &indices[0]);
	MPI_Address(value.coord, &indices[1]);
	MPI_Address(&value.jgr, &indices[2]);
	MPI_Address(&value.jph, &indices[3]);
	MPI_Address(&value.buff, &indices[4]);
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		T2_loop{
			sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		}
		EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		tmpVector[pIDX].angle[0] = t1;
		tmpVector[pIDX].angle[1] = ph;
		tmpVector[pIDX].angle[2] = t2;
		tmpVector[pIDX].coord[0] = px+lxs+1;
		tmpVector[pIDX].coord[1] = py+1;
		tmpVector[pIDX].coord[2] = pz+1;
		tmpVector[pIDX].jgr = grain_f[pIDX];
		tmpVector[pIDX].jph = phase_f[pIDX];
		tmpVector[pIDX].buff = 0;
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteTextureMPI()*/

