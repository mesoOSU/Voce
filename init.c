#include <stdio.h>
#include <math.h>
#include <fftw_mpi.h>
#include <assert.h>
#include "evp.h"

#define CC(i,j,k,l) (C_ijkl[i-1][j-1][k-1][l-1])

static void ElastStiffnessMatrix(real *ElConst, ten4th cijkl)
{
	real c11, c12, c44;
	real C_ijkl[3][3][3][3] = {0.0};
	int i,j,k,l;

	c11 = *ElConst;
	c12 = *(ElConst+1);
	c44 = *(ElConst+2);

	CC(1,1,1,1) = c11;
	CC(2,2,2,2) = c11;
	CC(3,3,3,3) = c11;
	CC(1,1,2,2) = c12;
	CC(2,2,3,3) = c12;
	CC(3,3,1,1) = c12;
	CC(1,2,1,2) = c44;
	CC(2,3,2,3) = c44;
	CC(3,1,3,1) = c44;

	CC(2,2,1,1) = c12;
	CC(3,3,2,2) = c12;
	CC(1,1,3,3) = c12;
	CC(2,1,2,1) = c44;
	CC(1,2,2,1) = c44;
	CC(2,1,1,2) = c44;
	CC(3,2,3,2) = c44;
	CC(3,2,2,3) = c44;
	CC(2,3,3,2) = c44;
	CC(1,3,1,3) = c44;
	CC(1,3,3,1) = c44;
	CC(3,1,1,3) = c44;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
					cijkl[i][j][k][l] = CC(i+1,j+1,k+1,l+1);

	return;
}/*end ElastStiffnessMatrix()*/

static void PlasticInit(char *s, int phase)
{
	FILE *fp;
	char buffer[80] = {0};
	int i, j, mx, m_count;
	int matcheck;
	voigt5 aux5;
	ten2nd aux33;
	voigt66 aux66;
	ten4th aux3333;
	
    


	fp = fopen(s, "r");

	fgets(buffer, 80, fp);

	fgets(buffer, 80, fp);
	sscanf(buffer,"%s",XtalSys);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);

	fgets(buffer, 80, fp);
	sscanf(buffer,"%lf %lf %lf",&XtalAxis[0], &XtalAxis[1], &XtalAxis[2]);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);

	fgets(buffer, 80, fp);
	sscanf(buffer,"%d", &N_modes_max);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);

	fgets(buffer, 80, fp);
	sscanf(buffer,"%d", &N_modes);
	AllocMem(iMode,N_modes,int);
//	AllocMem(bn, N_modes, real**);
//	AllocMem(bb, N_modes, real**);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);
	for(i=0;i<N_modes;i++){
		fgets(buffer, 80, fp);
		sscanf(buffer,"%d",&iMode[i]);
	}


	fgets(buffer, 80, fp);

	m_count = 0;
	for(i=0;i<N_modes_max;i++){
		fgets(buffer, 80, fp);
		sscanf(buffer, "%*s %*s %d %*s",&mx);
		matcheck = 0;
		j = 0;
		while(j<N_modes){
			if(mx==iMode[j]){
				matcheck = 1;
				m_count++;
				break;
			}
			j++;
		}
		if(matcheck==1){

			fgets(buffer,80,fp);
			fgets(buffer,80,fp);
			sscanf(buffer,"%d\t\t%lf\t\t%lf\t\t%lf\t\t%d", &nsmx, &nrsx, &gamd0x, &twshx, &isectwx);
			assert(nsmx<=NSYSMX);
			nSYS[phase] = nsmx;
	//		AllocMem(bn[m_count-1],nsmx,real*);
	//		AllocMem(bb[m_count-1],nsmx,real*);
//			for(j=0;j<nsmx;j++){
//				AllocMem(bn[m_count-1][j],3,real);
//				AllocMem(bb[m_count-1][j],3,real);
//			}

			fgets(buffer,80,fp);
			fgets(buffer,80,fp);
			sscanf(buffer,"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf", &tau0xf, &tau0xb, &tau1x, &thet0, &thet1);
			fgets(buffer,80,fp);
			fgets(buffer,80,fp);
			sscanf(buffer,"%lf\t\t%lf", &hselfx, &hlatex2);
			fgets(buffer,80,fp);
		/*	if(thet0<thet1){
				PError("Initial hardening rate is lower than final hardening rate !!",230);
			} */
			if(tau1x<1.E-6){	// linear hardening, independent of tau0
				tau1x =  1.E-6;	// avoid division by zero
				thet0 = thet1;
			}
			for(j=0;j<nsmx;j++){
				// define strain rate sensitivity and crss
				nSRS[phase][j] = nrsx;
				gam0[phase][j] = gamd0x;
// 				tau[phase][j][0] = tau0xf;
// 				tau[phase][j][1] = tau0xb;
// 				tau[phase][j][2] = tau1x;
// 				thet[phase][j][0] = thet0;
				//thet[phase][j][1] = thet1;
				if(j<nsmx/2) {
            	tau[phase][j][0] = tau0xf;
				tau[phase][j][1] = tau0xb;
 				tau[phase][j][2] = tau1x;
 				thet[phase][j][0] = thet0;
				thet[phase][j][1] = thet1;
                                    hlatex[j] = 1.1;
                                } 
                                if(j>=nsmx/2) {
                tau[phase][j][0] = tau0xf;
				tau[phase][j][1] = tau0xb;
 				tau[phase][j][2] = tau1x;
 				thet[phase][j][0] = thet0;
				thet[phase][j][1] = thet1;
                                    hlatex[j] = 1.1;
                                }

				fgets(buffer,80,fp);
//				sscanf(buffer,"%lf %lf %lf\t\t%lf %lf %lf",
//						&bn[m_count-1][j][0], &bn[m_count-1][j][1], &bn[m_count-1][j][2],
//						&bb[m_count-1][j][0], &bb[m_count-1][j][1], &bb[m_count-1][j][2]);
				sscanf(buffer,"%lf %lf %lf\t\t%lf %lf %lf",
						&bn[phase][j][0], &bn[phase][j][1], &bn[phase][j][2],
						&bb[phase][j][0], &bb[phase][j][1], &bb[phase][j][2]);
			}
		//	fgets(buffer,80,fp);
		}
		else{
			for(j=0;j<24;j++)
				fgets(buffer,80,fp);
		}
	}

	fclose(fp);

	/* normalize bn and bb */
	real norm;
	for(i=0;i<N_modes;i++){
		if(rank==0){
                    
			printf("Slip mode#%d:\n", i);
		}
		for(j=0;j<nsmx;j++){
			norm = sqrt(pow(bn[i][j][0],2.0) + pow(bn[i][j][1],2.0) + pow(bn[i][j][2],2.0));
			for(mx=0;mx<3;mx++){
				bn[i][j][mx] /= norm;
			}
			
			if(rank==0)printf("bn[%d][%d] = (%lf, %lf, %lf)\n",i,j,bn[i][j][0],bn[i][j][1],bn[i][j][2]);
			norm = sqrt(pow(bb[i][j][0],2.0) + pow(bb[i][j][1],2.0) + pow(bb[i][j][2],2.0));
			for(mx=0;mx<3;mx++){
				bb[i][j][mx] /= norm;
			}
			if(rank==0)printf("bb[%d][%d] = (%lf, %lf, %lf)\n",i,j,bb[i][j][0],bb[i][j][1],bb[i][j][2]);

			norm = 0.0;	// used to check orthogonality
			for(mx=0;mx<3;mx++){
				norm += bn[i][j][mx]*bb[i][j][mx];
			}
			if(norm>1E-3){
				PError("Slip systems is NOT orthogonal !!",7);
			}
                        //if(rank==0)printf("%lf\n",hlatex[j]);
			// Schmid tensors
			T2_loop{
				aux33[mi][mj] = (bn[i][j][mi]*bb[i][j][mj]+bn[i][j][mj]*bb[i][j][mi])/2.0;
			}
			chg_basis5(aux5,aux33,aux66,aux3333,2);
			for(mx=0;mx<5;mx++){
				Schm_xt[phase][j][mx] = aux5[mx];
			}
		}
	}

	// initialize self and latent hardening coeff. for each system of the phase
	// notice here we assue latent hardeing is the same for all the off diangonal
	// elements, which is of course not always the case.
	
	//cp
	//levi-Civita
	/*for (i=0;i<3;i++){
	    for (j=0;j<3;j++){
	        for(k=0;k<3;k++){
	            x=j-i;y=k-j;z=i-k;
	          s1=(x > 0) ? 1 : ((x < 0) ? -1 : 0);
	           s2=(y > 0) ? 1 : ((y < 0) ? -1 : 0);
	            s3=(z > 0) ? 1 : ((z < 0) ? -1 : 0);
	            l_c[i][j][k]=s1+s2+s3;
	        } 
	    }
	}*/
	//cp
/*	for(i=0;i<nSYS[phase];i++){
		for(j=0;j<nSYS[phase];j++){
                /*    if(i==j){
				Hard[phase][i][j] = hselfx;
			} 
			else if((bn[phase][i][0] == bn[phase][j][0]) && (bn[phase][i][1] == bn[phase][j][1]) && (bn[phase][i][2] == bn[phase][j][2])){
                            //if((i<nsmx/2) && (j>=nsmx/2)) {
                        Hard[phase][i][j] = hselfx;
                        //}
                        }
                        
                   else if((bb[phase][i][0] == bb[phase][j][0]) && (bb[phase][i][1] == bb[phase][j][1]) && (bb[phase][i][2] == bb[phase][j][2])){
                            //if((i<nsmx/2) && (j>=nsmx/2)) {
                        Hard[phase][i][j] = 1.8;
                        //}
                        }
                       else if((bb[phase][i][0]*bb[phase][j][0]) + (bb[phase][i][1]*bb[phase][j][1]) + (bb[phase][i][2]*bb[phase][j][2]) == 0.0){
                            //if((i<nsmx/2) && (j>=nsmx/2)) {
                        Hard[phase][i][j] = 0.5;
                        //}
                        }
                        
                        
                        //else{
                           // Hard[phase][i][j] = hlatex[j];
                        //}
                        
                        else{  
                        Hard[phase][i][j] = hlatex[j];
                          
                        }
                    
                        if(rank==0)printf("%lf\n",Hard[phase][i][j]);  
                        
        }
        }*/
        
        	Hard[phase][0][0] = self[phase];
			Hard[phase][0][1] = coplane[phase];
			Hard[phase][0][2] = coplane[phase];
		Hard[phase][0][3] = hirth[phase];
			Hard[phase][0][4] = glissile[phase];
			Hard[phase][0][5] = lomer[phase];
			Hard[phase][0][6] = cross[phase];
			Hard[phase][0][7] = glissile[phase];
			Hard[phase][0][8] = glissile[phase];
			Hard[phase][0][9] = hirth[phase];
			Hard[phase][0][10] = lomer[phase];
			Hard[phase][0][11] = glissile[phase];

			Hard[phase][1][1] = self[phase];
			Hard[phase][1][2] = coplane[phase];
			Hard[phase][1][3] = glissile[phase];
			Hard[phase][1][4] = cross[phase];
			Hard[phase][1][5] = glissile[phase];
			Hard[phase][1][6] = glissile[phase];
			Hard[phase][1][7] = hirth[phase];
			Hard[phase][1][8] = lomer[phase];
			Hard[phase][1][9] = lomer[phase];
			Hard[phase][1][10] = hirth[phase];
			Hard[phase][1][11] = glissile[phase];

			Hard[phase][2][2] = self[phase];
			Hard[phase][2][3] = lomer[phase];
			Hard[phase][2][4] = glissile[phase];
			Hard[phase][2][5] = hirth[phase];
			Hard[phase][2][6] = glissile[phase];
			Hard[phase][2][7] = lomer[phase];
			Hard[phase][2][8] = hirth[phase];
			Hard[phase][2][9] = glissile[phase];
			Hard[phase][2][10] = glissile[phase];
			Hard[phase][2][11] = cross[phase];

			Hard[phase][3][3] = self[phase];
			Hard[phase][3][4] = coplane[phase];
			Hard[phase][3][5] = coplane[phase];
			Hard[phase][3][6] = hirth[phase];
			Hard[phase][3][7] = lomer[phase];
			Hard[phase][3][8] = glissile[phase];
			Hard[phase][3][9] = cross[phase];
			Hard[phase][3][10] = glissile[phase];
			Hard[phase][3][11] = glissile[phase];

			Hard[phase][4][4] = self[phase];
			Hard[phase][4][5] = coplane[phase];
			Hard[phase][4][6] = lomer[phase];
			Hard[phase][4][7] = hirth[phase];
			Hard[phase][4][8] = glissile[phase];
			Hard[phase][4][9] = glissile[phase];
			Hard[phase][4][10] = hirth[phase];
			Hard[phase][4][11] = lomer[phase];

			Hard[phase][5][5] = self[phase];
			Hard[phase][5][6] = glissile[phase];
			Hard[phase][5][7] = glissile[phase];
			Hard[phase][5][8] = cross[phase];
			Hard[phase][5][9] = glissile[phase];
			Hard[phase][5][10] = lomer[phase];
			Hard[phase][5][11] = hirth[phase];

			Hard[phase][6][6] = self[phase];
			Hard[phase][6][7] = coplane[phase];
			Hard[phase][6][8] = coplane[phase];
			Hard[phase][6][9] = hirth[phase];
			Hard[phase][6][10] = glissile[phase];
			Hard[phase][6][11] = lomer[phase];

			Hard[phase][7][7] = self[phase];
			Hard[phase][7][8] = coplane[phase];
			Hard[phase][7][9] = glissile[phase];
			Hard[phase][7][10] = cross[phase];
			Hard[phase][7][11] = glissile[phase];

			Hard[phase][8][8] = self[phase];
			Hard[phase][8][9] = lomer[phase];
			Hard[phase][8][10] = glissile[phase];
			Hard[phase][8][11] = hirth[phase];

			Hard[phase][9][9] = self[phase];
			Hard[phase][9][10] = coplane[phase];
			Hard[phase][9][11] = coplane[phase];

			Hard[phase][10][10] = self[phase];
			Hard[phase][10][11] = coplane[phase];

			Hard[phase][11][11] = self[phase];

		/* symmetric matrix */
		for(i=1;i<12;i++){
			for(j=0;j<i;j++){
					Hard[phase][i][j] = 	Hard[phase][j][i];
                        
                       
                        
        }
        }

// 	for(i=0;i<nSYS[phase];i++){
// 		for(j=0;j<nSYS[phase];j++){
//                            if(i!=j){
// 				
// 			   }
// 		}
// 	}
	return;
}/*end PlasticInit()*/

static void InitMicroStruct(char *s)
{
	ten4th caux3333;
	voigt aux6;
	ten2nd aux33;
	voigt66 caux66;
	voigt66 c066_local = {0.0};
	real ph, th, om;
	int jgr, jph;
	int ii,jj,kk;
	int nph1, nph1_all;
	FILE *fp;
	char buffer[80] = {0};
	int i, idx;
	int EmptySteps = rank*Nxyz/NumPE;
	ten2nd sa2xt;	// transform matrix (sample -> xtal)

	fp = fopen(s,"r");
	/*empty reading to go to the corresponding slabbed region*/
	for(i=0;i<EmptySteps;i++){
		fgets(buffer,80,fp);
	}

	nph1 = 0;
	local_loop{
		fgets(buffer,80,fp);
		sscanf(buffer,"%lf %lf %lf %d %d %d %d %d", &ph, &th, &om, &ii, &jj, &kk, &jgr, &jph);
		idx = ((ii-1-lxs)*CellDim[1]+jj-1)*CellDim[2]+kk-1;
		if(idx!=pIDX){
			PError("Error in reading microstructure data file (inconsistent index system)!!", 1139);
		}

		if(jph==1) nph1++;
		grain_f[pIDX] = jgr;
		phase_f[pIDX] = jph;

		if(!Type_phases[jph-1]){	// NOT gas!!
			ph *= PI/180.;	th *= PI/180.;	om *= PI/180.;
			EulerToTransMatrix(&ph, &th, &om, sa2xt, 2);
			T2_loop{
				TranMat_xt2sa[pIDX][mi][mj] = sa2xt[mj][mi];
			}
			Ten4thTransform(caux3333,sa2xt,Cijkl[jph-1],2);	// Cijkl[jph] is in xtal ref. Do inverse transform
			chg_basis(aux6, aux33, caux66,caux3333,4);
			C6_loop{
				C_gr[pIDX][mi][mj] = caux66[mi][mj];
				c066_local[mi][mj] += caux66[mi][mj];
			}
		}
	}
	fclose(fp);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&nph1, &nph1_all, 1, MPI_INT,
			MPI_SUM, MPI_COMM_WORLD);
	Wgt_ph1 = (1.0*(real)nph1_all)/Nxyz;

	C6_loop{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&c066_local[mi][mj], &C066[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		C066[mi][mj] *= WGT;
	//	printf("C066[%d][%d]=%e\n",mi,mj,C066[mi][mj]);
	}

	C6_loop{
		S066[mi][mj] = C066[mi][mj];
	}

	LU_inv_66(S066);

	chg_basis(aux6,aux33,C066,C0,3);
	chg_basis(aux6,aux33,S066,S0,3);

//	/* test LU_inv_66 */
//	voigt66 s0 = {0.0};
//	s0[0][1] = 1.0;
//	s0[1][0] = 2.0; s0[1][2] = 2.0;
//	s0[2][1] = 3.0; s0[2][3] = 1.0;
//	s0[3][2] = 1.0; s0[3][4] = 2.0;
//	s0[4][3] = 3.0; s0[4][5] = 1.0;
//	s0[5][4] = 2.0;
//	printf("pyz: Before inverse:\n");
//	PrintVoigt66(s0);
//	LU_inv_66(s0);
//	printf("pyz: After inverse:\n");
//	PrintVoigt66(s0);


	return;
}/*end InitMicroStruct()*/


static void InitFields(void)
{
	int jph;
	voigt aux6;
	ten2nd aux33;
	voigt66 cgaux66;
	ten4th cg;

	// Macro strain
	T2_loop{
		DisGradAvg[mi][mj] = Udot[mi][mj]*TimeStep;
		dDisGradAvg[mi][mj] = 0.0;
		dDisGradAvg_acum[mi][mj] = 0.0;
	}

	local_loop{
		//jgr = grain_f[pIDX];	// there could be thermo strain associated with grains. Here we ignore it for the time being.
		// strain, strain rate, disp. grad.
		T2_loop{
			Edot[pIDX][mi][mj] = 0.0;
			Eps[pIDX][mi][mj] = 0.0;
			DisGrad[pIDX][mi][mj] = 0.0;
		}

		// inital guess for the stress, assume elasticity
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			C6_loop{
				cgaux66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			chg_basis(aux6,aux33,cgaux66,cg,3);
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
				fluctuation[pIDX][mi][mj] = 0.0;
				T2p_loop{
					Sig[pIDX][mi][mj] += cg[mi][mj][mip][mjp]*DisGradAvg[mip][mjp];
				}
				if(CREEP_FLAG==1){
					Sig[pIDX][mi][mj] += Scauchy[mi][mj];
				}
			}
			//PrintTensorNorm(Sig[pIDX]);
		}
		else{	// gas
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
			}
		}
	}

	e_vm = d_vm*TimeStep;

	if(rank==0){
		fprintf(fp_vm,"EVM\tEVMP\tDVM\tDVMP\tSVM\tSVMSVM_1\tIteration_tot\n");
		fprintf(fp_err,"IT\tERRE\tERRS\tSVM\n");
	}

	return;
}/*end InitFields()*/

static void M3Inverse(real a[3][3], real ia[3][3])
{
	int i, j;
	real deta = 0.0;

	deta = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) \
		   -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) \
		   +a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
	ia[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
	ia[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
	ia[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
	ia[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
	ia[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
	ia[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
	ia[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
	ia[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
	ia[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			ia[i][j] /=deta;

	return;
}/*end M3Inverse()*/

static void CalcGAMMA(void)
{
	int i,j,k,l;
	real nvector[3];
	real iomega[3][3], omega[3][3];

	local_loop{
		C4_loop{
			GAMMA[pIDX][mi][mj][mk][ml] = 0.0;
		}
	}

	local_loop{
		nvector[0] = g[pIDX].x;
		nvector[1] = g[pIDX].y;
		nvector[2] = g[pIDX].z;
		if(sIDX==0){		//k=0, boundary condition
			T2_loop{
				T2p_loop{
					GAMMA[pIDX][mi][mj][mip][mjp] = 0.0;
				}
			}
		}
		else{
			// Calculate Omega
			for(i=0;i<3;i++){
				for(k=0;k<3;k++){
					iomega[i][k] = 0.0;
					for(j=0;j<3;j++)
						for(l=0;l<3;l++)
							iomega[i][k] += C0[i][j][k][l]*nvector[j]*nvector[l];
				}
			}

			M3Inverse(iomega, omega);

			T2_loop{
				T2p_loop{
					GAMMA[pIDX][mi][mj][mip][mjp] = -1.0*omega[mi][mip]*nvector[mj]*nvector[mjp];
				}
			}
		}
	}

	return;
}/*end CalcGAMMA()*/

void SetupJob(void)
{
	int i;
	voigt aux6;
	ten2nd aux33;
	voigt66 aux66;
	ten4th aux3333;
	ten2nd dev_Udot_s;

	if(rank==0){
		printf("\n========================================\n");
		printf("------Start setting up the job......------\n");
	}

	// total number of grid points
	Nxyz = CellDim[0]*CellDim[1]*CellDim[2];
	WGT = 1./Nxyz;

	assert(N_phases<=NPHMAX);

	/************************
	  fft_mpi initialization
	 ************************/
	plan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, CellDim[0], CellDim[1], CellDim[2],
				FFTW_FORWARD, FFTW_ESTIMATE);
	iplan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, CellDim[0], CellDim[1], CellDim[2],
				FFTW_BACKWARD, FFTW_ESTIMATE);
	/* slab decomposition
	   all PEs work coorporatively
	   to perform FFT*/
	fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);
	fft_data = (fftw_complex*)fftw_malloc(lsize*sizeof(fftw_complex));
	fft_work = (fftw_complex*)fftw_malloc(lsize*sizeof(fftw_complex));
	if(rank==0) printf("\nfft_mpi initialized\n");

	/************************
	 k-space setup
	 ************************/
	AllocMem(g, lsize, Vec3R);
	AllocMem(g2, lsize, real);
	AllocMem(GAMMA, lsize, ten4th);
	real tmp[3];
	int m;
	local_loop{
		tmp[0] = (real)(px+lxs) * (real)(1.0/CellDim[0]);
		tmp[1] = (real)(py) * (real)(1.0/CellDim[1]);
		tmp[2] = (real)(pz) * (real)(1.0/CellDim[2]);
		for(m=0; m<3; m++){
			tmp[m] -= round(tmp[m]);
			tmp[m] *= TWO_PI;
		}

		VSet(g[pIDX], tmp[0], tmp[1], tmp[2]);
		g2[pIDX] = VSecNorm(g[pIDX]);

		if((px+lxs+py+pz) == 0){
			g2[pIDX] = 1.0;
			VZero(g[pIDX]);
		}
		VScale(g[pIDX], 1./g2[pIDX]);
	}
	if(rank==0) printf("\nk-space initialized\n");

	/*******************************
	  Material properties initialization
	  ******************************/
	chg_basis(aux6, aux33, aux66, aux3333, 0);	// calculate B[6]
	if(rank==0){
		PrintB();
	}
	AllocMem(nSRS,N_phases,real*);
	AllocMem(nSYS,N_phases,int);
	for(i=0;i<N_phases;i++){
		nSYS[i] = 0;
		AllocMem(nSRS[i],NSYSMX,real);
	}
	AllocMem(C_gr,lsize,voigt66);
	AllocMem(TranMat_xt2sa,lsize,ten2nd);
	AllocMem(Schm_gr,lsize,voigtch);
	if(!Type_phases[0]){
		ElastStiffnessMatrix(ElastConst_PhaseI, Cijkl[0]);
		PlasticInit(Slip_PhaseI, 0);
	}
	if(N_phases==2&&(!Type_phases[1])){
		ElastStiffnessMatrix(ElastConst_PhaseII, Cijkl[1]);
		PlasticInit(Slip_PhaseII, 1);
	}
	AllocMem(grain_f,lsize,int);
	AllocMem(phase_f,lsize,int);
	InitMicroStruct(initial_ms);
//	PrintSchmidXt();

	/************************
	  Boundary conditions
	 ************************/
	// strain rate
	if((VelGrad_BC_Flag[0]+VelGrad_BC_Flag[1]+VelGrad_BC_Flag[2])==2){
		PError("Cannot enforce only two deviatoric components (Check input VelGrad_BC_Flag)",116);
	}
	VoigtToFull(VelGrad_BC, Udot);
	if(rank==0){
		printf("Applied velocity gradient:\n");
		PrintTensor(Udot);
	}
	SymAntDecompose(Udot, Udot_s, Udot_a);
	chg_basis(D_bar6,Udot_s,aux66,aux3333,2);
	for(i=0;i<5;i++) D_bar5[i] = D_bar6[i];
	chg_basis5(D_bar5,dev_Udot_s,aux66,aux3333,1);
	d_vm = 0.0;
	T2_loop{
		d_vm += pow(dev_Udot_s[mi][mj],2.0);
	}
	d_vm = sqrt(2./3.*d_vm);

	// Cauchy stress
	for(i=0;i<6;i++){	// check BC validity
		if((Stress_BC_Flag[i]*VelGrad_BC_Flag[i]!=0)||((Stress_BC_Flag[i]+VelGrad_BC_Flag[i])!=1))
			PError("Check boundary conditions on strain rate and stress!!",114);
	}
  CREEP_FLAG = 0;
	for(i=0;i<6;i++){	// creep?
		CREEP_FLAG += Stress_BC_Flag[i];
	}
	CREEP_FLAG /= 6;
	VoigtToFull(Stress_BC, Scauchy);
	if(rank==0){
		printf("Applied (Cauchy) stress:\n");
		PrintTensor(Scauchy);
	}

	// Initialize CRSS and accum shear for grains
	AllocMem(gamacum,lsize,real);
	AllocMem(SVM,lsize,real);
	AllocMem(EVM,lsize,real);
	AllocMem(rho_rss,lsize,real);
	AllocMem(rho_crss,lsize,real);
	AllocMem(work_p,lsize,real);
        AllocMem(gammatotal,lsize,real*);
         AllocMem(tau_x,lsize,real*);
        //AllocMem(cum_gama,lsize,real*);
	AllocMem(gamdot,lsize,real*);
	AllocMem(displacement_fluct,lsize,real*);
	AllocMem(xkin,lsize,real*);
	AllocMem(trial_tau,lsize,real**);
	AllocMem(crss,lsize,real**);
	local_loop{
                AllocMem(gammatotal[pIDX],NSYSMX,real);
                AllocMem(tau_x[pIDX],NSYSMX,real);
		AllocMem(gamdot[pIDX],NSYSMX,real);
		AllocMem(displacement_fluct[pIDX],3,real);
                //AllocMem(cum_gama[pIDX],NSYSMX,real);
		AllocMem(xkin[pIDX],NSYSMX,real);
		AllocMem(trial_tau[pIDX],NSYSMX,real*);
		AllocMem(crss[pIDX],NSYSMX,real*);
		for(i=0;i<NSYSMX;i++){
			AllocMem(trial_tau[pIDX][i],2,real);
			AllocMem(crss[pIDX][i],2,real);
		}
	}
	int jph;
	local_loop{
work_p[pIDX]=0.0;
SVM[pIDX]=0.0;
EVM[pIDX]=0.0;
rho_rss[pIDX]=0.0;
rho_crss[pIDX]=0.0;
//             for(i=0;i<nSYS[jph-1];i++){
// 		gamacum[pIDX][i] = 0.0;
//             }
		jph = phase_f[pIDX];
		if(Type_phases[jph-1]==0){
			for(i=0;i<nSYS[jph-1];i++){
				crss[pIDX][i][0] = tau[jph-1][i][0];
				crss[pIDX][i][1] = tau[jph-1][i][1];
				trial_tau[pIDX][i][0] = tau[jph-1][i][0];
				trial_tau[pIDX][i][1] = tau[jph-1][i][1];
				xkin[pIDX][i] = 0.0;
			}
		}
	}



	/*******************************
	  global arrays for stress/strain...
	  arrays allocation and initialization
	  ******************************/
	AllocMem(Sig, lsize, ten2nd);
	AllocMem(fluctuation,lsize,ten2nd);
	AllocMem(Eps, lsize, ten2nd);
	AllocMem(Edot, lsize, ten2nd);
	AllocMem(DisGrad, lsize, ten2nd);
	AllocMem(VelGrad, lsize, ten2nd);



	/*******************************
	  Initialize the stress/strain,rates...
	  ******************************/
	InitFields();

	/*******************************
	  Green's operator in k-space
	  ******************************/
	CalcGAMMA();	// calculate Green operator in k-space

	if(rank==0){
		printf("\n-----------Job set up-----------------\n");
		printf("========================================\n");
	}
	return;
}/*end SetupJob()*/

void DestroyJob(void)
{
	int i;

	if(rank==0){
		printf("\n========================================\n");
		printf("--------Start destroying job......----------\n");
	}

	/*******************************
	  Material properties initialization
	  ******************************/
	for(i=0;i<N_phases;i++)
		free(nSRS[i]);
	free(nSRS);
	free(nSYS);
	free(C_gr);
	free(TranMat_xt2sa);
	free(grain_f);
	free(phase_f);
	free(Schm_gr);

	free(iMode);

//	for(i=0;i<N_modes;i++){
//		for(j=0;j<nsmx;j++){
//			free(bn[i][j]);
//			free(bb[i][j]);
//		}
//		free(bn[i]);
//		free(bb[i]);
//	}
//	free(bb);
//	free(bn);

	/************************
	  crss, hardening,...
	 ************************/
	free(gamacum);
	free(work_p);
	free(SVM);
	free(EVM);
	free(rho_crss);
	free(rho_rss);
        free(gammatotal);
        free(tau_x);
        //free(cum_gama);
	local_loop{
		for(i=0;i<NSYSMX;i++){
			free(trial_tau[pIDX][i]);
			free(crss[pIDX][i]);
		}
		free(gamdot[pIDX]);
                free(gammatotal[pIDX]);
                free(displacement_fluct[pIDX]);
                free(tau_x);
                //free(cum_gama[pIDX]);
		free(xkin[pIDX]);
		free(trial_tau[pIDX]);
		free(crss[pIDX]);
	}
	free(gamdot);
	free(displacement_fluct);
	free(xkin);
	free(trial_tau);
	free(crss);

	/************************
	  fft_mpi finalization
	 ************************/
	fftwnd_mpi_destroy_plan(plan);
	fftwnd_mpi_destroy_plan(iplan);
	fftw_free(fft_data);
	fftw_free(fft_work);
	if(rank==0) printf("\nfft_mpi finalized\n");

	/************************
	 free k-space
	 ************************/
	free(g);
	free(g2);
	free(GAMMA);
	if(rank==0) printf("\nk-space finalized\n");


	if(rank==0){
		printf("\n------------Job destroyed-----------\n");
		printf("========================================\n");
	}

	/************************
	 free global arrays (stress/strain)
	 ************************/
	free(Sig);
	free(fluctuation);
	free(Eps);
	free(Edot);
	free(DisGrad);
	free(VelGrad);
	return;
}/*end DestroyJob()*/

