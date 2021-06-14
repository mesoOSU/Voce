#include <stdio.h>
#include <mpi.h>
#include "evp.h"

void PrintTensor(ten2nd m)
{
	int i,j;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
			printf("%e\t",m[i][j]);
		printf("\n");
	}
	return;
}/*end PrintTensor()*/

void PrintTensorNorm(ten2nd m)
{
	real norm = 0.;
	T2_loop{
		norm += m[mi][mj]*m[mi][mj];
	}
	norm = sqrt(norm);
	printf("pyz: Tensor norm = %e\n",norm);

}/*end PrintTensorNorm()*/

void PrintVoigt(voigt m)
{
	int i;
	for(i=0;i<6;i++){
		printf("%.3e\t",m[i]);
	}
	printf("\n");
	return;
}/*end PrintVoigt()*/

void PrintVoigt66(voigt66 m)
{
	int i,j;
	for(i=0;i<6;i++){
		for(j=0;j<6;j++)
			printf("%e\t",m[i][j]);
		printf("\n");
	}
	return;
}/*end PrintVoigt66()*/

void PrintB(void)
{
	int m;
	for(m=0;m<6;m++){
		printf("B[%d]:\n",m+1);
		PrintTensor(B[m]);
	}

	return;
}/*end PrintB()*/

void WriteEtaMPI(char *s)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;

	AllocMem(tmpVector, lsize, real);
	
	local_loop{
		tmpVector[pIDX] = (real)grain_f[pIDX];
	}

	sprintf(fname, "%s.iout", s);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, rank*lsize*sizeof(double), MPI_DOUBLE,
			MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteEtaMPI()*/

void PrintSchmidXt(void)
{
	int i,j;
	for(j=0;j<N_modes;j++){
		for(i=0;i<nsmx;i++){
			printf("pyz: Schmid vector of mode#%d:",i);
			printf("%.3f %.3f %.3f %.3f %.3f\n",
					Schm_xt[j][i][0], Schm_xt[j][i][1], Schm_xt[j][i][2], Schm_xt[j][i][3], Schm_xt[j][i][4]);
		}
	}
	return;
}/*end PrintSchmidXt()*/

static void GenerateLamellae(char *s, real vol, real t1[],real ph[],real t2[])
{
	FILE *fp;
	char fname[80];
	int i,j,k;
	int jph,jgr;
	int N=64;
	real half_d; 

	half_d = sqrt(2)/2.*N*(1.-sqrt(1-vol));
	sprintf(fname,"%s.ms",s);
	fp=fopen(fname,"w");
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				if((k>=N-sqrt(2)*half_d-j)&&(k<N+sqrt(2)*half_d-j)){
					jph = 2;
					jgr = 1;
					fprintf(fp,"%2f\t%2f\t%2f\t%d\t%d\t%d\t%d\t%d\n",
							t1[jph-1],ph[jph-1],t2[jph-1],i+1,j+1,k+1,jgr,jph);
				}
				else{
					jph = 1;
					jgr = 0;
					fprintf(fp,"%2f\t%2f\t%2f\t%d\t%d\t%d\t%d\t%d\n",
							t1[jph-1],ph[jph-1],t2[jph-1],i+1,j+1,k+1,jgr,jph);
				}
			}
		}
	}
	fclose(fp);

	return;
}/*end GenerateLamellae()*/

void Ti_alpha_beta(void)
{
	/* calculate the Euler angle for seven samples of
	   single alpha/beta colony */

	ten2nd e_xt;	// xtal directions in computational basis (sample axes)
					// for hcp, we use [2 -1 -1 0] (x-axis), [0 1 -1 0] (y-axis), and [0 0 0 1] (z-axis)
	ten2nd tran_m;	// sample -> xtal(alpha)
	ten2nd tran_BurgerOR;	// xtal(alph) -> xtal(beta)
	ten2nd tran_tot;
	real t1, ph, t2;
	real et1[2],eph[2],et2[2];
	int jph;
	ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;
	real theta;

	/* transformation matrix */
	theta = atan(1./sqrt(2));
	e_xt[0][0]=cos(theta)/sqrt(2); e_xt[0][1]=-1.*sin(theta)/sqrt(2); e_xt[0][2]=1./sqrt(2);
	e_xt[1][0]=-1.*cos(theta)/sqrt(2); e_xt[1][1]=sin(theta)/sqrt(2); e_xt[1][2]=1./sqrt(2);
	e_xt[2][0]=-1.*sin(theta); e_xt[2][1]=-1.*cos(theta); e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_BurgerOR[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}


	/* prismatic a1 */
	e_xt[0][0]=0.0; e_xt[0][1]=1./sqrt(2); e_xt[0][2]=-1./sqrt(2);
	e_xt[1][0]=0.0; e_xt[1][1]=1./sqrt(2); e_xt[1][2]=1./sqrt(2);
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Prismatic a1 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Prismatic a1 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Prism_1_64x64x64",0.12,et1,eph,et2);
	}

	/* prismatic a2 */
	theta = 75.*PI/180.;
	e_xt[0][0]=0.0; e_xt[0][1]=cos(theta); e_xt[0][2]=-1.*sin(theta);
	e_xt[1][0]=0.0; e_xt[1][1]=sin(theta); e_xt[1][2]=cos(theta);
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Prismatic a2 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Prismatic a2 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Prism_2_64x64x64",0.12,et1,eph,et2);
	}

	/* prismatic a3 */
	theta = 75.*PI/180.;
	e_xt[0][0]=0.0; e_xt[0][2]=-1.*cos(theta); e_xt[0][1]=sin(theta);
	e_xt[1][0]=0.0; e_xt[1][2]=sin(theta); e_xt[1][1]=cos(theta);
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Prismatic a3 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Prismatic a3 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Prism_3_64x64x64",0.12,et1,eph,et2);
	}

	/* basal a1 */
	e_xt[0][0]=0.0; e_xt[0][1]=1.0/sqrt(2); e_xt[0][2]=-1.0/sqrt(2);
	e_xt[1][0]=1.0; e_xt[1][1]=0.; e_xt[1][2]=0.;
	e_xt[2][0]=0.0; e_xt[2][1]=1.0/sqrt(2); e_xt[2][2]=1.0/sqrt(2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Basal a1 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Basal a1 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Basal_1_64x64x64",0.12,et1,eph,et2);
	}

	/* do a rotation about xtal c axis to obtain proper coord.
	   of basal a2 and a3 */
	/* basal a2 */
	e_xt[0][0]=0.866025; e_xt[0][1]=-0.353553; e_xt[0][2]=0.353553;
	e_xt[1][0]=0.5; e_xt[1][1]=0.612372; e_xt[1][2]=-0.612372;
	e_xt[2][0]=0.0; e_xt[2][1]=1.0/sqrt(2); e_xt[2][2]=1.0/sqrt(2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Basal a2 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Basal a2 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Basal_2_64x64x64",0.12,et1,eph,et2);
	}

	/* basal a3 */
	e_xt[0][0]=-0.866025; e_xt[0][1]=-0.353553; e_xt[0][2]=0.353553;
	e_xt[1][0]=0.5; e_xt[1][1]=-0.612372; e_xt[1][2]=0.612372;
	e_xt[2][0]=0.0; e_xt[2][1]=1.0/sqrt(2); e_xt[2][2]=1.0/sqrt(2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Basal a3 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Basal a3 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Basal_3_64x64x64",0.12,et1,eph,et2);
	}

	/* pyramidal c+a */
	e_xt[0][0]=0.; e_xt[0][1]=1.; e_xt[0][2]=0.;
	e_xt[1][0]=-1.; e_xt[1][1]=0.; e_xt[1][2]=0.;
	e_xt[2][0]=0.; e_xt[2][1]=0.; e_xt[2][2]=1.;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Pyramidal c+a of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(rank==0){
		printf("Pyramidal c+a of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(rank==0){
		GenerateLamellae("Pyramidal_64x64x64",0.12,et1,eph,et2);
	}

	return;
}/*end Ti_alpha_beta()*/
