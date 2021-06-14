#ifndef EVP_H
#define EVP_H
#include <stdlib.h>
#include <mpi.h>
#include <fftw_mpi.h>
#include "evpFlags.h"
#include "V3math.h"

/**********************
  types and structures
  *********************/
typedef double real;
#define MPI_real MPI_DOUBLE
//typedef float real;
//#define MPI_real MPI_FLOAT
typedef struct{
	real x;
	real y;
	real z;
}Vec3R;
typedef real ten2nd[3][3];
typedef real ten4th[3][3][3][3];
typedef real voigt[6];
typedef real voigt5[5];
typedef real voigtch[NSYSMX][5];
typedef real voigt66[6][6];


/***************************
  loops, allocators, indices
 ***************************/
#define sIDX (((px+lxs)*CellDim[1] + py)*CellDim[2] + pz)	/* index for space array */
#define pIDX ((px*CellDim[1] + py)*CellDim[2] + pz)	/* index for local array */
#define AllocMem(a, n, t) a = (t*)malloc((n)*sizeof(t))
#define B(i,j,p) (B[p-1][i-1][j-1])
#define Np_loop for(ip=0;ip<N_phases;ip++)
#define local_loop for(px=0;px<lnx;px++)for(py=0;py<CellDim[1];py++)for(pz=0;pz<CellDim[2];pz++)
#define T2_loop for(mi=0;mi<3;mi++)for(mj=0;mj<3;mj++)
#define T2p_loop for(mip=0;mip<3;mip++)for(mjp=0;mjp<3;mjp++)
#define C6_loop for(mi=0;mi<6;mi++)for(mj=0;mj<6;mj++)
#define C4_loop for(mi=0;mi<3;mi++)for(mj=0;mj<3;mj++)for(mk=0;mk<3;mk++)for(ml=0;ml<3;ml++)
#define PError(s,i) perror(#s),\
	MPI_Abort(MPI_COMM_WORLD, i)
extern int ip;		// dummy indices for global usage
extern int px,py,pz;
extern int mi, mj, mk, ml, mip, mjp;	// loop over the 3x3/3x3x3x3 tensors (stress/strain, etc)
#ifdef EVP_GLOBALS_ONCE
int ip;
int px,py,pz;
int mi, mj, mk, ml, mip, mjp;
#endif


/***************************
  units & constants
 ***************************/
#define PI 3.14159265359
#define TWO_PI 6.28318530718


/***************************
  MPI parameters
 ***************************/
extern int rank, NumPE;
extern int lnx, lxs, lnyt, lyst, lsize;
extern MPI_Status status;
#ifdef EVP_GLOBALS_ONCE
int rank = 0, NumPE = 0;
int lnx = 0, lxs = 0, lnyt = 0, lyst = 0, lsize = 0;
MPI_Status status;
#endif


/***************************
  fft, k-space
 ***************************/
extern fftwnd_mpi_plan plan, iplan;
extern fftw_complex *fft_data, *fft_work;
extern Vec3R *g;
extern real *g2;
extern ten4th *GAMMA;	// Green operator in Fourier space, a function of Ref. stiffness and frequency
#ifdef EVP_GLOBALS_ONCE
fftwnd_mpi_plan plan, iplan;
fftw_complex *fft_data = 0, *fft_work = 0;
Vec3R *g = 0;
real *g2 = 0;
ten4th *GAMMA = 0;
#endif


/***************************
 Computational grid, ICs, BCs
 Simulation controllers
 ***************************/
extern int CellDim[3];
extern int Nxyz;
extern real TimeStep;
extern int N_steps;
extern real Err;
extern int IterMax;
extern int PrintControl[2];
extern int Update_Flag;
extern int Hard_Flag;
extern int Tex_Flag;
extern int VelGrad_BC_Flag[6];// pyz: we currently assume symmetric vel. gradient applied
extern real VelGrad_BC[6];
extern ten2nd Udot;		// full matrix form of VelGrad_BC
extern ten2nd Udot_s, Udot_a;	// sym and ant part of Udot
extern voigt D_bar6;	// voigt notation of Udot_s
extern voigt5 D_bar5;	// devi. part of D_bar6
extern ten2nd DisGradAvg;	// macro displ. gradent
extern ten2nd dDisGradAvg;
extern ten2nd dDisGradAvg_acum;
extern int Stress_BC_Flag[6];
extern int CREEP_FLAG;  // For fully constant stress-controlled boundary
extern real Stress_BC[6];
extern ten2nd Scauchy;	// full matrix form of Stress_BC
#ifdef EVP_GLOBALS_ONCE
int CellDim[3] = {0};
int Nxyz = 0;
real TimeStep = 0.0;
int N_steps = 0;
real Err = 0.0;
int IterMax = 0;
int PrintControl[2] = {0};
int Update_Flag = 0;
int Hard_Flag = 0;
int Tex_Flag = 0;
int VelGrad_BC_Flag[6] = {0};
real VelGrad_BC[6] = {0.};
ten2nd Udot = {0.0};
ten2nd Udot_s = {0.};
ten2nd Udot_a = {0.};
voigt D_bar6 = {0.};
voigt5 D_bar5 = {0.};
ten2nd DisGradAvg = {0.};
ten2nd dDisGradAvg = {0.};
ten2nd dDisGradAvg_acum = {0.0};
int Stress_BC_Flag[6] = {0};
int CREEP_FLAG=0;  // For fully constant stress-controlled boundary
real Stress_BC[6] = {0.};
ten2nd Scauchy = {0.};
#endif


/***************************
	measures and errors
 ***************************/
extern real d_vm;	//	von Mises (vm) strain rate
extern real s_vm;	//  von Mises stress
extern real e_vm;	//  von Mises strain
extern real s_vm1;	// von Mises stress in phase-1
extern real Err_e;	// error in strain
extern real Err_s;	// error in stress
#ifdef EVP_GLOBALS_ONCE
real d_vm = 0.;
real s_vm = 0.;
real e_vm = 0.;
real s_vm1 = 0.;
real Err_e = 0.;
real Err_s = 0.;
#endif


/***************************
  global files
 ***************************/
extern FILE *fp_vm;
extern FILE *fp_err;
#ifdef EVP_GLOBALS_ONCE
FILE *fp_vm = 0;
FILE *fp_err = 0;
#endif

/**********************
  Material types and phases
  *********************/
/*
   N_phases -- # of phases under consideration
   Slip_Phase* -- Path of the input file containing plasticity (slip) of phase *
   */
extern int N_phases;	// <= 2 for current version
extern int Type_phases[2];
extern real ElastConst_PhaseI[4];
extern real ElastConst_PhaseII[4];
extern char Slip_PhaseI[100];	// currently assume only 2 phases at most
extern char Slip_PhaseII[100];
extern char initial_ms[100]; // file name of inital microstructure
extern ten4th Cijkl[NPHMAX];
extern real **nSRS;		// strain rate sensitvity of each slip/twin system for each phase
extern int *nSYS;		// number of slip+twin systems for each phase
extern int N_modes_max, N_modes;
extern char XtalSys[10];
extern real XtalAxis[3];
extern int *iMode;
//extern real ***bn, ***bb;
extern real bn[NPHMAX][NSYSMX][3];
extern real bb[NPHMAX][NSYSMX][3];
extern voigt5 Schm_xt[NPHMAX][NSYSMX];	// Schmid tensors in xtal axes.<-> # of slip systems
extern voigtch *Schm_gr;	// Schmid tensors in sample axes <-> # of grid points
extern ten2nd B[6];		// used in chg_basis
extern real self[NPHMAX],coplane[NPHMAX],cross[NPHMAX],glissile[NPHMAX],hirth[NPHMAX],lomer[NPHMAX];
#ifdef EVP_GLOBALS_ONCE
int N_phases = 0;
int Type_phases[2] = {0};
real ElastConst_PhaseI[4] = {0.0};
real ElastConst_PhaseII[4] = {0.0};
char Slip_PhaseI[100] = {0};
char Slip_PhaseII[100] = {0};
char initial_ms[100] = {0};
ten4th Cijkl[NPHMAX] = {0.0};
real **nSRS = 0;
int *nSYS = 0;
int N_modes_max, N_modes;
char XtalSys[10] = {0};
real XtalAxis[3] = {0.0};
int *iMode = 0;
//real ***bn, ***bb;
real bn[NPHMAX][NSYSMX][3];
real bb[NPHMAX][NSYSMX][3];
voigt5 Schm_xt[NPHMAX][NSYSMX] = {0.0};
voigtch *Schm_gr;
ten2nd B[6] = {{0.0}};		// used in chg_basis
real self[NPHMAX],coplane[NPHMAX],cross[NPHMAX],glissile[NPHMAX],hirth[NPHMAX],lomer[NPHMAX];
#endif


/**********************
	Grains
  *********************/
extern real WGT;	// weight of each grid point, i.e. 1/Nxyz
extern real Wgt_ph1;	// weight of phase#1 (solid)
extern ten4th C0;	// avg Cijkl of grain essemble
extern ten4th S0;	// avg Sijkl of grain essemble
extern voigt66 C066;	// 6x6 form of C0
extern voigt66 S066;	// 6x6 form of S0
extern voigt66 *C_gr;	// Field containing the inhomogeneous Cijkl of grain essemble
extern int *grain_f;	// Field to identify grain #
extern int *phase_f;	// Field to identify phase #
extern ten2nd *TranMat_xt2sa;	// Field to store the transformation matrices (xtal -> sample)
#ifdef EVP_GLOBALS_ONCE
real WGT = 0.0;
real Wgt_ph1 = 0.0;
ten4th C0 = {0.0};
ten4th S0 = {0.0};
voigt66 C066 = {0.};
voigt66 S066 = {0.0};
voigt66 *C_gr = 0;
int *grain_f = 0;
int *phase_f = 0;
ten2nd *TranMat_xt2sa = 0;
#endif


/**********************
  hardening
  *********************/
extern int nsmx, isectwx;
extern real nrsx, gamd0x, twshx, tau0xf, tau0xb, tau1x, thet0, thet1, hselfx,hlatex2;
extern real ***trial_tau;	// trial crss
extern real ***crss;
extern real **gamdot;	// shear rate fields
extern real **xkin;
extern real **gammatotal;
extern real **tau_x;
extern real *gamacum;	// accumulated shear
extern real *SVM;
extern real *EVM;
extern real *work_p;
extern real *rho_crss;
extern real *rho_rss;
extern real tau[NPHMAX][NSYSMX][3];
extern real thet[NPHMAX][NSYSMX][2];
extern real gam0[NPHMAX][NSYSMX];
extern real hlatex[NSYSMX];
extern real Hard[NPHMAX][NSYSMX][NSYSMX];	// strain dependent hardening rates
#ifdef EVP_GLOBALS_ONCE
int nsmx, isectwx;
real nrsx, gamd0x, twshx, tau0xf, tau0xb, tau1x, thet0, thet1, hselfx,hlatex2;
real ***trial_tau = 0;	// trial crss
real ***crss = 0;
real **gamdot =0;	// shear rate fields
real **xkin = 0;
real **gammatotal = 0;
real **tau_x =0;
real *gamacum=0;	// accumulated shear
real *SVM=0;
real *EVM=0;
real *rho_rss=0;
real *rho_crss=0;
real *work_p=0;
real tau[NPHMAX][NSYSMX][3] = {0.0};
real thet[NPHMAX][NSYSMX][2] = {0.0};
real gam0[NPHMAX][NSYSMX] = {0.0};
real hlatex[NSYSMX] = {0.0};
real Hard[NPHMAX][NSYSMX][NSYSMX] = {0.0};
#endif


/**********************
  stress/strain fields
  macro stress/strain
  *********************/
extern ten2nd *Sig;	// stress in Eq. (2). After iteration, it should be the stress at time t+dt.
extern ten2nd *fluctuation;
extern ten2nd *Eps;		// total strain in Eq. (2). At beginning of iteration, it is the known quantity for time t.
extern ten2nd *Edot;	// plastic strain rate in Eq. (2). Note this is in correspondance to the current stress (Sig).
extern ten2nd SigAvg, EpsAvg, EdotAvg;
extern ten2nd SigDevAvg, VelGradAvg;
extern ten2nd SigAvg1;
extern ten2nd *DisGrad;		// displacement gradient field
extern ten2nd *VelGrad;		// velocity gradient field
extern real **displacement_fluct;
#ifdef EVP_GLOBALS_ONCE
ten2nd *Sig = 0;
ten2nd *fluctuation=0;
ten2nd *Eps = 0;
ten2nd *Edot = 0;
ten2nd SigAvg = {0.};
ten2nd SigAvg1 = {0.};
ten2nd SigDevAvg = {0.0};
ten2nd EpsAvg = {0.};
ten2nd EdotAvg = {0.0};
ten2nd VelGradAvg = {0.};
ten2nd *DisGrad = 0;
ten2nd *VelGrad = 0;
real **displacement_fluct = 0;
#endif


/**********************
  io.c
  *********************/
int GetNameList(char**);
void PrintNameList(void);
void OpenFiles(void);
void CloseFiles(void);
void WriteEdotMPI(char *s, int step);
void WriteEpsMPI(char *s, int step);
void WriteElsMPI(char *s, int step);
void WriteSVMMPI(char *s, int step);
void WriteEVMMPI(char *s, int step);
void WriteRhoRMPI(char *s, int step);
void WriteRhoCRMPI(char *s, int step);
void WriteWorkPMPI(char *s, int step);
void WriteWorkEMPI(char *s, int step);
void WriteSigMPI(char *s, int step);
void WriteTextureMPI(char *s, int step);
void WriteSLIPMPI(char *s, int step);
void WritetauxMPI(char *s, int step);
void WriteNewPositionMPI(char *s, int step);
void WriteCRSSMPI(char *s, int step);
void WriteTWINMPI(char *s, int step);


/**********************
  init.c
  *********************/
void SetupJob(void);
void DestroyJob(void);

/**********************
  kinematics.c
  *********************/
void VoigtToFull(real*, ten2nd);
void chg_basis(voigt V6, ten2nd T2, voigt66 C2, ten4th T4, int opt);
void chg_basis5(voigt5 V6, ten2nd T2, voigt66 C2, ten4th T4, int opt);
void SymAntDecompose(ten2nd t, ten2nd s, ten2nd a);
void Ten4thTransform(ten4th c1, ten2nd a, ten4th c2, int opt);
void EulerToTransMatrix(real *t1, real *ph, real *t2, ten2nd a, int opt);
void LU_inv_66(voigt66 c);
void update_schmid(void);


/**********************
  evolution.c
  *********************/
void Evolution(void);


/**********************
  constitutive.c
  *********************/
void StrainRate_eval(voigt,voigt,voigt66,int,int);
void get_trialtau(int, int);
void update_orient(void);
void harden(void);


/**********************
  testf.c
  *********************/
void PrintTensor(ten2nd);
void PrintTensorNorm(ten2nd);
void PrintB(void);
void PrintVoigt(voigt);
void WriteEtaMPI(char *s);	// record microstructure
void PrintVoigt66(voigt66 m);
void PrintSchmidXt(void);
void Ti_alpha_beta(void);

#endif	// end evp.h

