#include "evp.h"


void StrainRate_eval(voigt stress, voigt edot, voigt66 d_edot, int idx, int jph)
{
	int i,j,k;
	int isign;
	real rss[NSYSMX];	// resolved shear stress
	real rss1[NSYSMX];
	real rss2[NSYSMX];
	voigt5 sc[NSYSMX];
	real tau_t[NSYSMX][2];	// crss
	real nsr[NSYSMX];	// strain rate sensitivity
	real xkinaux[NSYSMX];

	for(i=0;i<nSYS[jph-1];i++){
		nsr[i] = nSRS[jph-1][i];	// SRS could be different for different slip systems
		tau_t[i][0] = trial_tau[idx][i][0];
		tau_t[i][1] = trial_tau[idx][i][1];
		xkinaux[i] = xkin[idx][i];
		for(j=0;j<5;j++) sc[i][j] = Schm_gr[idx][i][j];
	}

	// calculate resolved shear stress and shear rate
	for(i=0;i<nSYS[jph-1];i++){
		rss[i] = sc[i][0]*stress[0]+sc[i][1]*stress[1]+sc[i][2]*stress[2]+
			sc[i][3]*stress[3]+sc[i][4]*stress[4];
		if((rss[i]-xkinaux[i])<0){
			isign = 1;
		}
		else{
			isign = 0;
		}
		rss[i] = (rss[i]-xkinaux[i])/tau_t[i][isign];
                
              /*  if((i >= nSYS[jph-1]/2) && (rss[i] < 0)) {
                    rss[i] = 0.0;
                } */
                
		rss1[i] = gam0[jph-1][i]*nsr[i]*fabs(pow(rss[i],nsr[i]-1.))/tau_t[i][isign];
		rss2[i] = gam0[jph-1][i]*fabs(pow(rss[i],nsr[i]))*((real)(2*(rss[i]>0)-1));
	
		// shear rate for each system
		gamdot[idx][i] = rss2[i];
		tau_x[idx][i]  = rss1[i];
	}

	// plastic strain rate
	for(i=0;i<5;i++){
		edot[i] = 0.;
		for(j=0;j<nSYS[jph-1];j++){
			edot[i] += sc[j][i]*rss2[j];
		}
	}
	edot[5] = 0.0;

	for(i=0;i<5;i++){
		for(j=0;j<5;j++){
			d_edot[i][j] = 0.0;
			for(k=0;k<nSYS[jph-1];k++){
				d_edot[i][j] += sc[k][i]*sc[k][j]*rss1[k];
			}
		}
	}

	for(i=0;i<6;i++){
		d_edot[i][5] = 0.0;
		d_edot[5][i] = 0.0;
	}
	
	return;
}/*end StrainRate_eval()*/

void get_trialtau(int idx, int jph)
{
	int i,j;
	real gamtot, deltgam;
	real dtau;
	real tau0, tau1, thet0, thet1;
	real voce, fact, exp_ini, exp_del;
	real TINY;

	gamtot = gamacum[idx];
	deltgam = 0.0;
	for(i=0;i<nSYS[jph-1];i++)
		deltgam += fabs(gamdot[idx][i])*TimeStep;

	for(i=0;i<nSYS[jph-1];i++){
		dtau = 0.0;
		for(j=0;j<nSYS[jph-1];j++)
			dtau += Hard[jph-1][i][j]*fabs(gamdot[idx][j])*TimeStep;
		tau0 = tau[jph-1][i][0];
		tau1 = tau[jph-1][i][2];
		thet0 = thet[jph-1][i][0];
		thet1 = thet[jph-1][i][1];
		TINY = 1.E-4*tau0;

		voce = 0.0;
		if(fabs(thet0)>TINY){
			voce = thet1*deltgam;
			if(fabs(tau1)>TINY){
				fact = fabs(thet0/tau1);
				exp_ini = exp(-1.0*gamtot*fact);
				exp_del = exp(-1.0*deltgam*fact);
				voce += -1.*(fact*tau1-thet1)/fact*exp_ini*(exp_del-1.0)-
					thet1/fact*exp_ini*(exp_del*((gamtot+deltgam)*fact+1.0)-(gamtot*fact+1.0));
			}
		}

		trial_tau[idx][i][0] = crss[idx][i][0]+dtau*voce/deltgam;
		trial_tau[idx][i][1] = crss[idx][i][1]+dtau*voce/deltgam;
	}

	return;
}/*end get_trialtau()*/

static void Orientation(ten2nd tran, ten2nd rot)
{
	int i;
	real w[3];
	real w_norm, w_tan;
	ten2nd omega,omega2;
	ten2nd tran_new;

	w[0] = rot[2][1];
	w[1] = rot[0][2];
	w[2] = rot[1][0];

	w_norm = 0.0;
	for(i=0;i<3;i++)
		w_norm += w[i]*w[i];
	w_norm = sqrt(w_norm);
	w_tan = tan(w_norm/2.0);
	if(fabs(w_norm)<1E-6)
		w_norm = 1.0;
	for(i=0;i<3;i++)
		w[i] *= w_tan/w_norm;

	// construct "normalized" omega matrix
	omega[0][0]=0.; omega[0][1]=-1.*w[2]; omega[0][2]=w[1];
	omega[1][0]=w[2]; omega[1][1]=0.; omega[1][2]=-1.*w[0];
	omega[2][0]=-1.*w[1]; omega[2][1]=w[0]; omega[2][2]=0.;
	T2_loop{
		omega2[mi][mj] = 0.0;
		for(i=0;i<3;i++){
			omega2[mi][mj] += omega[mi][i]*omega[i][mj];
		}
	}

	w_tan = w_tan*w_tan;
	T2_loop{
		rot[mi][mj] = (real)(mi==mj) +
			2.*(omega[mi][mj]+omega2[mi][mj])/(1.+w_tan);
	}

	// transf. matrix at t+dt
	T2_loop{
		tran_new[mi][mj] = 0.0;
		for(i=0;i<3;i++)
			tran_new[mi][mj] += rot[mi][i]*tran[i][mj];
	}

	// record new transf. matrix
	T2_loop{
		tran[mi][mj] = tran_new[mi][mj];
	}

	return;
}/*end Orientation()*/

void update_orient(void)
{
	int is, i,j;
	int jph;
	ten2nd tranmat;	// xtal -> sample (grain)
	ten2nd RotLoc;
	ten2nd Lp, RotSlip;
	ten2nd Rot;	// total rotation
	real bb_sa[3], bn_sa[3];
	real rsl_bar, rlc_bar;
	real rsl_bar_tot, rlc_bar_tot;


	rsl_bar = 0.0;
	rlc_bar = 0.0;
	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			// local rotate rate
			T2_loop{
				RotLoc[mi][mj] = (VelGrad[pIDX][mi][mj]-VelGrad[pIDX][mj][mi])/2.;
				Lp[mi][mj] = 0.0;
				tranmat[mi][mj] = TranMat_xt2sa[pIDX][mi][mj];
			}

			// slip ration rate
			for(is=0;is<nSYS[jph-1];is++){
				for(i=0;i<3;i++){
					bb_sa[i] = 0.0;
					bn_sa[i] = 0.0;
					for(j=0;j<3;j++){
						bb_sa[i] += tranmat[i][j]*bb[jph-1][is][j];
						bn_sa[i] += tranmat[i][j]*bn[jph-1][is][j];
					}
				}
				T2_loop{
					Lp[mi][mj] += bb_sa[mi]*bn_sa[mj]*gamdot[pIDX][is];
				}
			}
			T2_loop{
				RotSlip[mi][mj] = (Lp[mi][mj]-Lp[mj][mi])/2.;
			}

			// avg rotation rate
			rsl_bar += sqrt(RotSlip[2][1]*RotSlip[2][1] +
					RotSlip[0][2]*RotSlip[0][2] + RotSlip[1][0]*RotSlip[1][0])*WGT;
			rlc_bar += sqrt(RotLoc[2][1]*RotLoc[2][1] +
					RotLoc[0][2]*RotLoc[0][2] + RotLoc[1][0]*RotLoc[1][0])*WGT;

			// total rotation
			T2_loop{
				// Udot_a is the applied rotation (anti-symm of Udot)
				// RotSlip does NOT change crystal axes orientation
				Rot[mi][mj] = (Udot_a[mi][mj]+RotLoc[mi][mj]-RotSlip[mi][mj])*TimeStep;
			}

			// reorientation
			Orientation(tranmat,Rot);

			// update Transformation matrix
			T2_loop{
				TranMat_xt2sa[pIDX][mi][mj] = tranmat[mi][mj];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&rsl_bar, &rsl_bar_tot, 1, MPI_real,
			MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rlc_bar, &rlc_bar_tot, 1, MPI_real,
			MPI_SUM, MPI_COMM_WORLD);

	if(rank==0){
		printf("Average plastic rotation = %e\n",rsl_bar_tot);
		printf("Average local rotation = %e\n",rlc_bar_tot);
	}

	return;
}/*end update_orient()*/

void harden(void)
{
	int i,j;
	int jph;
	real gamtot, deltgam;
        real deltagamma[NSYSMX];
	real dtau;
	real tau0, tau1, thet0, thet1;
	real voce, fact, exp_ini, exp_del;
	real TINY;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			gamtot = gamacum[pIDX];
                        
                        for(i=0;i<nSYS[jph-1];i++){
				deltagamma[i] = (gamdot[pIDX][i])*TimeStep;
                                gammatotal[pIDX][i]  = gammatotal[pIDX][i] + deltagamma[i];
                        }


			deltgam = 0.0;
			for(i=0;i<nSYS[jph-1];i++)
				deltgam += fabs(gamdot[pIDX][i])*TimeStep;

			for(i=0;i<nSYS[jph-1];i++){
				dtau = 0.;
				for(j=0;j<nSYS[jph-1];j++){
					dtau += Hard[jph-1][i][j]*fabs(gamdot[pIDX][j])*TimeStep;
				}
				tau0 = tau[jph-1][i][0];
				tau1 = tau[jph-1][i][2];
				thet0 = thet[jph-1][i][0];
				thet1 = thet[jph-1][i][1];
				TINY = 1.E-4*tau0;

				voce = 0.0;
				if(fabs(thet0)>TINY){
					voce = thet1*deltgam;
					if(fabs(tau1)>TINY){
						fact = fabs(thet0/tau1);
						exp_ini = exp(-1.0*gamtot*fact);
						exp_del = exp(-1.0*deltgam*fact);
						voce += -1.*(fact*tau1-thet1)/fact*exp_ini*(exp_del-1.0)-
							thet1/fact*exp_ini*(exp_del*((gamtot+deltgam)*fact+1.0)-(gamtot*fact+1.0));
					}
				}
				crss[pIDX][i][0] += dtau*voce/deltgam;
				crss[pIDX][i][1] += dtau*voce/deltgam;

				trial_tau[pIDX][i][0] = crss[pIDX][i][0];
				trial_tau[pIDX][i][1] = crss[pIDX][i][1];
			}
			gamacum[pIDX] = gamtot+deltgam;
		}
	}


	return;
}/*end harden()*/

