#include "evp.h"

static void UpdateStress(int istep, real *Err_e_local, real *Err_s_local)
{
	int jph;
	voigt66 sg66;
	ten2nd xlambda_aux, sig_aux, eps_aux, strain_aux,SigDev,eDev;
	voigt xlambda6, sig6, eps6, strain6;
	voigt edotp6;
	ten2nd edotp_aux,tot_e;
	voigt66 d_edotp66;
	voigt tot_eps;	// total strain
	voigt sig6_old;
	ten4th aux3333;
	voigt66 aux66;
	real signorm, epsnorm,trace_svm,trace_evm;
	real erroral, erral;
	real conv_NR, conv_istep_e, conv_istep_s;

	voigt res;	// residual R to be nullified
	voigt66	jacob_inv;	// Jacobian of R
	int itmaxal, iterl;
	int i,j,k;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			C6_loop{
				sg66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			LU_inv_66(sg66);
			
			T2_loop{
				/* the iteration starts with the current stress field */
				xlambda_aux[mi][mj] = Sig[pIDX][mi][mj];
				sig_aux[mi][mj] = Sig[pIDX][mi][mj];
				/* strain at time t. This is used in Eq. 4 to update total strain
				 from calculated strain rate*/
				eps_aux[mi][mj] = Eps[pIDX][mi][mj];
				/* DisGrad stores the updated displacement gradient obtained from Eq. 15.
				   So here strain_aux/strain6 stores the updated strain*/
				strain_aux[mi][mj] = (DisGrad[pIDX][mi][mj]+DisGrad[pIDX][mj][mi])/2.0;
			}
			chg_basis(xlambda6,xlambda_aux,aux66,aux3333,2);
			chg_basis(sig6,sig_aux,aux66,aux3333,2);
			chg_basis(eps6,eps_aux,aux66,aux3333,2);
			chg_basis(strain6,strain_aux,aux66,aux3333,2);

			signorm = 0.0;
			epsnorm = 0.0;
			T2_loop{
				signorm += xlambda_aux[mi][mj]*xlambda_aux[mi][mj];
				epsnorm += strain_aux[mi][mj]*strain_aux[mi][mj];
			}
			signorm = sqrt(signorm);
			epsnorm = sqrt(epsnorm);

			erroral = 1E-7;
			itmaxal = 100;
			iterl = 0;
			erral = 10*erroral;

			while((iterl<itmaxal)&&(fabs(erral)>fabs(erroral))){
				iterl++;

				for(i=0;i<6;i++) sig6_old[i] = sig6[i];
				StrainRate_eval(sig6,edotp6, d_edotp66,pIDX, jph);

				// total strain, Eq. 4
				for(i=0;i<6;i++){
					tot_eps[i] = eps6[i] + edotp6[i]*TimeStep;
					for(j=0;j<6;j++){
						tot_eps[i] += sg66[i][j]*sig6[j];
					}
				}

				// calculate the residual R, Eq. 16
				for(i=0;i<6;i++){
					res[i] = sig6[i] - xlambda6[i];
					for(j=0;j<6;j++){
						res[i] += C066[i][j]*(tot_eps[j]-strain6[j]);
					}
				}
				// calculate the Jacobian of R
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						jacob_inv[i][j] = (real)(i==j);
						for(k=0;k<6;k++){
							// Eq. 18
							jacob_inv[i][j] += C066[i][k]*(sg66[k][j]+d_edotp66[k][j]*TimeStep);
						}
					}
				}
				LU_inv_66(jacob_inv);

				// Newton-Raphson update
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						sig6[i] -= jacob_inv[i][j]*res[j];
					}
				}

				// convergence check
				conv_NR= 0.0;
				conv_istep_e= 0.0;
				conv_istep_s = 0.0;
				for(i=0;i<6;i++){
					conv_NR += (sig6[i]-sig6_old[i])*(sig6[i]-sig6_old[i]);
					conv_istep_e += (tot_eps[i]-strain6[i])*(tot_eps[i]-strain6[i]);
					conv_istep_s += (sig6[i]-xlambda6[i])*(sig6[i]-xlambda6[i]);
				}
				erral = sqrt(conv_NR)/signorm;
				conv_istep_e=sqrt(conv_istep_e);
				conv_istep_s=sqrt(conv_istep_s);

				// update crss
			//	if((Hard_Flag==1)&&(istep>2))
					/* because the shear rate that is about to use
					   is calcualted (in StrainRate_eval()) based on
					   trial stress field, we use the get_trialtau()
					   subroutine which is a trial versio nof harden() */
				//	get_trialtau(pIDX,jph);
			}	// end of while() loop

			chg_basis(sig6,sig_aux,aux66,aux3333,1);
			chg_basis(edotp6,edotp_aux,aux66,aux3333,1);
				chg_basis(tot_eps,tot_e,aux66,aux3333,1);
			// update stress and strain rate fields
			
			T2_loop{
				Sig[pIDX][mi][mj] = sig_aux[mi][mj];
				Edot[pIDX][mi][mj] = edotp_aux[mi][mj];
			}
				trace_svm = sig_aux[0][0]+sig_aux[1][1]+sig_aux[2][2];
	T2_loop{
		SigDev[mi][mj] = (sig_aux[mi][mj]+sig_aux[mj][mi])/2. - (real)(mi==mj)*trace_svm/3.0;
	}


	SVM[pIDX] = 0.0;
	T2_loop{
		SVM[pIDX] += SigDev[mi][mj]*SigDev[mi][mj];
	}
	SVM[pIDX] = sqrt(3./2.*SVM[pIDX]);	// for stress, it is 3/2
	
				trace_evm = tot_e[0][0]+tot_e[1][1]+tot_e[2][2];
	T2_loop{
		eDev[mi][mj] = (tot_e[mi][mj]+tot_e[mj][mi])/2. - (real)(mi==mj)*trace_evm/3.0;
	}
	
		EVM[pIDX] = 0.0;
	T2_loop{
		EVM[pIDX] += eDev[mi][mj]*eDev[mi][mj];
	}
	EVM[pIDX] = sqrt(2./3.*EVM[pIDX]);
			*Err_s_local += conv_istep_s;
			*Err_e_local += conv_istep_e;
		}
		else{
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
				Edot[pIDX][mi][mj] = 0.0;
			}
		}
	}

	return;
}/*end UpdateStress()*/

static void get_smacro(void)
{
	real local_sigavg, local_sigavg1;
	voigt sav6;
	voigt5 sav5;
	voigt66 aux66;
	ten4th aux3333;
	int i, ii,jj, k, kk,ll;

	int IJV[2][6] = {{0,1,2,1,0,0,},{0,1,2,2,2,1}};

	// overal stress
	T2_loop{
		local_sigavg = 0.0;
		local_sigavg1 = 0.0;
		local_loop{
			local_sigavg += Sig[pIDX][mi][mj];
			if(phase_f[pIDX]==1)
				local_sigavg1 += Sig[pIDX][mi][mj];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg, &SigAvg[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&local_sigavg1, &SigAvg1[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		SigAvg[mi][mj] *= WGT;
		SigAvg1[mi][mj] *= WGT/Wgt_ph1;
	}

	for(i=0;i<6;i++){
		ii=IJV[0][i];
		jj=IJV[1][i];
		dDisGradAvg[ii][jj] = 0.0;
		if(VelGrad_BC_Flag[i]==0){	// the component controlled by stress, implying a disp. variation
			for(k=0;k<6;k++){
				kk=IJV[0][k];
				ll=IJV[1][k];
				dDisGradAvg[ii][jj] += S0[ii][jj][kk][ll]*Stress_BC_Flag[k]*
					(Scauchy[kk][ll]-SigAvg[kk][ll]);
			}
		}
	}
	T2_loop{
		dDisGradAvg_acum[mi][mj] += dDisGradAvg[mi][mj];
	}
	if(rank==0){
		printf("dDisGradAvg(1,1),(2,2) = %e,%e\n",dDisGradAvg[0][0],dDisGradAvg[1][1]);
	}
	chg_basis(sav6,SigAvg,aux66,aux3333,2);
	for(i=0;i<5;i++) sav5[i] = sav6[i];
	chg_basis5(sav5,SigDevAvg,aux66,aux3333,1);

	s_vm = 0.0;
	T2_loop{
		s_vm += SigDevAvg[mi][mj]*SigDevAvg[mi][mj];

	}
	s_vm = sqrt(3./2.*s_vm);	// for stress, it is 3/2

	chg_basis(sav6,SigAvg1,aux66,aux3333,2);
	for(i=0;i<5;i++) sav5[i] = sav6[i];
	chg_basis5(sav5,SigDevAvg,aux66,aux3333,1);

	s_vm1 = 0.0;
	T2_loop{
		s_vm1 += SigDevAvg[mi][mj]*SigDevAvg[mi][mj];
	
	}
	s_vm1 = sqrt(3./2.*s_vm1);
				
	return;
}/*end get_smacro()*/

static real VonMises(ten2nd t)
{
	/* Calculate the von Mises equivalent of a non-symmetric,
	   non-traceless (disp. gradent or velocity gradent) tenor */
	real trace;
	ten2nd dt;
	real vm;

	trace = t[0][0]+t[1][1]+t[2][2];
	T2_loop{
		dt[mi][mj] = (t[mi][mj]+t[mj][mi])/2. - (real)(mi==mj)*trace/3.0;
	}
	vm = 0.0;
	T2_loop{
		vm += dt[mi][mj]*dt[mi][mj];
	}

	return (sqrt(2./3.*vm));	// NOTE: for strain, it's 2/3, NOT 3/2(for stress)

}/* end VonMises()*/


void Evolution(void)
{
	int istep;		// step # of deformation test
	int iter;		// step # of iteration (within a given istep)
	ten2nd *kSig_r, *kSig_i;
	ten2nd DisGradAvg_t = {0.0};	// store the actual macro disp. grad. at time t
	ten2nd DisGradAvg_actual = {0.0};
	ten2nd sym_du_r, sym_du_i;
   	voigt str;
	ten2nd sig;
		ten4th aux3333;
	voigt66 aux66;
	real rss[NSYSMX];
	voigt5 sc[NSYSMX];
	int i,j,jph;
	ten2nd EpsAvg_local, EdotAvg_local;
	real evmp, dvmp;
	real Err_e_local, Err_s_local;
    int iter_total=0;
    real tmp[3];
    real tmp_re[3];
    real tmp_im[3];
    real tmp1;
	AllocMem(kSig_r,lsize,ten2nd);
	AllocMem(kSig_i,lsize,ten2nd);

	if(rank==0){
		printf("\n\n======================================\n");
		printf("-------------Simulation starts----------------\n");
	}

	for(istep=1;istep<=N_steps;istep++){
		if(rank==0){
			printf("\n****************************************\n");

			printf("STEP = %d\n",istep);
			if(N_steps!=1){
				fprintf(fp_err,"STEP = %d\n",istep);
			}
		}

		local_loop{
			T2_loop{
				// store displacement gradient of current time into velgrad
				VelGrad[pIDX][mi][mj] = DisGrad[pIDX][mi][mj];
				// enforce macro deformation
				if(CREEP_FLAG==1){
					DisGrad[pIDX][mi][mj] = DisGradAvg_actual[mi][mj];
				}
        else{
				  DisGrad[pIDX][mi][mj] += Udot[mi][mj]*TimeStep;
        }
			}
		}

		if(istep==1 || Update_Flag==1){
			update_schmid();
		}

		/* "dDisGradAvg" is the prescribed strain E_pq in E. 23.
		   It is updated iteratively during solving the
		   micromechanical problem (the fowllowing while() loop).
		   So the intital guess can simply be set to zero, since
		   in get_smacro() it will be updated based on updated
		   stress field */
		T2_loop{
			dDisGradAvg[mi][mj] = 0.0;
			dDisGradAvg_acum[mi][mj] = 0.0;
		}

		/***********************************
		  iteration to update stress at t+dt
		  ***********************************/
		iter = 0;
		Err_e = 2.*Err;
		Err_s = 2.*Err;
		
		while((iter<IterMax)&&(MAX(fabs(Err_s),fabs(Err_e))>fabs(Err))){
			iter++;
			iter_total++;
			if(rank==0){
				printf("\nITER = %d\n",iter);
				printf("Forward FFT of stress field\n");
			}
			// k-space stress field
			T2_loop{
				local_loop{
					fft_data[pIDX].re = Sig[pIDX][mi][mj];
					fft_data[pIDX].im = 0.0;
				}
				MPI_Barrier(MPI_COMM_WORLD);
				fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
				local_loop{
					kSig_r[pIDX][mi][mj] = fft_data[pIDX].re;
					kSig_i[pIDX][mi][mj] = fft_data[pIDX].im;
				}
			}

			local_loop{
				T2_loop{
					sym_du_r[mi][mj] = 0.0;;
					sym_du_i[mi][mj] = 0.0;;
					T2p_loop{
						sym_du_r[mi][mj] += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_r[pIDX][mip][mjp];
						sym_du_i[mi][mj] += GAMMA[pIDX][mi][mj][mip][mjp]*kSig_i[pIDX][mip][mjp];
					}
				}
				T2_loop{
					kSig_r[pIDX][mi][mj] = sym_du_r[mi][mj];
					kSig_i[pIDX][mi][mj] = sym_du_i[mi][mj];
				}

			}

			if(rank==0){
				printf("Inverse FFT to get strain field\n");
			}
			// update strain field in real space
			T2_loop{
				local_loop{
					fft_data[pIDX].re = kSig_r[pIDX][mi][mj];
					fft_data[pIDX].im = kSig_i[pIDX][mi][mj];
				}
				MPI_Barrier(MPI_COMM_WORLD);
				fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
				local_loop{
					DisGrad[pIDX][mi][mj] += dDisGradAvg[mi][mj] + fft_data[pIDX].re/Nxyz;
					fluctuation[pIDX][mi][mj] +=  fft_data[pIDX].re/Nxyz;
				}
			}

			if(rank==0){
				printf("Update stress field\n");
			}
			// update stress, which requires Newton-Raphson method
			// NOTE: N-R is run locally. So PEs could run for different iteration steps
			Err_e_local = 0.0;
			Err_s_local = 0.0;
			UpdateStress(istep, &Err_e_local, &Err_s_local);	
			// collect errors
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&Err_e_local, &Err_e, 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&Err_s_local, &Err_s, 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			Err_e *= WGT;
			Err_s *= WGT;

			//if(rank==0){
			//	printf("ERRE = %e\n",Err_e);
			//	printf("ERRS = %e\n",Err_s);
			//}

			get_smacro();

      if(fabs(e_vm)>1E-5){
			  Err_e /= e_vm;	
      }
			Err_s /= s_vm;
			if(rank==0){
			//	printf("Strain field error = %f\n", Err_e);
			//	printf("Stress field error = %e\n", Err_s);
				fprintf(fp_err,"%d\t%e\t%e\t%e\n",iter,Err_e,Err_s,s_vm);
        fflush(fp_err);
			}

		}

		local_loop{
			T2_loop{
				VelGrad[pIDX][mi][mj] = (DisGrad[pIDX][mi][mj] - VelGrad[pIDX][mi][mj])/TimeStep;
			}
		}
		T2_loop{
			// dDisGradAvg_acum is updated in get_smacro()
			DisGradAvg_actual[mi][mj] = DisGradAvg[mi][mj] + dDisGradAvg_acum[mi][mj];
			VelGradAvg[mi][mj] = (DisGradAvg_actual[mi][mj]-DisGradAvg_t[mi][mj])/TimeStep;
		}

		if(rank==0){
			printf("DisGradAvg(1,1),DisGradAvg(2,2),DisGradAvg(3,3)\n");
			printf("%e,%e,%e\n",DisGradAvg_actual[0][0], DisGradAvg_actual[1][1], DisGradAvg_actual[2][2]);
			printf("DisGradAvg(1,1)/DisGradAvg(3,3)\n");
			printf("%e\n",DisGradAvg_actual[0][0]/DisGradAvg_actual[2][2]);
			printf("SigAvg1(1,1),SigAvg1(2,2),SigAvg1(3,3)\n");
			printf("%e,%e,%e\n", SigAvg1[0][0],SigAvg1[1][1],SigAvg1[2][2]);
		}
		e_vm = VonMises(DisGradAvg_actual);
		d_vm = VonMises(VelGradAvg);
		T2_loop{
			DisGradAvg_t[mi][mj] = DisGradAvg_actual[mi][mj];
		}

		// Initial guess of macro disp. gradient at t+dt, always elastic
		T2_loop{
			DisGradAvg[mi][mj] = DisGradAvg_t[mi][mj] + Udot[mi][mj]*TimeStep;
		}

		// update strain field
		local_loop{
			T2_loop{
				Eps[pIDX][mi][mj] += Edot[pIDX][mi][mj]*TimeStep;	// Edot was updated in UpdateStress()
			}
		}

		// Plastic VM
		T2_loop{
			EpsAvg_local[mi][mj] = 0.0;
			EdotAvg_local[mi][mj] = 0.0;
		}
		local_loop{
			T2_loop{
				EpsAvg_local[mi][mj] += Eps[pIDX][mi][mj];
				EdotAvg_local[mi][mj] += Edot[pIDX][mi][mj];
			}
		}
		T2_loop{
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&EpsAvg_local[mi][mj], &EpsAvg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&EdotAvg_local[mi][mj], &EdotAvg[mi][mj], 1, MPI_real,
					MPI_SUM, MPI_COMM_WORLD);
			EpsAvg[mi][mj] *= WGT;
			EdotAvg[mi][mj] *= WGT;
		}

		evmp = 0.0;
		dvmp = 0.0;
		T2_loop{
			evmp += EpsAvg[mi][mj]*EpsAvg[mi][mj];
			dvmp += EdotAvg[mi][mj]*EdotAvg[mi][mj];
		}
		evmp = sqrt(2./3.*evmp);
		dvmp = sqrt(2./3.*dvmp);

		if((Update_Flag==1)&&(istep>1)){
			// grain reorientation
			update_orient();
		}

		if((Hard_Flag==1)&&(istep>2)){
			// hardening
			// update "crss" and "gamacum"
			harden();
						local_loop{
				jph = phase_f[pIDX];
				T2_loop{
			    sig[mi][mj]=Sig[pIDX][mi][mj];
			}
			rho_rss[pIDX]=0.0;
			rho_crss[pIDX]=0.0;
			work_p[pIDX]=0.0;
			 	chg_basis(str,sig,aux66,aux3333,2);
                        for(i=0;i<nSYS[jph-1];i++){
                            
                            	for(j=0;j<5;j++) sc[i][j] = Schm_gr[pIDX][i][j];
					rss[i] = sc[i][0]*str[0]+sc[i][1]*str[1]+sc[i][2]*str[2]+
			sc[i][3]*str[3]+sc[i][4]*str[4];
				rho_rss[pIDX]+=pow(rss[i]/(0.35*47.565*0.256),2);
					rho_crss[pIDX]+=pow(crss[pIDX][i][0]/(0.35*47.565*0.256),2);
			//		fract[pIDX][i]=abs(rss[pIDX][i]/crss[pIDX]);
			//work_p[pIDX]+=fabs((gamtot+deltgam)*rss[i]);
				work_p[pIDX]+=fabs((gammatotal[pIDX][i])*rss[i]);
			}
		}
		}
		
		    for(mi=0;mi<3;mi++){
				    tmp_re[mi]=0.0;
				    tmp_im[mi]=0.0;
				    for(mj=0;mj<3;mj++){
				
					local_loop{
						fft_data[pIDX].re = fluctuation[pIDX][mi][mj];
						fft_data[pIDX].im = 0.0;
					}
					MPI_Barrier(MPI_COMM_WORLD);
					fftwnd_mpi(plan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
				local_loop{
                     tmp[0] = g[pIDX].x;
                      tmp[1] = g[pIDX].y;
                      tmp[2] = g[pIDX].z;  
                     tmp1 = fft_data[pIDX].re;
           tmp_re[mi] += (-1)*fft_data[pIDX].im*tmp[mj];
           tmp_im[mi] += tmp1*tmp[mj];
}
}
            fft_data[pIDX].re = tmp_re[mi];
            fft_data[pIDX].im = tmp_im[mi];
			MPI_Barrier(MPI_COMM_WORLD);
            fftwnd_mpi(iplan, 1, fft_data, fft_work, FFTW_NORMAL_ORDER);
            local_loop{
            displacement_fluct[pIDX][mi] = - fft_data[pIDX].re/Nxyz;
				}
							
				}

		if(rank==0){
			fprintf(fp_vm,"%e\t%e\t%e\t%e\t%e\t%e\t%d\n",
					e_vm,evmp,d_vm,dvmp,s_vm,s_vm1,iter_total);
      fflush(fp_vm);
		}

		// record stress/strain/strain rate fields
		if((PrintControl[0]==1)&&(istep/PrintControl[1]*PrintControl[1]==istep)){
			WriteEpsMPI("eps",istep);
			WriteSigMPI("sig",istep);
	    	WriteElsMPI("els",istep);
		    WriteEdotMPI("edot",istep);
                        if(Tex_Flag==1){
		WriteSLIPMPI("slip_sys",istep);
		WritetauxMPI("taux_sys",istep);
		WriteCRSSMPI("crss_sys",istep);
                //WriteTWINMPI("twin_sys",istep);
        WriteTextureMPI("tex",istep);
        WriteRhoRMPI("rho_rss",istep);
        WriteRhoCRMPI("rho_crss",istep);
        WriteWorkPMPI("work_p",istep);
      //  WriteWorkEMPI("work_e",istep);
        WriteSVMMPI("SVM",istep);
        WriteEVMMPI("EVM",istep);
        WriteNewPositionMPI("new_position",istep);

		}
		
		
	}
        }// end of mechanical test

	if(Tex_Flag==1){
		
WriteSLIPMPI("slip_final",0);
WritetauxMPI("taux_final",0);
WriteCRSSMPI("crss_final",0);
WriteTextureMPI("tex_final",0);
WriteSVMMPI("SVM_final",0);
WriteEVMMPI("EVM_final",0);
WriteRhoRMPI("rho_rss_final",0);
WriteRhoCRMPI("rho_crss_final",0);
WriteWorkPMPI("work_p_final",0);
WriteWorkEMPI("work_e_final",0);
                //WriteTWINMPI("twin_sys",0);
WriteEpsMPI("eps_final",0);
WriteSigMPI("sig_final",0);
WriteNewPositionMPI("new_position_final",0);
WriteEdotMPI("edot_final",0);
WriteElsMPI("els_final",0);
	}

	if(rank==0){
		printf("-------------Simulation ends----------------\n");
		printf("======================================\n\n");
	}

	free(kSig_r);
	free(kSig_i);
	return;
}/*end Evolution()*/
