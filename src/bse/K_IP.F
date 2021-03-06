!
!        Copyright (C) 2000-2020 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS AM
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine K_IP(iq,Ken,Xk,Dip,W_bss)
 !
 use pars,          ONLY:SP,IP,pi,cZERO,cONE,cI
 use IO_m,          ONLY:io_DIP
 use frequency,     ONLY:w_samp
 use LIVE_t,        ONLY:live_timing
 use drivers,       ONLY:l_rpa_IP,l_bs_fxc
 use com,           ONLY:msg
 use units,         ONLY:HA2EV
 use vec_operate,   ONLY:v_norm,sort,degeneration_finder
 use R_lattice,     ONLY:bz_samp,d3k_factor,bare_qpg,q0_def_norm,qindx_X
 use stderr,        ONLY:intc
 use wave_func,     ONLY:WF
 use interfaces,    ONLY:WF_load,WF_free
 use electrons,     ONLY:levels,spin_occ,spin,nel,n_sp_pol
 use X_m,           ONLY:Epsilon_ii,Joint_DOS,alpha_dim,eps_2_alpha,&
&                        X_drude_term,global_gauge,BS_E_sorted,   &
&                        skip_cond_sum_rule,l_drude
 use DIPOLES,       ONLY:DIPOLE_t,DIP_alloc,DIP_rotated
 use BS_solvers,    ONLY:BSS_Vnl_included,BSS_add_antiresonant,&
&                        Co_factor,diam_term_exact,BSS_Wd
 use BS,            ONLY:BS_T_grp,BS_K_coupling,BS_nT_grps,BS_perturbative_SOC,L_kind,BS_bands,&
 &                       BS_dip_size,l_BS_jdos,l_BS_esort,l_BS_trace,BS_res_ares_dim,BS_K_dim,BS_H_dim,&
 &                       l_BS_abs,l_BS_kerr,l_BS_magnons,l_BS_dichroism,l_BS_photolum,l_BS_optics,&
 &                       BS_mat_res_ares_dim
 use parallel_int,  ONLY:PP_redux_wait
 use parallel_m,    ONLY:PAR_BS_nT_col_grps,PAR_COM_eh_INDEX,PAR_IND_T_groups,&
&                        PAR_COM_Xk_ibz_INDEX,myid
 use collision_el,  ONLY:elemental_collision,elemental_collision_free,elemental_collision_alloc
 use parser_m,      ONLY:parser
 use timing_m,      ONLY:timing
 use vec_operate,   ONLY:rot_v_xy_plane
 use BS,            ONLY:l_BS_anomalous_Hall
 use BS_solvers,    ONLY:A_drude,B_Hall
 use PHOTOLUM,      ONLY:PL,PL_prefactor,PL_weights
 use X_m,           ONLY:X_magnons,Epsilon_ij,X_dichroism
 use units,         ONLY:SPEED_OF_LIGHT
 !
#include<memory.h>
 !
 type(levels) ::Ken
 type(bz_samp)::Xk
 type(DIPOLE_t)::Dip
 type(w_samp) ::W_bss
 integer      ::iq
 !
 ! Work space
 !
 integer           ::ik_bz,ok_bz,ik,ok,isymm,osymm,iGo,iv,ic,i_sp_pol_c,i_sp_pol_v,i_T_g,i_T,&
 &                   n_SOC,i_pert_SOC,io_err,i_res_ares,i_E_sort
 real(SP)          ::f_
 logical           ::dip_v,dip_r
 complex(SP)       ::drude_GreenF(W_bss%n_freqs),res_abs,E_plus_W0,E_plus_W(2/n_sp_pol),&
&                    para_term_w0,diam_term_used,Z_,E_
 !
 integer  :: i_VAR,i_grp,i_cv,nVAR,EDIM,n_deg_grp
 real(SP), allocatable     :: vtmpR(:)
 complex(SP), allocatable  :: BS_E_sorted_tmp(:,:)
 integer, allocatable      :: sort_indx(:),first_el(:),n_of_el(:)
 !
 real(SP)          ::f_PL,res_PL
 complex(SP)       ::DIP_Splus,DIP_Smins,P_weighted(3)
 complex(SP)       ::res_kerr,res_magn(2),res_dich
 complex(SP)       ::DIP_projected_dich(BS_dip_size)
 !
 character(6)             :: DIP_kind
 complex(SP)              :: DIP_projected(BS_dip_size),DIP_expanded(3)
 real(SP)                 :: field_dir(3,BS_dip_size)
 type(elemental_collision):: BSE_scatt
 !
 !
 call section('=','Independent Particles properties @q'//trim(intc(iq)))
 !
 if (l_bs_fxc) return
 !
 !===================================================
 !
 call timing('T_space IP',OPR='start')
 !
 ! Polarizability ?
 !==================
 if (trim(alpha_dim)/='adim') then
   call msg('r', 'Optical renormalization       ',eps_2_alpha,"[a.u.]")
   call msg('rn','Polarizability dimension      ',trim(alpha_dim))
 endif
 !
 if(trim(global_gauge)=="velocity") call parser('NoCondSumRule',skip_cond_sum_rule)
 !
 !...Drude
 call X_Drude(1,Ken,Xk,W_bss,BSS_Wd,drude_GreenF)
 !
 ! Oscillators I/O & Eps0
 !========================
 !
 if((.not.BS_perturbative_SOC).or.n_sp_pol==2) then
   n_SOC=1
   Co_factor=real(spin_occ)/(2._SP*pi)**3*d3k_factor*4._SP*pi
 endif
 !
 !
 if( BS_perturbative_SOC.and.n_sp_pol==1) then
   n_SOC=2
   Co_factor=    1._SP     /(2._SP*pi)**3*d3k_factor*4._SP*pi
 endif
 !
 diam_term_exact=cmplx(nel*real(Xk%nbz,SP)/spin_occ,0._SP,SP)
 !
 para_term_w0=cZERO
 !
 ! The residuals of the BSE are <wf_k|e^{iqr}|wf_k+q>
 ! These are trivial to compute at finite q, instead they are 
 ! obtained via the dipoles in the q->0 limit and stored to disk
 !
 if (iq==1) then
   call DIPOLE_IO(Xk,Ken,Dip,'read ',io_err,'K')
 else
   call WF_load(WF,1,maxval(qindx_X(:,:,2)),BS_bands,(/1,Xk%nibz/),title='-Oscillators/R space/K_IP')
   call elemental_collision_free(BSE_scatt)
   call elemental_collision_alloc(BSE_scatt,NG=1,TITLE="BSE") 
 endif
 !
 BSS_Vnl_included=Dip%Vnl_included
 !
 if (iq==1) then
   !
   ! Directions for the external field
   !
   field_dir(:,1)=Dip%q0/v_norm(Dip%q0)
   !
   if(l_BS_kerr) call rot_v_xy_plane(-pi/2._SP,field_dir(:,1),v_out=field_dir(:,2))
   !
   if (l_BS_trace) then
     field_dir(:,1)=(/1._SP,0._SP,0._SP/)
     field_dir(:,2)=(/0._SP,1._SP,0._SP/)
     field_dir(:,3)=(/0._SP,0._SP,1._SP/)
   endif
   !
 endif
 !
 if (trim(global_gauge)=='length'  ) then
   DIP_kind ="DIP_iR"
   dip_r=.true.
   dip_v=.false.
   field_dir=field_dir*q0_def_norm
 endif
 if (trim(global_gauge)=='velocity') then
   DIP_kind ="DIP_v"
   dip_r=.false.
   dip_v=.true.
 endif
 !
 if(l_BS_esort) then
   EDIM=BS_K_dim(1)
   nVAR=2
   if(l_BS_abs) nVAR=nVAR+1
   if(l_BS_kerr)      nVAR=nVAR+1
   if(l_BS_dichroism) nVAR=nVAR+1
   if(l_BS_magnons)   nVAR=nVAR+2
   YAMBO_ALLOC(BS_E_sorted_tmp,(EDIM,nVAR))
   BS_E_sorted_tmp=cZERO
 endif
 !
 if (l_rpa_IP.or.iq>1) call live_timing('IP properties',PAR_BS_nT_col_grps)
 !
 do i_T_g=1,BS_nT_grps
   !
   if (.not.PAR_IND_T_groups%element_1D(i_T_g)) cycle
   !
   i_res_ares=BS_T_grp(i_T_g)%i_res_ares
   !
   do i_T=1,BS_T_grp(i_T_g)%size
     !
     ! a.Indexes
     !==========
     !
     ik_bz = BS_T_grp(i_T_g)%table(i_T,1)
     iv    = BS_T_grp(i_T_g)%table(i_T,2)
     ic    = BS_T_grp(i_T_g)%table(i_T,3)
     i_sp_pol_c = BS_T_grp(i_T_g)%table(i_T,4)
     i_sp_pol_v = BS_T_grp(i_T_g)%table(i_T,5)
     !
     if(iq/=1) then
       !
       ok_bz=qindx_X(iq,ik_bz,1)
       iGo  =qindx_X(iq,ik_bz,2)
       !
       ik   =Xk%sstar(ik_bz,1) 
       isymm=Xk%sstar(ik_bz,2) 
       !
       ok   =Xk%sstar(ok_bz,1)
       osymm=Xk%sstar(ok_bz,2)
       !
       BSE_scatt%is=(/ic,ik,isymm,i_sp_pol_c/)
       BSE_scatt%os=(/iv,ok,osymm,i_sp_pol_v/)
       BSE_scatt%qs=(/iGo,iq,1/)
       !
     endif
     !
     ! b.Energy
     !=========
     !
     if(allocated(BS_T_grp(i_T_g)%W)) then
       E_plus_W0=BS_T_grp(i_T_g)%E(i_T)-cI*BS_T_grp(i_T_g)%W(i_T)
     else
       E_plus_W0=BS_T_grp(i_T_g)%E(i_T)
     endif
     !
     if (.not.BS_perturbative_SOC) E_plus_W=E_plus_W0
     !
     if (     BS_perturbative_SOC) E_plus_W(:)=E_plus_W0+BS_T_grp(i_T_g)%E_SOC_corr(:,i_T)
     !
     if(l_BS_esort.and.i_res_ares==1) then
       i_E_sort=sum(BS_T_grp(BS_T_grp(i_T_g)%i_T_ref:i_T_g-1)%size)+i_T !+(i_res_ares-1)*BS_K_dim(1)
       BS_E_sorted_tmp(i_E_sort,1)=E_plus_W0
       nVAR=2
     endif
     !
     ! c1.Compute "dipoles"
     !=====================
     !
     ! *length*   : DIP_projected(c,v,k) is     q*<v|r|c> while I need     q*<c|r|v> = conjg(    q*<v|r|c>)
     ! *velocity* : DIP_projected(c,v,k) is q_ver*<v|p|c> while I need q_ver*<c|r|v> = conjg(q_ver*<v|r|c>)
     !              see definition of conjg_fac few lines before for the -1
     !
     ! Absorption, PL & KERR
     !
     if(l_BS_optics) then
       if(iq==1) then
         DIP_expanded=-cI*DIP_rotated(ic,iv,ik_bz,i_sp_pol_c,DIP_kind,Xk)
         DIP_projected(1)=dot_product(field_dir(:,1),DIP_expanded)
         if(BS_dip_size>1) DIP_projected(2)=dot_product(field_dir(:,2),DIP_expanded)
         if(BS_dip_size>2) DIP_projected(3)=dot_product(field_dir(:,3),DIP_expanded)
       else
         if(trim(global_gauge)=="length")   call scatter_Bamp(BSE_scatt)
         if(trim(global_gauge)=="velocity") call scatter_Bamp_using_the_gradient(BSE_scatt,"o")
         DIP_projected(1)=-conjg(BSE_scatt%rhotw(1))
       endif
     endif
     !
     if (l_BS_dichroism) then
       if(iq==1) then   ! EM: "iq==1":  kpt=Gamma (molecules)
         DIP_expanded=DIP_rotated(ic,iv,ik_bz,i_sp_pol_c,"DIP_L",Xk)
         DIP_projected_dich(1)=dot_product(field_dir(:,1),DIP_expanded)
         if(BS_dip_size>1) DIP_projected_dich(2)=dot_product(field_dir(:,2),DIP_expanded)
         if(BS_dip_size>2) DIP_projected_dich(3)=dot_product(field_dir(:,3),DIP_expanded)
         if (trim(global_gauge)=='length'  ) DIP_projected_dich=DIP_projected_dich/q0_def_norm
       else
         call error(' finite q not implemented for dichroism')
       endif
     endif
     !
     ! Magnons
     !
     if (l_BS_magnons) then
       !
       ! notice:
       ! In optics q=0 requires a special treatment due to the non analytic behaviour
       ! This is not needed for magnons.
       ! However to follow the same structure in K_IP.F fow WFs load / DIPOLE_IO here I keep it
       ! The implementation with scatter can be used to double check the spin dipoles
       !
       if(iq==1) then
         DIP_expanded   = DIP_rotated(ic,iv,ik_bz,i_sp_pol_c,"DIP_S",Xk)
       else
         call scatter_Bamp_spin(BSE_scatt,'x')
         DIP_expanded(1)=-conjg(BSE_scatt%rhotw(1))
         call scatter_Bamp_spin(BSE_scatt,'y')
         DIP_expanded(2)=-conjg(BSE_scatt%rhotw(1))
         !The following is not used at present
         call scatter_Bamp_spin(BSE_scatt,'z')
         DIP_expanded(3)=-conjg(BSE_scatt%rhotw(1))
       endif
       !
       ! DIP_Smins(c,v,k) is  <v|S-|c>(i_sp_pol_c), to get <c|S-|v>(i_sp_pol_v) = conjg(<v|S+|c>(i_sp_pol_c))
       ! DIP_Splus(c,v,k) is  <v|S+|c>(i_sp_pol_c), to get <c|S+|v>(i_sp_pol_v) = conjg(<v|S-|c>(i_sp_pol_c))
       !
       DIP_Splus = 0.5_SP*(DIP_expanded(1)+cI*DIP_expanded(2))  ! S+ for c\dn> transitions: <v\up|c\dn>
       DIP_Smins = 0.5_SP*(DIP_expanded(1)-cI*DIP_expanded(2))  ! S- for c\up> transitions: <v\dn|c\up>
       !
     endif
     !
     ! c.2 Store Dipoles
     !====================
     !
     ! Absorption, KERR, Dichroism, PL
     if(l_BS_optics) then
       BS_T_grp(i_T_g)%dipoles_opt(:BS_dip_size,i_T,1)=conjg(DIP_projected(:BS_dip_size))
       if (BS_K_coupling.and..not.BS_res_ares_dim==2) &
       &  BS_T_grp(i_T_g)%dipoles_opt(:BS_dip_size,i_T,2)=DIP_projected(:BS_dip_size)
     endif
     !
     ! Dichroism
     if(l_BS_dichroism) then
       BS_T_grp(i_T_g)%dipoles_dic(:BS_dip_size,i_T,1)=conjg(DIP_projected_dich(:BS_dip_size))
       if (BS_K_coupling.and..not.BS_res_ares_dim==2) &
       & BS_T_grp(i_T_g)%dipoles_dic(:BS_dip_size,i_T,2)=DIP_projected_dich(:BS_dip_size)
     endif
     !
     ! Magnons
     if (l_BS_magnons) then                                      ! This is the right R residual associated to (eps_c\dn-eps_v\up)
         BS_T_grp(i_T_g)%dipoles_mag(1,i_T,1)=conjg(DIP_Splus)   ! S- for |c\dn> trans: <c\dn|S-|v\up>(i_sp_pol_v) = conjg(<v\up|S+|c\dn>(i_sp_pol_c))
         BS_T_grp(i_T_g)%dipoles_mag(2,i_T,1)=conjg(DIP_Smins)   ! S+ for |c\up> trans: <c\dn|S+|v\up>(i_sp_pol_v) = conjg(<v\dn|S-|c\up>(i_sp_pol_c))
                                                                 !      or the right R residual associated to (eps_c\dn-eps_v\up)
                                                                 ! S- for |v\dn> trans: <v\dn|S-|c\up>(i_sp_pol_c) = conjg(<c\up|S+|v\dn>(i_sp_pol_v))
                                                                 ! S+ for |v\up> trans: <v\up|S+|c\dn>(i_sp_pol_c) = conjg(<c\dn|S-|v\up>(i_sp_pol_v))
         if (BS_K_coupling.and..not.BS_res_ares_dim==2) then     ! This is the right A residual associated to (eps_v\up-eps_c\dn)
           BS_T_grp(i_T_g)%dipoles_mag(1,i_T,2)=DIP_Smins        ! S- for |v\dn> trans: <v\dn|S-|c\up>(i_sp_pol_c)
           BS_T_grp(i_T_g)%dipoles_mag(2,i_T,2)=DIP_Splus        ! S+ for |v\up> trans: <v\up|S+|c\dn>(i_sp_pol_c)
         endif
     endif
     !
     ! c.3 Residuals from "Dipoles"
     !==============================
     !
     Z_=cONE
     E_=cONE
     f_=BS_T_grp(i_T_g)%f_RES(i_T)
     if (allocated(BS_T_grp(i_T_g)%Z))                                      Z_=BS_T_grp(i_T_g)%Z(i_T)
     if (allocated(BS_T_grp(i_T_g)%W) .and. trim(global_gauge)=="velocity") E_=E_plus_W0/conjg(E_plus_W0)
     !
     ! Residuals DIP_cv * DIP_vc
     !
     if(l_BS_abs) then
       res_abs=BS_T_grp(i_T_g)%dipoles_opt(1,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_opt(1,i_T,1))
       if(l_BS_trace) then
         res_abs=res_abs+BS_T_grp(i_T_g)%dipoles_opt(2,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_opt(2,i_T,1))
         res_abs=res_abs+BS_T_grp(i_T_g)%dipoles_opt(3,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_opt(3,i_T,1))
         res_abs=res_abs/3._SP
       endif
       res_abs=res_abs*f_*Z_*E_
       if(l_BS_esort.and.i_res_ares==1) then
         BS_E_sorted_tmp(i_E_sort,nVAR+1:nVAR+1)=res_abs
         nVAR=nVAR+1
       endif
     endif
     !
     if (l_BS_kerr) then
       res_kerr=BS_T_grp(i_T_g)%dipoles_opt(1,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_opt(2,i_T,1))
       res_kerr=res_kerr*Z_*E_*f_
       if(l_BS_esort.and.i_res_ares==1) then
         BS_E_sorted_tmp(i_E_sort,nVAR+1:nVAR+1)=res_kerr
         nVAR=nVAR+1
       endif
     endif
     !
     if (l_BS_photolum) then
       !
       ! Note that res_PL is always positive as f_PL>0
       !
       f_PL=BS_T_grp(i_T_g)%f_PL(i_T,1)
       P_weighted(:)=BS_T_grp(i_T_g)%dipoles_opt(:,i_T,1)*PL_weights(:)
       res_PL=real(dot_product(P_weighted,P_weighted)*f_PL*Z_*E_)
     endif
     !
     ! Dichroism
     !
     if(l_BS_dichroism) then
       res_dich=BS_T_grp(i_T_g)%dipoles_opt(1,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_dic(1,i_T,1)) ! EM:  x*L = DIP_x*DIP_orbital;
       if(l_BS_trace) then
         res_dich=res_dich+BS_T_grp(i_T_g)%dipoles_opt(2,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_dic(2,i_T,1))
         res_dich=res_dich+BS_T_grp(i_T_g)%dipoles_opt(3,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_dic(3,i_T,1))
         res_dich=res_dich/3._SP
       endif
       res_dich=res_dich*f_*Z_/(q0_def_norm)                    !       f_ = occupations
       if(l_BS_esort.and.i_res_ares==1) then
         BS_E_sorted_tmp(i_E_sort,nVAR+1:nVAR+1)=res_dich
         nVAR=nVAR+1
       endif
     endif
     !
     ! Magnons
     !
     ! Res (RES): (Smp) <v|S-|c>(i_sp_pol_c)<c|S+|v>(i_sp_pol_c) = dip_R_magn(1)*conjg(dip_R_magn(1))
     !            (Spm) <v|S+|c>(i_sp_pol_c)<c|S-|v>(i_sp_pol_c) = dip_R_magn(2)*conjg(dip_R_magn(2))
     !
     if (l_BS_magnons) then
       res_magn(:)=BS_T_grp(i_T_g)%dipoles_mag(:,i_T,1)*conjg(BS_T_grp(i_T_g)%dipoles_mag(:,i_T,1))
       res_magn=res_magn*f_*Z_
       if(l_BS_esort.and.i_res_ares==1) then
         BS_E_sorted_tmp(i_E_sort,nVAR+1:nVAR+2)=res_magn
         nVAR=nVAR+2
       endif
     endif
     !
     if ( abs(real(E_plus_W0)) < Dip%Energy_treshold .and. iq==1 .and. dip_r ) then
       res_abs =0._SP
       res_kerr=0._SP
       res_PL  =0._SP
       res_dich=0._SP
     endif
     !
     ! d.Response functions and pert SOC loop     
     !=========================================
     !
     do i_pert_SOC=1,n_SOC
       !
       if (l_BS_jdos) then
         Joint_DOS(:,3) =Joint_DOS(:,3) -1._SP/(W_bss%p(:)-E_plus_W(i_pert_SOC))
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) &
         &  Joint_DOS(:,3) =Joint_DOS(:,3) +    1._SP/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC)))
       endif
       !
       if (l_BS_abs) then
         Epsilon_ii(:,3)=Epsilon_ii(:,3)-res_abs/(W_bss%p(:)-E_plus_W(i_pert_SOC))
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) &
         &  Epsilon_ii(:,3)=Epsilon_ii(:,3)+conjg(res_abs)/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC)))
       endif
       !
       if(abs(E_plus_W(i_pert_SOC))>0._SP) then
         if(trim(global_gauge)=='length')   para_term_w0=para_term_w0+res_abs*E_plus_W(i_pert_SOC)/bare_qpg(iq,1)**2
         if(trim(global_gauge)=='velocity') para_term_w0=para_term_w0+res_abs/E_plus_W(i_pert_SOC)
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) then
           if(trim(global_gauge)=='length')   para_term_w0=para_term_w0+conjg(res_abs*E_plus_W(i_pert_SOC))/bare_qpg(iq,1)**2
           if(trim(global_gauge)=='velocity') para_term_w0=para_term_w0+conjg(res_abs/E_plus_W(i_pert_SOC))
         endif
       endif
       !
       if (l_BS_kerr) then
         !
         Epsilon_ij(:,3)=Epsilon_ij(:,3)-  res_kerr/(W_bss%p(:)-E_plus_W(i_pert_SOC))
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) then
           Epsilon_ij(:,3)=Epsilon_ij(:,3)+conjg(res_kerr)/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC)))
         endif
         !
         ! Anomalous Hal term, Eq.(1.12) PRB 48, 11705 (1993)
         if(trim(global_gauge)=='length')   B_Hall(1)=B_Hall(1)+res_kerr/   bare_qpg(1,1)    **2
         if(trim(global_gauge)=='velocity') B_Hall(1)=B_Hall(1)+res_kerr/E_plus_W(i_pert_SOC)**2
         !
       endif
       !
       if (l_BS_photolum) then
         PL(:,3)=PL(:,3) +PL_prefactor*res_PL*aimag(-1._SP/(W_bss%p(:)-E_plus_W(i_pert_SOC)))/pi
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) then
           PL(:,3)=PL(:,3)-PL_prefactor*res_PL*aimag(-1._SP/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC))))/pi
         endif
       endif
       !
       if (l_BS_dichroism) then
         X_dichroism(:,3)=X_dichroism(:,3)-res_dich/(W_bss%p(:)-E_plus_W(i_pert_SOC))  ! EM: w(:) = omega;  
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) then
           X_dichroism(:,3)=X_dichroism(:,3)+res_dich/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC))) ! DS to check congj factor
         endif
       endif
       !
       if (l_BS_magnons) then
         X_magnons(:,1,3)=X_magnons(:,1,3)-res_magn(1)/(W_bss%p(:)-E_plus_W(i_pert_SOC))
         X_magnons(:,2,3)=X_magnons(:,2,3)-res_magn(2)/(W_bss%p(:)-E_plus_W(i_pert_SOC))
         if ((BSS_add_antiresonant.or.BS_K_coupling).and..not.BS_res_ares_dim==2) then
           X_magnons(:,1,3)=X_magnons(:,1,3)+conjg(res_magn(2))/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC)))
           X_magnons(:,2,3)=X_magnons(:,2,3)+conjg(res_magn(1))/(W_bss%p(:)+conjg(E_plus_W(i_pert_SOC)))
         endif
       endif
       !
     enddo
     !
   enddo
   !
   if (l_rpa_IP.or.iq>1) call live_timing(steps=1)
   !
 enddo
 !
 if (l_rpa_IP.or.iq>1) call live_timing()
 !
 ! ALL 2 ALL
 !==============
 !
 if (l_BS_abs) then
   call PP_redux_wait( Epsilon_ii(:,3) ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( Epsilon_ii(:,3) ,COMM=PAR_COM_eh_INDEX%COMM )
   !
   call PP_redux_wait( para_term_w0 ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( para_term_w0 ,COMM=PAR_COM_eh_INDEX%COMM )
 endif
 !
 if (l_BS_jdos) then
   call PP_redux_wait( Joint_DOS(:,3),COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( Joint_DOS(:,3),COMM=PAR_COM_eh_INDEX%COMM )
 endif
 !
 if(l_BS_esort) then
   call PP_redux_wait( BS_E_sorted_tmp ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( BS_E_sorted_tmp ,COMM=PAR_COM_eh_INDEX%COMM )
   YAMBO_ALLOC(vtmpR,(EDIM)); YAMBO_ALLOC(sort_indx,(EDIM))
   vtmpR=real(BS_E_sorted_tmp(:,1),SP)
   call sort(arrin=vtmpR,indx=sort_indx)
   EDIM=min(BS_K_dim(1),10000)
   allocate(first_el(EDIM)); allocate(n_of_el(EDIM))
   ! Threshold is 1meV
   call degeneration_finder(vtmpR(:EDIM),EDIM,first_el,n_of_el,n_deg_grp,0.001_SP/HA2EV,Include_single_values=.TRUE.)
   EDIM=min(BS_K_dim(1),1000)
   do i_VAR=1,nVAR
     BS_E_sorted(:,i_VAR,1)=cZERO
     do i_grp=1,n_deg_grp
       if(i_grp>EDIM) cycle
       if(i_VAR==2) then
         BS_E_sorted(i_grp,i_VAR,1)=real(n_of_el(i_grp),SP)
         cycle
       endif
       do i_cv=first_el(i_grp),first_el(i_grp)+n_of_el(i_grp)-1
         BS_E_sorted(i_grp,i_VAR,1)=BS_E_sorted(i_grp,i_VAR,1)+BS_E_sorted_tmp(sort_indx(i_cv),i_VAR)
       enddo
       BS_E_sorted(i_grp,i_VAR,1)=BS_E_sorted(i_grp,i_VAR,1)/real(n_of_el(i_grp),SP)
     enddo
   enddo
   deallocate(first_el);   deallocate(n_of_el)
   YAMBO_FREE(vtmpR);      YAMBO_FREE(sort_indx)
   YAMBO_FREE(BS_E_sorted_tmp)
 endif
 !
 if(l_BS_kerr) then
   call PP_redux_wait( Epsilon_ij(:,3) ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( Epsilon_ij(:,3) ,COMM=PAR_COM_eh_INDEX%COMM )
   !
   call PP_redux_wait( B_Hall(1) ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( B_Hall(1) ,COMM=PAR_COM_eh_INDEX%COMM )
 endif
 !
 if (l_BS_photolum) then
   call PP_redux_wait(PL(:,3)        ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait(PL(:,3)        ,COMM=PAR_COM_eh_INDEX%COMM )
 endif
 !
 if(l_BS_dichroism) then
   call PP_redux_wait( X_dichroism(:,3) ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( X_dichroism(:,3) ,COMM=PAR_COM_eh_INDEX%COMM )
 endif
 !
 if(l_BS_magnons) then
   call PP_redux_wait( X_magnons(:,:,3) ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait( X_magnons(:,:,3) ,COMM=PAR_COM_eh_INDEX%COMM )
 endif
 !
 if (trim(L_kind)=="full") then
   do i_T_g=1,BS_nT_grps
     if (.not.allocated(BS_T_grp(i_T_g)%dipoles_opt) ) then
       YAMBO_ALLOC(BS_T_grp(i_T_g)%dipoles_opt,(BS_dip_size,BS_T_grp(i_T_g)%size,2/BS_res_ares_dim))
       BS_T_grp(i_T_g)%dipoles_opt=cZERO
     endif
     call PP_redux_wait(BS_T_grp(i_T_g)%dipoles_opt, COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BS_T_grp(i_T_g)%dipoles_opt, COMM=PAR_COM_eh_INDEX%COMM )
   enddo
 endif
 !
 ! Factors and more
 !====================
 !
 ! DOS
 !
 if (l_BS_jdos) Joint_DOS(:,3)=Joint_DOS(:,3)*Co_factor !/(4._SP*pi) This factor is needed for the correct Jdos
 !                                                          ! However it is commented for now since all test-suite
 !                                                          ! references needs to be updated before insterting it
 !
 if (l_BS_abs) then
   !
   if (trim(global_gauge)=='velocity') then
     if (     skip_cond_sum_rule) diam_term_used=diam_term_exact
     if (.not.skip_cond_sum_rule) diam_term_used=para_term_w0
   endif
   !
   if (trim(global_gauge)=='length'  ) Epsilon_ii(:,3)=Co_factor* Epsilon_ii(:,3)                / bare_qpg(iq,1)**2
   if (trim(global_gauge)=='velocity') Epsilon_ii(:,3)=Co_factor*(Epsilon_ii(:,3)-diam_term_used)/W_bss%p(:)**2
   !
   ! Drude term
   !
   if (l_drude           ) X_drude_term(:)=       -drude_GreenF(:)       *4._SP*pi/ bare_qpg(iq,1)**2
   if (skip_cond_sum_rule) X_drude_term(:)=(para_term_w0-diam_term_exact)*Co_factor/W_bss%p(:)**2
   !
   Epsilon_ii(:,3)=Epsilon_ii(:,3)+X_drude_term(:)
   !
 endif
 !
 if (l_BS_abs) then
   !
   ! Eq.(1.12) PRB 48, 11705 (1993)
   !
   A_drude(1)=para_term_w0-diam_term_exact
   !
 endif
 !
 if (l_BS_kerr) then
   !
   B_Hall(1)=B_Hall(1)*Co_factor
   if (BSS_add_antiresonant.or.BS_K_coupling) B_Hall(1)=B_Hall(1)-conjg(B_Hall(1))
   !
   if(trim(global_gauge)=='length'   ) Epsilon_ij(:,3)= Epsilon_ij(:,3)*Co_factor/bare_qpg(iq,1)**2
   if(trim(global_gauge)=='velocity' ) Epsilon_ij(:,3)= Epsilon_ij(:,3)*Co_factor/ W_bss%p(:)   **2
   !
   ! I add the term describing the Anomalous Hall effect which is
   ! missing in the length gauge (this should be zero for systems with a gap)
   if (l_BS_anomalous_Hall.and.trim(global_gauge)=='length') Epsilon_ij(:,3)=Epsilon_ij(:,3)+B_Hall(1)/W_bss%p(:)
   !
 endif
 !
 if (l_BS_dichroism) then
  !  beta_ij=-((i*hbar*me)/(omega*qe^2))*G_ij;
  !  DeltaN=((8*pi*N*omega/(3*c))*Tr(beta_ij);
  !  N = molecular density = ??
  !  SPEED_OF_LIGHT =  c in a.u.
  X_dichroism(:,3)=X_dichroism(:,3)*(-1/SPEED_OF_LIGHT)*Co_factor/(4._SP*pi)   !!  this is (omega/c)*beta, with beta=-((i*hbar*me)/(omega*qe^2))*G_ij
 endif
 !
 if (l_BS_magnons) X_magnons(:,:,3)=X_magnons(:,:,3)*Co_factor/(4._SP*pi)
 !
 ! CLEAN
 !
 if (iq==1.and.io_DIP) then
   call DIP_alloc('DIP_iR')
   call DIP_alloc('DIP_P')
   call DIP_alloc('DIP_v')
   call DIP_alloc('DIP_S')
   call DIP_alloc('DIP_spin')
   call DIP_alloc('DIP_orbital')
   call DIP_alloc('P_square') 
   call DIP_alloc('DIP_P_spinor')
 endif
 !
 if(iq>1) then
   call elemental_collision_free(BSE_scatt)
   call WF_free(WF)
 endif
 !
 if (l_rpa_IP) then
   !
   ! Initialize & write the output file
   !
   call K_output_file(iq,"open IP")
   call K_output_file(iq,"IP")
   !
 endif
 !
 call timing('T_space IP',OPR='stop')
 !
end subroutine
