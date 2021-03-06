!
!        Copyright (C) 2000-2020 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MG, DS, DV
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
subroutine XC_libxc_driver(en,Xk,FUNCTIONAL,ORDER,EXX_FRACTION,EXX_SCREENING)
 !
 ! After libxc is distributed together with Yambo (rev 1180)  
 ! this routine evaluates the xc functional in yambo 
 ! (03/05/12: LDA OK, GGA OK => EXTENSIVE TESTS NEEDED) 
 ! (21/11/13: for HYB_GGA the GGA part calculated => TO TEST) 
 !
 use pars,           ONLY:SP,DP,pi,RZERO
 use R_lattice,      ONLY:bz_samp
 use D_lattice,      ONLY:DL_vol,alat
#if defined _RT
 use drivers,        ONLY:l_real_time,l_nl_optics
 use hamiltonian,    ONLY:rho_n,magn_n
#endif
 use electrons,      ONLY:levels,n_spin,n_spinor
 use FFT_m,          ONLY:fft_size
 use wrapper,        ONLY:V_dot_V
 use xc_functionals, ONLY:E_xc,E_xc_val,magn,F_xc,V_xc,xc_setup_hyb,&
&                         XC_LDA_C_KP,xc_fnctl_t,xc_setup_fnctl,    &
&                         xc_spin,xc_spin_map1, &
&                         XC_HYB_GGA_XC_GAUPBE,XC_FACTOR
 use xc_f90_lib_m,   ONLY:XC_FAMILY_LDA,XC_FAMILY_GGA,XC_FAMILY_HYB_GGA
 use libxc_funcs_m,  ONLY:XC_GGA_X_PBE, XC_GGA_C_PBE
 use interfaces,     ONLY:el_density_and_current
 use pseudo,         ONLY:pp_rho_nlcc,pp_has_nlcc,PP_nlcc_free
 !
#include<memory.h>
 !
 type(levels) ::en
 type(bz_samp)::Xk
 integer      ::FUNCTIONAL,FUNCTIONALTMP,ORDER
 real(SP), optional :: EXX_FRACTION, EXX_SCREENING
 !
 ! Work Space
 !
 integer  :: i1,i2,is_xc,ic,is1,is2,ixc,ifft
 real(SP) :: mod_magn
 real(SP) :: local_exx_fraction 
 real(DP) :: exsr
 !
 ! XC routines arrays
 !
 real(SP) :: spin_rho_SP(fft_size,n_spin),rho(fft_size)
 real(DP) :: spin_rho_DP(n_spin,fft_size),rspts(fft_size)
 real(SP),allocatable:: drho(:,:,:)
 real(DP),allocatable:: sigma(:,:)
 complex(SP), allocatable:: drho_cmplx(:,:)
 complex(SP), allocatable:: spin_rho_SP_cmplx(:)
 !
 ! note that multidimensional fxc & vxc is allowed only in xcspol
 !
 real(DP),allocatable :: tmp_exc(:)
 real(DP),allocatable :: v1rho(:,:),v1sigma(:,:)
 real(DP),allocatable :: vx1(:,:),vx2(:,:)
 real(DP),allocatable :: v2rho2(:,:),v2rhosigma(:,:),v2sigma2(:,:)
 logical              :: EvaluateDGradient
 !
 ! Functional and Functional infos
 !
 type(xc_fnctl_t) :: fnctl(2) 
 !
 !=============================== 
 ! Convert and Initialize
 !===============================
 !
 exsr=0._DP
 FUNCTIONALTMP=0
!Set the functional as a PBE to get Exc, next calculate the GAU-PBE Exc as
!Exc=Exc(PBE)-alfa*Ex(PBE,SR) next add the HF part.
 !
 if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) then 
   FUNCTIONALTMP=XC_HYB_GGA_XC_GAUPBE*XC_FACTOR
   FUNCTIONAL=XC_GGA_X_PBE*XC_FACTOR+XC_GGA_C_PBE  !101130   This is semi-local PBE  
 endif
   call xc_setup_fnctl(fnctl,FUNCTIONAL,n_spin)
   if (FUNCTIONALTMP==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) FUNCTIONAL=XC_HYB_GGA_XC_GAUPBE*XC_FACTOR
 !
   local_exx_fraction = RZERO
   if (present(EXX_FRACTION).and.fnctl(1)%family==XC_FAMILY_HYB_GGA) then 
     call xc_setup_hyb(fnctl,EXX_FRACTION,EXX_SCREENING)
     local_exx_fraction = fnctl(1)%exx_fraction
   end if
   if (FUNCTIONALTMP==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) local_exx_fraction=EXX_FRACTION
 !
   EvaluateDGradient = any((fnctl(:)%family==XC_FAMILY_GGA).or.any(fnctl(:)%family==XC_FAMILY_HYB_GGA))
 !
 !
 !==========================
 ! Electronic/spin density 
 !==========================
 !
#if defined _RT
 if (l_real_time.or.l_nl_optics) then
   rho=rho_n
   if(n_spin> 1) magn=magn_n
 else
#endif
   call el_density_and_current(en,Xk,rho=rho)
   if (n_spin> 1) call el_magnetization(en,Xk,magn)
#if defined _RT
 endif
#endif
 !
 ! non-linear core corrections (NLCC)
 ! to be improved for yambo_rt, otherwise this is done at each time step
 !
 call XC_nlcc_setup()
 if (pp_has_nlcc) rho(:)=rho(:)+real(pp_rho_nlcc(:),SP)
 call PP_nlcc_free()
 !
 if(n_spin==1) spin_rho_SP(:,1)=rho(:)
 if(n_spin> 1) then
   if (n_spinor==1) then
!$omp parallel do default(shared), private(ifft)
     do ifft=1,fft_size
       spin_rho_SP(ifft,1)=0.5*(rho(ifft)+magn(ifft,3))
       spin_rho_SP(ifft,2)=0.5*(rho(ifft)-magn(ifft,3))
     enddo
!$omp end parallel do
   else
!$omp parallel do default(shared), private(ifft,mod_magn)
     do ifft=1,fft_size
       mod_magn=sqrt(magn(ifft,1)**2+magn(ifft,2)**2+magn(ifft,3)**2)
       spin_rho_SP(ifft,1)=0.5*(rho(ifft)+mod_magn)
       spin_rho_SP(ifft,2)=0.5*(rho(ifft)-mod_magn)
     enddo
!$omp end parallel do
   endif
 endif
 !
 forall(i1=1:fft_size,i2=1:n_spin)  spin_rho_DP(i2,i1)=spin_rho_SP(i1,i2)*real(fft_size,DP)/real(DL_vol,DP)
 !
 !=============================
 ! Electronic density gradient
 !=============================
 if (ORDER==0) then
   YAMBO_ALLOC(tmp_exc,(fft_size))
   tmp_exc=0._SP
 endif
 if (ORDER==1.or.ORDER==2) then
   YAMBO_ALLOC(v1rho,(n_spin,fft_size))
   v1rho=0._DP
 endif
 if (ORDER==2) then
   YAMBO_ALLOC(v2rho2,(xc_spin(1),fft_size))
   v2rho2=0._DP
 endif
 !
 if (EvaluateDGradient) then 
   !
   YAMBO_ALLOC(drho,(fft_size,n_spin,3))
   YAMBO_ALLOC(sigma,(xc_spin(1),fft_size))
   YAMBO_ALLOC(drho_cmplx,(fft_size,3))
   YAMBO_ALLOC(spin_rho_SP_cmplx,(fft_size))
   if (ORDER==1.or.ORDER==2) then
     YAMBO_ALLOC(v1sigma,(xc_spin(1),fft_size))
     v1sigma=0._DP
   endif
   if (ORDER==2) then
     YAMBO_ALLOC(v2rhosigma,(xc_spin(1)*n_spin,fft_size))
     YAMBO_ALLOC(v2sigma2,(xc_spin(1)*n_spin,fft_size))
     v2rhosigma=0._DP
     v2sigma2=0._DP
   endif
   if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) then
     YAMBO_ALLOC(vx1,(n_spin,fft_size))
     YAMBO_ALLOC(vx2,(xc_spin(1),fft_size))
     vx1=0._DP
     vx2=0._DP
   endif 
   !
   do is1 = 1,n_spin
     spin_rho_SP_cmplx(:)=cmplx(spin_rho_SP(:,is1),0.0_SP,kind=SP)
     !
     call eval_Gradient(spin_rho_SP_cmplx,drho_cmplx,1,"density")
     !
     drho(:,is1,:)=real(drho_cmplx)
     forall (i1=1:fft_size) drho(i1,is1,:) = drho(i1,is1,:)*fft_size/DL_vol*2._SP*pi/alat(:)
   end do
   !
   do i1 = 1,fft_size
     do is1 = 1,n_spin
       do is2 = is1,n_spin
         is_xc = xc_spin_map1(is1,is2)
         sigma(is_xc,i1) = real(V_dot_V(3,drho(i1,is1,:),drho(i1,is2,:)),DP)
       end do
     end do
   end do
   !
   YAMBO_FREE(spin_rho_SP_cmplx)
   YAMBO_FREE(drho_cmplx)
   !
 endif
 !
 !===========================
 ! Evaluate the xc functional
 !===========================
 !
 if (order==0) E_xc=0._SP
 if (order==1) V_xc=0._SP
 if (order==2) F_xc=0._SP
 !
 if (FUNCTIONAL == XC_LDA_C_KP) then ! Special case, XC_LDA_C_KP not available from libxc
   select case(ORDER)
     case(0)
       forall (i1=1:fft_size) rspts(i1)=(3._DP/4._DP/pi/spin_rho_SP(i1,1)/DL_vol)**(1._DP/3._DP)
       call xc_rpa_kp(rspts,v1rho)
       E_xc = real(tmp_exc,SP) 
       E_xc_val = V_dot_V(fft_size,rho,E_xc(:))
     case(1)
       call error(' Vxc not coded for LDA C KP')
     case(2)
       call error(' Fxc not coded for LDA C KP')
   end select
   return
 end if
 !
 do ixc = 1,2
   if (fnctl(ixc)%id == 0) cycle
   if (fnctl(ixc)%family==XC_FAMILY_LDA) then 
     select case(ORDER)
       case(0)
         call xc_f90_lda_exc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), tmp_exc(1))
         E_xc = E_xc + real(tmp_exc,SP)    
       case(1)
         call xc_f90_lda_vxc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), v1rho(1,1))
         forall(i1=1:fft_size,i2=1:n_spin) V_xc(i1,i2) = V_xc(i1,i2) + real(v1rho(i2,i1),SP)
       case(2)
         call xc_f90_lda_vxc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), v1rho(1,1))
         call xc_f90_lda_fxc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), v2rho2(1,1))
         ! Yambo internal for Kxc
         call XC_eval_lda_kernel(v1rho,v2rho2)
     end select
   elseif (fnctl(ixc)%family==XC_FAMILY_GGA.or.fnctl(ixc)%family==XC_FAMILY_HYB_GGA) then
        !GAUPBE short range x-energy  and potentials part
        if(FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) call XC_eval_pbegaux(spin_rho_DP(1,1), sigma(1,1), exsr,vx1(1,1),vx2(1,1))
     select case(ORDER)
       case(0)
         call xc_f90_gga_exc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), sigma(1,1), tmp_exc(1))
         E_xc = E_xc + real(tmp_exc,SP)    
       case(1)
         call xc_f90_gga_vxc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), sigma(1,1), v1rho(1,1), v1sigma(1,1))
         if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR.and.ixc==1) then
           !GAUPBE short energy part of the x potential,factor 2 in vx2 for
           !consistency qith QE
           v1rho=v1rho-EXX_FRACTION*vx1
           v1sigma=v1sigma-(EXX_FRACTION*vx2)/2._DP
         endif
         ! Yambo internal for V_xc
         call XC_eval_gga_potential(v1rho,v1sigma,drho)
       case(2)
         if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) call error(' Fxc not coded for GAU-PBEP')
         call xc_f90_gga_vxc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), sigma(1,1), v1rho(1,1), v1sigma(1,1))
         call xc_f90_gga_fxc(fnctl(ixc)%conf, fft_size, spin_rho_DP(1,1), sigma(1,1), v2rho2(1,1), v2rhosigma(1,1), v2sigma2(1,1))
         ! Yambo internal for Kxc
         call XC_eval_gga_kernel(v1rho,v1sigma,v2rho2,v2rhosigma,v2sigma2)
     end select
   end if
 enddo
 !
 if (order==0) then 
   E_xc_val = V_dot_V(fft_size,rho,E_xc(:))
   if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR)  E_xc_val= E_xc_val-EXX_FRACTION*exsr
 endif
 !
 if (EvaluateDGradient) then
   YAMBO_FREE(sigma)
   YAMBO_FREE(drho)
 endif
 if (ORDER==1) then
   YAMBO_FREE(tmp_exc)
 endif
 if (EvaluateDGradient.and.(ORDER==1)) then
   YAMBO_FREE(v1rho)
   YAMBO_FREE(v1sigma)
 endif
 if (EvaluateDGradient.and.(ORDER==2)) then
   YAMBO_FREE(v2rho2)
   YAMBO_FREE(v2rhosigma)
   YAMBO_FREE(v2sigma2)
 endif
 if (EvaluateDGradient.and.(FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR)) then
   YAMBO_FREE(vx1)
   YAMBO_FREE(vx2)
 endif
 !
 if (present(EXX_FRACTION)) EXX_FRACTION = local_exx_fraction 
 !
end subroutine XC_libxc_driver
