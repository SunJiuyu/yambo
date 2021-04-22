!
!        Copyright (C) 2000-2019 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM, DS
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
subroutine QP_hyb_xDFT(E,Xk,k,q)
 !
 ! Screened Local V_xc in hybrid functional
 !
 use pars,          ONLY:SP,pi,DP,cZERO
 use com,           ONLY:error,msg
 use electrons,     ONLY:levels,n_spin,n_spinor,spin
 use QP_m,          ONLY:QP_Vxc,QP_scr_Vxc,QP_n_states,QP_table,QP_hyb_screen,QP_ng_Sex,QP_ng_Vxc,QP_nb,QP_ng_Sx
 use R_lattice,     ONLY:bz_samp,qindx_S
 use FFT_m,         ONLY:fft_size,fft_dim,fftw_plan,fft_g_table
 use wave_func,     ONLY:WF
 use interfaces,    ONLY:PARALLEL_global_indexes,PARALLEL_WF_distribute,PARALLEL_WF_index,WF_load,QP_state_extract_print 
 use xc_functionals,ONLY:E_xc,V_xc,xc_string,XC_potential_driver,magn,XC_potential_driver,Vxc_separate ,Vxc_x, Vxc_c ! Jiuyu separately count x and c
 use global_XC,     ONLY:WF_kind,WF_xc_functional,WF_exx_fraction,WF_exx_screening
 use wrapper_omp,   ONLY:Vstar_dot_V_omp
 use parallel_m,    ONLY:PP_redux_wait,PAR_IND_WF_linear,PAR_IND_Xk_ibz,PAR_IND_QP,&
&                        PAR_IND_QP_ID,PAR_IND_Q,PAR_IND_G_b,&
&                        PAR_IND_Q_ID,n_WF_bands_to_load
 use pseudo,        ONLY:pp_is_uspp,PP_nlcc_free
#if defined _RT || defined _SC
 use drivers,       ONLY:l_real_time,l_sc_run
 use hamiltonian,   ONLY:H_ref_bands
#endif
#if defined _RT
 use drivers,       ONLY:l_real_time
 use parallel_m,    ONLY:PAR_IND_B_mat_ordered,PAR_IND_WF_k 
 use matrix_operate,ONLY:UP_matrix_index,matrix_index
#endif
 use timing_m,      ONLY:timing
 !
#include<memory.h>
 !
 type(levels) ::E       
 type(bz_samp)::Xk,k,q
 !
 ! Work Space
 !
 integer     ::i1,ifft,iwf,iwfp,ib,ibp,ik,i_sp_pol,ig
 real(SP)    ::E_xc_
 complex(SP) ::c_ws(fft_size,n_spinor)
 complex(SP),allocatable ::  xcx_ws(:,:),xcc_ws(:,:)
 complex(DP),allocatable ::  QP_hyb_screenAUX(:)
 real(SP),allocatable :: QP_hyb_screenR(:)
 !
 ! Spin support (collinear and not collinear)
 !
 complex(SP), allocatable :: V_xc_mat(:,:,:)
 !
 logical :: hyb2D
 !
 call XC_nlcc_setup()
 !
   call PARALLEL_global_indexes(E,k,q,"Self_Energy")
   !
   ! Final states (k-q,n')...
   call PARALLEL_WF_distribute(K_index=PAR_IND_Xk_ibz,B_index=PAR_IND_G_b,CLEAN_UP=.TRUE.)
   ! Initial states (k,n \in QP)...
   call PARALLEL_WF_distribute(QP_index=PAR_IND_QP)
   call PARALLEL_WF_index() 
   call WF_load(WF,QP_ng_Vxc,maxval(qindx_S(:,:,2)),(/1,E%nbm/),(/1,k%nibz/),title='-HF/Rho')   
 !
 call timing('XCo_local',OPR="start")
 !
 ! Allocation 
 !
 hyb2D = .true.
 Vxc_separate = .true.
   if(Vxc_separate) then
   allocate(Vxc_x(fft_size,n_spin))
   allocate(Vxc_c(fft_size,n_spin))
!   allocate(xcx_ws(fft_size,n_spinor))
!   allocate(xcc_ws(fft_size,n_spinor))
   endif 
 if (.not.allocated(V_xc)) then
   YAMBO_ALLOC(V_xc,(fft_size,n_spin))
 endif
 if (.not.allocated(E_xc)) then
   YAMBO_ALLOC(E_xc,(fft_size))
 endif
 if (.not.allocated(magn)) then
   YAMBO_ALLOC(magn,(fft_size,3))
 endif
 if(.not.allocated(QP_scr_Vxc))   YAMBO_ALLOC(QP_scr_Vxc,(QP_n_states))
 QP_scr_Vxc=cZERO
 !
 if(hyb2D)  then
  allocate(QP_hyb_screenR(fft_size))
 endif
 !
 if(hyb2D) then
   allocate(QP_hyb_screenAUX(fft_size))
   QP_hyb_screenR=0.0
   QP_hyb_screenAUX=(0.0,0.0)
   do ig=1,QP_ng_Sex
     QP_hyb_screenR(ig) = 1._SP-real(QP_hyb_screen(ig,1))
     QP_hyb_screenAUX(fft_g_table(ig,1))=cmplx(real(QP_hyb_screenR(ig),DP), 0._DP)
   enddo
    ! do FT to R-space
#if defined _FFTW
     call dfftw_destroy_plan(fftw_plan)
     fftw_plan = 0
     call fft_3d(QP_hyb_screenAUX/real(fft_size,DP),fft_dim,1,fftw_plan)
#else
     call fft_3d(QP_hyb_screenAUX/real(fft_size,DP),fft_dim,1)
#endif
     QP_hyb_screenR(:) = real(QP_hyb_screenAUX(:),SP)
!     do ifft=1,fft_size
!       call msg('r','QP_hyb_screenR: ',QP_hyb_screenR(ifft))
!     enddo
     deallocate( QP_hyb_screenAUX )
     deallocate( QP_hyb_screen )
 endif
 !
 ! Precompiler options in case a non-DFT potential is required
 !
     ! evaluate E_xc as a double check against DFT
     call XC_potential_driver(E,Xk,WF_KIND,WF_xc_functional,0,&
&         exx_fraction=WF_exx_fraction,exx_screening=WF_exx_screening)
     ! evaluate v_xc
     call XC_potential_driver(E,Xk,WF_KIND,WF_xc_functional,1,&
&                             exx_fraction=WF_exx_fraction,exx_screening=WF_exx_screening)
   !
 !
 if (n_spin==1) then
   ! 
   do i1=1,QP_n_states
     !
     ib =QP_table(i1,1)
     ibp=QP_table(i1,2)
     ik =QP_table(i1,3)
     !
       !
       if (.not.PAR_IND_WF_linear%element_2D(ib,ik)) cycle
       !
     ! <iwf|V_xc|iwfp>
     !
     iwf =WF%index(ib,ik,1)
     iwfp=WF%index(ibp,ik,1)
     ! 
     !
!$omp parallel do default(shared), private(ifft)
     do ifft=1,fft_size
!       if(Vxc_separate) then
!        V_xc(ifft,1) = Vxc_x(ifft,1)*QP_hyb_screenR(ifft)+Vxc_c(ifft,1)
 !        xcx_ws(ifft,1)=Vxc_x(ifft,1)*WF%c(ifft,1,iwfp)
 !        xcc_ws(ifft,1)=Vxc_c(ifft,1)*WF%c(ifft,1,iwfp)
!       endif
       c_ws(ifft,1)=(Vxc_x(ifft,1)*QP_hyb_screenR(ifft)+Vxc_c(ifft,1))*WF%c(ifft,1,iwfp)
     enddo
!$omp end parallel do
     !
     QP_scr_Vxc(i1)=Vstar_dot_V_omp(fft_size*n_spinor,WF%c(:,:,iwf),c_ws)
     !
     !
     if (pp_is_uspp) then
       call PP_vloc_augment(QP_Vxc(i1),V_xc,ik,ib,ibp)
     endif
     !
   enddo
   !
   call PP_redux_wait(QP_scr_Vxc)
   !
   call LOCAL_clean( )
   !
   call timing('XCo_local',OPR="stop")
   !
   call PP_nlcc_free()
   call WF_free(WF)
   !
   return
   !
 endif

 if(hyb2D)  deallocate(QP_hyb_screenR)
 !
#if defined _RT
 nb=H_ref_bands(2)-H_ref_bands(1)+1
#endif
 !
 if (n_spin==2) then
   !
   YAMBO_ALLOC(V_xc_mat,(fft_size,n_spin,n_spin))
   !
   call Build_V_xc_mat(V_xc,V_xc_mat)
   !
   do i1=1,QP_n_states
     !
     ib =QP_table(i1,1)
     ibp=QP_table(i1,2)
     ik =QP_table(i1,3)
     i_sp_pol=spin(QP_table(i1,:))
     !
#if defined _RT
     if (l_real_time) then
       if (.not.PAR_IND_WF_k%element_1D(ik)) cycle
       if (.not.PAR_IND_B_mat_ordered%element_1D( UP_matrix_index(ib-H_ref_bands(1)+1,ibp-H_ref_bands(1)+1)-1 )) cycle
     else
#endif
       if (.not.PAR_IND_WF_linear%element_2D(ib,ik)) cycle
#if defined _RT
     endif
#endif
     !
     if (n_spinor==1) then
       !
       iwf =WF%index(ib,ik,i_sp_pol)
       iwfp=WF%index(ibp,ik,i_sp_pol)
       !
       c_ws(:,1)=V_xc_mat(:,i_sp_pol,i_sp_pol)*WF%c(:,1,iwfp)
       !
       QP_Vxc(i1)=Vstar_dot_V_omp(fft_size*n_spinor,WF%c(:,:,iwf),c_ws)
       !
       if (pp_is_uspp) then
         call error('[PPs] USPP: nspin=2 and n_spinor=1 not implemented')
       endif
       !
     else
       !
       iwf=WF%index(ib,ik,1)
       iwfp=WF%index(ibp,ik,1)
       !
       c_ws(:,1)=V_xc_mat(:,1,1)*WF%c(:,1,iwfp)+V_xc_mat(:,1,2)*WF%c(:,2,iwfp)
       c_ws(:,2)=V_xc_mat(:,2,1)*WF%c(:,2,iwfp)+V_xc_mat(:,2,2)*WF%c(:,2,iwfp)
       !
       QP_Vxc(i1)=Vstar_dot_V_omp(fft_size*n_spinor,WF%c(:,:,iwf),c_ws)
       !
       if (pp_is_uspp) then
         call error('[PPs] USPP: nspin=2 and n_spinor>1 not implemented')
       endif
       !
     endif
     !
   enddo
   !
   YAMBO_FREE(V_xc_mat)
   !
   call LOCAL_clean( )
   !
   call PP_redux_wait(QP_Vxc)
   !
 endif
 !
 call timing('XCo_local',OPR="stop")
 !
 return
 !
 contains
   !
   subroutine LOCAL_clean()
     !
#if defined _SC || defined _RT
     if (.not.(l_sc_run.or.l_real_time)) then
#endif
       YAMBO_FREE(V_xc)
       YAMBO_FREE(E_xc)
       YAMBO_FREE(magn)
#if defined _SC || defined _RT
     endif
#endif
   end subroutine
   !
end subroutine
