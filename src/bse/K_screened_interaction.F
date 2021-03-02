!
!        Copyright (C) 2000-2020 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS
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
subroutine K_screened_interaction(X,Xw,q,io_X_err)
 !
 use pars,         ONLY:SP,pi,cZERO,cI,DP
 use stderr,       ONLY:intc
 use units,        ONLY:HBAR_eVfs,BO2ANG,HBAR_SI
 use frequency,    ONLY:w_samp
 use drivers,      ONLY:l_td_hf
 use com,          ONLY:warning,msg
 use collision_el, ONLY:elemental_collision,elemental_collision_free
 use R_lattice,    ONLY:bz_samp,RIM_qpg,RIM_anisotropy,bare_qpg
 use D_lattice,    ONLY:i_time_rev,i_space_inv,DL_vol
 use electrons,    ONLY:levels,spin_occ,spin,n_spin,nel
 use X_m,          ONLY:X_alloc,X_t,X_mat
#ifdef _CUDA
 use X_m,          ONLY:X_mat_d
#endif
 use BS,           ONLY:BS_n_g_W,BS_W_is_diagonal,BS_W,   &
& BS_K_is_HYB,BS_K_is_SXX,BS_K_is_MDL,BS_K_eMread,  &
& BS_hybrd_fracH,BS_hybrd_type,BS_hybrd_frac, BS_ALDA_is_G
 use IO_m,         ONLY:io_control,OP_RD_CL,REP,VERIFY,NONE,OP_RD,RD,RD_CL,RD_CL_IF_END
 use deviceXlib_m, ONLY:dev_memcpy
 use fft_m,        ONLY: fft_size,fft_dim,fftw_plan,fft_g_table
 !
#include<dev_defs.h>
#include<memory.h>
 !
 type(bz_samp) ::q
 type(X_t)     ::X
 type(w_samp)  ::Xw
 !
 integer,     intent(out) :: io_X_err
 !
 type(elemental_collision), target ::isc
#ifdef _CUDA
 complex(SP), pointer, device :: gamp_d(:,:)
 real(SP)   :: rtmp
 integer    :: nibz
#endif
 !
 integer    :: BS_W_size2,i_g2,i_g3,iq_W
 logical    :: W_is_em1s,W_is_pp

 !
 ! I/O
 !
 integer           :: XID
 integer, external :: io_X

 logical              ::  BS_eM_is_manual=.false.
 logical              ::  BS_hyb_is_global,BS_W_is_model=.false.
 real(SP),allocatable  :: model_epsilon(:,:)
 complex(DP),allocatable :: BS_hybrd_fracR(:)
 real(SP)             :: ele_density,q_ThmsFrm2,m2wp2,BS_K_eM

 ele_density = nel/DL_vol/(BO2ANG**3)  ! density of electron, unit: A^{-3}
 q_ThmsFrm2 = 4.0*(3.0*ele_density/pi)**(1.0/3.0)/BO2ANG  !Thomas-Fermi vector: q_tf^2 = 
 m2wp2 = ele_density*9.10938/8.854*1.602176462**2  ! free electron plasma freqency: wp^2 = density*e^2/e0/m_e
 ! 

! open(81,file='model_inv_eps.dat')
 if (l_td_hf) then
   !
   io_X_err=0
   !
 else
   !
   call section('+','Screened interaction header I/O')
   !===================================================
   !
   ! X%ng are the G's in the X DB while BS_n_g_W the G's I want to read.
   ! Now to read only BS_n_g_W components I need to redefine X%ng
   !

   if(BS_K_is_SXX) then
    call msg('r','SXX kernel')
!    allocate(inv_epsilon(1,1))
    if(BS_K_eMread.gt.0.99) then
      BS_eM_is_manual = .true.
      BS_K_eM = BS_K_eMread
      call msg('r','We use manual epsilon: ', BS_K_eMread)
    else
      BS_K_eMread = 0.0
      BS_K_eM = 0.0
      call msg('r','We use RPA epsilon: ', BS_K_eMread) 
    endif
   endif

   if(BS_K_is_MDL) then
    BS_W_is_model = .true.
    call msg('r','Model dielectric function kernel')
      allocate(model_epsilon(BS_n_g_W,1))
    if(BS_K_eMread.gt.0.99) then
      BS_eM_is_manual = .true.
      BS_K_eM   = BS_K_eMread
      call msg('r','We use manual epsilon: ', BS_K_eMread)
    else
      BS_K_eMread = 0.0
      BS_K_eM = 0.0
      call msg('r','We use RPA epsilon: ', BS_K_eMread) 
    endif      
   endif

   if(BS_K_is_HYB) then
!     allocate(inv_epsilon(1,1))
     if(BS_hybrd_fracH.gt.1E-4) then
      BS_eM_is_manual = .true.
!      inv_epsilon(1,1) = BS_hybrd_frac
      BS_K_eM = 1._SP/BS_hybrd_fracH
      BS_K_eMread = BS_K_eM
      call msg('r','We use manual epsilon: ', BS_K_eM)
     elseif(BS_hybrd_fracH.lt.-0.1) then
      BS_K_eMread = 0.0
      BS_K_eM = 0.0
      call msg('r','We use RPA epsilon: ', BS_K_eM) 
     endif

     select case(BS_hybrd_type)
    case(1)
 !   if(BS_n_g_W)
      call msg('r','We are doing SXX + (1-a)*ALDA ',BS_hybrd_type)
      BS_hyb_is_global=.true.
    case(2)
      call msg('r','We are doing SXX + (1-a)*ALDA_x + ALDA_c ',BS_hybrd_type)
      BS_hyb_is_global=.true.
    case(3)
      call msg('r','We are doing SXX + (1-a)*ALDA_x + ALDA_c, together ',BS_hybrd_type)
      BS_hyb_is_global=.true.
    case(4)
      call msg('r','We are doing Range-separated DDH in G-space, as well as ALDA',BS_hybrd_type)
      BS_ALDA_is_G = .true.
      BS_hyb_is_global=.false.
      BS_W_is_model=.true.
      allocate(BS_hybrd_frac(BS_n_g_W))
      allocate(model_epsilon(BS_n_g_W,1))
    case(5)
      call msg('r','We are doing Range-separated DDH in R-space, as well as ALDA',BS_hybrd_type)
      BS_ALDA_is_G = .false.
      BS_hyb_is_global=.false.
      BS_W_is_model=.true.
      allocate(BS_hybrd_frac(fft_size))
      allocate(BS_hybrd_fracR(fft_size))
 call msg('r','we are ok after alloc', fft_size)
      allocate(model_epsilon(BS_n_g_W,1))
    case(6)  ! rpa epsilon for 2D test, further 2D model
      call msg('r','We are doing Range-separated DDH in R-space, epsilon read from RPA',BS_hybrd_type)
      BS_ALDA_is_G = .false.
      BS_hyb_is_global=.false.
      BS_W_is_model=.false.
      allocate(BS_hybrd_frac(fft_size))
      allocate(BS_hybrd_fracR(fft_size))
 call msg('r','we are ok after alloc', fft_size)
!      allocate(model_epsilon(BS_n_g_W,1))
    end select

   endif



   if(BS_K_eMread .lt. 1.0) then

     if(BS_K_eMread.gt.1E-4) then
       io_X_err=-1
         call msg('r','Wrong number of eMread: should > 1.0 or =0.0: ', BS_K_eMread)
         return
     else  
     X%ng=BS_n_g_W
     call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=XID)
     io_X_err=io_X(X,Xw,XID) 
     !
     call msg('r', 'io_X_err: ', io_X_err)
     if(io_X_err<0) return

     endif
   endif
   !
   ! X%ng is set to BS_n_g_W to VERIFY if there are enough bands.
   ! From now on it is used to keep track of the actual sixe of W in the DB
   !
   X%ng=X%ng_db
   !
 endif
 !
 call elemental_collision_free(isc)
 !
 X%ng=BS_n_g_W
 if ((.not.l_td_hf).and.(BS_K_eMread.lt.1E-4)) then
   !
   ! PP/Epsm1S DBs
   !
   W_is_em1s=X%whoami==2
   W_is_pp  =X%whoami==4
   !
   if (W_is_em1s) call X_alloc('X',(/BS_n_g_W,BS_n_g_W,1/))
   if (W_is_pp)   call X_alloc('X',(/BS_n_g_W,BS_n_g_W,2/))
   !
 endif
 !
 BS_W_size2=BS_n_g_W
 if (BS_W_is_diagonal) BS_W_size2=1
 !
 ! When TR is present but not the SI X_mat indexes need to be exchanged 
 ! when the TR is applied
 !
 if (i_space_inv==0.and.i_time_rev==1.and..not.BS_W_is_diagonal) then
   YAMBO_ALLOC(BS_W,(BS_n_g_W,BS_W_size2,2*q%nibz))
 else
   YAMBO_ALLOC(BS_W,(BS_n_g_W,BS_W_size2,q%nibz))
 endif
 !
 YAMBO_ALLOC(isc%gamp,(BS_n_g_W,BS_W_size2))
#ifdef _CUDA
   YAMBO_ALLOC(isc%gamp_d,(BS_n_g_W,BS_W_size2))
#endif
 !
 if ((.not.l_td_hf).and.(BS_K_eMread.lt.1E-4)) then
   call io_control(ACTION=OP_RD,COM=NONE,SEC=(/1/),ID=XID)
   io_X_err=io_X(X,Xw,XID)
   if (io_X_err<0) then
     call warning('Reading of io_X heared failed depite previous check.')
     return
   endif
 endif
 !
 do iq_W=1,q%nibz
   !
   isc%qs(2)=iq_W
   call DEV_SUBNAME(scatter_Gamp)(isc,'c')
   !
   if (.not.l_td_hf) then
     !
    if(.not.BS_eM_is_manual) then
     !
     call io_control(ACTION=RD_CL_IF_END,COM=NONE,SEC=(/2*iq_W,2*iq_W+1/),ID=XID)
     io_X_err=io_X(X,Xw,XID)
     if (io_X_err<0) then
       call warning('Reading of io_X for iq='//trim(intc(iq_W))//' failed depite previous check.')
       return
     endif
     !
     !forall(i_g2=1:BS_n_g_W) X_mat(i_g2,i_g2,1)=X_mat(i_g2,i_g2,1)+1._SP
     if(iq_W.eq.1) then
      BS_K_eM = 1._SP/(real(X_mat(1,1,1))+1._SP)
      if(BS_K_is_HYB) BS_hybrd_fracH = 1._SP/BS_K_eM
     endif
    endif
     !
#ifdef _CUDA
     ! device <- host transfer
     X_mat_d=X_mat
     !
     ! CUDA Fortran workarounds
     gamp_d => isc%gamp_d
     nibz=q%nibz
     !
     if (BS_W_is_diagonal) then
       !$cuf kernel do(1) <<<*,*>>>
       do i_g2=1, BS_n_g_W
          BS_W(i_g2,1,iq_W)=real(X_mat_d(i_g2,i_g2,1),SP)*gamp_d(i_g2,1)
       enddo
     else if (.not.BS_W_is_diagonal) then
       !$cuf kernel do(2) <<<*,*>>>
       do i_g3=1,BS_n_g_W
         do i_g2=1,BS_n_g_W
           BS_W(i_g2,i_g3,iq_W)=X_mat_d(i_g2,i_g3,1)*gamp_d(i_g2,i_g3)
           if (i_space_inv==0.and.i_time_rev==1) BS_W(i_g2,i_g3,nibz+iq_W)=X_mat_d(i_g3,i_g2,1)*gamp_d(i_g2,i_g3)
         enddo
       enddo
     endif
#else
     if (BS_W_is_diagonal) then
      if(BS_K_is_SXX.or.(BS_K_is_HYB.and.BS_hyb_is_global)) then
       forall(i_g2=1:BS_n_g_W) BS_W(i_g2,1,iq_W)=1._SP/BS_K_eM*isc%gamp(i_g2,1)
      ! test hybrid2D
      if(BS_hybrd_type==6 .and.(.not.BS_eM_is_manual)) then
      ! W part
       forall(i_g2=1:BS_n_g_W) X_mat(i_g2,i_g2,1)=X_mat(i_g2,i_g2,1)+1._SP
       forall(i_g2=1:BS_n_g_W) BS_W(i_g2,1,iq_W)=real(X_mat(i_g2,i_g2,1),SP)*isc%gamp(i_g2,1) 
       ! fxc part
       if(iq_W==1) then

         BS_hybrd_fracR=(0._DP, 0._DP)

         do i_g2=1,BS_n_g_W
! call msg('r','fft_g_table: ',  fft_g_table(i_g2,1))
          BS_hybrd_frac(i_g2)=1._SP-real(X_mat(i_g2,i_g2,1),SP)
          BS_hybrd_fracR(fft_g_table(i_g2,1))=cmplx(real(BS_hybrd_frac(i_g2),DP), 0._DP)
         enddo
          ! do FT to R-space
#if defined _FFTW
    call dfftw_destroy_plan(fftw_plan)
    fftw_plan = 0
    call fft_3d(BS_hybrd_fracR/real(fft_size,DP),fft_dim,1,fftw_plan)
#else
    call fft_3d(BS_hybrd_fracR/real(fft_size,DP),fft_dim,1)
#endif
!
         BS_hybrd_frac(:) = real(BS_hybrd_fracR(:),SP)
         do i_g2=1,fft_size
           call msg('r','BS_hybrd_frac: ',BS_hybrd_frac(i_g2))
!            call msg('r','BS_hybrd_frac/fft_size: ',BS_hybrd_frac(i_g2)/fft_size)
         enddo
         deallocate( BS_hybrd_fracR )
       endif
      !       
      elseif(BS_W_is_model) then
       ! 
       do i_g2=1,BS_n_g_W
        !
        model_epsilon(i_g2,1) = 1._SP+1._SP/(1._SP/(BS_K_eM-1._SP)+ &
&   1.563*(real(bare_qpg(iq_W,i_g2))/BO2ANG)**2/q_ThmsFrm2   +   &
&   (real(bare_qpg(iq_W,i_g2))/BO2ANG)**4*1.054571726**2*1.0E-1/4._SP/m2wp2    )
        !
        BS_W(i_g2,1,iq_W)=1._SP/model_epsilon(i_g2,1)*isc%gamp(i_g2,1)
        !
        ! begin fractional of hybrid functional
        if(iq_W==1 .and. BS_K_is_HYB) BS_hybrd_frac(i_g2)=1._SP-1._SP/model_epsilon(i_g2,1)
       enddo ! ig < ng_W

        if(iq_W==1 .and. BS_hybrd_type==5) then
          BS_hybrd_fracR=(0._DP, 0._DP)

          do i_g2=1,BS_n_g_W
! call msg('r','fft_g_table: ',  fft_g_table(i_g2,1))
            BS_hybrd_fracR(fft_g_table(i_g2,1))=cmplx(real(BS_hybrd_frac(i_g2),DP), 0._DP)
          enddo
          !
!          do i_g2=BS_n_g_W+1,fft_size
!            BS_hybrd_fracR(fft_g_table(i_g2,1))=(1._DP, 0._DP)
!          enddo ! ig > ng_W
          !
          ! do FT to R-space
#if defined _FFTW
    call dfftw_destroy_plan(fftw_plan)
    fftw_plan = 0
    call fft_3d(BS_hybrd_fracR/real(fft_size,DP),fft_dim,1,fftw_plan)
#else
    call fft_3d(BS_hybrd_fracR/real(fft_size,DP),fft_dim,1)
#endif
!
          BS_hybrd_frac(:) = real(BS_hybrd_fracR(:),SP)
          do i_g2=1,fft_size
            call msg('r','BS_hybrd_frac: ',BS_hybrd_frac(i_g2))
!            call msg('r','BS_hybrd_frac/fft_size: ',BS_hybrd_frac(i_g2)/fft_size)
          enddo
          deallocate( BS_hybrd_fracR )
        endif
        ! end of fractional of hybrid functionals
       
      endif
     else if (.not.BS_W_is_diagonal) then
       call warning("off-diagonal W is not supported currently!!")
       return
       do i_g3=1,BS_n_g_W
         do i_g2=1,BS_n_g_W
           BS_W(i_g2,i_g3,iq_W)=X_mat(i_g2,i_g3,1)*isc%gamp(i_g2,i_g3)
           if (i_space_inv==0.and.i_time_rev==1) BS_W(i_g2,i_g3,q%nibz+iq_W)=X_mat(i_g3,i_g2,1)*isc%gamp(i_g2,i_g3)
         enddo
       enddo
     endif
#endif
     !
   else
     !
     call dev_memcpy(BS_W(:,1,iq_W),DEV_VARNAME(isc%gamp)(:,1))
     !
   endif
   ! 
 enddo
 !
 ! Anisotropy correction. Note that isc%gamp(1,1)=RIM_anisotropy while
 ! the \delta(G,G') term must be multiplied by the standard RIM_qpg(1,1)
 !
#ifdef _CUDA
 if (RIM_anisotropy/=0.) then
    rtmp=BS_W(1,1,1)
    rtmp=rtmp+RIM_qpg(1,1,1)/2._SP-RIM_anisotropy/2._SP
    BS_W(1,1,1)=rtmp
 endif
#else
 if (RIM_anisotropy/=0.) BS_W(1,1,1)=BS_W(1,1,1)+RIM_qpg(1,1,1)/2._SP-RIM_anisotropy/2._SP
#endif
 !
 call elemental_collision_free(isc)
 !
 if (.not.l_td_hf) call X_alloc('X')
 X%ng=X%ng_db
 !
end subroutine
