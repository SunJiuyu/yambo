!
! License-Identifier: GPL
!
! Copyright (C) 2018 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine K_screened_interaction_hyb(X,Xw,q,io_X_err)
 !
 use pars,         ONLY:SP,pi,cZERO,cI
 use stderr,       ONLY:intc
 use units,        ONLY:HBAR_eVfs,BO2ANG,HBAR_SI
 use com,          ONLY:exp_user,msg
 use frequency,    ONLY:w_samp
 use drivers,      ONLY:l_td_hf,l_tddft,l_rim_w
 use collision_el, ONLY:elemental_collision,elemental_collision_free
 use R_lattice,    ONLY:bz_samp,RIM_qpg,RIM_anisotropy,RIM_W,RIM_W_d,RIM_W_is_diagonal,&
&                       RIM_W_ng
 use D_lattice,    ONLY:i_time_rev,i_space_inv
 use R_lattice,    ONLY:bz_samp,RIM_qpg,RIM_anisotropy,bare_qpg
 use D_lattice,    ONLY:i_time_rev,i_space_inv,DL_vol
 use electrons,    ONLY:levels,spin_occ,spin,n_spin,nel
 use X_m,          ONLY:X_ALLOC_elemental,X_t,X_mat,X_mat_d
 use BS,           ONLY:BS_n_g_W,BS_W_is_diagonal,BS_W, &
 & BS_K_is_HYB,BS_K_is_SXX,BS_K_is_MDL,BS_K_eMread,BS_hybrd_frac,BS_hybrd_type
 use IO_int,       ONLY:io_control
 use IO_m,         ONLY:manage_action,OP_RD_CL,REP,VERIFY,NONE,&
&                       OP_RD,RD,RD_CL,RD_CL_IF_END
 use deviceXlib_m, ONLY:dev_memcpy
 use cuda_m,       ONLY:have_cuda
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
 complex(SP), pointer DEV_ATTR :: gamp_d(:,:)
 real(SP)   :: rtmp
 integer    :: nibz
 !
 integer    :: BS_W_size2,i_g2,i_g3,i_g2_mat,i_g3_mat,iq_W
 logical    :: W_is_em1s,W_is_pp,W_is_unscreened,l_RIM_W_g
 !
 real(SP),   allocatable          :: RIM_vpW(:,:)
 real(SP),   allocatable DEV_ATTR :: RIM_vpW_d(:,:)
 !
 ! I/O
 !
 integer           :: XID,IO_ACT
 integer, external :: io_X
 !
 ! 2021/02/11 to be modified for functionals with screening
 !
 ! Jiuyu setting model BSE and parameters, to be optimized to dielectric module
 logical              ::  BS_eM_is_manual=.false.
 real(SP),allocatable  :: model_epsilon(:,:)
 real(SP)             :: ele_density,q_ThmsFrm2,m2wp2,BS_K_eM

 ele_density = nel/DL_vol/(BO2ANG**3)  ! density of electron, unit: A^{-3}
 q_ThmsFrm2 = 4.0*(3.0*ele_density/pi)**(1.0/3.0)/BO2ANG  !Thomas-Fermi vector: q_tf^2 = 
 m2wp2 = ele_density*9.10938/8.854*1.602176462**2  ! free electron plasma freqency: wp^2 = density*e^2/e0/m_e
!
 W_is_unscreened=l_td_hf !.or.l_tddft: here we comment l_tddft because it is in screen_hyb subroutine. We do not consider standard runs.
 !
 if (W_is_unscreened) then
   !
   io_X_err=0
   !
 else
   !
   if(.not.l_rim_w) call section('+','Screened interaction header I/O')
   if(     l_rim_w) call section('=','Screened interaction header I/O')
   !===================================================
   !
   if (X%ng<BS_n_g_W) then
     if (     exp_user) call warning("Using screening with less G than Weh interaction")
     if (.not.exp_user) call   error("Using screening with less G than Weh interaction")
   endif
   !
   call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=XID)
   io_X_err=io_X(X,Xw,XID) 
   !
   if(io_X_err<0) return
   !
   ! Jiuyu setting SXX and hybrid
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
     if(BS_hybrd_frac.gt.1E-4) then
      BS_eM_is_manual = .true.
!      inv_epsilon(1,1) = BS_hybrd_frac
      BS_K_eM = 1._SP/BS_hybrd_frac
      BS_K_eMread = BS_K_eM
      call msg('r','We use manual epsilon: ', BS_K_eM)
     elseif(BS_hybrd_frac.lt.-0.1) then
      BS_K_eMread = 0.0
      BS_K_eM = 0.0
      call msg('r','We use RPA epsilon: ', BS_K_eM) 
     endif

     select case(BS_hybrd_type)
    case(1)
 !   if(BS_n_g_W)
      call msg('r','We are doing SXX + (1-a)*ALDA ',BS_hybrd_type)
    case(2)
      call msg('r','We are doing SXX + (1-a)*ALDA_x + ALDA_c ',BS_hybrd_type)
    case(3)
      call msg('r','We are doing SXX + (1-a)*ALDA_x + ALDA_c, together ',BS_hybrd_type)
!    call msg('r','We are doing SXX, empirical epsilon ',BS_SXX_type)
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

 endif
 !
 call elemental_collision_free(isc)
 !
 if ((.not.W_is_unscreened).and.(BS_K_eMread.lt.1E-4)) then
   !
   ! PP/Epsm1S DBs
   !
   W_is_em1s=X%whoami==2
   W_is_pp  =X%whoami==4
   !
   if (W_is_em1s) call X_ALLOC_elemental('X',(/X%ng,X%ng,1/))
   if (W_is_pp)   call X_ALLOC_elemental('X',(/X%ng,X%ng,2/))
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
 if (have_cuda) then
   YAMBO_ALLOC(isc%gamp_d,(BS_n_g_W,BS_W_size2))
   YAMBO_ALLOC(RIM_vpW_d,(RIM_W_ng,RIM_W_ng))
 else
  YAMBO_ALLOC(RIM_vpW,(RIM_W_ng,RIM_W_ng))
 endif
 !
 if ((.not.W_is_unscreened).and.(BS_K_eMread.lt.1E-4)) then
   call io_control(ACTION=OP_RD,COM=NONE,SEC=(/1/),ID=XID)
   io_X_err=io_X(X,Xw,XID)
   if (io_X_err<0) then
     call warning('Reading of io_X heared failed despite previous check.')
     return
   endif
 endif
 !
 do iq_W=1,q%nibz
   !
   isc%qs(2)=iq_W
   call DEV_SUB(scatter_Gamp)(isc,'c')
   !
   if (.not.W_is_unscreened) then
     !
     if(.not.BS_eM_is_manual) then
     !
     IO_ACT=manage_action(RD_CL_IF_END,iq_W,1,q%nibz)
     call io_control(ACTION=IO_ACT,COM=NONE,SEC=(/2*iq_W,2*iq_W+1/),ID=XID)
     io_X_err=io_X(X,Xw,XID)
     if (io_X_err<0) then
       call warning('Reading of io_X for iq='//trim(intc(iq_W))//' failed despite previous check.')
       return
     endif
     !
     !forall(i_g2=1:BS_n_g_W) X_mat(i_g2,i_g2,1)=X_mat(i_g2,i_g2,1)+1._SP
     if(iq_W.eq.1) then
      BS_K_eM = 1._SP/(real(X_mat(1,1,1))+1._SP)
      if(BS_K_is_HYB)  BS_hybrd_frac = 1._SP/BS_K_eM
     endif
    endif
     !
#if defined _CUDA
     ! device <- host transfer
     X_mat_d=X_mat
     !
     ! CUDA Fortran workarounds
     gamp_d => isc%gamp_d
     nibz=q%nibz
#endif
     !
     if(l_rim_w) then
       !$cuf kernel do(2) <<<*,*>>>
       do i_g3=1,RIM_W_ng
         do i_g2=1,RIM_W_ng
#if defined _CUDA
           RIM_vpW_d(i_g2,i_g3)= RIM_W_d(1,iq_W,i_g2,i_g3)
           if(i_g2==i_g3) RIM_vpW_d(i_g2,i_g2)=RIM_vpW_d(i_g2,i_g2)+2._SP*gamp_d(i_g2,i_g2)
#else
           RIM_vpW(i_g2,i_g3) = RIM_W(1,iq_W,i_g2,i_g3)
           if(i_g2==i_g3) RIM_vpW(i_g2,i_g2)=RIM_vpW(i_g2,i_g2)+2._SP*isc%gamp(i_g2,i_g2)
#endif
         enddo
       enddo
       !
     end if
     !
     if (BS_W_is_diagonal) then
      if(BS_K_is_SXX.or.BS_K_is_HYB) then
       forall(i_g2=1:BS_n_g_W) BS_W(i_g2,1,iq_W)=1._SP/BS_K_eM*isc%gamp(i_g2,1)
      elseif(BS_K_is_MDL) then
        !$cuf kernel do(1) <<<*,*>>>
        do i_g2=1, BS_n_g_W
          i_g2_mat=min(i_g2,X%ng)
#if defined _CUDA
          BS_W(i_g2,1,iq_W)=real(X_mat_d(i_g2_mat,i_g2_mat,1),SP)*gamp_d(i_g2,1)
#else
          model_epsilon(i_g2,1) = 1._SP+1._SP/(1._SP/(BS_K_eM-1._SP)+ &
&   1.563*(real(bare_qpg(iq_W,i_g2))/BO2ANG)**2/q_ThmsFrm2   +   &
&   (real(bare_qpg(iq_W,i_g2))/BO2ANG)**4*1.054571726**2*1.0E-1/4._SP/m2wp2    )
          BS_W(i_g2,1,iq_W)=1._SP/model_epsilon(i_g2,1)*isc%gamp(i_g2,1)
#endif
        enddo
      endif
     else if (.not.BS_W_is_diagonal) then
       call warning("off-diagonal W is not supported currently!!")
       return
       !$cuf kernel do(1) <<<*,*>>>
       do i_g3=1,BS_n_g_W
         i_g3_mat=min(i_g3,X%ng)
         do i_g2=1,BS_n_g_W
           i_g2_mat=min(i_g2,X%ng)
           !
           !RIM_W support
           l_RIM_W_g=(l_rim_w.and.i_g2<=RIM_W_ng.and.i_g3<=RIM_W_ng)
           if (RIM_W_is_diagonal.and.l_RIM_W_g) l_RIM_W_g=(i_g2==i_g3)
           !
           if (l_RIM_W_g) then
             !
#if defined _CUDA
             BS_W(i_g2,i_g3,iq_W) = cmplx(RIM_vpW_d(i_g2_mat,i_g3_mat)/2._SP, &
&                             aimag(X_mat_d(i_g2_mat,i_g3_mat,1))*real(gamp_d(i_g2_mat,i_g3_mat)),kind=SP)
             if (i_space_inv==0.and.i_time_rev==1) BS_W(i_g2_mat,i_g3_mat,nibz+iq_W)= cmplx(RIM_vpW_d(i_g3_mat,i_g2_mat)/2._SP, &
&                             aimag(X_mat_d(i_g3_mat,i_g2_mat,1))*real(gamp_d(i_g2_mat,i_g3_mat)),kind=SP)
#else
             BS_W(i_g2,i_g3,iq_W) = cmplx(RIM_vpW(i_g2_mat,i_g3_mat)/2._SP, &
&                             aimag(X_mat(i_g2_mat,i_g3_mat,1))*real(isc%gamp(i_g2_mat,i_g3_mat)),kind=SP)
             if (i_space_inv==0.and.i_time_rev==1) BS_W(i_g2_mat,i_g3_mat,q%nibz+iq_W)= cmplx(RIM_vpW(i_g3_mat,i_g2_mat)/2._SP, &
&                             aimag(X_mat(i_g3_mat,i_g2_mat,1))*real(isc%gamp(i_g2_mat,i_g3_mat)),kind=SP)
#endif
           else

#if defined _CUDA
             BS_W(i_g2,i_g3,iq_W)=X_mat_d(i_g2_mat,i_g3_mat,1)*gamp_d(i_g2,i_g3)
             if (i_space_inv==0.and.i_time_rev==1) BS_W(i_g2,i_g3,nibz+iq_W)=X_mat_d(i_g3_mat,i_g2_mat,1)*gamp_d(i_g2,i_g3)
#else
             BS_W(i_g2,i_g3,iq_W)=X_mat(i_g2_mat,i_g3_mat,1)*isc%gamp(i_g2,i_g3)
             if (i_space_inv==0.and.i_time_rev==1) BS_W(i_g2,i_g3,q%nibz+iq_W)=X_mat(i_g3_mat,i_g2_mat,1)*isc%gamp(i_g2,i_g3)
#endif
           end if
         enddo
       enddo
     endif
     !
   else
     !
     call dev_memcpy(BS_W(:,1,iq_W),DEV_VAR(isc%gamp)(:,1))
     !
   endif
   ! 
 enddo
 !
 ! Anisotropy correction. Note that isc%gamp(1,1)=RIM_anisotropy while
 ! the \delta(G,G') term must be multiplied by the standard RIM_qpg(1,1)
 !
 if (RIM_anisotropy/=0.) then
#ifdef _CUDA
   rtmp=BS_W(1,1,1)
   rtmp=rtmp+RIM_qpg(1,1,1)/2._SP-RIM_anisotropy/2._SP
   BS_W(1,1,1)=rtmp
#else
   BS_W(1,1,1)=BS_W(1,1,1)+RIM_qpg(1,1,1)/2._SP-RIM_anisotropy/2._SP
#endif
 endif
 !
 call elemental_collision_free(isc)
 !
 if (have_cuda) then
   YAMBO_FREE(isc%gamp_d)
   YAMBO_FREE(RIM_vpW_d)
   YAMBO_FREE(RIM_W_d)
 else
   YAMBO_FREE(RIM_vpW)
 endif
 !
 if (.not.W_is_unscreened) call X_ALLOC_elemental('X')
 !
end subroutine
