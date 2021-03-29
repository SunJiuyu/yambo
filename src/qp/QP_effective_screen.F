!
!        Copyright (C) 2000-2019 the YAMBO team
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
subroutine QP_effective_screen(X,Xw,q,io_X_err)
 !
 use pars,         ONLY:SP
 use stderr,       ONLY:intc
 use frequency,    ONLY:w_samp
 use drivers,      ONLY:l_td_hf
 use com,          ONLY:warning,msg,error
 use collision_el, ONLY:elemental_collision,elemental_collision_free
 use R_lattice,    ONLY:bz_samp,RIM_qpg,RIM_anisotropy,bare_qpg,g_vec
 use D_lattice,    ONLY:i_time_rev,i_space_inv
 use X_m,          ONLY:X_alloc,X_t,X_mat
 use QP_m,         ONLY:QP_hyb_screen,QP_ng_Sx
#ifdef _CUDA
 use X_m,          ONLY:X_mat_d
#endif
 use IO_m,         ONLY:io_control,OP_RD_CL,REP,VERIFY,NONE,OP_RD,RD,RD_CL,RD_CL_IF_END
 use deviceXlib_m, ONLY:dev_memcpy
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
 !
 integer    :: QP_W_size2,i_g2,i_g3,iq_W
 logical    :: QP_W_em1s,QP_W_pp
 !
 ! I/O
 !
 integer           :: ID
 integer, external :: io_X
 ! 
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
   allocate(QP_hyb_screen(QP_ng_Sx,q%nbz))
   QP_hyb_screen=(0.0,0.0)
   X%ng=QP_ng_Sx
   call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=ID)
   io_X_err=io_X(X,Xw,ID) 
   !
   if(io_X_err<0) call error('Incomplete and/or broken PPA/Static diel. fun. database')
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
 X%ng=QP_ng_Sx
 if (.not.l_td_hf) then
   !
   ! PP/Epsm1S DBs
   !
   QP_W_em1s=X%whoami==2
   QP_W_pp  =X%whoami==4
   !
   if (QP_W_em1s) call X_alloc('X',(/QP_ng_Sx,QP_ng_Sx,1/))
   if (QP_W_pp)   call X_alloc('X',(/QP_ng_Sx,QP_ng_Sx,2/))
   !
 endif
 !
 QP_W_size2=QP_ng_Sx
 !if (BS_W_is_diagonal) BS_W_size2=1
 !
 ! When TR is present but not the SI X_mat indexes need to be exchanged 
 ! when the TR is applied
 !
! if (i_space_inv==0.and.i_time_rev==1) then
!   YAMBO_ALLOC(QP_hyb_screen,(QP_ng_Sx,2*q%nibz))
! else
!   YAMBO_ALLOC(QP_hyb_screen,(QP_ng_Sx,q%nibz))
! endif
 !
 YAMBO_ALLOC(isc%gamp,(QP_ng_Sx,QP_W_size2))

 !
 if (.not.l_td_hf) then
   call io_control(ACTION=OP_RD,COM=NONE,SEC=(/1/),ID=ID)
   io_X_err=io_X(X,Xw,ID)
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
     !
     call io_control(ACTION=RD_CL_IF_END,COM=NONE,SEC=(/2*iq_W,2*iq_W+1/),ID=ID)
     io_X_err=io_X(X,Xw,ID)
     if (io_X_err<0) then
       call warning('Reading of io_X for iq='//trim(intc(iq_W))//' failed depite previous check.')
       return
     endif
     !
     forall(i_g2=1:QP_ng_Sx) X_mat(i_g2,i_g2,1)=X_mat(i_g2,i_g2,1)+1._SP
        if(iq_W.eq.1) call msg('r',"inv epsilon0", real(X_mat(1,1,1)))
     !
     ! effecitve sum-over W_bar
     do i_g3=1,QP_ng_Sx
      QP_hyb_screen(i_g3,iq_W) = X_mat(i_g3,i_g3,1)
       do i_g2=1,QP_ng_Sx
         if(i_g2.ne.i_g3 ) then !G// diag only
         QP_hyb_screen(i_g3,iq_W)=QP_hyb_screen(i_g3,iq_W) + X_mat(i_g2,i_g3,1)/bare_qpg(iq_W,i_g2)*bare_qpg(iq_W,i_g3)
!         if (i_space_inv==0.and.i_time_rev==1) QP_hyb_screen(i_g3,q%nibz+iq_W)=QP_hyb_screen(i_g3,q%nibz+iq_W) +X_mat(i_g3,i_g2,1)*isc%gamp(i_g2,i_g3)
	endif
       enddo
!       QP_hyb_screen(i_g3,iq_W)=QP_hyb_screen(i_g3,iq_W)/isc%gamp(i_g3,1)
       if (i_space_inv==0.and.i_time_rev==1) QP_hyb_screen(i_g3,iq_W+q%nibz)=QP_hyb_screen(i_g3,iq_W)
     enddo
     !
   ! 
 enddo
 !
 ! Anisotropy correction. Note that isc%gamp(1,1)=RIM_anisotropy while
 ! the \delta(G,G') term must be multiplied by the standard RIM_qpg(1,1)
 !
! if (RIM_anisotropy/=0.) QP_hyb_screen(1,1,1)=QP_hyb_screen(1,1,1)+RIM_qpg(1,1,1)/2._SP-RIM_anisotropy/2._SP
 !
 call elemental_collision_free(isc)
 !
 if (.not.l_td_hf) call X_alloc('X')
 X%ng=X%ng_db
 !
end subroutine
