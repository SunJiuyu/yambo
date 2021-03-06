!
!        Copyright (C) 2000-2020 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM MG
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
subroutine K_components_folded_in_serial_arrays(iq)
 !
 ! Here I fill the kernel (coupling included) in an entire BIG matrix to be
 ! used for diagonalization and/or invertion
 !
 !      | (K_r)     (cI*K_c)    |  
 !  K = |                       |
 !      | (-cI*K_c^*)  (-K_r^*) |
 !
 !
 use electrons,      ONLY:n_sp_pol
 use parallel_int,   ONLY:PP_redux_wait
 use parallel_m,     ONLY:PAR_COM_eh_INDEX,PAR_IND_T_groups,PAR_COM_Xk_ibz_INDEX
 use BS_solvers,     ONLY:BSS_eh_E_SOC_corr,&
&                         BSS_dipoles_opt,BSS_eh_table,BSS_eh_E,BSS_eh_f,BSS_eh_W,BSS_eh_Z,BSS_alloc
 use BS,             ONLY:BS_K_dim,BS_nT_grps,BS_T_grp,BS_K_coupling,BS_perturbative_SOC,&
&                         BS_res_ares_dim,l_BS_optics,l_BS_magnons,l_BS_dichroism,l_BS_photolum
 use MAGNONS,        ONLY:BSS_dipoles_magn,BSS_MAGN_alloc
 use DICHROISM,      ONLY:BSS_dipoles_dich,BSS_DICH_alloc
 use PHOTOLUM,       ONLY:BSS_dipoles_PL,BSS_PL_alloc,BSS_PL_f,PL_weights
 use timing_m,       ONLY:timing
 !
 implicit none
 !
 integer      :: iq
 !
 ! Work Space
 !
 integer     :: i_Tr,H_pos,i_g,i_g_start,i_k_bz,i_v,i_c,i_sp_c,i_sp_v,&
 &              BS_res_ares_dim_solver,i_res_ares,i_res_ares_solver,ares_fac
 !
 !========================
 ! BSS COMPONENTS SECTION
 !========================
 !
 if (allocated(BSS_eh_table)) return
   !
   call timing('BSE components folding',OPR='start')
   !
   ! ALLOCATION
   !============
   !
   call BSS_alloc( )
   call BSS_MAGN_alloc( )
   call BSS_DICH_alloc( )
   call BSS_PL_alloc( )
   !
   BS_res_ares_dim_solver=1
   if(BS_res_ares_dim==1.and.BS_K_coupling) BS_res_ares_dim_solver=2
   !
   ! FILLING BSS_eh_* and BSS_dipoles
   !==================================
   !
   i_g_start=1
   i_res_ares=1
   do i_g=1,BS_nT_grps
     !
     if (.not.PAR_IND_T_groups%element_1D(i_g)) cycle
     !
     if(BS_res_ares_dim==2) then
       i_res_ares=BS_T_grp(i_g)%i_res_ares
       i_g_start =BS_T_grp(i_g)%i_T_ref
     endif
     !
     do i_Tr=1,BS_T_grp(i_g)%size
       !
       i_k_bz=BS_T_grp(i_g)%table(i_Tr,1)
       i_v   =BS_T_grp(i_g)%table(i_Tr,2)
       i_c   =BS_T_grp(i_g)%table(i_Tr,3)
       i_sp_c=BS_T_grp(i_g)%table(i_Tr,4)
       i_sp_v=BS_T_grp(i_g)%table(i_Tr,5)
       !
       !  This loop is for the case
       !  R/A symm    RES/RET/CPL   BSS_f/BSS_E/BSS_dip   res_ares_dim    BS_mats_to_solve     BS_mat_size    
       !    S          CPL             BS_H_dim                 1              1                 BS_H_dim
       !  where the antiresonant transitions are not in the groups
       !
       do i_res_ares_solver=1,BS_res_ares_dim_solver
         !
         H_pos=sum(BS_T_grp(i_g_start:I_g-1)%size)+i_Tr+(max(i_res_ares_solver,i_res_ares)-1)*BS_K_dim(1)
         ares_fac=(-1)**(i_res_ares_solver-1)
         !
         BSS_eh_table(H_pos,:)=(/i_k_bz,i_v,i_c,i_sp_c,i_sp_v/)
         !
         BSS_eh_E(H_pos)=BS_T_grp(i_g)%E(i_Tr)*ares_fac
         BSS_eh_f(H_pos)=BS_T_grp(i_g)%f(i_Tr)*ares_fac
         if (allocated(BS_T_grp(i_g)%Z)) BSS_eh_Z(H_pos)=BS_T_grp(i_g)%Z(i_Tr)*ares_fac
         if (allocated(BS_T_grp(i_g)%W)) BSS_eh_W(H_pos)=BS_T_grp(i_g)%W(i_Tr)*ares_fac
         !
         if (BS_perturbative_SOC) BSS_eh_E_SOC_corr(:,H_pos)=BS_T_grp(i_g)%E_SOC_corr(:,i_Tr)*ares_fac
         !
         if (l_BS_photolum) then
           BSS_PL_f(H_pos)=BS_T_grp(i_g)%f_PL(i_Tr,i_res_ares_solver)
           if (BS_res_ares_dim==1.and.BS_res_ares_dim_solver==1) BSS_PL_f(H_pos+BS_K_dim(1))=BS_T_grp(i_g)%f_PL(i_Tr,2)
         endif
         !
         if (l_BS_optics) BSS_dipoles_opt(:,H_pos)=BS_T_grp(i_g)%dipoles_opt(:,i_Tr,i_res_ares_solver)
         !
         if (l_BS_magnons) BSS_dipoles_magn(:,H_pos)=BS_T_grp(i_g)%dipoles_mag(:,i_Tr,i_res_ares_solver)
         !
         if (l_BS_dichroism) BSS_dipoles_dich(:,H_pos)=BS_T_grp(i_g)%dipoles_dic(:,i_c,i_res_ares_solver)
         !
         if (l_BS_photolum) BSS_dipoles_PL(:,H_pos)=BS_T_grp(i_g)%dipoles_opt(:,i_Tr,i_res_ares_solver)*PL_weights(:)
         !
       enddo
       !
     enddo
     !
   enddo
   !
   ! ALL 2 ALL
   !===========
   !
   ! BSS_eh_* and BSS_dipoles are distributed ONLY among k-points and eh transitions.
   !
   call PP_redux_wait(BSS_eh_E,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait(BSS_eh_E,COMM=PAR_COM_eh_INDEX%COMM)
   call PP_redux_wait(BSS_eh_f,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait(BSS_eh_f,COMM=PAR_COM_eh_INDEX%COMM)
   call PP_redux_wait(BSS_eh_table,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
   call PP_redux_wait(BSS_eh_table,COMM=PAR_COM_eh_INDEX%COMM)
   if (allocated(BSS_eh_W)) then
     call PP_redux_wait(BSS_eh_W,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_eh_W,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   if (allocated(BSS_eh_Z)) then
     call PP_redux_wait(BSS_eh_Z,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_eh_Z,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   if (l_BS_optics) then
     call PP_redux_wait(BSS_dipoles_opt,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_dipoles_opt,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   if (l_BS_magnons) then
     call PP_redux_wait(BSS_dipoles_magn,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_dipoles_magn,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   if (l_BS_dichroism) then
     call PP_redux_wait(BSS_dipoles_dich,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_dipoles_dich,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   if (l_BS_photolum) then
     call PP_redux_wait(BSS_dipoles_PL,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_dipoles_PL,COMM=PAR_COM_eh_INDEX%COMM)
     call PP_redux_wait(BSS_PL_f,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_PL_f,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   if ( BS_perturbative_SOC ) then
     call PP_redux_wait(BSS_eh_E_SOC_corr,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
     call PP_redux_wait(BSS_eh_E_SOC_corr,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   call timing('BSE components folding',OPR='stop')
   !
 !
end subroutine K_components_folded_in_serial_arrays
