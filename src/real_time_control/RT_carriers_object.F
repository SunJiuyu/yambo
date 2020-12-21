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
subroutine RT_carriers_object(RT,E,K,DB_id,WHAT)
 !
 ! This subroutine contains the RT objects, included the RT table
 ! In the future the RT code can be based on those, at least partially
 !
 use pars,           ONLY:SP,schlen
 use iso_c_binding,  ONLY:C_LOC,C_F_POINTER
 use drivers,        ONLY:l_elel_scatt,l_elph_scatt,l_elphoton_scatt
 use electrons,      ONLY:n_sp_pol,spin_occ,levels
 use R_lattice,      ONLY:bz_samp
 use RT_control,     ONLY:CARR_RT_IO_t
 use real_time,      ONLY:RT_t,REF_lifetime,E_E_Hole_lifetime,E_E_Electron_lifetime,   &
&                         E_P_Hole_lifetime,E_P_Electron_lifetime,E_P_Ho_abs_lifetime, &
&                         E_P_El_abs_lifetime,RTibz,&
&                         RT_nk,RT_delta_occupations,RT_bands,&
&                         RT_obj_alloc,RT_obj_free,RT_levels,&
&                         E_PHOT_Hole_lifetime,E_PHOT_Electron_lifetime
 use hamiltonian,    ONLY:E_reference
 use IO_m,           ONLY:def_variable_bulk,io_variable_bulk,NF90_unlimited
 !
 implicit none
 !
 ! WHAT=1 allocate and define, 2 deallocate
 !
 type(RT_t),      target, intent(inout):: RT
 type(levels),  intent(inout),optional :: E
 type(bz_samp), intent(in)   ,optional :: K
 integer,       intent(in)   ,optional :: DB_id
 character(*),  intent(in)   ,optional :: WHAT
 !
 integer  :: i_RT,i_b,i_k,iE_fg,i_k_RT,i_spin
 real(SP) :: Hole_lifetime,Elec_lifetime,df_occ,f_occ
 !
 integer           :: VAR_SZ(2),IPOS(2),TMP_SIZE(2)
 real(SP), pointer :: TMP_VAR(:,:)
 character(schlen) :: dim_names(2)
 !
 select case(WHAT)
   !
   case('allocate')
     !
     ! Allocate and define
     !
     RT%nk=RT_nk
     RT%nb=RT_bands
     RT%nstates=RT_nk*(RT_bands(2)-RT_bands(1)+1)
     !
     call RT_obj_alloc(RT)
     !
     RT%E_bare=0._SP
     RT%f_bare=0._SP
     if (allocated(REF_lifetime)) RT%GAMMA_bare=0._SP
     !
     do i_k=1,K%nibz
       do i_k_RT=RTibz%k_range(i_k,1),RTibz%k_range(i_k,2)
         RT%kpt(i_k_RT,:)    = RTibz%pt(i_k_RT,:)
         RT%k_weight(i_k_RT) = RTibz%weights( i_k_RT)
       enddo
     enddo
     !
     call build_RT_table()
     !
     do i_RT=1,RT%nstates
       !
       i_b    =RT%table(i_RT,1)
       i_k    =RT%table(i_RT,2)
       i_k_RT =RT%table(i_RT,3)
       i_spin =RT%table(i_RT,4)
       !
       if(i_k_RT/=RTibz%k_range(i_k,1)) then
         iE_fg=RTibz%E_map(i_k_RT)
         RT%E_bare(i_RT)=E_reference%FG%E(i_b,iE_fg,i_spin)
         if(.not.allocated(E%fo)) RT%f_bare(i_RT)=E%FG%f(i_b,iE_fg,i_spin)
         if(     allocated(E%fo)) RT%f_bare(i_RT)=E%FG%f(i_b,iE_fg,i_spin)-E%FG%df(i_b,iE_fg,i_spin)
       else
         RT%E_bare(i_RT)=E_reference%E(i_b,i_k,i_spin)
         if(.not.allocated(E%fo)) RT%f_bare(i_RT)=E%f(i_b,i_k,i_spin)
         if(     allocated(E%fo)) RT%f_bare(i_RT)=E%fo(i_b,i_k,i_spin)
       endif
       !
     enddo
     !
   case('fill','dump')
     !
     do i_RT=1,RT%nstates
       !
       Hole_lifetime=0._SP
       Elec_lifetime=0._SP
       !
       i_b    =RT%table(i_RT,1)
       i_k    =RT%table(i_RT,2)
       i_k_RT =RT%table(i_RT,3)
       i_spin =RT%table(i_RT,4)
       !
       ! Occupations
       !
       if (what=="fill") then
         RT%delta_f(i_RT)=RT_delta_occupations(i_b,i_k_RT)
       else if (what=="dump".and.allocated(RT_delta_occupations)) then
         RT_delta_occupations(i_b,i_k_RT)=RT%delta_f(i_RT)
       endif
       !
       ! Energies
       !
       if (what=="fill") then
         RT%delta_E(i_RT)=RT_levels%E(i_b,i_k_RT,i_spin)-RT%E_bare(i_RT)
       else if (what=="dump".and.allocated(RT_levels%E)) then
         RT_levels%E(i_b,i_k_RT,i_spin)=RT%delta_E(i_RT)+RT%E_bare(i_RT)
       endif
       !
       if(l_elphoton_scatt) then
         if (what=="fill") then
           RT%dGAMMA_E_PHOT_Hole(i_RT)=E_PHOT_Hole_lifetime(i_b,i_k_RT)
           RT%dGAMMA_E_PHOT_Electron(i_RT)=E_PHOT_Electron_lifetime(i_b,i_k_RT)
         else if (allocated(E_PHOT_Hole_lifetime)) then
           E_PHOT_Hole_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_PHOT_Hole(i_RT)
           E_PHOT_Electron_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_PHOT_Electron(i_RT)
         endif
         if (allocated(E_PHOT_Hole_lifetime)) then
           Hole_lifetime=Hole_lifetime+E_PHOT_Hole_lifetime(i_b,i_k_RT)
           Elec_lifetime=Elec_lifetime+E_PHOT_Electron_lifetime(i_b,i_k_RT)
         endif
       endif
       !
       if(l_elel_scatt) then
         if (what=="fill") then
           RT%dGAMMA_E_E_Hole(i_RT)=E_E_Hole_lifetime(i_b,i_k_RT)
           RT%dGAMMA_E_E_Electron(i_RT)=E_E_Electron_lifetime(i_b,i_k_RT)
         else if (allocated(E_E_Hole_lifetime)) then
           E_E_Hole_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_E_Hole(i_RT)
           E_E_Electron_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_E_Electron(i_RT)
         endif
         if (allocated(E_E_Hole_lifetime)) then
           Hole_lifetime=Hole_lifetime+E_E_Hole_lifetime(i_b,i_k_RT)
           Elec_lifetime=Elec_lifetime+E_E_Electron_lifetime(i_b,i_k_RT)
         endif
       endif
       !
       if(l_elph_scatt) then
         if (what=="fill") then
           RT%dGAMMA_E_P_Hole(i_RT)=E_P_Hole_lifetime(i_b,i_k_RT)
           RT%dGAMMA_E_P_Electron(i_RT)=E_P_Electron_lifetime(i_b,i_k_RT)
           RT%dGAMMA_E_P_Ho_abs(i_RT)=E_P_Ho_abs_lifetime(i_b,i_k_RT)
           RT%dGAMMA_E_P_El_abs(i_RT)=E_P_El_abs_lifetime(i_b,i_k_RT)
         else if (allocated(E_P_Hole_lifetime)) then
           E_P_Hole_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_P_Hole(i_RT) 
           E_P_Electron_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_P_Electron(i_RT) 
           E_P_Ho_abs_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_P_Ho_abs(i_RT) 
           E_P_El_abs_lifetime(i_b,i_k_RT)=RT%dGAMMA_E_P_El_abs(i_RT) 
         endif
         if (allocated(E_P_Hole_lifetime)) then
           Hole_lifetime=Hole_lifetime+E_P_Hole_lifetime(i_b,i_k_RT)    +E_P_Ho_abs_lifetime(i_b,i_k_RT)
           Elec_lifetime=Elec_lifetime+E_P_Electron_lifetime(i_b,i_k_RT)+E_P_El_abs_lifetime(i_b,i_k_RT)
         endif
       endif
       !
       if (what=="fill") then
         !
         f_occ =RT%f_bare(i_RT)/spin_occ
         df_occ=RT%delta_f(i_RT)/spin_occ
         if(l_elel_scatt.or.l_elph_scatt) then
           RT%delta_GAMMA(i_RT)=  Hole_lifetime+Elec_lifetime
           Hole_lifetime=Hole_lifetime+REF_lifetime(i_b,i_k_RT,1)+REF_lifetime(i_b,i_k_RT,3)+&
&                                      REF_lifetime(i_b,i_k_RT,5)+REF_lifetime(i_b,i_k_RT,7)
           Elec_lifetime=Elec_lifetime+REF_lifetime(i_b,i_k_RT,2)+REF_lifetime(i_b,i_k_RT,4)+&
&                                      REF_lifetime(i_b,i_k_RT,6)+REF_lifetime(i_b,i_k_RT,8)
           !
           RT%NEQ_Lifetime(i_RT)=  (Hole_lifetime+Elec_lifetime)
           if(RT%E_bare(i_RT)> 0._SP) then
             if ((      f_occ)>0._SP) RT%NEQ_Lifetime(i_RT)=RT%NEQ_Lifetime(i_RT)-Hole_lifetime/       f_occ
             if(  abs(df_occ)>0._SP)  RT%NEQ_Lifetime(i_RT)=RT%NEQ_Lifetime(i_RT)-Hole_lifetime/(     df_occ)
           else if(RT%E_bare(i_RT)<=0._SP) then
             if((1._SP-f_occ)>0._SP) RT%NEQ_Lifetime(i_RT)=RT%NEQ_Lifetime(i_RT)-Elec_lifetime/(1._SP-f_occ)
             if(  abs(df_occ)>0._SP) RT%NEQ_Lifetime(i_RT)=RT%NEQ_Lifetime(i_RT)-Elec_lifetime/(    -df_occ)
           endif
           !
         endif
         !
       endif
       !
     enddo
     !
   case('fill_ref')
     !
     ! Reference Lifetimes
     !
     if(.not.allocated(REF_lifetime)) return
     !
     do i_RT=1,RT%nstates
       !
       Hole_lifetime=0._SP
       Elec_lifetime=0._SP
       !
       i_b    =RT%table(i_RT,1)
       i_k    =RT%table(i_RT,2)
       i_k_RT =RT%table(i_RT,3)
       i_spin =RT%table(i_RT,4)
       !
       RT%GAMMA_bare(i_RT,:)=REF_lifetime(i_b,i_k_RT,:)
       !
     enddo
     !
   case('ref_io')
     !
     call def_variable_bulk(DB_id,"RT_carriers_E_bare",7,(/RT%nstates/),SP)
     call  io_variable_bulk(DB_id,7,R1=RT%E_bare)
     call def_variable_bulk(DB_id,"RT_carriers_f_bare",8,(/RT%nstates/),SP)
     call  io_variable_bulk(DB_id,8,R1=RT%f_bare)
     !
     if (allocated(RT%GAMMA_bare)) then
       call def_variable_bulk(DB_id,"RT_carriers_GAMMA_bare",9,(/RT%nstates,8/),SP)
       call  io_variable_bulk(DB_id,9,R2=RT%GAMMA_bare)
     endif
     !
   case('def_io_variables')
     !
     VAR_SZ=(/RT%nstates,NF90_unlimited/)
     dim_names  =[character(schlen) :: 'RT_nstates', 'time']
     !
     ! Energies
     !
     call def_variable_bulk(DB_id,"RT_carriers_delta_E",10,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
     !
     ! Occupations
     !
     call def_variable_bulk(DB_id,"RT_carriers_delta_f",11,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
     !
     ! Lifetimes
     !
     if (l_elel_scatt.or.l_elph_scatt.or.l_elphoton_scatt) then
       call def_variable_bulk(DB_id,"RT_carriers_delta_GAMMA",12,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
       call def_variable_bulk(DB_id,"RT_NEQ_Lifetime",        13,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
     endif
     if (l_elel_scatt) then
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_E_Hole",    14,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_E_Electron",15,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
     endif
     if (l_elphoton_scatt) then
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_PHOT_Hole",    16,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_PHOT_Electron",17,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
     endif
     if (l_elph_scatt) then
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_P_Hole",    18,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_P_Electron",19,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_P_Ho_abs",  20,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
       call def_variable_bulk(DB_id,"RT_dGAMMA_E_P_El_abs",  21,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
     endif
     !
   case('do_io_variables')
     !
     IPOS=(/1,CARR_RT_IO_t%N/)
     TMP_SIZE=(/RT%nstates,1/)
     !
     ! Energies
     !
     call C_F_POINTER(C_LOC(RT%delta_E), TMP_VAR, TMP_SIZE)
     call io_variable_bulk(DB_id,10,R2=TMP_VAR,IPOS=IPOS)
     !
     ! Occupations
     !
     call C_F_POINTER(C_LOC(RT%delta_f), TMP_VAR, TMP_SIZE)
     call io_variable_bulk(DB_id,11,R2=TMP_VAR,IPOS=IPOS)
     !
     ! Lifetimes
     !
     if (l_elel_scatt.or.l_elph_scatt.or.l_elphoton_scatt) then
       !
       call C_F_POINTER(C_LOC(RT%delta_GAMMA), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,12,R2=TMP_VAR,IPOS=IPOS)
       !
       call C_F_POINTER(C_LOC(RT%NEQ_Lifetime), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,13,R2=TMP_VAR,IPOS=IPOS)
       !
     endif
     !
     if (l_elel_scatt) then
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_E_Hole), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,14,R2=TMP_VAR,IPOS=IPOS)
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_E_Electron), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,15,R2=TMP_VAR,IPOS=IPOS)
       !
     endif
     !
     if (l_elphoton_scatt) then
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_PHOT_Hole), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,16,R2=TMP_VAR,IPOS=IPOS)
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_PHOT_Electron), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,17,R2=TMP_VAR,IPOS=IPOS)
       !
     endif
     !
     if (l_elph_scatt) then
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_P_Hole), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,18,R2=TMP_VAR,IPOS=IPOS)
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_P_Electron), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,19,R2=TMP_VAR,IPOS=IPOS)
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_P_Ho_abs), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,20,R2=TMP_VAR,IPOS=IPOS)
       !
       call C_F_POINTER(C_LOC(RT%dGAMMA_E_P_El_abs), TMP_VAR, TMP_SIZE)
       call io_variable_bulk(DB_id,21,R2=TMP_VAR,IPOS=IPOS)
       !
     endif
     !
   case('free')
     !
     ! Free
     !
     call RT_obj_free(RT)
     !
 end select
 !
 contains
   !
   subroutine build_RT_table()
     !
     i_RT=0
     !
     do i_k_RT=1,RT%nk
       do i_b=RT%nb(1),RT%nb(2)
         do i_spin=1,n_sp_pol
           !
           i_RT=i_RT+1
           RT%table(i_RT,1)=i_b
           RT%table(i_RT,2)=RTibz%k_map(i_k_RT)
           RT%table(i_RT,3)=i_k_RT
           RT%table(i_RT,4)=i_spin
           !
         enddo
       enddo
     enddo
     !
   end subroutine
   !
end subroutine RT_carriers_object