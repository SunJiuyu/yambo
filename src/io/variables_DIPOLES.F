!
!        Copyright (C) 2000-2020 the YAMBO team
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
integer function variables_DIPOLES(Dip,ID)
 !
 use pars,          ONLY:SP
 use units,         ONLY:HA2EV
 use DIPOLES,       ONLY:DIPOLE_t,DIP_desc
 use IO_m,          ONLY:io_variable_elemental,def_variable_elemental,io_status,io_elemental
 use global_XC,     ONLY:Dipole_WF_xc_string,loaded_WF_xc_string
 use descriptors,   ONLY:IO_desc_reset
 !
 implicit none
 !
 type(DIPOLE_t) :: Dip
 integer        :: ID
 !
 ! Workspace
 !
 call IO_desc_reset(DIP_desc)
 !
 call def_variable_elemental(ID,"PARS",11,SP,0,par_io_kind="independent")
 !
 call io_variable_elemental(ID,DB_I1=Dip%db_ib,VAR="DIP band range",I1=Dip%ib,CHECK=.true.,OP=(/">=","<="/),&
&                           DESCRIPTOR=DIP_desc)
 call io_variable_elemental(ID,DB_I1=Dip%db_ib_lim,VAR="DIP band range limits",I1=Dip%ib_lim,CHECK=.true.,OP=(/"<=",">="/),&
&                           DESCRIPTOR=DIP_desc)
 call io_variable_elemental(ID,UNIT=HA2EV,&
&     VAR="DIP e/h energy range ",TERMINATOR="[eV]",R1=Dip%ehe,CHECK=.true.,OP=(/">=","<="/),DESCRIPTOR=DIP_desc)
 call io_variable_elemental(ID,VAR="RL vectors in the sum",I0=Dip%ng,WARN=.true.,OP=(/"<="/),&
&                           DESCRIPTOR=DIP_desc)
 call io_variable_elemental(ID,VAR="[r,Vnl] included",L0=Dip%Vnl_included,CHECK=.true.,OP=(/"=="/),&
&                           DESCRIPTOR=DIP_desc)
 call io_variable_elemental(ID,VAR="Bands ordered",L0=Dip%bands_ordered,CHECK=.true.,OP=(/">="/),&
&                           DESCRIPTOR=DIP_desc)
 call io_variable_elemental(ID,VAR="Direct v evaluation",L0=Dip%force_v_g_space,CHECK=.true.,OP=(/"=="/),&
&                           DESCRIPTOR=DIP_desc)
 !
 ! Character variables
 !
 call def_variable_elemental(ID,"APPROACH",1,0,0,par_io_kind="independent")
 call  io_variable_elemental(ID,VAR="Approach used",CH0=Dip%approach,CHECK=.true.,OP=(/"=="/),&
&                           DESCRIPTOR=DIP_desc)
 !
 call def_variable_elemental(ID,"KINDS",1,0,0,par_io_kind="independent")
 call  io_variable_elemental(ID,VAR="Dipoles computed",CH0=Dip%computed,CHECK=.true.,OP=(/"<="/),&
&                           DESCRIPTOR=DIP_desc)
 !
 ! Wavefunctions xc 
 !
 call def_variable_elemental(ID,'WAVE_FUNC_XC',1,0,0)
 call  io_variable_elemental(ID,DB_CH0=Dipole_WF_xc_string,CH0=loaded_WF_xc_string,VAR=' Wavefunctions',CHECK=.true.,OP=(/"=="/),&
&                           DESCRIPTOR=DIP_desc)
 call io_elemental(ID,VAR="",VAR_SZ=0)
 !
 variables_DIPOLES=io_status(ID)
 !
end function variables_DIPOLES
