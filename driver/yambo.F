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
integer function yambo(np,pid,lnstr,iinf,iind,iod,icd,ijs,instr,inf,ind,od,com_dir,js)
 !
 !Y88b    /   e           e    e      888~~'    ,88~-_   
 ! Y88b  /   d8b         d8b  d8b     888   |  d888   '  
 !  Y88b/   /Y88b       d888bdY88b    888 _/  88888    | 
 !   Y8Y   /  Y88b     / Y88Y Y888b   888  '  88888    | 
 !    Y   /____Y88b   /   YY   Y888b  888   |  Y888   /  
 !   /   /      Y88b /          Y888b 888__/    `88_-~   
 !
 use drivers,       ONLY:l_HF_and_locXC,l_sc_run,l_eval_collisions,l_real_time,&
&                        l_optics,l_chi,l_em1s,l_ppa,l_bss,l_bse,l_em1d,&
&                        l_gw0,l_life,l_acfdt,l_dipoles,l_nl_optics,l_elel_corr,l_phel_corr
#if defined _SCALAPACK
 use drivers,       ONLY:l_slk_test
#endif
#if defined _ELPH
 use drivers,       ONLY:l_elph_Hamiltonian
#endif
 use X_m,           ONLY:X_t
 use DIPOLES,       ONLY:DIPOLE_t
 use parallel_m,    ONLY:myid,ncpu
 use it_tools,      ONLY:it_reset
 use it_m,          ONLY:infile
 use QP_m,          ONLY:QP_ng_Sx,QP_ng_SH
 use wave_func,     ONLY:wf_ng
 use electrons,     ONLY:levels,E_reset
 use R_lattice,     ONLY:bz_samp,bz_samp_reset
 use frequency,     ONLY:w_samp
 !
 implicit none
 !
 type(levels) ::en,Xen,Ken,DIPen
 type(bz_samp)::q,k,Xk
 ! 
 ! 4 response function & energies types: x s d p
 !
 type(X_t)    ::X(4) 
 type(w_samp) ::Xw(4)
 !
 type(DIPOLE_t) :: Dip
 !
 integer,          intent(in)    :: lnstr,iind,iod,ijs,np,pid,icd
 integer,          intent(inout) :: iinf
 character(lnstr), intent(in)    :: instr
 character(iinf),  intent(in)    :: inf
 character(iind),  intent(in)    :: ind
 character(iod),   intent(in)    :: od
 character(ijs),   intent(in)    :: js
 character(icd),   intent(in)    :: com_dir
 !
 ! Local Variables
 !
 integer           ::i_err,i_X_kind
 logical           ::l_sta_X,l_dyn_X,driver_now
 !
 ! Externals
 !
 integer, external ::X_dielectric_matrix,INIT
 !
 ! Presets
 !
 ncpu    = np
 myid    = pid
 infile  = inf
 yambo   = 0
 iinf    = 0 
 !
 ! ... Internal Defaults
 !
 call SET_defaults(instr,ind,od,com_dir)
 !
 ! ... C_driver -> YPP information transfer
 !
 call C_driver_transfer()
 !
 call it_reset(1)
#if defined _MEM_CHECK 
 call MEM_defaults( )
#endif
 call SET_job_strings_and_dirs(js)
 !
 ! not resetted in INIT call
 !
 call bz_samp_reset(Xk)  
 call E_reset(Xen) 
 call E_reset(Ken) 
 call E_reset(DIPen) 
 !
 ! INPUT FILE PARSING
 !====================
 !
 iinf=INIT(en,q,k,X,Xw,Dip,instr,lnstr,.false.,.false.)
 !
 ! iinf=0 => No error message from driver.c 
 !
 if (iinf/=0) then
   yambo=iinf
   return
 endif
 !
 ! PARALLEL and I/O Setup 
 !=======================
 !
 call PARALLEL_and_IO_Setup(en,k)
 !
 ! SETUPS ...
 !============
 !
 ! Basical 
 !
 call setup(en,Xen,Ken,k,Xk)
 !
 ! Update the default variables
 !
 i_err=INIT(en,q,k,X,Xw,Dip,instr,lnstr,.true.,.false.)
 !
 ! BZ sampling
 !==============
 !
 call bz_samp_indexes(k,Xk,q) 
 !
 ! COLOUMB DRIVER
 !================
 !
 call col_driver(maxval((/wf_ng,QP_ng_Sx,QP_ng_SH/)),q)
 !
 ! MAGN and DEN
 !=============
 !
 !call eval_Mtot(en,k,q)
 !
 ! reset the section numbering
 !
 call section('*','') 
 !
 ! QP and RT CORRECTIONS
 !=======================
 !
 call External_corrections(en,Xen,Ken,DIPen,Xk,k,X)
 !
 ! Logicals for \eps^-1
 !======================
 !
 l_dyn_X = l_em1d .or. (l_optics.and.l_chi) .or. (l_life.and.l_elel_corr) .or. l_ppa
 l_sta_X = l_em1s .or. (l_optics.and.l_bse)
 !
 i_X_kind=1
 if ( l_real_time .or. l_sc_run) i_X_kind=1
 if ( l_sta_X .and. .not.l_ppa ) i_X_kind=2
 if ( l_dyn_X .and. .not.l_ppa ) i_X_kind=3
 if (                    l_ppa ) i_X_kind=4
 !
 ! DIPOLES
 !=========
 !
 driver_now= l_sta_X.or.l_dyn_X.or.l_dipoles.or.l_real_time.or.l_sc_run.or.l_nl_optics
 !
 if (driver_now) call DIPOLE_driver(DIPen, Xk, q, Dip)
 if (driver_now) call mem_manager_report
 !
 ! OPTICS
 !========
 !
 driver_now=l_optics.and.l_chi
 !
 if (driver_now) call OPTICS_driver(Xen,Xk,q,Xw(3),X(3),Dip)
 if (driver_now) call mem_manager_report
 !
 ! INVERSE DIELECTRIC FUNCTIONS
 !===============================
 !
 driver_now=l_em1s.or.(l_em1d.and..not.l_life).or.l_ppa
 !
 if (driver_now) i_err=X_dielectric_matrix(Xen,Xk,q,X(i_X_kind),Xw(i_X_kind),Dip,.false.)
 if (driver_now) call mem_manager_report
 !
 ! HARTREE-FOCK AND VXC
 !======================
 !
 driver_now=l_HF_and_locXC.and..not.any((/l_sc_run,l_eval_collisions,l_real_time/))
 !
 if (driver_now) call XCo_driver(en,k,Xk,q)
 if (driver_now) call mem_manager_report
 !
 ! EXTENDED COLLISIONS 
 !=====================
 !
 driver_now=l_eval_collisions
 !
#if defined _SC || defined _RT
 if (driver_now) call COLLISIONS_eval(en,X,k,Xk,q,Xw(:2),Dip)
 if (driver_now) call mem_manager_report
#endif
 !
 ! SELF-ENERGIES
 !===============
 !
 driver_now=l_gw0.or.l_life
 !
 if (driver_now.and..not.l_phel_corr) call QP_driver(X(i_X_kind),Xen,Xk,en,k,q,Xw(i_X_kind),Dip)
#if defined _PHEL
 if (driver_now.and.     l_phel_corr) call PHEL_driver(en,k,q)
#endif
 if (driver_now) call mem_manager_report
 !
 ! BETHE-SALPETER
 !================
 !
 driver_now=(l_optics.and.l_bse).or.l_bss
 !
 if (driver_now) call K_driver(Xen,Ken,Xk,q,X(i_X_kind),Xw(i_X_kind),Dip)
 if (driver_now) call mem_manager_report
 !
 ! ACFDT
 !=======
 !
 driver_now=l_acfdt
 !
 if (driver_now) call acfdt_tot_energy(Xen,Xk,En,k,q,X(1),Dip)
 if (driver_now) call mem_manager_report
 !
#if defined _SC
 !
 ! SC
 !====
 !
 driver_now=l_sc_run
 !
 if (driver_now) call SC_driver(X(:2),Xw(:2),Xk,en,k,q,Dip)
 if (driver_now) call mem_manager_report
 !
#endif
 !
#if defined _RT
 !
 ! RT
 !====
 !
 driver_now=l_real_time
 !
 if (driver_now) call RT_driver(en,X(:2),k,q,Dip)
 if (driver_now) call mem_manager_report
 !
#endif
 !
#if defined _NL
 !
 ! RT
 !====
 !
 driver_now=l_nl_optics
 !
 if (driver_now) call NL_driver(en,X(:2),k,q,Dip)
 if (driver_now) call mem_manager_report
 !
#endif
 !
#if defined _ELPH
 !
 ! Electron-Phonon Hamiltonian
 !=============================
 !
 driver_now=l_elph_Hamiltonian
 !
 if (driver_now) call ELPH_Hamiltonian(en,k,q)
 if (driver_now) call mem_manager_report
 !
#endif
 !
#if defined _SCALAPACK 
 !
 ! SLK TESTER
 !============
 !
 driver_now=l_slk_test
 !
 if (driver_now) call SLK_test(en,k,q)
 if (driver_now) call mem_manager_report
 !
#endif
 ! 
 ! CLOSING
 !=========
 !
 call CLOSE_the_run(.TRUE.,.TRUE.,.FALSE.)
 !
 i_err=INIT(en,q,k,X,Xw,Dip,instr,lnstr,.false.,.true.)
 !
 call it_reset(-1)
 !
 contains
  !
  subroutine mem_manager_report
#if defined _MEM_CHECK 
    call MEM_report("check")
#endif
  end subroutine mem_manager_report
  !
end function yambo
