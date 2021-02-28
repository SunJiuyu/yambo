!
!        Copyright (C) 2000-2020 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AF DS
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
function K_ALDA_kernel_resonant(iq, BS_n_g_exch, BS_T_grp_ip, i_Tp, BS_T_grp_ik, i_Tk) result(H_alda)
  !==============================
  !
  ! :::ALDA G FXC    (resonant):::
  !=============================
  !
  use pars,         ONLY: SP,rZERO,pi
  use BS,           ONLY: BS_T_group
  use wrapper_omp,  ONLY: Vstar_dot_V_omp
  use wrapper,      ONLY: M_by_V,Vstar_dot_V
  use R_lattice,    ONLY: bare_qpg, bare_qpg_d
  use TDDFT,        ONLY: FXC,  FXC_n_g_corr
  use com,          ONLY: warning,msg
#include<memory.h>
!  implicit none
  !
  integer,                  intent(in) :: iq,BS_n_g_exch,i_Tp,i_Tk
  type(BS_T_group), target, intent(in) :: BS_T_grp_ip,BS_T_grp_ik
  complex(SP) :: H_alda, H_x
  complex(SP), allocatable :: M1(:), M2(:),M_times_FXC(:)
  integer :: i_g1, i_g2
  !
  if( FXC_n_g_corr .ne. BS_n_g_exch) then
    call warning("input FXC_G is not equal to FXC_n_g_corr in yambo.in")
  endif
  !
#ifdef _CUDA
  complex(SP), pointer, device :: BS_T_grp_ip_O_x_d(:,:)
  complex(SP), pointer, device :: BS_T_grp_ik_O_x_d(:,:)
  integer :: ig
  real(SP):: H_x_r,H_x_i
#endif

#ifdef _CUDA
  !
  call warning("ONLY hartree term will be calculated, i.e. exhcnage kernel")
  BS_T_grp_ip_O_x_d => BS_T_grp_ip%O_x_d
  BS_T_grp_ik_O_x_d => BS_T_grp_ik%O_x_d
  !
  H_x_r=rZERO
  H_x_i=rZERO
  !$cuf kernel do(1) <<<*,*>>>
  do ig = 1, BS_n_g_exch
    H_x=conjg(BS_T_grp_ip_O_x_d(ig,i_Tp))*BS_T_grp_ik_O_x_d(ig,i_Tk)/bare_qpg_d(ig,iq)**2
    H_x_r=H_x_r+real(H_x,SP)
    H_x_i=H_x_i+aimag(H_x)
  enddo
  H_x=cmplx(H_x_r,H_x_i,kind=SP)
  !
#else
  YAMBO_ALLOC(M1,(FXC_n_g_corr))
  YAMBO_ALLOC(M2,(FXC_n_g_corr))
  YAMBO_ALLOC(M_times_FXC,(FXC_n_g_corr))
!  if (i_g1 == 1 ) then
    M1(1) = BS_T_grp_ik%dipoles_R(i_Tk)
    M2(1) = BS_T_grp_ip%dipoles_R(i_Tp)
!  endif
  if(real(M1(1)*FXC(1,1,1)*conjg(M2(1))) .gt. 1e-6) then
    call msg('r', 'Dipole wrong M1: ', real(M1(1)))
    call msg('r', 'Dipole wrong M2: ', real(M2(1)))
  endif


  if(FXC_n_g_corr.gt.1) then
    do i_g1=2,FXC_n_g_corr  
      M1(i_g1) = BS_T_grp_ik%O_x(i_g1,i_Tk)!*bare_qpg(iq,i_g1)
      M2(i_g1) = BS_T_grp_ip%O_x(i_g1,i_Tp)!*bare_qpg(iq,i_g1)
    enddo   
  endif
  !
  !H_alda=Vstar_dot_V_omp(FXC_n_g_corr, M1(:),M2(:))
  H_alda=cmplx(0._SP,0._SP)
  !
  call M_by_V('T',FXC_n_g_corr,FXC(:,:,1),M1,M_times_FXC)
  !
  H_alda=Vstar_dot_V(FXC_n_g_corr,M2,M_times_FXC)/4._SP/pi

  YAMBO_FREE(M1)
  YAMBO_FREE(M2)
  YAMBO_FREE(M_times_FXC)
  !
  ! regular double loop for FXC, negelectable difference from M_by_V if G is large enough
!  do i_g1 = 1,FXC_n_g_corr
!    do i_g2 = 1,FXC_n_g_corr
!      H_alda= H_alda + M1(i_g1)*FXC(i_g1,i_g2,1)*conjg(M2(i_g2))/4._SP/pi
      !if(i_g1==i_g2)  H_alda = H_alda*(1.0-X_store(1,i_g2))
      !H_c = H_c+ H_alda
!       if(i_g1==23 .and. i_g2==1) call msg('s',"FXC: ",real(FXC(!i_g1,i_g2,1)))         
!    enddo
!  enddo
  !
#endif
  !
end function K_ALDA_kernel_resonant
!
!
!
!
function K_ALDA_kernel_coupling(iq, BS_n_g_exch, BS_T_grp_ip, i_Tp, BS_T_grp_ik, i_Tk) result(H_x)
  !==============================
  !
  ! :::EXCHANGE    (coupling):::
  !=============================
  !
  use pars,         ONLY: SP,rZERO
  use BS,           ONLY: BS_T_group
  use wrapper_omp,  ONLY: V_dot_V_omp
  use R_lattice,    ONLY: bare_qpg, minus_G, bare_qpg_d, minus_G_d
  implicit none
  !
  integer,                  intent(in) :: iq, BS_n_g_exch,i_Tp,i_Tk
  type(BS_T_group), target, intent(in) :: BS_T_grp_ip,BS_T_grp_ik
  complex(SP) :: H_x
  !
#ifdef _CUDA
  complex(SP), pointer, device :: BS_T_grp_ip_O_x_d(:,:)
  complex(SP), pointer, device :: BS_T_grp_ik_O_x_d(:,:)
  integer :: ig
  real(SP):: H_x_r,H_x_i
#endif
         
#ifdef _CUDA
  !
  BS_T_grp_ip_O_x_d => BS_T_grp_ip%O_x_d
  BS_T_grp_ik_O_x_d => BS_T_grp_ik%O_x_d
  !
  H_x_r=rZERO
  H_x_i=rZERO
  !$cuf kernel do(1) <<<*,*>>>
  do ig = 1, BS_n_g_exch
    H_x = BS_T_grp_ip_O_x_d( minus_G_d(ig),i_Tp) *BS_T_grp_ik_O_x_d(ig,i_Tk)/bare_qpg_d(ig,iq)**2
    H_x_r=H_x_r+real(H_x,SP)
    H_x_i=H_x_i+aimag(H_x)
  enddo
  H_x=cmplx(H_x_r,H_x_i,SP)
#else
  H_x=V_dot_V_omp(BS_n_g_exch, BS_T_grp_ip%O_x( minus_G(1:BS_n_g_exch),i_Tp),&
                               BS_T_grp_ik%O_x(:,i_Tk)/bare_qpg(iq,:BS_n_g_exch)**2)
#endif
 !
end function K_ALDA_kernel_coupling
