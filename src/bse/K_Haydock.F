!        Copyright (C) 2000-2020 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM MG DS
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
subroutine K_Haydock(iq,W,mode)
 !
 use pars,           ONLY:SP,schlen,cI,cZERO,cONE,pi
 use functions,      ONLY:NAN
 use X_m,            ONLY:global_gauge
 use wrapper_omp,    ONLY:V_by_V_plus_V_omp
 use frequency,      ONLY:w_samp
 use com,            ONLY:msg,isec
 use BS_solvers,     ONLY:BSS_desc,Haydock_threshold,io_BSS_Haydock,&
&                         Haydock_vector,BS_HAYVEC_alloc,BS_HAYVEC_free,Haydock_terminator,Haydock_v0_mod, &
&                         Haydock_iterIO,Haydock_iterMAX,HAYVEC_size
 use BS,             ONLY:BS_K_Coupling,BS_Blocks_symmetrize_K,BS_T_grp,BS_nT_grps,BS_K_dim,&
&                         l_BS_ares_asymm,BS_H_dim,BS_mat_res_ares_dim
 use LIVE_t,         ONLY:LIVE_message
 use stderr,         ONLY:tty_size,intc,real2ch
 use IO_m,           ONLY:io_control,OP_RD_CL,OP_WR_CL,VERIFY,REP
 use parser_m,       ONLY:parser
 use parallel_int,   ONLY:PP_bcast
 use parallel_m,     ONLY:PAR_IND_T_groups,PAR_IND_T_Haydock,PAR_COM_T_Haydock,PAR_COM_Xk_ibz_INDEX, &
                          PAR_COM_eh_INDEX
 use timing_m,       ONLY:timing
 use descriptors,    ONLY:IO_desc_add
 !
#include<memory.h>
 !
 type(w_samp), intent(in) :: W
 integer,      intent(in) :: iq
 character(*), intent(in) :: mode
 !
 ! Work Space
 !
 integer            ::it,it_on_disk,i_g,ID,io_err,i_desc_acc,i_desc_iter
 type(Haydock_vector),allocatable :: Vn(:),Vnm1(:),Vnp1(:),V0(:),V_tmp(:)
 complex(SP), pointer :: Vn_f(:)   => null()
 complex(SP), pointer :: Vnm1_f(:) => null()
 complex(SP), pointer :: Vnp1_f(:) => null()
 complex(SP), allocatable :: Cf(:) 
 !
 type(Haydock_vector),allocatable :: W0(:)
 !
 real(SP)           ::reached_threshold
 real(SP)           ::Bf(Haydock_iterMAX+1),Af(Haydock_iterMAX)
 integer            ::ItParity,Cf_size,Haydock_n_converged,Haydock_required_converged
 complex(SP)        ::v0_mod_
 character(schlen)  ::dumb_ch
 logical            ::impose_pseudo_herm_scheme,use_left_and_right,velocity_correction,tmp_converged
 logical, external  ::STOP_now
 character(1)       ::sec
 character(11)      ::scheme
 !
 if (BS_mat_res_ares_dim==2) call error(" Haydock not coded for case  with two BSE matricies")
 !
 call K_output_file(iq,"open Haydock "//trim(mode))
 !
 call PARALLEL_Haydock_VEC_COMMs('assign')
 !
 call timing('Haydock Solver',OPR='start')
 !
 if (isec(2)/=0) sec="="
 if (isec(2)==0) sec="+"
 !
 Haydock_n_converged=0
 Haydock_required_converged=2
 !
 call BS_Blocks_symmetrize_K()
 !
 ! Use terminator?
 if (.not.Haydock_terminator) call parser('BSHayTer',Haydock_terminator)
 ! Developer option
 call parser('ImposePseudoHerm',impose_pseudo_herm_scheme)
 !
 ! Which approach
 use_left_and_right=.false.
 scheme="hermitian"
 if(BS_K_Coupling.or.impose_pseudo_herm_scheme) scheme="pseudo-herm"
 !if(BS_K_Coupling.or.impose_pseudo_herm_scheme.or.allocated(BS_T_grp(1)%Z)) scheme="pseudo-herm"
 use_left_and_right=trim(mode)=="kerr".or.trim(mode)=="dichroism"
 velocity_correction= trim(global_gauge)=="velocity" .and. (trim(mode)=="abs" .or. trim(mode)=="kerr")
 !
 call section(sec,'Haydock Solver for '//trim(mode)//' @q'//trim(intc(iq))//", scheme "//trim(scheme))
 !
 if (.not.l_BS_ares_asymm) HAYVEC_size = BS_K_dim(1)
 if (     l_BS_ares_asymm) then
   if(.not.BS_K_coupling) HAYVEC_size = BS_K_dim(1)
   if(     BS_K_coupling) HAYVEC_size = BS_H_dim
 endif
 !
 ! Allocate Haydock vectors:
 !
 allocate(V0(BS_nT_grps))
 call BS_HAYVEC_alloc(V0)
 !
 if(velocity_correction) then
   allocate(V_tmp(BS_nT_grps))
   call BS_HAYVEC_alloc(V_tmp)
 endif
 !
 if(use_left_and_right) then
   allocate(W0(BS_nT_grps))
   call BS_HAYVEC_alloc(W0)
 endif
 !
 ItParity=1
 it_on_disk=0
 reached_threshold=100._SP
 !
 Cf_size=2
 if (trim(scheme)=="pseudo-herm".or.use_left_and_right) Cf_size=Haydock_iterMAX+1
 YAMBO_ALLOC(Cf,(Cf_size))
 Cf=cZERO
 !
 allocate(Vn(BS_nT_grps))
 allocate(Vnp1(BS_nT_grps))
 allocate(Vnm1(BS_nT_grps))
 call BS_HAYVEC_alloc(Vn)
 call BS_HAYVEC_alloc(Vnm1)
 call BS_HAYVEC_alloc(Vnp1)
 !
 ! Haydock Coefficients and restart vector from disk
 !=================================================
 call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/1/),ID=ID)
 io_err=io_BSS_Haydock(ID,iq,1,reached_threshold,mode,Af(1),Bf(:2),Cf(:2),Vnm1_f,Vn_f,Vnp1_f)
 !
 if (io_err>0) then
   !
   YAMBO_ALLOC_P(Vn_f,(HAYVEC_size))
   YAMBO_ALLOC_P(Vnm1_f,(HAYVEC_size))
   if (scheme=="pseudo-herm") then
     YAMBO_ALLOC_P(Vnp1_f,(HAYVEC_size))
   endif
   !
   it_on_disk=io_err
   if (scheme=="pseudo-herm") Cf_size=it_on_disk
   call io_control(ACTION=OP_RD_CL,SEC=(/1,2/),ID=ID)
   io_err=io_BSS_Haydock(ID,iq,it_on_disk,reached_threshold,mode,Af(:it_on_disk),&
&                   Bf(:it_on_disk+1),Cf(:Cf_size),Vnm1_f,Vn_f,Vnp1_f)
   !
   call msg('s','[Haydock] Restart DB with accuracy ',reached_threshold)
   !  
   call PP_bcast(Af,0,COMM=PAR_COM_Xk_ibz_INDEX%COMM)
   call PP_bcast(Af,0,COMM=PAR_COM_eh_INDEX%COMM)
   call PP_bcast(Bf,0,COMM=PAR_COM_Xk_ibz_INDEX%COMM)
   call PP_bcast(Bf,0,COMM=PAR_COM_eh_INDEX%COMM)
   if (scheme=="pseudo-herm".or.use_left_and_right) then
     call PP_bcast(Cf,0,COMM=PAR_COM_Xk_ibz_INDEX%COMM)
     call PP_bcast(Cf,0,COMM=PAR_COM_eh_INDEX%COMM)
   endif
   !
   call HAYVEC_scatter(Vn,Vn_f,HAYVEC_size)
   call HAYVEC_scatter(Vnm1,Vnm1_f,HAYVEC_size)
   if (scheme=="pseudo-herm") call HAYVEC_scatter(Vnp1,Vnp1_f,HAYVEC_size)
   !
   YAMBO_FREE_P(Vnm1_f)
   YAMBO_FREE_P(Vn_f)
   if (scheme=="pseudo-herm") then
     YAMBO_FREE_P(Vnp1_f)
   endif
   !
 end if
 !
 ! Update Descriptions
 !=====================
 call IO_desc_add(BSS_desc,kind="s",str=" ",S=" ")
 !
 call IO_desc_add(BSS_desc,kind="r",str='[Haydock] Accuracy (requested)',term='o/o',R=(/Haydock_threshold/))
 call IO_desc_add(BSS_desc,kind="r",str='[Haydock] Accuracy (reached)',term='o/o',R=(/reached_threshold/))
 i_desc_acc=BSS_desc%n
 call IO_desc_add(BSS_desc,kind="i",str='[Haydock] Iterations',I=(/it_on_disk/))
 i_desc_iter=BSS_desc%n
 !
 call msg('r','Accuracy (requested)      ',Haydock_threshold,"[o/o]")
 !
 ! Initialization of Haydock vectors and tridiagonal matrix elements
 !
 if (it_on_disk==0) then
   !
   do i_g=1,BS_nT_grps
     !
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     !
     Vn(i_g)%fragment  = cZERO
     Vnm1(i_g)%fragment= cZERO
     Vnp1(i_g)%fragment= cZERO
     !
     V0(i_g)%fragment  = cZERO
     if(use_left_and_right) W0(i_g)%fragment  = cZERO
     if(velocity_correction) V_tmp(i_g)%fragment = cZERO
     !
   enddo
   !
   Af= 0._SP
   Bf= 0._SP
   !
 end if
 !
 do i_g=1,BS_nT_grps
   if (.not.PAR_IND_T_groups%element_1D(i_g)) cycle
   if (trim(mode)=="abs")       V0(i_g)%fragment = BS_T_grp(i_g)%dipoles_opt(1,:,1)
   if (trim(mode)=="kerr")      V0(i_g)%fragment = BS_T_grp(i_g)%dipoles_opt(2,:,1)
   if (trim(mode)=="magnons")   V0(i_g)%fragment = BS_T_grp(i_g)%dipoles_mag(1,:,1)
   V0(i_g)%fragment = V0(i_g)%fragment*sqrt(cmplx(BS_T_grp(i_g)%f_RES,kind=SP))
   if(velocity_correction) V_tmp(i_g)%fragment = V0(i_g)%fragment/BS_T_grp(i_g)%E
 enddo
 !
 do i_g = 1,BS_nT_grps
   if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
#if defined _MPI
   if(     velocity_correction) call PP_bcast(V_tmp(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
   if(.not.velocity_correction) call PP_bcast(   V0(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
#endif
 enddo
 !
 if(velocity_correction) call K_multiply_by_V(V_tmp,V0,1)
 !
 if(use_left_and_right) then
   !
   do i_g=1,BS_nT_grps
     if (.not.PAR_IND_T_groups%element_1D(i_g)) cycle
     if (trim(mode)=="kerr") W0(i_g)%fragment = BS_T_grp(i_g)%dipoles_opt(1,:,1)
     W0(i_g)%fragment = W0(i_g)%fragment*(sqrt(cmplx(BS_T_grp(i_g)%f_RES,kind=SP)))
     if(velocity_correction) V_tmp(i_g)%fragment = W0(i_g)%fragment/BS_T_grp(i_g)%E
   enddo
   !
   do i_g = 1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
#if defined _MPI
     if(     velocity_correction) call PP_bcast(V_tmp(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
     if(.not.velocity_correction) call PP_bcast(   W0(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
#endif
   enddo
   !
   if(velocity_correction) call K_multiply_by_V(V_tmp,W0,1)
   !
 endif
 !
 if(velocity_correction) then
   call BS_HAYVEC_free(V_tmp)
   deallocate(V_tmp)
 endif
 !
 if (it_on_disk==0) then
   !
   if (scheme=="pseudo-herm") then       ! [***] Initialization non Hermitian case
     !
     ! Before starting iterate one needs
     !
     ! | Vn > =     |V0>/(<V0|F(H|V0>))
     ! | Vnp1 > = H |Vn>
     !
     Cf= cZERO
     !
     call K_multiply_by_V(V0,Vnp1,-1)
     call K_dot_product(v0_mod_,V0,Vnp1,1)
     !
     ! For the coupling case if A = -R* (l_BS_ares_asymm=.true.), V0 = (u0, -iu0*)
     ! At each iteration however the relation changes V1= (u1, iu1*),
     ! this is why the ItParity is needed 
     !
     if(.not.l_BS_ares_asymm) Haydock_v0_mod=sqrt(2*real(v0_mod_))
     if(     l_BS_ares_asymm) Haydock_v0_mod=sqrt(  real(v0_mod_))
     !
     do i_g=1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       V0(i_g)%fragment   = V0(i_g)%fragment/Haydock_v0_mod
       if (use_left_and_right) W0(i_g)%fragment = W0(i_g)%fragment/Haydock_v0_mod
       Vn(i_g)%fragment   = V0(i_g)%fragment
       Vnp1(i_g)%fragment = Vnp1(i_g)%fragment/Haydock_v0_mod
     enddo
     !
     ! <w0|u0>=<V0|F|V0>  since <w0|=<V0|F
     ! Is it correct that also V0 was renormalized ??
     if (.not.use_left_and_right) call K_dot_product(v0_mod_,V0,V0,1)
     if (     use_left_and_right) call K_dot_product(v0_mod_,W0,V0,1)
     !
     if(.not.l_BS_ares_asymm) Cf(1)=2*cI*aimag(v0_mod_)
     if(     l_BS_ares_asymm) Cf(1)=           v0_mod_
     !
   else                           ! [***] Initialization Hermitian case
     !
     ! it=1:|Vn>=|q->0>/(<q->0|q->0>)^(1/2)
     !
     call K_dot_product(v0_mod_,V0,V0,1)
     Haydock_v0_mod=sqrt(real(v0_mod_,SP))
     !
     do i_g=1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       V0(i_g)%fragment  = V0(i_g)%fragment/Haydock_v0_mod
       if (use_left_and_right) W0(i_g)%fragment = W0(i_g)%fragment/Haydock_v0_mod
       Vn(i_g)%fragment  = V0(i_g)%fragment
     enddo
     !
     if (use_left_and_right) then
       call K_dot_product(v0_mod_,W0,V0,1)
       Cf(1)= v0_mod_
     endif
     !
   endif
   !
 endif
 !
 ! DEBUG <
 ! write(*,*) "Haydock_v0_mod= ", Haydock_v0_mod
 ! DEBUG >
 !
 ! [***] Loop ...
 !
 do it=it_on_disk+1,Haydock_iterMAX
   !
   if (trim(scheme)=="hermitian") then  ! [***] Hermitian case:
     !
     !|Vn+1> = H |Vn>
     !
     call  K_multiply_by_V(Vn,Vnp1,1)
     !
     ! A(n) = <Vn|Vn+1>
     !
     call K_dot_product(v0_mod_,Vn,Vnp1,1)
     Af(it)=real(v0_mod_,SP)
     !
     !|Vn+1> = |Vn+1> - A(n)|Vn> - B(n)|Vn-1>
     !
     do i_g = 1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       call V_by_V_plus_V_omp(BS_T_grp(i_g)%size,-Af(it)*cONE,Vn(i_g)%fragment,Vnp1(i_g)%fragment)
       call V_by_V_plus_V_omp(BS_T_grp(i_g)%size,-Bf(it)*cONE,Vnm1(i_g)%fragment,Vnp1(i_g)%fragment)
     end do
     !
     ! B(n+1)= <Vn+1|Vn+1>^(1/2)
     !
     call K_dot_product(v0_mod_,Vnp1,Vnp1,1)
     Bf(it+1)=sqrt(real(v0_mod_,SP))
     if( Bf(it+1)==(Bf(it+1)+1._SP) ) call error('Bf=NaN likely because some eigenvalue of the BSE is negative.')
     !
     !|Vn-1> = |Vn>
     !
     do i_g=1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       Vnm1(i_g)%fragment  =Vn(i_g)%fragment
     enddo
     !
     !|Vn> =|Vn+1> / B(n+1)
     !
     do i_g=1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       Vn(i_g)%fragment  =Vnp1(i_g)%fragment/Bf(it+1)
     enddo
     !
     ! DEBUG <
     !! check orthonotmality
     !call K_dot_product(v0_mod_,V0,Vn,1)
     !write(100,*) "Ortho check = ",v0_mod_
     ! DEBUG
     !
     if (use_left_and_right) then
       Cf_size=it+1
       call K_dot_product(v0_mod_,W0,Vn,1)
       Cf(it+1)= v0_mod_
     endif
     !
   else                           ! [***] Non-Hermitian case:
     !
     ! A(n) = <Vn+1|F|Vn+1> (that is <Vn|HFH|Vn>),
     !
     if (.not.l_BS_ares_asymm) then
       ! A(n) =  0   by symmetry
       Af(it)=0._SP
       ItParity=(-1)**(it+1)
     else
       call K_dot_product(v0_mod_,Vnp1,Vnp1,1)
       Af(it)=real(v0_mod_,SP)
     endif
     !
     !|Vn+1> = |Vn+1> - A(n)|Vn> - B(n)|Vn-1>
     !
     if (l_BS_ares_asymm) then
       do i_g = 1,BS_nT_grps
         if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
         call V_by_V_plus_V_omp(BS_T_grp(i_g)%size,-Af(it)*cONE,Vn(i_g)%fragment,Vnp1(i_g)%fragment)
       enddo
     endif
     !
     do i_g = 1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       call V_by_V_plus_V_omp(BS_T_grp(i_g)%size,-Bf(it)*cONE,Vnm1(i_g)%fragment,Vnp1(i_g)%fragment)
     end do
     !
     ! |Vn-1> = |Vn>
     ! |Vn>   = |Vn+1>
     !
     do i_g=1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       Vnm1(i_g)%fragment = Vn(i_g)%fragment
       Vn(i_g)%fragment   = Vnp1(i_g)%fragment
     enddo
     !
     !|Vn+1> = H |Vn>
     !
     call K_multiply_by_V(Vn,Vnp1,ItParity)
     !
     ! B(n+1)= <Vn|F|Vn+1>^(1/2) (that is <Vn|FH|Vn>^(1/2))
     !
     call K_dot_product(v0_mod_,Vn,Vnp1,1) 
     !       = (2*Re(<Vn|Vnp1>))^(1/2) by symmetry,
     !          where the dot_product is just on eh pair space
     if (.not.l_BS_ares_asymm) Bf(it+1)=sqrt(2*real(v0_mod_,SP))
     if (     l_BS_ares_asymm) Bf(it+1)=sqrt(  real(v0_mod_,SP))
     ! Simple test for INF
     if( Bf(it+1)==(Bf(it+1)+1._SP) ) call error('Bf=INF likely because some eigenvalue of the BSE is negative.')
     ! Tests for NaN
     if (NAN(Bf(it+1)))  call error(' Bf=NaN likely because some eigenvalue of BSE is negative.')
     !
     !|Vn> =|Vn+1> / B(n+1)
     !
     do i_g=1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       Vnp1(i_g)%fragment  =Vnp1(i_g)%fragment/Bf(it+1)
       Vn(i_g)%fragment    =Vn(i_g)%fragment/Bf(it+1)
     enddo
     !
     ! Cf(n) = <w0|q_n> = <V0|F|Vn>
     ! F is needed since <w0| = <V0| F
     !
     Cf_size=it+1
     !
     if (.not.use_left_and_right) call K_dot_product(v0_mod_,V0,Vn,1)
     if (     use_left_and_right) call K_dot_product(v0_mod_,W0,Vn,1)
     !
     if (.not.l_BS_ares_asymm) then
       ! Using symmetry arguments
       ! Cf=(1-ItParity)Re(<dipoles|Vn>) + i*(1+ItParity)Im(<dipoles|Vn>)
       ! with all product evaluated just in the eh pairs subspace
       Cf(it+1)=   ((1._SP+real(ItParity,SP))*real(v0_mod_,SP)) +&
       &        cI*((1._SP-real(ItParity,SP))*aimag(v0_mod_))
     else
       Cf(it+1)=   v0_mod_
     endif
     !
   endif
   !
   !
   if (it>2) then
     !
     call K_Haydock_response(iq,it,Cf_size-1,mode,scheme,W,Af(:it),Bf(:it+1),Cf(:Cf_size-1),reached_threshold)
     !
     BSS_desc%rval(1,i_desc_acc) =reached_threshold
     BSS_desc%ival(1,i_desc_iter)=it
     !
     dumb_ch='[Haydock] Iteration '//trim(intc(it))//' Accuracy '//&
&            trim(real2ch(reached_threshold))//'|'//trim(real2ch(Haydock_threshold))
     if (tty_size>0.and.(it_on_disk==0.or. it_on_disk< it-1)) call LIVE_message(dumb_ch,"r","","%s")
     if (tty_size>0.and.(it_on_disk> 0.and.it_on_disk==it-1)) call LIVE_message(dumb_ch,"n","","%s")
     if (tty_size<0) call LIVE_message(dumb_ch,"n","","%s")
     !
   else
     !
     dumb_ch='[Haydock] Iteration '//trim(intc(it))
     call LIVE_message(dumb_ch,"n","","%s")
     !
   endif
   !
   !
   ! DEBUG<
   !write(*,*) it,"Af= ",Af(it)," ; Bf= ",Bf(it)," ; Cf= ",real(Cf(Cf_size-1)),aimag(Cf(Cf_size-1)),reached_threshold
   ! DEBUG>
   !
   tmp_converged=(reached_threshold<abs(Haydock_threshold)) .or. STOP_now(.FALSE.)
   !
   if(     tmp_converged) Haydock_n_converged=Haydock_n_converged+1
   if(.not.tmp_converged) Haydock_n_converged=0
   !
   if ( mod(it,max(1,Haydock_iterIO/10))==0) call K_output_file(iq,"Haydock "//trim(mode))
   !
   if ( mod(it,Haydock_iterIO)==0 .or. Haydock_n_converged==Haydock_required_converged) then
     !
     ! Haydock Coefficients and restart vectors to disk
     !--------------------------------------------------
     !
     YAMBO_ALLOC_P(Vn_f,(HAYVEC_size))
     YAMBO_ALLOC_P(Vnm1_f,(HAYVEC_size))
     if (trim(scheme)=="pseudo-herm") then
       YAMBO_ALLOC_P(Vnp1_f,(HAYVEC_size))
     endif
     !
     call HAYVEC_gather(Vnm1,Vnm1_f,HAYVEC_size)
     call HAYVEC_gather(Vn,Vn_f,HAYVEC_size)
     if (trim(scheme)=="pseudo-herm") then
       call HAYVEC_gather(Vnp1,Vnp1_f,HAYVEC_size) 
     endif
     !
     call io_control(ACTION=OP_WR_CL,SEC=(/1,2/),ID=ID)
     io_err=io_BSS_Haydock(ID,iq,it,reached_threshold,mode,Af(:it),Bf(:it+1),Cf(:Cf_size),Vnm1_f,Vn_f,Vnp1_f)
     !
     YAMBO_FREE_P(Vn_f)
     YAMBO_FREE_P(Vnm1_f)
     if (trim(scheme)=="pseudo-herm") then
       YAMBO_FREE_P(Vnp1_f)
     endif
     !
   endif
   !
   if (Haydock_n_converged==Haydock_required_converged) then
     !
     call K_output_file(iq,"Haydock final "//trim(mode))
     !
     call msg('r','Accuracy (reached)        ',reached_threshold,"[o/o]")
     call msg('r','Iterations                ',it)
     exit
     !
   endif
   !
 enddo
 !
 ! CLEAN
 !=======
 !
 call BS_HAYVEC_free(Vn)
 call BS_HAYVEC_free(Vnp1)
 call BS_HAYVEC_free(Vnm1)
 deallocate(Vn,Vnp1,Vnm1)
 !
 call BS_HAYVEC_free(V0)
 deallocate(V0)
 if (use_left_and_right)  then
   call BS_HAYVEC_free(W0)
   deallocate(W0)
 endif
 if (use_left_and_right.or.trim(scheme)=="pseudo-herm")  then
   YAMBO_FREE(Cf)
 end if
 !
 call PARALLEL_Haydock_VEC_COMMs('reset')
 !
 call timing('Haydock Solver',OPR='stop')
 !
end subroutine
!
!
subroutine HAYVEC_gather(V,fullV,fullV_size)
 !
 ! Gather the Haydock vector fragments into one 
 ! large vector of full dimension  
 ! on the master for I/O 
 ! 
 use pars,          ONLY:SP,cZERO
 use BS,            ONLY:BS_T_grp,BS_nT_grps
 use BS_solvers,    ONLY:Haydock_vector
 use parallel_int,  ONLY:PP_redux_wait
 use parallel_m,    ONLY:PAR_COM_Xk_ibz_INDEX,PAR_IND_T_groups, &
                         PAR_COM_eh_INDEX
 !
 implicit none
 !
 type(Haydock_vector),intent(in) :: V(BS_nT_grps)
 integer,    intent(in)  :: fullV_size
 complex(SP),intent(out) :: fullV(fullV_size)
 ! ws: 
 integer                 :: istart_g,istop_g,i_g
 !
 fullV = cZERO
 istop_g = 0
 !
 ! In this loop I need T_groups indexes since I do not want duplicates 
 !
 do i_g=1,BS_nT_grps
   istart_g = istop_g + 1 
   istop_g  = istart_g + BS_T_grp(i_g)%size - 1
   if (.not.PAR_IND_T_groups%element_1D(i_g)) cycle      
   fullV(istart_g:istop_g) = V(i_g)%fragment(1:BS_T_grp(i_g)%size)
 enddo
 !
 call PP_redux_wait( fullV ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
 call PP_redux_wait( fullV ,COMM=PAR_COM_eh_INDEX%COMM )
 !
end subroutine HAYVEC_gather  
!
!
subroutine HAYVEC_scatter(V,fullV,fullV_size)
 !
 ! Scatter the full dimension vector on the master (read from disk/restart)
 ! into Haydock vector fragments 
 ! 
 use pars,          ONLY:SP
 use BS,            ONLY:BS_T_grp,BS_nT_grps
 use BS_solvers,    ONLY:Haydock_vector
 use parallel_int,  ONLY:PP_bcast
 use parallel_m,    ONLY:PAR_COM_Xk_ibz_INDEX,PAR_COM_eh_INDEX,PAR_IND_T_Haydock
 !
 implicit none
 !
 integer,    intent(in) :: fullV_size
 complex(SP),intent(in) :: fullV(fullV_size)
 !
 type(Haydock_vector),intent(inout) :: V(BS_nT_grps)
 ! ws: 
 integer                 :: istart_g,istop_g,i_g
 !
 ! In this loop I need T_groups Haydock since I do need duplicates 
 !
 call PP_bcast(fullV,0,COMM=PAR_COM_Xk_ibz_INDEX%COMM)  
 call PP_bcast(fullV,0,COMM=PAR_COM_eh_INDEX%COMM)  
 istop_g = 0
 !
 do i_g=1,BS_nT_grps
   istart_g = istop_g + 1 
   istop_g  = istart_g + BS_T_grp(i_g)%size - 1
   if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
   V(i_g)%fragment(1:BS_T_grp(i_g)%size) = fullV(istart_g:istop_g)  
 enddo
 !
end subroutine HAYVEC_scatter
