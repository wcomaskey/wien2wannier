!!! wien2wannier/SRC_w2w/modules_rc.F
!!!
!!!    Modules for wien2wannier.  This file contains modules that
!!!    depend on real/complex compilation.
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2016 Elias Assmann

!!/=== Defined here: =============================
!!
!! xa3: init_xa3(), BKX(:),BKY(:),BKZ(:), GX(:,:),GY(:,:),GZ(:,:),
!!      XK(:),YK(:),ZK(:), vecsz(:)
!!
!! Procedure modules:
!!
!!    almgen_m, l2Mmn_m, l2Amn_m, planew_m, read_vec_m
!!
!!\===============================================

module xa3
  use const, only: R8
#ifdef _COMPLEX_
  use const, only: C16
#endif

  implicit none
  private; save
  public :: init_xa3, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, vecsz, a

#ifndef _COMPLEX_
  real(R8), &
#else
  complex(C16), &
#endif
  allocatable :: a(:,:,:)

  integer,  allocatable :: GX(:,:),GY(:,:),GZ(:,:), vecsz(:)
  real(R8), allocatable :: BKX(:), BKY(:), BKZ(:)
  real(R8), allocatable :: XK(:),  YK(:),  ZK(:)

contains
  subroutine init_xa3(nb, nmat, num_kpts)
    integer, intent(in) :: NB, Nmat, num_kpts

    allocate( a(Nmat, Nb, num_kpts), vecsz(num_kpts) )

    allocate( BKX(Nmat),         BKY(Nmat),         BKZ(Nmat),          &
         &     XK(     num_kpts), YK(     num_kpts), ZK(     num_kpts), &
         &     GX(Nmat,num_kpts), GY(Nmat,num_kpts), GZ(Nmat,num_kpts))

    GX = 0; GY = 0; GZ = 0
  end subroutine init_xa3
end module xa3


!---------------------------  Procedure modules  ---------------------------
module     almgen_m; contains
subroutine almgen(Alm, stru, jatom, lfirst, NB, kkk)
  use param,     only: lomax, NRF, Lmax2
  use w2w,       only: Nmat, iBlock
  use xa3,       only: vecsz, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, A
  use xa,        only: fj, dfj, bk, phs
  use lolog,     only: nlo,nlov,nlon,ilo,lapw,n_rad
  use atspdt,    only: P, DP
  use loabc,     only: alo
  use const,     only: R8, C16, TAU
  use structmod, only: struct_t

  !! procedure includes
  use harmon_m
  use Ylm_m

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(in)  :: jatom, lfirst, NB, kkk
  complex(C16),   intent(out) :: Alm(Nb, Nrf, (LMAX2+1)**2, stru%mult(jatom))

  real(R8) :: h_al(iblock), h_bl(iblock)

  complex(C16) :: PHSHEL,CFAC
  complex(C16) :: YL    ((LMAX2+1)**2)
  complex(C16) :: h_yl  ((LMAX2+1)**2, iblock)
  complex(C16) :: h_alyl((LMAX2+1)**2, iblock)
  complex(C16) :: h_blyl((LMAX2+1)**2, iblock)
  complex(C16) :: ALMt  ((LMAX2+1)**2, Nb, Nrf, stru%mult(JATOM))

  ! Precomputed conjugate spherical harmonics for each plane wave
  complex(C16), allocatable :: YL_all(:,:)
  integer  :: Npw

  real(R8) :: fac, arg, arg1, arg2, arg3
  integer  :: N, num, i,ii,i3, latom, mu, index, l, max, m,m1, irf,mirf, ibb, lda,ldb,ldc, jlo,jneq

  FAC=2*TAU*stru%RMT(JATOM)**2/sqrt(stru%vol)

  N=vecsz(kkk)
  ALMt=0
  do I=1,N
     BKX(I)=XK(kkk)+GX(I,kkk)                    ! plane wave Q+k
     BKY(I)=YK(kkk)+GY(I,kkk)
     BKZ(I)=ZK(kkk)+GZ(I,kkk)
  end do

  ! calculate Bessel functions j_l(|Q+k|*RMT)
  call harmon(stru, jatom, N, BKX, BKY, BKZ ,Lmax2, FJ, DFJ)

  ! Precompute conjugate spherical harmonics for all LAPW plane waves.
  ! YLM depends only on the direction of Q+k (via conv_rec), not on
  ! atom position, so this can be done once before the mu loop.
  Npw = N-(nlo+nlon+nlov)
  allocate(YL_all((LMAX2+1)**2, Npw))
  do I=1,Npw
     BK = matmul(stru%conv_rec, (/BKX(I), BKY(I), BKZ(I)/))
     call YLM(BK, LMAX2, YL)
     YL_all(:,I) = conjg(YL)
  end do

  do MU=1, stru%mult(jatom)
     latom=lfirst-1+mu

     ! matching of plane wave and atomic functions
     do ii=1,Npw,iblock
        i3=0
        do i=ii,min(ii+iblock-1,Npw)
           i3=i3+1

           ! Phase factor exp(i(Q+k).R_j) depends on atom position
           ARG1=BKX(I) * stru%POS(1,LATOM)*TAU
           ARG2=BKY(I) * stru%POS(2,LATOM)*TAU
           ARG3=BKZ(I) * stru%POS(3,LATOM)*TAU

           PHSHEL=exp((0,1)*(ARG1+ARG2+ARG3)) ! exp(i(Q+k)R_j)
           ! Combine precomputed conjg(YLM) with phase factor
           h_yl(:,i3) = YL_all(:,I) * phshel
        end do

        INDEX=0
        do L=0,LMAX2
           i3=0
           do i=ii,min(ii+iblock-1,Npw)
              i3=i3+1
              if (lapw(l)) then
                 h_AL(i3)=DFJ(L,I)*P(L,2)-FJ(L,I)*DP(L,2)
                 h_BL(i3)=FJ(L,I)*DP(L,1)-DFJ(L,I)*P(L,1)
              else
                 h_AL(i3)=FJ(L,I)/P(L,1)/stru%RMT(JATOM)**2
                 h_BL(i3) = 0.d0
              end if
           end do

           MAX=2*L+1
           do M=1,MAX
              INDEX=INDEX+1
              i3=0
              do i=ii,min(ii+iblock-1,Npw)
                 i3=i3+1
                 h_alyl(index,i3)=h_AL(i3)*h_YL(INDEX,i3)
                 h_blyl(index,i3)=h_BL(i3)*h_YL(INDEX,i3)
              end do
           end do
        end do

        ibb = min(iblock, Npw - ii + 1)
#ifndef _COMPLEX_
        lda=2*(LMAX2+1)**2
#else
        lda=(LMAX2+1)**2
#endif
        ldc=lda
        ldb=nmat
#ifndef _COMPLEX_
        call dgemm('N', 'N', 2*index, Nb, ibb,           &
             &     1._R8, h_alyl, lda, a(ii,1,kkk), ldb, &
             &     1._R8, almt(1,1,1,mu), ldc)
        call dgemm('N', 'N', 2*index, Nb, ibb,           &
             &     1._R8, h_blyl, lda, a(ii,1,kkk), ldb, &
             &     1._R8, almt(1,1,2,mu), ldc)
#else
        call zgemm('N', 'N', index, Nb, ibb,                  &
             &     (1._R8, 0), h_alyl, lda, a(ii,1,kkk), ldb, &
             &     (1._R8, 0), almt(1,1,1,mu), ldc)
        call zgemm('N', 'N', index, Nb, ibb,                  &
             &     (1._R8, 0), h_blyl, lda, a(ii,1,kkk), ldb, &
             &     (1._R8, 0), almt(1,1,2,mu), ldc)
#endif
     end do

     if (nlo.ne.0) then
        i=n-(nlo+nlon)
        do L=0,LoMAX
           mirf=n_rad(l)
           do jlo=1,ilo(l)
              do jneq=1,stru%mult(jatom)
                 do M1=-l,l
                    i=i+1
                    BK = matmul(stru%conv_rec, (/BKX(I), BKY(I), BKZ(I)/))

                    call YLM (BK,LOMAX,YL)
                    ARG = dot_product((/BKX(I), BKY(I), BKZ(I)/), &
                         &            stru%pos(:,latom))
                    PHSHEL=exp((0,1)*(ARG*TAU))
                    do NUM=1,NB
                       PHS(NUM)=PHSHEL*A(I,NUM,kkk)
                    enddo
                    do M=-l,l
                       index=l*(l+1)+m+1
                       do NUM=1,NB
                          do irf=1,mirf
                             ALMt(index,num,irf,mu)=ALMt(INDEX,num,irf,mu)+ &
                                  ALo(l,jlo,irf)*dconjg(YL(INDEX))*PHS(NUM)
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end if

     INDEX=0
     do L=0,LMAX2
        mirf=n_rad(l)
        MAX=2*L+1
        CFAC=(0,1)**L
        do M=1,MAX
           INDEX=INDEX+1
           do irf=1,mirf
              do NUM=1,NB
                 ALM(num,irf,index,mu)=ALMt(INDEX,NUM,irf,mu)*FAC*CFAC
              end do
           end do
        end do
     end do
  end do
  deallocate(YL_all)
end subroutine almgen
end module     almgen_m


module     l2Mmn_m; contains
subroutine l2Mmn(stru, inwf, num_kpts, NNTOT)
  use const,     only: R8, C16, TAU
  use param,     only: unit_out, unit_vsp, unit_vector, NRF, Lmax2
  use xa,        only: bk
  use xa3,       only: xk, yk, zk
  use bessel,    only: ri_mat
  use Amn_Mmn,   only: overlap
  use pairs,     only: kp, kpb, bqx, bqy, bqz
  use lolog,     only: n_rad
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t

  !! procedure includes
  use almgen_m
  use gaunt_cache, only: get_gaunt
  use radint_m
  use atpar_m
  use Ylm_m

  implicit none

  type(struct_t), intent(in) :: stru
  type(inwf_t),   intent(in) :: inwf
  integer,        intent(in) :: num_kpts, NNtot

  integer  :: pair_index, r_index, iscf, jatom, lfirst, itap, itape, jtape
  integer  :: k1, k2, kkk, Nb
  real(R8) :: tALM, tMeas1, tMeas2, t1, tt0, tt1, fac, bx, by, bz
  real(R8) :: KX1,KY1,KZ1     ! k-points in u.c. coordinates k and k+b

  complex(C16) :: YLB((LMAX2+1)**2) ! spherical harmonics expansion of b
  complex(C16) :: PHSHEL, tmp, tmp1

  complex(C16), allocatable :: alm(:,:,:,:),blm(:,:,:,:)

  real(R8) :: BM, tt2,tt3, t2, arg1, arg2, arg3
  integer  :: index1,index2,indexj, irf1,irf2, lj, m1,m2,mj, mu
  integer  :: latom, l_index, l1, l2

  ! --- Batched zgemm variables (Opt 3: zgeru → zgemm) ---
  ! Column map: col_a(irf, index) → flat column 1..Na for active (irf, L,M)
  integer  :: col_map(NRF, (LMAX2+1)**2)
  integer  :: Na          ! number of active columns
  integer  :: ca, cb, ll  ! column indices for ALM/BLM
  ! Coefficient matrix: S_base(col_b, col_a) accumulates ri_mat*YLB*gaunt
  ! (without PHSHEL, which is factored out per mu)
  complex(C16), allocatable :: S_base(:,:)
  ! Gathered 2D views and workspace for zgemm
  complex(C16), allocatable :: ALM_2d(:,:), BLM_2d(:,:), W(:,:)

  integer :: k1_prog_itvl
  k1_prog_itvl = min(max(num_kpts/10, 1), 100)

  Nb = inwf%bmax - inwf%bmin + 1

  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST=1
  atoms: do jatom = 1, stru%Nneq
     write(unit_out, "(/, '===== atom', I5, ' /', I5, ' =====' /)") &
          jatom, stru%Nneq

     allocate(Alm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)),  &
          &   Blm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)))

     talm=0.
     tmeas1=0.
     tmeas2=0.
     call cpu_time(t1)

     if(jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     ITAP=30+JATOM
     itape=unit_vector
     jtape=unit_vsp
     rewind(itape)
     call atpar(stru, JATOM, itape, jtape)
     fac=2*TAU * stru%RMT(jatom)**2 / sqrt(stru%vol)
     rewind(itape)

     ! Build column map: enumerate active (irf, L, M) combinations
     ! col_map(irf, index) = flat column number 1..Na, or 0 if inactive
     col_map = 0
     Na = 0
     do ll = 0, LMAX2
        do m1 = -ll, ll
           index1 = ll*(ll+1) + m1 + 1
           do irf1 = 1, n_rad(ll)
              Na = Na + 1
              col_map(irf1, index1) = Na
           end do
        end do
     end do
     write(unit_out, "(' Batched zgemm: Na=', I4, ' active columns')") Na

     ! Allocate work arrays for batched zgemm
     allocate(S_base(Na, Na))
     allocate(ALM_2d(Nb, Na), BLM_2d(Nb, Na), W(Nb, Na))

     pair_index=0
     k1loop: do k1=1,num_kpts      !   loop over k1
        pair_index=pair_index+1
        kkk=KP(pair_index)
        KX1=XK(kkk)
        KY1=YK(kkk)
        KZ1=ZK(kkk)

        call almgen(Alm, stru, jatom, lfirst, Nb, kkk)

        k2loop: do k2=1,NNTOT
           if (k2 > 1) pair_index=pair_index+1
           kkk=KPB(pair_index)
           call cpu_time(tt0)
           call almgen(Blm, stru, jatom, lfirst, Nb, kkk)
           BLM=conjg(BLM)

           call cpu_time(tt1)

           ! calculate b=k2-k1 add BQ if going around BZ
           BX=XK(kkk)-KX1+BQX(pair_index)
           BY=YK(kkk)-KY1+BQY(pair_index)
           BZ=ZK(kkk)-KZ1+BQZ(pair_index)

           ! transform to cartesian coordinates
           BK = matmul(stru%conv_rec, (/BX, BY, BZ/))
           BM = sqrt(sum(BK**2))

           ! compute radial intergrals <R(r)|j_(|b|*r)|R'(r)>
           call radint(stru, jatom, inwf%LJmax, bm)
           ! computer Y_lm(b)
           call YLM (BK, inwf%LJmax, YLB)

           indexj=0
           do LJ = 0, inwf%LJmax
              do MJ=-LJ,LJ
                 indexj=indexj+1
                 YLB(indexj)=conjg(YLB(indexj))*(0,1)**LJ
              enddo
           enddo
           call cpu_time(tt2)

           ! --- Opt 3: Build S_base coefficient matrix (without PHSHEL) ---
           ! S_base(col_b, col_a) = sum over (LJ,MJ) of
           !   ri_mat(R_index,L_index) * YLB(indexj) * gaunt(L2,LJ,L1,M2,MJ,M1)
           ! where col_a = col_map(irf1, index1), col_b = col_map(irf2, index2).
           ! PHSHEL is factored out (varies per mu only).
           S_base = (0.0_R8, 0.0_R8)
           L_index=0
           do L1=0,LMAX2
              do L2=0,LMAX2
                 do LJ = abs(L1-L2), min(L1+L2, inwf%LJmax), 2
                    if (mod((L1+L2+LJ),2) == 1) cycle
                    if ((L1+L2-LJ) < 0         &
                         .or. (L1-L2+LJ) < 0   &
                         .or. (-L1+L2+LJ) < 0) &
                         & cycle
                    L_index=L_index+1
                    do M1=-L1,L1
                       index1=L1*(L1+1)+M1+1
                       do MJ = max(-LJ, -L2-M1), min(LJ, L2-M1)
                          M2=M1+MJ
                          index2=L2*(L2+1)+M2+1
                          indexj=LJ*(LJ+1)+MJ+1
                          tmp=YLB(indexj)*get_gaunt(L2,LJ,L1,M2,MJ,M1)
                          R_index=0
                          do irf1=1,n_rad(L1)
                             do irf2=1,n_rad(L2)
                                R_index=R_index+1
                                ca = col_map(irf1, index1)
                                cb = col_map(irf2, index2)
                                S_base(cb, ca) = S_base(cb, ca) &
                                     + ri_mat(R_index,L_index) * tmp
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo

           ! --- Per mu: scale by PHSHEL, gather ALM/BLM, two zgemms ---
           muloop: do mu=1, stru%mult(jatom)
              latom=lfirst-1+mu
              arg1=Bx * stru%pos(1,latom)*TAU
              arg2=By * stru%pos(2,latom)*TAU
              arg3=Bz * stru%pos(3,latom)*TAU
              PHSHEL=exp((0,1)*(ARG1+ARG2+ARG3))*2*TAU

              ! Gather ALM and BLM into contiguous 2D arrays (Nb x Na)
              do ll = 0, LMAX2
                 do m1 = -ll, ll
                    index1 = ll*(ll+1) + m1 + 1
                    do irf1 = 1, n_rad(ll)
                       ca = col_map(irf1, index1)
                       ALM_2d(:, ca) = ALM(:, irf1, index1, mu)
                       BLM_2d(:, ca) = BLM(:, irf1, index1, mu)
                    end do
                 end do
              end do

              ! W(Nb, Na) = BLM_2d(Nb, Na) * (PHSHEL * S_base(Na, Na))
              ! overlap += W * ALM_2d^T
              ! zgemm 1: W = BLM_2d * S_base, scaled by PHSHEL (alpha)
              call zgemm('N', 'N', Nb, Na, Na, PHSHEL, &
                   BLM_2d, Nb, S_base, Na, &
                   (0.0_R8, 0.0_R8), W, Nb)
              ! zgemm 2: overlap += W * ALM_2d^T
              call zgemm('N', 'T', Nb, Nb, Na, (1.0_R8, 0.0_R8), &
                   W, Nb, ALM_2d, Nb, &
                   (1.0_R8, 0.0_R8), overlap(1,1,pair_index), &
                   size(overlap,1))
           enddo muloop
           call cpu_time(tt3)
           talm=talm+tt1-tt0
           tmeas1=tmeas1+tt2-tt1
           tmeas2=tmeas2+tt3-tt2
        end do k2loop

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
     end do k1loop
     deallocate(ALM, BLM, S_base, ALM_2d, BLM_2d, W)

     call cpu_time(t2)
     talm=talm*(NNTOT+1)/NNTOT
     write(unit_out,*) 'CPU time used for atom ',JATOM,' ='&
          ,t2-t1,talm,tmeas1,tmeas2

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
  end do atoms
end subroutine l2MMN
end module     l2Mmn_m


module     l2Amn_m; contains
subroutine l2Amn(stru, inwf, num_kpts)
  use param,     only: Lmax2, NRF, lomax, unit_vsp, unit_out, unit_vector
  use w2w,       only: Nmat, iblock, unit_amn
  use const,     only: R8, C16, TAU
  use xa,        only: BK, FJ, DFJ, PHS, BKrot, BKRloc
  use xa3,       only: A, XK,YK,ZK, GX,GY,GZ, BKX,BKY,BKZ, vecsz
  use lolog,     only: nlo,nlov,nlon,loor,ilo,lapw
  use loabc,     only: alo
  use atspdt,    only: P, DP
  use const,     only: R8, C16
  use radfu,     only: RF1, RF2
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t

  !! procedure includes
  use rint13_m
  use harmon_m
  use atpar_m
  use Ylm_m

  implicit none

  type(struct_t), intent(in) :: stru
  type(inwf_t),   intent(in) :: inwf
  integer,        intent(in) :: num_kpts

  complex(C16) :: Yl((Lmax2+1)**2)        ! list of spherical harmics
  complex(C16) :: projection(inwf%bmax-inwf%bmin+1, inwf%Nproj, num_kpts)
  complex(C16) :: phshel,cfac

  real(R8)     :: h_al(iblock), h_bl(iblock)
  complex(C16) :: h_yl  ( (Lmax2+1)**2, iblock)
  complex(C16) :: h_alyl( (Lmax2+1)**2, iblock)
  complex(C16) :: h_blyl( (Lmax2+1)**2, iblock)
  complex(C16) :: Alm   ( (Lmax2+1)**2, inwf%bmax-inwf%bmin+1, stru%Nat, Nrf)

  integer  :: i,ii,i3,l, m,mu,m1, n, ipr, ib, iY
  integer  :: jatom,jlo,jneq, idx, ibb,iscf,irf, lda,ldb,ldc, latom,lfirst
  integer  :: kkk, max, itap,itape,jtape, Nb
  real(R8) :: fac, arg, pi12lo(0:Lmax2)

  Nb = inwf%bmax-inwf%bmin+1
  projection = 0

  rewind(unit_vsp)
  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST=1
  atoms: do JATOM = 1, stru%Nneq
     if(jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     ITAP=30+JATOM
     itape=unit_vector
     jtape=unit_vsp
     rewind(itape)
     call atpar(stru, jatom, itape, jtape)
     rewind(itape)
     write(unit_out, *)
     write(unit_out, '(7X,  "rotij", 10X, "transij")')
     write(unit_out, '(6X, 23("-"))')

     pi12lo = 0
     do l = 0,LOmax
        if (loor(l)) call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
             &                   rf1(:,l,3), rf2(:,l,3), pi12lo(l))
     end do

     kpoints: do KKK=1,NUM_KPTS
        N=vecsz(kkk)

        ! plane wave Q+k
        BKX(1:N) = XK(kkk) + GX(1:N, kkk)
        BKY(1:N) = YK(kkk) + GY(1:N, kkk)
        BKZ(1:N) = ZK(kkk) + GZ(1:N, kkk)

        ! calculate Bessel functions j_l(|Q+k|*RMT)
        call harmon(stru, jatom, N, Bkx, Bky, Bkz, Lmax2, FJ, DFJ)

        fac = 2*TAU * stru%RMT(jatom)**2 / sqrt(stru%vol)
        Alm = 0

        equiv1: do mu = 1, stru%mult(jatom) !  loop over equivalent atoms
           latom=lfirst-1+mu
           if (kkk == 1) then
              do i=1,3
                 write(unit_out, '(4x,3I3,4x,f12.4)') &
                      stru%rotij(i,:,latom), stru%transij(i,latom)
              enddo
              write(unit_out,*)
           endif

           ! matching of plane wave and atomic functions
           do ii=1,N-(nlo+nlon+nlov),iblock
              i3=0
              do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                 i3=i3+1

                 ! go to coord. of latom
                 bkrot = matmul(stru%rotij(:,:,latom), &
                      &         (/bkx(i), bky(i), bkz(i)/))

                 ! transform Q+k to cartesian coordinates
                 bk = matmul(stru%conv_rec, bkrot)

                 ! go to local (rotated) coords
                 bkrloc = matmul(stru%locrot(:,:,jatom), bk)

                 ! calculate Y_lm(Q+k)
                 call YLM (bkrloc,LMAX2,YL)
                 arg =  dot_product(bkrot, stru%pos(:, lfirst)) &
                      + dot_product((/bkx(i), bky(i), bkz(i)/), &
                      &             stru%transij(:,latom))
                 phshel=exp((0,1)*arg*TAU)

                 idx=0
                 do l=0,Lmax2
                    do m=1, 2*l+1
                       idx=idx+1
                       h_yl(idx,i3)=conjg(yl(idx))*phshel
                    enddo
                 enddo
              enddo

              IDX=0
              do L=0,LMAX2
                 i3=0
                 do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                    i3=i3+1
                    if (lapw(l)) then
                       h_AL(i3)=DFJ(L,I)*P(L,2)-FJ(L,I)*DP(L,2)
                       h_BL(i3)=FJ(L,I)*DP(L,1)-DFJ(L,I)*P(L,1)
                    else
                       h_AL(i3)=FJ(L,I)/P(L,1)/stru%RMT(jatom)**2
                       h_BL(i3) = 0
                    endif
                 enddo
                 MAX=2*L+1
                 do M=1,MAX
                    IDX=IDX+1
                    i3=0
                    do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                       i3=i3+1
                       h_alyl(idx,i3)=h_AL(i3)*h_YL(idx,i3)
                       h_blyl(idx,i3)=h_BL(i3)*h_YL(idx,i3)
                    enddo
                 enddo
              enddo
              ibb=min(iblock,N-(nlo+nlon+nlov)-ii+1)
#ifndef _COMPLEX_
              lda=2*(LMAX2+1)*(LMAX2+1)
#else
              lda=(LMAX2+1)*(LMAX2+1)
#endif
              ldc=lda
              ldb=nmat

#ifndef _COMPLEX_
              call dgemm('N','N',2*idx,nb,ibb,1.d0, &
                   h_alyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu,1),ldc)
              call dgemm('N','N',2*idx,nb,ibb,1.d0, &
                   h_blyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu,2),ldc)
#else
              call zgemm('N','N',idx,nb,ibb,(1.d0,0.d0), &
                   h_alyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu,1),ldc)
              call zgemm('N','N',idx,nb,ibb,(1.d0,0.d0), &
                   h_blyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu,2),ldc)
#endif
           enddo
!********  LOCAL ORBITALS *****************************
           lo: if (nlo.ne.0) then
              i=n-(nlo+nlon)
              do L=0,LoMAX
                 do jlo=1,ilo(l)
                    do jneq=1,stru%mult(jatom)
                       do M1=-l,l
                          i=i+1
                          BK(1)=BKX(I)
                          BK(2)=BKY(I)
                          BK(3)=BKZ(I)
                          bkrot = matmul(stru%rotij(:,:,latom), BK)
                          BK    = matmul(stru%conv_rec, BKROT)

                          bkrloc = matmul(stru%locrot(:,:,jatom), BK)
                          call Ylm (Bkrloc, LOmax, YL)
                          arg =  dot_product(bkrot(:), stru%pos(:, lfirst)) &
                               + dot_product((/bkx(i), bky(i), bkz(i)/),    &
                               &             stru%transij(:,latom))

                          phshel=exp((0,1)*TAU*arg)
                          PHS(:) = phshel * A(i,:,kkk)

                          do m=-l,l
                             idx=l*(l+1)+m+1
                             do ib=1,Nb
                                Alm(idx, ib, mu, :)=Alm(idx, ib, mu, :) + &
                                     Alo(l, jlo, :) * conjg(Yl(idx))*phs(ib)
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end if lo
!*********************************************************
           idx=0
           do l=0,Lmax2
              cfac = (0,1)**l
              do m=1, 2*l+1
                 idx=idx+1
                 do ib=1,Nb
                    do irf=1,Nrf
                       Alm(idx,ib,mu,irf)=Alm(idx,ib,mu,irf)*fac*cfac
                    enddo
                 enddo
              enddo
           enddo
        enddo equiv1

        equiv2: do mu=1,stru%mult(jatom)
           latom = lfirst - 1 + mu

           do ib = 1,Nb
              do ipr = 1, inwf%Nproj
                 do iY = 1, inwf%projections(ipr)%NY
                    if (inwf%projections(ipr)%iat(iY) /= latom) &
                         cycle
                    l = inwf%projections(ipr)%l(iY)
                    m = inwf%projections(ipr)%m(iY)
                    ! idx is the compound index to Alm
                    idx = l**2 + m + l + 1

                    projection(ib, ipr, kkk) = projection(ib, ipr, kkk) &
                         & + conjg(Alm(idx, ib, mu, 1))                 &
                         &   * inwf%projections(ipr)%coeff(iY)          &
                         & + conjg(Alm(idx, ib, mu, 3))                 &
                         &   * inwf%projections(ipr)%coeff(iY)          &
                         &   * pi12lo(l)
                 enddo
              enddo
           enddo
        enddo equiv2
     enddo kpoints
  enddo atoms

  do kkk = 1,num_kpts
     do ipr = 1, inwf%Nproj
        do ib = 1,Nb
           write(unit_amn, '(2i4,1x,i5,1x,2e18.5)') &
                ib, ipr, kkk, projection(ib, ipr, KKK)
        enddo
     enddo
  enddo

  return
end subroutine l2Amn
end module     l2Amn_m


module     planew_m; contains
subroutine planew(stru, NB, num_kpts, Nntot, maxx, maxy, maxz)
  use param,     only: unit_out
  use w2w,       only: unit_mmn
  use xa3,       only: gx,gy,gz, a, vecsz
  use amn_mmn,   only: overlap
  use pairs,     only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1
  use lolog,     only: nlo,nlov,nlon
  use const,     only: R8, C16, TAU
  use structmod, only: struct_t

  implicit none

  type(struct_t), intent(in) :: stru
  integer,        intent(in) :: NB, num_kpts, NNTOT, maxx,maxy,maxz

  integer  :: index, pair_index, m2x,m2y,m2z, max, i,j,k, ii,jj,kk
  integer  :: i_x,i_y,i_z, k1,k2, kkk, lfirst, jatom, latom, mu
  integer  :: num1,num2, id, ix1,ix2,iy1,iy2,iz1,iz2
  integer  :: nvalid
  real(R8) :: t1,t2, tt1,tt2,tt3, t_b, t_prod, vv, bx,by,bz, x, arg

  logical, allocatable, dimension(:,:) :: a_null, a_null_fin
  complex(C16)       sum(NB),SUMA(NB,NB),FF
  complex(c16), allocatable, dimension(:,:,:) :: FAC
  complex(C16)       PHSHEL
  complex(C16), allocatable, dimension(:,:,:) :: A_
  complex(C16), allocatable, dimension(:,:) :: B_
  ! Workspace for gathering valid columns for BLAS zgemm in Step 2
  complex(C16), allocatable, dimension(:,:) :: Agath, Bgath
  ! Sparse compressed index lists for Step 1 loop acceleration
  integer, allocatable :: idx_null(:,:), idx_fin(:,:)
  integer, allocatable :: idx_null_ix(:,:), idx_null_iy(:,:), idx_null_iz(:,:)
  integer, allocatable :: idx_fin_ix(:,:), idx_fin_iy(:,:), idx_fin_iz(:,:)
  integer, allocatable :: n_null(:), n_fin(:)
  integer :: ni, nj
#ifdef HAVE_FFTW
  include 'fftw3.f'
  integer  :: Nfx, Nfy, Nfz, Nfft_total
  integer*8 :: plan_fwd, plan_bwd
  complex(C16), allocatable :: FAC_hat(:,:,:)
  complex(C16), allocatable :: work3d(:,:,:)
  real(R8) :: inv_Nfft
  integer  :: band, wrap_ix, wrap_iy, wrap_iz
  logical  :: use_fft
#endif

  allocate(FAC(-(2*maxx+1):(2*maxx+1),-(2*maxy+1):(2*maxy+1),-(2*maxz+1):(2*maxz+1)) )
  allocate(A_(NB,0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
  allocate(B_(NB,0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1) )
  allocate(a_null(0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
  allocate(a_null_fin(0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
  allocate(Agath(NB, (2*maxx+1)*(2*maxy+1)*(2*maxz+1)) )
  allocate(Bgath(NB, (2*maxx+1)*(2*maxy+1)*(2*maxz+1)) )
  allocate(idx_null((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_fin((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_null_ix((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_null_iy((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_null_iz((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_fin_ix((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_fin_iy((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(idx_fin_iz((2*maxx+1)*(2*maxy+1)*(2*maxz+1), NUM_KPTS))
  allocate(n_null(NUM_KPTS), n_fin(NUM_KPTS))

  m2x=2*maxx+1
  m2y=2*maxy+1
  m2z=2*maxz+1
  max=m2x*m2y*m2z-1
  write(unit_out,*)'MAX Q:',maxx,maxy,maxz
  t_b=0.d0
  t_prod=0.d0

#ifdef HAVE_FFTW
  ! FFT grid: 2*m2x per dim guarantees no circular aliasing
  ! (kernel support [-(m2x-1), m2x-1] has width 2*m2x-1 <= 2*m2x = Nfx)
  Nfx = 2 * m2x
  Nfy = 2 * m2y
  Nfz = 2 * m2z
  Nfft_total = Nfx * Nfy * Nfz
  inv_Nfft = 1.0_R8 / Nfft_total
  allocate(FAC_hat(Nfx, Nfy, Nfz), work3d(Nfx, Nfy, Nfz))
  ! Create FFTW plans (ESTIMATE: no destructive benchmarking)
  call dfftw_plan_dft_3d(plan_fwd, Nfx, Nfy, Nfz, &
       work3d, work3d, FFTW_FORWARD, FFTW_ESTIMATE)
  call dfftw_plan_dft_3d(plan_bwd, Nfx, Nfy, Nfz, &
       work3d, work3d, FFTW_BACKWARD, FFTW_ESTIMATE)
  write(unit_out,*) 'FFT grid:', Nfx, Nfy, Nfz, ' total=', Nfft_total
#endif

  call cpu_time(t1)
  A_=0
  a_null=.false.
  a_null_fin=.false.
  pair_index=0
  do k1=1,num_kpts
     pair_index=pair_index+1
     kkk=KP(pair_index)
     do I=1,vecsz(kkk)-(nlo+nlon+nlov)
        I_X=GX(I,kkk)+maxx
        I_Y=GY(I,kkk)+maxy
        I_Z=GZ(I,kkk)+maxz
        index=(I_Z*m2y+I_Y)*m2x+I_X
        a_null(index,kkk)=.true.

        A_(:,index,kkk)=A(I,:,kkk)
     end do
     do k2=1,NNTOT
        if (k2 > 1) pair_index=pair_index+1
        kk=KPB(pair_index)
        do I=1,vecsz(kkk)
           I_X=GX(I,kk)+maxx-BQX(pair_index)
           I_Y=GY(I,kk)+maxy-BQY(pair_index)
           I_Z=GZ(I,kk)+maxz-BQZ(pair_index)
           index=(I_Z*m2y+I_Y)*m2x+I_X
           a_null_fin(index,kkk)=.true.
        end do
     end do
  end do

  ! Build compressed index lists from a_null and a_null_fin.
  ! Reduces Step 1 from O(max^2) to O(n_null * n_fin).
  n_null = 0
  n_fin  = 0
  do kkk=1,num_kpts
     do J=0,max
        if (a_null(J, kkk)) then
           n_null(kkk) = n_null(kkk) + 1
           ni = n_null(kkk)
           idx_null(ni, kkk) = J
           JJ = J/m2x
           idx_null_iz(ni, kkk) = JJ/m2y
           idx_null_iy(ni, kkk) = JJ - idx_null_iz(ni,kkk)*m2y
           idx_null_ix(ni, kkk) = J - JJ*m2x
        end if
        if (a_null_fin(J, kkk)) then
           n_fin(kkk) = n_fin(kkk) + 1
           nj = n_fin(kkk)
           idx_fin(nj, kkk) = J
           JJ = J/m2x
           idx_fin_iz(nj, kkk) = JJ/m2y
           idx_fin_iy(nj, kkk) = JJ - idx_fin_iz(nj,kkk)*m2y
           idx_fin_ix(nj, kkk) = J - JJ*m2x
        end if
     end do
  end do

  FAC=0
  LFIRST=1
  do JATOM=1,stru%Nneq
     if (jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     VV=2*TAU * stru%RMT(jatom)**3/3
     do K=-m2z,m2z
        do J=-m2y,m2y
           do I=-m2x,m2x
              BX = dot_product(stru%conv_rec(1, :), (/ i, j, k /))
              BY = dot_product(stru%conv_rec(2, :), (/ i, j, k /))
              BZ = dot_product(stru%conv_rec(3, :), (/ i, j, k /))

              X = stru%RMT(JATOM) * sqrt(BX**2 + BY**2 + BZ**2)
              phshel=0
              do mu=1,stru%mult(jatom)
                 latom=lfirst+mu-1
                 arg =  dot_product((/i,j,k/), stru%pos(:,LATOM))
                 phshel = phshel + exp(-(0,1)*arg*TAU)
              end do

              if (i==0 .and. j==0 .and. k==0) then
                 fac(i,j,k) = fac(i,j,k) - stru%mult(jatom) * VV/stru%vol
              else
                 fac(i,j,k) = fac(i,j,k) &
                      - 3*phshel * (sin(x) - x*cos(x)) / x**3 * VV/stru%vol
              end if
           end do
        end do
     end do
  end do
  FAC(0,0,0)=FAC(0,0,0)+1

#ifdef HAVE_FFTW
  ! Pre-transform FAC kernel into Fourier space.
  ! Place FAC(dx,dy,dz) into padded grid with circular wraparound:
  !   non-negative d → position d+1 (Fortran 1-based)
  !   negative d     → position Nf+d+1 (wraps to end of array)
  ! Populate only [-(m2x-1), m2x-1] — the actual convolution support.
  FAC_hat = (0.0_R8, 0.0_R8)
  do k = -(m2z-1), m2z-1
     wrap_iz = merge(k + 1, Nfz + k + 1, k >= 0)
     do j = -(m2y-1), m2y-1
        wrap_iy = merge(j + 1, Nfy + j + 1, j >= 0)
        do i = -(m2x-1), m2x-1
           wrap_ix = merge(i + 1, Nfx + i + 1, i >= 0)
           FAC_hat(wrap_ix, wrap_iy, wrap_iz) = FAC(i, j, k)
        end do
     end do
  end do
  call dfftw_execute_dft(plan_fwd, FAC_hat, FAC_hat)
#endif

  call cpu_time(t2)
  write(unit_out,*)'t_init=',t2-t1

  ! read in k-point and basis vectors (original loop structure)
  pair_index=0
  do k1=1,num_kpts
     pair_index=pair_index+1
     kkk=KP(pair_index)
     call cpu_time(tt1)
     ! Step 1: Compute B_ = A_ conv FAC (interstitial overlap kernel).
     b_=0
#ifdef HAVE_FFTW
     ! Runtime crossover: use FFT when sparse convolution is expensive.
     ! The constant 16 accounts for FFT overhead (3 transforms + gather/scatter)
     ! vs the simple multiply-add of the sparse loop.
     use_fft = (int(n_null(kkk),8) * int(n_fin(kkk),8)) > &
               (int(Nfft_total,8) * 16_8)
     if (k1 == 1) write(unit_out,'(A,I6,A,I6,A,I10,A,L1)') &
          ' FFT crossover: n_null=', n_null(kkk), ' n_fin=', n_fin(kkk), &
          ' product=', n_null(kkk)*n_fin(kkk), ' use_fft=', use_fft
     if (use_fft) then
        ! FFT-based convolution: B = IFFT(FFT(A) * FAC_hat) per band
        do band = 1, NB
           ! (a) Scatter A_ for this band into zero-padded 3D grid
           work3d = (0.0_R8, 0.0_R8)
           do ni = 1, n_null(kkk)
              work3d(idx_null_ix(ni,kkk)+1, &
                     idx_null_iy(ni,kkk)+1, &
                     idx_null_iz(ni,kkk)+1) = A_(band, idx_null(ni,kkk), k1)
           end do
           ! (b) Forward FFT
           call dfftw_execute_dft(plan_fwd, work3d, work3d)
           ! (c) Pointwise multiply by pre-transformed kernel
           work3d = work3d * FAC_hat
           ! (d) Inverse FFT
           call dfftw_execute_dft(plan_bwd, work3d, work3d)
           ! (e) Gather B_ at output positions, with 1/N normalization
           !     (FFTW backward transform is unnormalized)
           do nj = 1, n_fin(kkk)
              B_(band, idx_fin(nj,kkk)) = &
                   work3d(idx_fin_ix(nj,kkk)+1, &
                          idx_fin_iy(nj,kkk)+1, &
                          idx_fin_iz(nj,kkk)+1) * inv_Nfft
           end do
        end do
     else
#endif
     ! Sparse convolution using compressed index lists.
     do nj=1,n_fin(kkk)
        J   = idx_fin(nj, kkk)
        ix2 = idx_fin_ix(nj, kkk)
        iy2 = idx_fin_iy(nj, kkk)
        iz2 = idx_fin_iz(nj, kkk)
        SUM=0
        do ni=1,n_null(kkk)
           I   = idx_null(ni, kkk)
           ix1 = idx_null_ix(ni, kkk)
           iy1 = idx_null_iy(ni, kkk)
           iz1 = idx_null_iz(ni, kkk)
           FF=FAC(ix2-ix1,iy2-iy1,iz2-iz1)
           !! Use explicit scalar loop to avoid gfortran
           !! vectorizer miscompilation (-ftree-vectorize,
           !! enabled at -O2 and above) on ARM targets.
           do num1=1,NB
              sum(num1)=sum(num1)+a_(num1,I,k1)*FF
           end do
        end do
        do num1=1,NB
           b_(num1,j)=sum(num1)
        end do
     end do
#ifdef HAVE_FFTW
     end if
#endif
     call cpu_time(tt2)
     do k2=1,NNTOT
        if (k2 > 1) pair_index=pair_index+1
        kkk=KPB(pair_index)

        id=-(BQZ(pair_index)*m2y+BQY(pair_index))*m2x-BQX(pair_index)-1

        ! Gather valid column pairs into contiguous arrays for zgemm.
        ! A_ and B_ are always complex(C16) in planew.
        nvalid = 0
        do J=0,max
           id=id+1
           if (a_null(j,kkk) .and. 0 <= id .and. id <= max) then
              nvalid = nvalid + 1
#ifndef _COMPLEX_
              Agath(:,nvalid) = a_(1:NB,j,kkk)
#else
              Agath(:,nvalid) = conjg(a_(1:NB,j,kkk))
#endif
              Bgath(:,nvalid) = b_(1:NB,id)
           endif
        end do

        ! SUMA = Agath * Bgath^T via zgemm
        if (nvalid > 0) then
           call zgemm('N', 'T', NB, NB, nvalid,           &
                &      (1._R8, 0._R8), Agath, NB, Bgath, NB, &
                &      (1._R8, 0._R8), overlap(1,1,pair_index), size(overlap,1))
        end if

        write(unit_mmn,235)KP(pair_index),KPB(pair_index),BQX1(pair_index),BQY1(pair_index), &
             BQZ1(pair_index)

        do NUM2=1,NB
           do NUM1=1,NB
              write(unit_mmn,236) conjg(overlap(num2,num1,pair_index))
           end do
        end do
     end do
     call cpu_time(tt3)
     t_b=t_b+tt2-tt1
     t_prod=t_prod+tt3-tt2
     if (k1 == 1) write(unit_out,*)'PLANEW TIME PER k-point:',t_b+t_prod
  end do
  write(unit_out,*)'t_b=',t_b
  write(unit_out,*)'t_prod=',t_prod


235 format(5i8)
236 format(32f18.12)

  !release the memory
#ifdef HAVE_FFTW
  call dfftw_destroy_plan(plan_fwd)
  call dfftw_destroy_plan(plan_bwd)
  deallocate(FAC_hat, work3d)
#endif
  deallocate(FAC)
  deallocate(A_)
  deallocate(B_)
  deallocate(a_null)
  deallocate(a_null_fin)
  deallocate(Agath)
  deallocate(Bgath)
  deallocate(idx_null, idx_fin)
  deallocate(idx_null_ix, idx_null_iy, idx_null_iz)
  deallocate(idx_fin_ix, idx_fin_iy, idx_fin_iz)
  deallocate(n_null, n_fin)
end subroutine PLANEW
end module planew_m


module     read_vec_m; contains
!!!    Read a ‘vector’ file
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2015 Elias Assmann
subroutine read_vec(NEmin, NEmax, Nneq, kkk, maxx, maxy, maxz, Efermi)
  use param,  only: unit_out, unit_vector
  use w2w,    only: unit_eig
  use const,  only: R8, Ryd_eV
  use xa3,    only: XK,YK,ZK, GX,GY,GZ, vecsz, A

  implicit none

  integer, intent(in)    :: NEmin, NEmax, Nneq
  integer, intent(out)   :: maxx, maxy, maxz
  integer, intent(inout) :: kkk
  real(R8),intent(in)    :: efermi

  integer :: i, j, n, NB, NE, num, E_size

  real(r8), allocatable :: E(:)

  character(10)    BNAME

  E_size = 0

  do I = 1, Nneq
     read(unit_vector)
     read(unit_vector)
  enddo

  maxx=0; maxy=0; maxz=0
  kpoint: do
     read(unit_vector,end=998) XK(kkk+1),YK(kkk+1),ZK(kkk+1),BNAME,N,NE
     KKK=KKK+1
     vecsz(kkk)=N
     read(unit_vector) (GX(I,kkk),GY(I,kkk),GZ(I,kkk), I=1,N)
     do I=1,N
        if (abs(GX(I,kkk)) > maxx) maxx=abs(GX(I,kkk))
        if (abs(GY(I,kkk)) > maxy) maxy=abs(GY(I,kkk))
        if (abs(GZ(I,kkk)) > maxz) maxz=abs(GZ(I,kkk))
     enddo

     ! Ensure E array is large enough for NE bands
     if (NE > E_size) then
        if (allocated(E)) deallocate(E)
        E_size = NE
        allocate(E(E_size))
     end if

     do J=1,NE
        read(unit_vector) num, E(num)
        if (num >= NEmin .and. num <= NEmax) then
           read(unit_vector)(A(I,NUM-NEMIN+1,kkk),I=1,N)
        else
           read(unit_vector)
        endif
     enddo

     do NUM=NEMIN,NEMAX
        NB=NUM-NEMIN+1
        write(unit_eig,"(2I12,F22.16)")NB,kkk, Ryd_eV*(E(NUM)-efermi)
     enddo
  end do kpoint

998 if (allocated(E)) deallocate(E)
  write(unit_out,*)'vector read in',kkk
end subroutine read_vec
end module read_vec_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
