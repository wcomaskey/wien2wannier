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
  public :: init_xa3, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, vecsz, a, Eb

#ifndef _COMPLEX_
  real(R8), &
#else
  complex(C16), &
#endif
  allocatable :: a(:,:,:)

  integer,  allocatable :: GX(:,:),GY(:,:),GZ(:,:), vecsz(:)
  !> Band energies (Ry, absolute), (band, k), filled by read_vec.  read_vec
  !> already reads them to write case.eig and then throws them away; keeping
  !> them costs Nb*nk reals and is what any diagnostic that needs to compare
  !> spectra at k and g.k requires.  Absolute, NOT E_F-relative: the .eig file
  !> is written E_F-relative, and mixing the two conventions is a mistake this
  !> tree has made before in the plotting layer.
  real(R8), allocatable :: Eb(:,:)
  real(R8), allocatable :: BKX(:), BKY(:), BKZ(:)
  real(R8), allocatable :: XK(:),  YK(:),  ZK(:)

contains
  subroutine init_xa3(nb, nmat, num_kpts)
    integer, intent(in) :: NB, Nmat, num_kpts

    allocate( a(Nmat, Nb, num_kpts), vecsz(num_kpts), Eb(Nb, num_kpts) )
    Eb = 0

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
  use xa3,       only: xk, yk, zk, A, gx, gy, gz, vecsz, Eb
  use bessel,    only: ri_mat
  use Amn_Mmn,   only: overlap
  use pairs,     only: kp, kpb, bqx, bqy, bqz
  use lolog,     only: n_rad, nlo, nlov, nlon
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t
  use mpi_w2w

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
  ! C1 Alm cache: on a uniform MP grid every k is hit once as the head of
  ! k1loop and nntot times as a neighbour, so almgen used to run (1+nntot)
  ! times per (atom, k).  Cache the RAW almgen output per k for the CURRENT
  ! atom (the atom loop is outermost; atpar state is atom-scoped) and copy
  ! (conjugating for the neighbour use) instead of recomputing.
  complex(C16), allocatable :: alm_cache(:,:,:,:,:)

  logical,      allocatable :: alm_built(:)

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
  logical  :: use_zgemm   ! runtime crossover: zgemm for Nb>=24, zgeru otherwise

  integer :: k1_prog_itvl
  k1_prog_itvl = min(max(num_kpts/10, 1), 100)

  Nb = inwf%bmax - inwf%bmin + 1

  ! Runtime crossover: batched zgemm wins when Nb is large enough
  ! for BLAS-3 efficiency to overcome gather/S_base overhead.
  ! Below Nb=24, the original zgeru rank-1 updates are faster.
  use_zgemm = (Nb >= 24)
  if (is_root) &
       write(unit_out, "(' Overlap method: ', A, ' (Nb=', I4, ')')") &
       merge('batched zgemm', 'zgeru (small)', use_zgemm), Nb

  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST=1

  atoms: do jatom = 1, stru%Nneq
     if (is_root) &
          write(unit_out, "(/, '===== atom', I5, ' /', I5, ' =====' /)") &
          jatom, stru%Nneq

     allocate(Alm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)),  &
          &   Blm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)))
     allocate(alm_cache(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom), num_kpts), &
          &   alm_built(num_kpts))
     alm_built = .false.
     if (is_root .and. jatom == 1) write(unit_out, &
          '(" Alm cache: ", F8.1, " MB per atom class")') &
          16.0_R8*size(alm_cache)/1024.0_R8**2

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

     if (use_zgemm) then
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
        if (is_root) &
             write(unit_out, "(' Batched zgemm: Na=', I4, ' active columns')") Na

        ! Allocate work arrays for batched zgemm
        allocate(S_base(Na, Na))
        allocate(ALM_2d(Nb, Na), BLM_2d(Nb, Na), W(Nb, Na))
     end if

     k1loop: do k1=1,num_kpts      !   loop over k1
#ifdef MPI
        if (k_owner(k1, num_kpts) /= myrank) cycle k1loop
#endif
        pair_index = (k1-1)*NNTOT + 1
        kkk=KP(pair_index)
        KX1=XK(kkk)
        KY1=YK(kkk)
        KZ1=ZK(kkk)

        if (.not. alm_built(kkk)) then
           call almgen(alm_cache(:,:,:,:,kkk), stru, jatom, lfirst, Nb, kkk)
           alm_built(kkk) = .true.
        end if
        Alm = alm_cache(:,:,:,:,kkk)

        k2loop: do k2=1,NNTOT
           pair_index = (k1-1)*NNTOT + k2
           kkk=KPB(pair_index)
           call cpu_time(tt0)
           if (.not. alm_built(kkk)) then
              call almgen(alm_cache(:,:,:,:,kkk), stru, jatom, lfirst, Nb, kkk)
              alm_built(kkk) = .true.
           end if
           BLM = conjg(alm_cache(:,:,:,:,kkk))

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
           ! Y_lm(b) in the GLOBAL cartesian frame: almgen builds the Alm
           ! for every equivalent atom in the global frame (no rotij/locrot
           ! -- those belong to l2Amn's separate local-frame builder), so
           ! the Gaunt contraction needs the global-frame Y_lm(b) too, and
           ! one YLB serves all atoms of the multiplet.
           call YLM (BK, inwf%LJmax, YLB)
           indexj=0
           do LJ = 0, inwf%LJmax
              do MJ=-LJ,LJ
                 indexj=indexj+1
                 YLB(indexj)=conjg(YLB(indexj))*(0,1)**LJ
              enddo
           enddo
           call cpu_time(tt2)

           muloop: do mu=1, stru%mult(jatom)
              latom=lfirst-1+mu

              arg1=Bx * stru%pos(1,latom)*TAU
              arg2=By * stru%pos(2,latom)*TAU
              arg3=Bz * stru%pos(3,latom)*TAU
              PHSHEL=exp((0,1)*(ARG1+ARG2+ARG3))*2*TAU

              if (use_zgemm) then
                 ! --- Opt 3: Batched zgemm path (efficient for Nb >= 24) ---
                 ! Build S_base on first mu; YLB is global-frame so it is
                 ! valid for every equivalent atom of the multiplet
                 if (mu == 1) then
                    ! S_base(col_b, col_a) = sum over (LJ,MJ) of
                    !   ri_mat * YLB * gaunt (without PHSHEL)
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
                                              + ri_mat(R_index,L_index)*tmp
                                      enddo
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 end if

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
                 call zgemm('N', 'N', Nb, Na, Na, PHSHEL, &
                      BLM_2d, Nb, S_base, Na, &
                      (0.0_R8, 0.0_R8), W, Nb)
                 call zgemm('N', 'T', Nb, Nb, Na, (1.0_R8, 0.0_R8), &
                      W, Nb, ALM_2d, Nb, &
                      (1.0_R8, 0.0_R8), overlap(1,1,pair_index), &
                      size(overlap,1))

              else
                 ! --- Original zgeru path (efficient for small Nb) ---
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
                                tmp=YLB(indexj)*PHSHEL &
                                     *get_gaunt(L2,LJ,L1,M2,MJ,M1)
                                R_index=0
                                do irf1=1,n_rad(L1)
                                   do irf2=1,n_rad(L2)
                                      R_index=R_index+1
                                      tmp1=ri_mat(R_index,L_index)*tmp
                                      call zgeru(NB,NB,tmp1, &
                                           BLM(1,irf2,index2,mu),1, &
                                           ALM(1,irf1,index1,mu),1, &
                                           overlap(1,1,pair_index), &
                                           size(overlap,1))
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
              end if
           enddo muloop
           call cpu_time(tt3)
           talm=talm+tt1-tt0
           tmeas1=tmeas1+tt2-tt1
           tmeas2=tmeas2+tt3-tt2
        end do k2loop

        if (is_root .and. mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
     end do k1loop
     if (use_zgemm) deallocate(S_base, ALM_2d, BLM_2d, W)

     deallocate(ALM, BLM, alm_cache, alm_built)

     call cpu_time(t2)
     talm=talm*(NNTOT+1)/NNTOT
     if (is_root) &
          write(unit_out,*) 'CPU time used for atom ',JATOM,' ='&
          ,t2-t1,talm,tmeas1,tmeas2
  end do atoms
end subroutine l2MMN
end module     l2Mmn_m


module     l2Amn_m; contains
subroutine l2Amn(stru, inwf, num_kpts, proj_out, write_files, proj_full_out)
  use param,     only: Lmax2, NRF, lomax, unit_vsp, unit_out, unit_vector
  use w2w,       only: Nmat, iblock, unit_amn
  use clio,      only: close_checked
  use const,     only: R8, C16, TAU
  use xa,        only: BK, FJ, DFJ, PHS, BKrot, BKRloc
  use xa3,       only: A, XK,YK,ZK, GX,GY,GZ, BKX,BKY,BKZ, vecsz
  use lolog,     only: nlo,nlov,nlon,loor,ilo,lapw
  use loabc,     only: alo
  use atspdt,    only: P, DP
  use radfu,     only: RF1, RF2
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t
  use mpi_w2w

  !! procedure includes
  use rint13_m
  use harmon_m
  use atpar_m
  use Ylm_m
  use sphbes_m

  implicit none


  type(struct_t), intent(in) :: stru
  type(inwf_t),   intent(in) :: inwf
  integer,        intent(in) :: num_kpts
  ! optional export of the reduced projection (= the amn array, layout
  ! (Nb, Nproj, Nk)) for the in-memory MODE=HYBR hook; only meaningful
  ! on root (the reduce is to root).
  complex(C16), intent(out), optional, allocatable :: proj_out(:,:,:)
  ! suppress the .amn text+sidecar write (SOC sigma-loop: the caller
  ! assembles per-sigma column blocks and writes the combined file)
  logical, intent(in), optional :: write_files
  ! optional SECOND export: the full-cell (tail-extended) projection,
  ! accumulated in the same pass.  When present, the tails machinery is
  ! forced on internally while `projection` (and the written amn) keeps
  ! its W2W_TAILS-controlled legacy content -- the HYBR one-shot gets
  ! MT amn for the SCDM gauge/hand-off AND full-cell amn for Loewdin P
  ! from a single plane-wave sweep.
  complex(C16), intent(out), optional, allocatable :: proj_full_out(:,:,:)

  complex(C16) :: Yl((Lmax2+1)**2)        ! list of spherical harmics
  complex(C16) :: projection(inwf%bmax-inwf%bmin+1, inwf%Nproj, num_kpts)
  complex(C16), allocatable :: projection_t(:,:,:)
  complex(C16) :: phshel,cfac
  logical      :: wfiles, want_full, env_tails

  real(R8)     :: h_al(iblock), h_bl(iblock)
  complex(C16) :: h_yl  ( (Lmax2+1)**2, iblock)
  complex(C16) :: h_alyl( (Lmax2+1)**2, iblock)
  complex(C16) :: h_blyl( (Lmax2+1)**2, iblock)
  complex(C16) :: Alm   ( (Lmax2+1)**2, inwf%bmax-inwf%bmin+1, stru%Nat, Nrf)

  integer  :: i,ii,i3,l, m,mu,m1, n, ipr, ib, iY
  integer  :: jatom,jlo,jneq, idx, ibb,iscf,irf, lda,ldb,ldc, latom,lfirst
  integer  :: kkk, max, itap,itape,jtape, Nb
  real(R8) :: fac, arg, pi12lo(0:Lmax2)

  ! full-cell projector tails (env W2W_TAILS=1): u_l is continued past
  ! R_MT with a log-derivative-matched decaying exponential; the tail's
  ! overlap with the interstitial plane-wave sum is added to the MT
  ! matching (same conjg(Ylm)*phase structure, Rayleigh prefactor
  ! 4pi/sqrt(vol) i^l instead of 4pi R^2/sqrt(vol)).  MT-only projectors
  ! systematically understate sp/interstitial character (GaAs valence
  ! P ~ 0.4), which breaks projectability-driven subspace selection.
  logical      :: do_tails
  character(len=8) :: envtl
  integer, parameter :: NQT = 800, NRT = 1200
  real(R8)     :: tail_tab(NQT, 0:Lmax2), qtab_max, tail_nrm(0:Lmax2)
  real(R8)     :: tkap, tuR, tduR, qq, rr, tdr, tsum, tw, fac2
  real(R8)     :: xj(Lmax2+1), h_ql(iblock), tql
  real(R8)     :: lo_p2a(0:Lmax2, Nrf), lo_q2a(0:Lmax2, Nrf), t_p1, t_p2, t_oo
  integer      :: irf_y
  complex(C16) :: h_tlyl((Lmax2+1)**2, iblock)
  complex(C16) :: Alm_t((Lmax2+1)**2, inwf%bmax-inwf%bmin+1, stru%Nat)
  complex(C16) :: tail_c
  integer      :: iq, irt

  Nb = inwf%bmax-inwf%bmin+1
  projection = 0

  call get_environment_variable('W2W_TAILS', envtl)
  env_tails = (trim(envtl) == '1')
  want_full = present(proj_full_out)
  ! tails machinery runs if EITHER the legacy env switch is set (the
  ! written amn becomes full-cell) or the caller wants the separate
  ! full-cell export (written amn stays MT)
  do_tails = env_tails .or. want_full
  if (want_full) then
     allocate(projection_t(Nb, inwf%Nproj, num_kpts))
     projection_t = 0
  end if
  qtab_max = 1.0_R8
  if (do_tails) then
     do kkk = 1, num_kpts
        do i = 1, vecsz(kkk)
           BK = matmul(stru%conv_rec, (/ XK(kkk)+GX(i,kkk), &
                YK(kkk)+GY(i,kkk), ZK(kkk)+GZ(i,kkk) /))
           qq = sqrt(dot_product(BK, BK))
           if (qq > qtab_max) qtab_max = qq
        end do
     end do
     qtab_max = qtab_max * 1.02_R8
     if (is_root) write(unit_out, '(A,F10.4)') &
          ' full-cell projector tails ON; q_max =', qtab_max
  end if

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
     if (is_root) then
        write(unit_out, *)
        write(unit_out, '(7X,  "rotij", 10X, "transij")')
        write(unit_out, '(6X, 23("-"))')
     end if

     pi12lo = 0
     do l = 0,LOmax
        if (loor(l)) call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
             &                   rf1(:,l,3), rf2(:,l,3), pi12lo(l))
     end do

     ! Orthogonalized LO-shell projector overlaps per (l, rf): needed for
     ! projection columns with rf > 1 (semicore shells, see inwfmod).
     ! <phi_perp|psi> = (A2*p2 + A_rf*q2)/sqrt(q2), q2 = ||LO||^2-<u|LO>^2
     lo_p2a = 0; lo_q2a = 0
     do l = 0, LOmax
        if (.not. loor(l)) cycle
        do irf = 3, Nrf
           call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
                &      rf1(:,l,irf), rf2(:,l,irf), t_p1)
           call rint13(stru, jatom, rf1(:,l,irf), rf2(:,l,irf), &
                &      rf1(:,l,irf), rf2(:,l,irf), t_oo)
           call rint13(stru, jatom, rf1(:,l,2), rf2(:,l,2), &
                &      rf1(:,l,irf), rf2(:,l,irf), t_p2)
           lo_p2a(l, irf) = t_p2
           lo_q2a(l, irf) = t_oo - t_p1**2
        end do
     end do

     tail_nrm = 1
     if (do_tails) then
        tail_tab = 0
        do l = 0, Lmax2
           tuR = P(l,1); tduR = DP(l,1)
           if (abs(tuR) < 1.0e-8_R8) cycle       ! node at R_MT: no tail
           tkap = -tduR/tuR
           ! floor at the free-atom decay scale (~sqrt(2*I_p) a.u.): the
           ! linearization-energy log-derivative is often non-decaying,
           ! and a slower tail overlaps neighbouring spheres unphysically
           if (tkap < 0.9_R8) tkap = 0.9_R8
           rr = stru%RMT(jatom)
           ! <u|u>_MT = 1 plus the analytic exponential-tail norm
           tsum = tuR**2 * (rr**2/(2*tkap) + rr/(2*tkap**2) + 1/(4*tkap**3))
           tail_nrm(l) = 1/sqrt(1 + tsum)
           ! t_l(q) = int_R^inf u(R) e^{-kap(r-R)} j_l(qr) r^2 dr, cached in q
           tdr = min(25.0_R8/tkap, 12.0_R8)/NRT
           do iq = 1, NQT
              qq = qtab_max*(iq-1)/real(NQT-1, R8)
              tsum = 0
              do irt = 0, NRT
                 rr = stru%RMT(jatom) + irt*tdr
                 tw = 1; if (irt == 0 .or. irt == NRT) tw = 0.5_R8
                 call sphbes(Lmax2, qq*rr, xj)
                 tsum = tsum + tw*exp(-tkap*(rr-stru%RMT(jatom)))*xj(l+1)*rr**2
              end do
              tail_tab(iq, l) = tuR * tsum * tdr
           end do
        end do
        if (is_root) write(unit_out, '(A,I3,A,6F9.4)') '  atom', jatom, &
             ' tail norms:', (tail_nrm(l), l=0,Lmax2)
     end if

     kpoints: do KKK=1,NUM_KPTS
#ifdef MPI
        if (k_owner(KKK, num_kpts) /= myrank) cycle kpoints
#endif
        N=vecsz(kkk)

        ! plane wave Q+k
        BKX(1:N) = XK(kkk) + GX(1:N, kkk)
        BKY(1:N) = YK(kkk) + GY(1:N, kkk)
        BKZ(1:N) = ZK(kkk) + GZ(1:N, kkk)

        ! calculate Bessel functions j_l(|Q+k|*RMT)
        call harmon(stru, jatom, N, Bkx, Bky, Bkz, Lmax2, FJ, DFJ)

        fac = 2*TAU * stru%RMT(jatom)**2 / sqrt(stru%vol)
        fac2 = 2*TAU / sqrt(stru%vol)
        Alm = 0
        if (do_tails) Alm_t = 0

        equiv1: do mu = 1, stru%mult(jatom) !  loop over equivalent atoms
           latom=lfirst-1+mu
           if (kkk == 1 .and. is_root) then
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
                 h_ql(i3) = sqrt(dot_product(bkrloc, bkrloc))

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
              if (do_tails) then
                 idx=0
                 do L=0,LMAX2
                    do M=1,2*L+1
                       idx=idx+1
                       i3=0
                       do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                          i3=i3+1
                          tql = h_ql(i3)*(NQT-1)/qtab_max + 1
                          iq = int(tql)
                          if (iq < 1) iq = 1
                          if (iq > NQT-1) iq = NQT-1
                          h_tlyl(idx,i3) = ( tail_tab(iq,L) + (tql-iq) * &
                               (tail_tab(iq+1,L)-tail_tab(iq,L)) ) * h_YL(idx,i3)
                       enddo
                    enddo
                 enddo
#ifndef _COMPLEX_
                 call dgemm('N','N',2*idx,nb,ibb,1.d0, &
                      h_tlyl,lda,a(ii,1,kkk),ldb,1.d0, &
                      Alm_t(1,1,mu),ldc)
#else
                 call zgemm('N','N',idx,nb,ibb,(1.d0,0.d0), &
                      h_tlyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                      Alm_t(1,1,mu),ldc)
#endif
              end if
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
                    if (do_tails) Alm_t(idx,ib,mu)=Alm_t(idx,ib,mu)*fac2*cfac
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

                    irf_y = 1
                    if (allocated(inwf%projections(ipr)%rf)) &
                         irf_y = inwf%projections(ipr)%rf(iY)
                    if (irf_y > 1) then
                       ! orthogonalized LO-shell column (semicore)
                       if (lo_q2a(l, irf_y) > 1.0e-8_R8) then
                          projection(ib, ipr, kkk) = &
                            projection(ib, ipr, kkk) &
                            & + ( conjg(Alm(idx, ib, mu, 2))*lo_p2a(l, irf_y) &
                            &     + conjg(Alm(idx, ib, mu, irf_y))            &
                            &       * lo_q2a(l, irf_y) )                      &
                            &   / sqrt(lo_q2a(l, irf_y))                      &
                            &   * inwf%projections(ipr)%coeff(iY)
                          ! LO-shell columns have no plane-wave tail: the
                          ! full-cell export takes the identical value
                          if (want_full) projection_t(ib, ipr, kkk) = &
                            projection_t(ib, ipr, kkk) &
                            & + ( conjg(Alm(idx, ib, mu, 2))*lo_p2a(l, irf_y) &
                            &     + conjg(Alm(idx, ib, mu, irf_y))            &
                            &       * lo_q2a(l, irf_y) )                      &
                            &   / sqrt(lo_q2a(l, irf_y))                      &
                            &   * inwf%projections(ipr)%coeff(iY)
                       end if
                    else
                       ! legacy amn content: full-cell iff W2W_TAILS=1
                       ! (tail_nrm(l)=1 exactly when tails are off, so the
                       ! merge keeps the MT path bit-identical)
                       tail_c = 0
                       if (env_tails) tail_c = conjg(Alm_t(idx, ib, mu))
                       projection(ib, ipr, kkk) = projection(ib, ipr, kkk) &
                            & + ( conjg(Alm(idx, ib, mu, 1))               &
                            &     + conjg(Alm(idx, ib, mu, 3))*pi12lo(l)   &
                            &     + tail_c )                               &
                            &   * inwf%projections(ipr)%coeff(iY)          &
                            &   * merge(tail_nrm(l), 1.0_R8, env_tails)
                       if (want_full) then
                          tail_c = conjg(Alm_t(idx, ib, mu))
                          projection_t(ib, ipr, kkk) = projection_t(ib, ipr, kkk) &
                               & + ( conjg(Alm(idx, ib, mu, 1))               &
                               &     + conjg(Alm(idx, ib, mu, 3))*pi12lo(l)   &
                               &     + tail_c )                               &
                               &   * inwf%projections(ipr)%coeff(iY)          &
                               &   * tail_nrm(l)
                       end if
                    end if
                 enddo
              enddo
           enddo
        enddo equiv2
     enddo kpoints
  enddo atoms

  ! Reduce projection from all ranks to rank 0
  ! promote to 64-bit BEFORE multiplying (see mpi_w2w.F)
  call mpi_w2w_reduce_sum_complex_inplace(projection, &
       int(Nb,8) * int(inwf%Nproj,8) * int(num_kpts,8))
  if (want_full) call mpi_w2w_reduce_sum_complex_inplace(projection_t, &
       int(Nb,8) * int(inwf%Nproj,8) * int(num_kpts,8))

  wfiles = .true.
  if (present(write_files)) wfiles = write_files

  ! Rank 0 writes AMN output
  if (is_root .and. wfiles) then
     ! es24.15: full double precision.  The old 2e18.5 kept only 5 significant
     ! digits, putting a 1e-5 floor under P=|amn|^2 masks and the SCDM gauge
     ! downstream.
     !
     ! Like the text .mmn, this is an INTERMEDIATE under HYBR -- the binary
     ! sidecar below is what in-tree readers use -- so it is skipped unless
     ! `text_io = true`.  When written, one WRITE per k replaces Nproj*Nb of
     ! them; FORMAT REVERSION keeps the record stream byte-identical.
     if (.not. (inwf%hybrid .and. .not. inwf%hyb_text_io)) then
        do kkk = 1,num_kpts
           write(unit_amn, '(2i4,1x,i5,1x,2es24.15)') &
                ((ib, ipr, kkk, projection(ib, ipr, KKK), ib=1,Nb), &
                  ipr=1,inwf%Nproj)
        enddo
     end if

     ! Binary sidecar (B4): projection(ib,ipr,kkk) is exactly the
     ! amn(nb,nproj,nk) layout hybrid_files_m reads.  Close the text
     ! unit first so the sidecar's mtime stays >= the text file's.
     block
       character(len=1024) :: amnfn
       integer :: u_sc, ios_sc
       inquire(unit=unit_amn, name=amnfn)
       call close_checked(unit_amn, 'the .amn text file')
       open(newunit=u_sc, file=trim(amnfn)//'.bin', form='unformatted', &
            access='stream', status='replace', action='write', iostat=ios_sc)
       if (ios_sc == 0) then
          write(u_sc) 'AMNB', 1, Nb, num_kpts, inwf%Nproj
          write(u_sc) projection
          call close_checked(u_sc, 'the .amn binary sidecar')
       end if
     end block
  end if

  if (present(proj_out) .and. is_root) proj_out = projection
  if (want_full .and. is_root) proj_full_out = projection_t

  return
end subroutine l2Amn
end module     l2Amn_m


module     autoprojection_m; contains

subroutine compute_projectability(stru, inwf, num_kpts, proj_char, proj_state, &
     &                             eweight, proj_max, proj_char_lo, proj_max_lo, &
     &                             lo_rf_out, accumulate, finalize, &
     &                             proj_state_lo, component_sum_max, &
     &                             reduce_proj_state, degtol)
  !! Compute projectability of all bands onto (atom, l) channels.
  !! Uses almgen to decompose Bloch states, then accumulates |Alm|^2.
  !! Output: proj_char(Nb, Nneq, 0:Lmax2) summed over k-points and m.
  !! Optional eweight(Nb, Nk): per-state weight applied to the proj_char
  !! (channel-ranking) aggregate ONLY -- an occupation-like erfc factor that
  !! suppresses the Rayleigh-expansion tails of high-energy NFE states,
  !! which otherwise pollute the (atom, l) ranking (P2).  proj_state stays
  !! UNWEIGHTED: the per-state masks need raw projectability.
  !! Optional: proj_state(Nb, Nk, Nneq, 0:Lmax2) per-state (band, k)
  !! projectability |<u_l Y_lm|psi_nk>|^2 = |Alm(irf=1) + pi12lo*Alm(irf=3)|^2
  !! summed over m and equivalent atoms.  This matches the projector
  !! convention of l2Amn (u_dot excluded), so per-state masks derived from
  !! it are consistent with the .amn projections.  Used for external
  !! per-state disentanglement analysis.
  use param,     only: Lmax2, NRF, lomax, unit_vsp, unit_out, unit_vector
  use w2w,       only: Nmat, iblock
  use const,     only: R8, C16, TAU, Ryd_eV
  use xa3,       only: vecsz, Eb
  use lolog,     only: loor, n_rad
  use radfu,     only: RF1, RF2
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t
  use clio,      only: croak
  use projection_support_m, only: state_support_max, state_support_max_trace
  use mpi_w2w

  !! procedure includes
  use rint13_m
  use atpar_m
  use almgen_m

  implicit none

  type(struct_t), intent(in)  :: stru
  type(inwf_t),   intent(in)  :: inwf
  integer,        intent(in)  :: num_kpts
  ! The accumulated outputs are intent(inout), not intent(out): with
  ! accumulate=.true. their values must survive from the previous call,
  ! which intent(out) does not guarantee (they become undefined on entry).
  real(R8),       intent(inout) :: proj_char(:,:,0:)
  ! proj_char(Nb, Nneq, 0:Lmax2)
  real(R8), intent(inout), optional :: proj_state(:,:,:,0:)
  ! proj_state(Nb, Nk, Nneq, 0:Lmax2)
  real(R8), intent(in),  optional :: eweight(:,:)
  ! eweight(Nb, Nk): P2 ranking weight (see header)
  real(R8), intent(inout), optional :: proj_max(:,0:)
  ! proj_max(Nneq, 0:Lmax2): max over all states (band, k) of the raw
  ! (unweighted) channel projectability.  A physical channel OWNS a
  ! state somewhere in the window (large maxP); a ghost/tail channel
  ! never does.  Used by generate_projections as a support gate (P3).
  real(R8), intent(inout), optional :: proj_char_lo(:,:,0:), proj_max_lo(:,0:)
  ! Orthogonalized-shell LO statistics: same aggregates for the projector
  ! phi_perp = (u_LO - <u|u_LO> u)/||.||, built per (atom, l) from the LO
  ! radial MOST orthogonal to u_l.  Semicore shells (e.g. Ga-3d) live in
  ! an LO radial and are invisible to the u-based projector (amplitude
  ! only pi12lo ~ 0.1); phi_perp sees them at ~1 while vanishing
  ! automatically for valence APW+lo radials parallel to u_l.
  integer,  intent(out), optional :: lo_rf_out(:,0:)
  ! lo_rf_out(Nneq, 0:Lmax2): chosen LO rf index per channel (0 = none)
  logical, intent(in), optional :: accumulate, finalize
  ! SOC spinor support: the caller invokes this once per spinor component
  ! (the two vector files), re-filling the module A array between calls.
  ! accumulate=.true. skips the entry-zeroing so the |Alm|^2 of the second
  ! component ADDS to the first (proj_char/proj_state/proj_char_lo are sums).
  ! proj_max/proj_max_lo do NOT compose across calls by themselves: the
  ! streaming per-k max() gives max over (n, k, s) of the PER-SPIN P_s(n,k),
  ! because stmp is re-zeroed per k.  The channel support the caller wants is
  !     proj_max(j,l) = max_{n,k} sum_s P_s(n,k; j,l)
  ! which needs the component-summed per-state array: pass
  ! component_sum_max=.true. (below).  On a single pass state_support_max
  ! over proj_state reproduces the streaming maximum bit-for-bit (same val,
  ! same (mu, m) order, same zero, and max is exact), so component_sum_max
  ! is a MEMORY choice (Nb*Nk*Nneq*(Lmax2+1) scratch), not a correctness
  ! switch.
  ! finalize=.false. defers the MPI reduce + per-band normalization so they
  ! run once, on the last (accumulating) call.  Both default to the scalar
  ! single-call behaviour (accumulate=F, finalize=T) -> no regression.
  real(R8), intent(inout), optional :: proj_state_lo(:,:,:,0:)
  ! proj_state_lo(Nb, Nk, Nneq, 0:Lmax2): unweighted per-state support of
  ! the orthogonalized LO-shell projector -- the LO twin of proj_state, fed
  ! val_lo under exactly the guards that feed stmp_lo.  Scratch: never
  ! MPI-reduced, never normalized.
  logical, intent(in), optional :: component_sum_max
  ! .true.: skip the streaming max() and derive proj_max/proj_max_lo at
  ! finalize from the accumulated proj_state/proj_state_lo, which must then
  ! be present (fails closed; no component-wise fallback).
  logical, intent(in), optional :: reduce_proj_state
  ! .false.: no full-array MPI allreduce of proj_state at finalize (it is
  ! scratch for component_sum_max; the maximum is reduced instead).
  ! Default .true. = existing behaviour.
  real(R8), intent(in), optional :: degtol
  ! Degeneracy tolerance (Ry) for the component_sum_max multiplet-trace
  ! metric (D5/D4-2): proj_max/proj_max_lo are then max_{k,c} of the
  ! per-cluster TRACE of proj_state/proj_state_lo, normalized by cluster
  ! size |c| (state_support_max_trace), c a contiguous run of bands within
  ! degtol of their neighbour in energy -- the fix for state_support_max's
  ! basis dependence inside an exactly degenerate multiplet (see
  ! projection_support.F).  Absent or <= 0: explicit opt-out, falls back to
  ! plain state_support_max (no clustering, no dependence on Eb).

  complex(C16), allocatable :: Alm(:,:,:,:)
  real(R8),     allocatable :: stmp(:,:)

  integer  :: jatom, kkk, lfirst, itape, jtape, iscf
  integer  :: l, m, idx, ib, mu, irf, Nb, mirf
  real(R8) :: pi12lo(0:Lmax2), pei_l, val
  logical  :: want_state, want_w, want_max, want_lo, do_acc, do_fin
  logical  :: want_state_lo, sum_max, red_state, use_trace
  real(R8) :: wfac
  ! orthogonalized-shell LO projector data (per atom, per l)
  integer  :: lo_rf(0:Lmax2), jlo_c, rf_c
  real(R8) :: lo_p1(0:Lmax2), lo_p2(0:Lmax2), lo_q2(0:Lmax2)
  real(R8) :: c_p1, c_p2, c_oo, c_ratio, best_ratio, val_lo
  real(R8), allocatable :: stmp_lo(:,:)
  ! multiplet-trace census (D5/D4-2), printed once on the root rank; a
  ! function of Eb/degtol alone, so identical on every rank without an
  ! MPI reduce (Eb is complete everywhere -- see read_vec)
  integer  :: ncl_min, ncl_max, nclus
  real(R8) :: gap_merged_max, gap_open_min

  Nb = inwf%bmax - inwf%bmin + 1
  do_acc = .false.; if (present(accumulate)) do_acc = accumulate
  do_fin = .true.;  if (present(finalize))   do_fin = finalize
  want_state = present(proj_state)
  want_w = present(eweight)
  want_max = present(proj_max)
  want_lo = present(proj_char_lo)
  want_state_lo = present(proj_state_lo)
  sum_max = .false.
  if (present(component_sum_max)) sum_max = component_sum_max
  use_trace = .false.
  if (present(degtol)) use_trace = degtol > 0.0_R8
  red_state = .true.
  if (present(reduce_proj_state)) red_state = reduce_proj_state
  if (sum_max) then
     if (want_max .and. .not. want_state) call croak('compute_projectability: ' // &
          'component_sum_max needs proj_state to derive proj_max')
     if (present(proj_max_lo) .and. .not. want_state_lo) call croak( &
          'compute_projectability: component_sum_max needs proj_state_lo ' // &
          'to derive proj_max_lo')
  end if
  if (.not. do_acc) then
     proj_char = 0.0_R8
     if (want_state) proj_state = 0.0_R8
     if (want_state_lo) proj_state_lo = 0.0_R8
     if (want_max) proj_max = 0.0_R8
     if (want_lo) then
        proj_char_lo = 0.0_R8
        if (present(proj_max_lo)) proj_max_lo = 0.0_R8
        if (present(lo_rf_out))   lo_rf_out = 0
     end if
  end if

  rewind(unit_vsp)
  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST = 1

  atoms: do jatom = 1, stru%Nneq
     if (jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     itape = unit_vector
     jtape = unit_vsp
     rewind(itape)
     call atpar(stru, jatom, itape, jtape)
     rewind(itape)

     ! Compute LO overlap integrals (same as l2Amn)
     pi12lo = 0
     do l = 0, LOmax
        if (loor(l)) call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
             &                   rf1(:,l,3), rf2(:,l,3), pi12lo(l))
     end do

     ! Orthogonalized-shell LO projector: per l pick the LO radial most
     ! orthogonal to u_l (min |<u|u_LO>|/||u_LO||); that radial carries a
     ! DISTINCT shell (semicore).  q2 = ||u_LO||^2 - <u|u_LO>^2 is the
     ! squared norm of the orthogonal component; q2 -> 0 for APW+lo
     ! radials parallel to u, killing the channel automatically.
     lo_rf = 0; lo_p1 = 0; lo_p2 = 0; lo_q2 = 0
     if (want_lo) then
        do l = 0, LOmax
           if (.not. loor(l)) cycle
           best_ratio = huge(1.0_R8)
           do jlo_c = 1, n_rad(l) - 2
              rf_c = 2 + jlo_c
              call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
                   &      rf1(:,l,rf_c), rf2(:,l,rf_c), c_p1)
              call rint13(stru, jatom, rf1(:,l,rf_c), rf2(:,l,rf_c), &
                   &      rf1(:,l,rf_c), rf2(:,l,rf_c), c_oo)
              if (c_oo < 1.0e-12_R8) cycle
              c_ratio = abs(c_p1)/sqrt(c_oo)
              if (c_ratio < best_ratio) then
                 best_ratio = c_ratio
                 lo_rf(l) = rf_c
                 lo_p1(l) = c_p1
                 lo_q2(l) = c_oo - c_p1**2
              end if
           end do
           if (lo_rf(l) > 0) then
              call rint13(stru, jatom, rf1(:,l,2), rf2(:,l,2), &
                   &      rf1(:,l,lo_rf(l)), rf2(:,l,lo_rf(l)), lo_p2(l))
              if (lo_q2(l) < 1.0e-8_R8) lo_rf(l) = 0
           end if
           if (present(lo_rf_out)) lo_rf_out(jatom, l) = lo_rf(l)
        end do
     end if

     ! Alm shape from almgen: (Nb, Nrf, (Lmax2+1)^2, mult(jatom))
     allocate(Alm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)))
     if (want_max .and. .not. allocated(stmp)) allocate(stmp(Nb, 0:Lmax2))
     if (want_lo .and. .not. allocated(stmp_lo)) allocate(stmp_lo(Nb, 0:Lmax2))

     kpoints: do kkk = 1, num_kpts
#ifdef MPI
        if (k_owner(kkk, num_kpts) /= myrank) cycle kpoints
#endif
        ! almgen computes the full Alm decomposition for this atom/k-point
        call almgen(Alm, stru, jatom, lfirst, Nb, kkk)

        if (want_max) stmp = 0.0_R8
        if (want_lo) stmp_lo = 0.0_R8

        ! Accumulate |Alm|^2 summed over equivalent atoms and m
        do mu = 1, stru%mult(jatom)
           idx = 0
           do l = 0, Lmax2
              mirf = n_rad(l)
              do m = 1, 2*l+1
                 idx = idx + 1
                 do ib = 1, Nb
                    ! irf=1: u_l contribution (always present)
                    val = real(Alm(ib, 1, idx, mu))**2 &
                        + aimag(Alm(ib, 1, idx, mu))**2

                    ! irf=2: u_dot contribution (always present for LAPW)
                    if (mirf >= 2) then
                       val = val + real(Alm(ib, 2, idx, mu))**2 &
                                 + aimag(Alm(ib, 2, idx, mu))**2
                    end if

                    ! irf=3: LO contribution (if local orbital exists)
                    if (mirf >= 3 .and. loor(l)) then
                       val = val + (real(Alm(ib, 3, idx, mu))**2 &
                                  + aimag(Alm(ib, 3, idx, mu))**2) &
                                 * pi12lo(l)**2
                    end if

                    wfac = 1.0_R8
                    if (want_w) wfac = eweight(ib, kkk)
                    proj_char(ib, jatom, l) = proj_char(ib, jatom, l) &
                         + val*wfac

                    ! Orthogonalized LO-shell projectability:
                    ! <phi_perp|psi> = (A2*p2 + A_r*q2)/sqrt(q2)
                    ! (cross-LO overlaps neglected)
                    if (want_lo) then
                       if (l <= LOmax) then
                          if (lo_rf(l) > 0) then
                             val_lo = abs(Alm(ib, 2, idx, mu)*lo_p2(l) &
                                  + Alm(ib, lo_rf(l), idx, mu)*lo_q2(l))**2 &
                                  / lo_q2(l)
                             wfac = 1.0_R8
                             if (want_w) wfac = eweight(ib, kkk)
                             proj_char_lo(ib, jatom, l) = &
                                  proj_char_lo(ib, jatom, l) + val_lo*wfac
                             stmp_lo(ib, l) = stmp_lo(ib, l) + val_lo
                             if (want_state_lo) proj_state_lo(ib, kkk, jatom, l) = &
                                  proj_state_lo(ib, kkk, jatom, l) + val_lo
                          end if
                       end if
                    end if

                    ! Per-state projectability with l2Amn's projector
                    ! convention: |<u_l Y_lm|psi>|^2, u_dot excluded,
                    ! LO entering via <u|u_LO> = pi12lo.
                    if (want_state .or. want_max) then
                       if (mirf >= 3 .and. loor(l)) then
                          val = abs(Alm(ib, 1, idx, mu) &
                               &    + pi12lo(l)*Alm(ib, 3, idx, mu))**2
                       else
                          val = abs(Alm(ib, 1, idx, mu))**2
                       end if
                       if (want_state) proj_state(ib, kkk, jatom, l) = &
                            proj_state(ib, kkk, jatom, l) + val
                       ! per-state channel total (summed over m, mu)
                       if (want_max) stmp(ib, l) = stmp(ib, l) + val
                    end if
                 end do
              end do
           end do
        end do

        ! streaming per-k maximum (per component; see component_sum_max)
        if (want_max .and. .not. sum_max) then
           do l = 0, Lmax2
              do ib = 1, Nb
                 proj_max(jatom, l) = max(proj_max(jatom, l), stmp(ib, l))
              end do
           end do
        end if
        if (want_lo .and. present(proj_max_lo) .and. .not. sum_max) then
           do l = 0, Lmax2
              do ib = 1, Nb
                 proj_max_lo(jatom, l) = &
                      max(proj_max_lo(jatom, l), stmp_lo(ib, l))
              end do
           end do
        end if
     end do kpoints

     deallocate(Alm)
  end do atoms

  ! Finalize once (deferred until the last spinor component in the SOC
  ! accumulate sequence): MPI reduce so all ranks agree, then normalize.
  if (do_fin) then
     ! component_sum_max: maxima from the component-summed per-state arrays
     ! over the rank-local k FIRST, then the allreduce_max below -- k
     ! ownership partitions the k set, so the max of the local maxima is
     ! the global maximum (unowned rows are still zero and cannot win).
     if (sum_max) then
        if (use_trace) then
           ! D5/D4-2: multiplet-trace metric.  ncl_min/ncl_max/nclus/the two
           ! gap stats depend on Eb and degtol only (not on proj_state), so
           ! they come out identical from either call below and identical on
           ! every rank -- Eb is complete everywhere -- with no MPI reduce.
           if (want_max) call state_support_max_trace(proj_state, &
                Eb(:,1:num_kpts), degtol, proj_max, ncl_min, ncl_max, &
                nclus, gap_merged_max, gap_open_min)
           if (want_lo .and. present(proj_max_lo)) &
                call state_support_max_trace(proj_state_lo, &
                Eb(:,1:num_kpts), degtol, proj_max_lo, ncl_min, ncl_max, &
                nclus, gap_merged_max, gap_open_min)
           if (is_root .and. (want_max .or. &
                (want_lo .and. present(proj_max_lo)))) write(unit_out, &
                '(A,ES9.2,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,ES9.2,A)') &
                ' AUTO SOC support: multiplet-trace metric, degtol = ', &
                degtol*Ryd_eV, ' eV; clusters min/max size = ', ncl_min, &
                '/', ncl_max, '; ', nclus, ' clusters over ', num_kpts, &
                ' k; largest gap merged = ', gap_merged_max*Ryd_eV, &
                ' eV, smallest gap left open = ', gap_open_min*Ryd_eV, ' eV'
        else
           if (want_max) call state_support_max(proj_state, proj_max)
           if (want_lo .and. present(proj_max_lo)) &
                call state_support_max(proj_state_lo, proj_max_lo)
        end if
     end if
     ! MPI allreduce projectability so all ranks have identical data
     ! (proj_state_lo is scratch and is never reduced)
     call mpi_w2w_allreduce_sum_real_inplace(proj_char, &
          Nb * stru%Nneq * (Lmax2+1))
     if (want_state .and. red_state) &
          call mpi_w2w_allreduce_sum_real_inplace(proj_state, &
          Nb * num_kpts * stru%Nneq * (Lmax2+1))
     if (want_max) call mpi_w2w_allreduce_max_real_inplace(proj_max, &
          stru%Nneq * (Lmax2+1))
     if (want_lo) then
        call mpi_w2w_allreduce_sum_real_inplace(proj_char_lo, &
             Nb * stru%Nneq * (Lmax2+1))
        if (present(proj_max_lo)) call mpi_w2w_allreduce_max_real_inplace( &
             proj_max_lo, stru%Nneq * (Lmax2+1))
     end if

     ! Normalize: divide by Nk * (2l+1) * mult to get average per band
     do jatom = 1, stru%Nneq
        do l = 0, Lmax2
           proj_char(:, jatom, l) = proj_char(:, jatom, l) &
                / (num_kpts * (2*l+1) * stru%mult(jatom))
           if (want_lo) proj_char_lo(:, jatom, l) = proj_char_lo(:, jatom, l) &
                / (num_kpts * (2*l+1) * stru%mult(jatom))
        end do
     end do
  end if

end subroutine compute_projectability


subroutine write_projstate(stru, inwf, Nb, Nk, proj_state)
  !! Write per-state projectabilities to 'projectability.dat' for
  !! external analysis (per-state disentanglement threshold derivation).
  !! Values are raw sum_{mu,m} |Alm|^2 per (band, k, atom, l) channel --
  !! NOT normalized and NOT a partition of unity (interstitial charge
  !! is not included in the muffin-tin Alm decomposition).
  use param,     only: Lmax2, unit_out
  use const,     only: R8
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t

  implicit none

  type(struct_t), intent(in) :: stru
  type(inwf_t),   intent(in) :: inwf
  integer,        intent(in) :: Nb, Nk
  real(R8),       intent(in) :: proj_state(:,:,:,0:)

  integer, parameter :: UNIT_PST = 97

  integer :: jatom, l, ib, kkk, nchan

  nchan = stru%Nneq * (Lmax2 + 1)

  open(UNIT_PST, FILE='projectability.dat', STATUS='replace', &
       FORM='formatted')

  write(UNIT_PST, '(A)') '# w2w per-state projectability'
  write(UNIT_PST, '(A)') '# sum_{mu,m} |<u_l Y_lm|psi_nk>|^2 per (band, k,&
       & atom, l); l2Amn projector convention (u_dot excluded),&
       & interstitial not included'
  write(UNIT_PST, '(A, 5(1X,I0))') '# bmin bmax Nb Nk Nchan:', &
       inwf%bmin, inwf%bmax, Nb, Nk, nchan
  write(UNIT_PST, '(A)', ADVANCE='no') '# channels (jatom,l):'
  do jatom = 1, stru%Nneq
     do l = 0, Lmax2
        write(UNIT_PST, '(1X,"(",I0,",",I0,")")', ADVANCE='no') jatom, l
     end do
  end do
  write(UNIT_PST, *)
  write(UNIT_PST, '(A)') '# rows: kkk ib P(chan_1) ... P(chan_Nchan)'

  do kkk = 1, Nk
     do ib = 1, Nb
        write(UNIT_PST, '(I6, 1X, I5, *(1X, ES16.8))') kkk, ib, &
             ((proj_state(ib, kkk, jatom, l), l = 0, Lmax2), &
             jatom = 1, stru%Nneq)
     end do
  end do

  close(UNIT_PST)
  write(unit_out, '(A)') ' Per-state projectability written:&
       & projectability.dat'

end subroutine write_projstate



! Otsu two-cluster split position for a descending positive sequence
! (computed in log space): returns ncut = size of the upper cluster.
! Parameter-free; robust to outliers at the top of the list.
subroutine otsu_cut(v, n, ncut)
  use const, only: R8
  implicit none
  integer, intent(in) :: n
  real(R8), intent(in) :: v(n)
  integer, intent(out) :: ncut
  real(R8) :: lv(n), m1, m2, w1, w2, best, crit
  integer :: i, j, npos
  npos = 0
  do i = 1, n
     if (v(i) <= 0.0_R8) exit
     npos = npos + 1
     lv(npos) = log(v(i))
  end do
  ncut = max(npos, 1)
  if (npos < 3) return
  best = -1.0_R8
  do i = 1, npos-1
     m1 = sum(lv(1:i))/i; m2 = sum(lv(i+1:npos))/(npos-i)
     w1 = real(i,R8)/npos; w2 = 1.0_R8 - w1
     crit = w1*w2*(m1-m2)**2          ! between-class variance
     if (crit > best) then
        best = crit; ncut = i
     end if
  end do
end subroutine otsu_cut

subroutine generate_projections(stru, proj_char, Nb, Nk, config, inwf, &
     &                          proj_max, proj_char_lo, proj_max_lo, lo_rf_in)
  !! Analyze projectability and auto-generate projection specifications.
  !! Ranks (atom, l) channels by aggregate weight, selects top channels,
  !! and populates inwf%projections.
  !! Optional proj_max(Nneq, 0:Lmax2): per-channel max of the raw
  !! per-state projectability.  Channels whose best state falls below
  !! SUPPORT_RATIO of the overall best channel's are gated out (P3):
  !! they own no state anywhere in the window (Rayleigh tails,
  !! hybridization admixture of a semicore shell cut at bmin) and make
  !! the SCDM anchor set singular if selected.
  use param,     only: Lmax2, unit_out
  use const,     only: R8, DPk
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t, auto_config_t
  use mpi_w2w,   only: is_root

  implicit none

  type(struct_t),      intent(in)    :: stru
  real(R8),            intent(in)    :: proj_char(:,:,0:)
  integer,             intent(in)    :: Nb, Nk
  type(auto_config_t), intent(in)    :: config
  type(inwf_t),        intent(inout) :: inwf
  real(R8), intent(in), optional     :: proj_max(:,0:)
  ! orthogonalized-shell LO channel statistics (see compute_projectability):
  ! adds one candidate channel per (atom, l) whose LO radial carries a
  ! shell distinct from u_l (semicore).  The perp-projector support
  ! vanishes for APW+lo radials parallel to u, so no spurious channels.
  real(R8), intent(in), optional     :: proj_char_lo(:,:,0:), proj_max_lo(:,0:)
  integer,  intent(in), optional     :: lo_rf_in(:,0:)

  ! P3 support gate: a channel must own a state somewhere in the window
  ! at >= SUPPORT_RATIO of the best channel's ownership.  Measured dead
  ! zone across Fe / SrVO3 / GaAs / MgB2: physical channels >= 0.56,
  ! tail channels <= 0.45 (relative).
  real(R8) :: SUPPORT_RATIO   ! from config%support_ratio (0 disables)

  ! Work arrays for ranking (dynamically sized)
  integer  :: max_channels
  real(R8), allocatable :: weight(:), chan_sup(:)
  integer,  allocatable :: chan_atom(:), chan_l(:), chan_rf(:)
  integer  :: nchan, nsel, nproj_total, best_gap_pos
  integer  :: jatom, l, m, i, j, w, lfirst, mu, latom
  real(R8) :: total_weight, cumul_weight, max_weight, swap_w
  real(R8) :: gap_ratio, ratio, supmax
  integer  :: swap_a, swap_l, swap_r, nelig, n_add, ngated
  logical  :: mask_flag(2*stru%Nneq*(Lmax2+1)), swap_f
  logical  :: want_lo

  real(R8) :: band_proj_min, band_proj

  SUPPORT_RATIO = config%support_ratio

  max_channels = 2 * stru%Nneq * (Lmax2 + 1)
  allocate(weight(max_channels), chan_atom(max_channels), &
       chan_l(max_channels), chan_rf(max_channels), chan_sup(max_channels))
  chan_rf = 1
  chan_sup = 0.0_R8
  want_lo = present(proj_char_lo) .and. present(proj_max_lo) &
       .and. present(lo_rf_in)

  ! Build list of all (atom, l) channels with their weights
  nchan = 0
  do jatom = 1, stru%Nneq
     do l = 0, Lmax2
        nchan = nchan + 1
        chan_atom(nchan) = jatom
        chan_l(nchan) = l
        ! Average projectability over all bands
        weight(nchan) = sum(proj_char(1:Nb, jatom, l)) / Nb
        if (present(proj_max)) chan_sup(nchan) = proj_max(jatom, l)
     end do
  end do
  ! Orthogonalized LO-shell channels (semicore): distinct candidates with
  ! their own weight/support; the ranking machinery treats them like any
  ! other channel (parameter-free -- parallel-to-u LOs have ~zero support
  ! and sink below the numerical noise floor)
  if (want_lo) then
     do jatom = 1, stru%Nneq
        do l = 0, Lmax2
           if (lo_rf_in(jatom, l) <= 0) cycle
           nchan = nchan + 1
           chan_atom(nchan) = jatom
           chan_l(nchan) = l
           chan_rf(nchan) = lo_rf_in(jatom, l)
           weight(nchan) = sum(proj_char_lo(1:Nb, jatom, l)) / Nb
           chan_sup(nchan) = proj_max_lo(jatom, l)
        end do
     end do
  end if

  ! P3 support gate: zero the weight of channels that own no state
  ! anywhere in the window; they sort to the bottom and are never
  ! eligible.
  ngated = 0
  mask_flag = .false.
  if (present(proj_max) .and. config%support_ratio < 0.0_R8) then
     ! ALGORITHMIC support gate (support_ratio < 0): the cut is the
     ! dominant cliff in the sorted per-channel support values -- the
     ! hardcoded 0.5 ratio encoded exactly this observed cliff.
     block
       real(R8) :: sup(nchan), sw, cliff_val, wnoise
       integer :: ii, jj, nsup
       ! Otsu over REAL candidates only: channels whose ranking weight is
       ! numerical noise (l channels with no bands in the window) carry
       ! near-zero support and would pull the split down to zero-vs-rest
       wnoise = 1.0e-6_R8 * maxval(weight(1:nchan))
       nsup = 0
       do ii = 1, nchan
          if (weight(ii) > wnoise) then
             nsup = nsup + 1
             sup(nsup) = chan_sup(ii)
          end if
       end do
       do ii = 2, nsup                       ! insertion sort, descending
          sw = sup(ii); jj = ii - 1
          do while (jj >= 1)
             if (sup(jj) >= sw) exit
             sup(jj+1) = sup(jj); jj = jj - 1
          end do
          sup(jj+1) = sw
       end do
       ! CLIFF-MODE RANKING KEY: energy-relevance (P2 weight) x state
       ! ownership (best per-state projectability, clamped at 1).  The P2
       ! weight tail is flat/noisy across d and f channels; the support
       ! factor breaks those ties physically.
       do ii = 1, nchan
          weight(ii) = weight(ii) * min(chan_sup(ii), 1.0_R8)
       end do
       block
         integer :: ncut_s
         call otsu_cut(sup, nsup, ncut_s)
         if (is_root) then
            write(unit_out, '(A,I3,A,I3)') '   support spectrum (nsup=', &
                 nsup, '), otsu ncut=', ncut_s
            write(unit_out, '(A,12F8.4)') '    sup: ', sup(1:min(nsup,12))
         end if
         if (ncut_s < nsup .and. ncut_s >= 1) then
            cliff_val = sup(min(ncut_s+1, nsup))
            ! cliff mode does NOT gate (the gauge wants breadth) -- the
            ! upper-support cluster defines the MASK channel set (this is
            ! the per-state-support cliff the hardcoded 0.45/0.56 numbers
            ! encoded)
            do ii = 1, nchan
               mask_flag(ii) = weight(ii) > wnoise .and. &
                    chan_sup(ii) > cliff_val
            end do
         else
            do ii = 1, nchan
               mask_flag(ii) = weight(ii) > wnoise
            end do
         end if
       end block
     end block
  else if (present(proj_max)) then
     supmax = maxval(proj_max)
     if (supmax > 0.0_R8) then
        do i = 1, nchan
           if (weight(i) > 0.0_R8 .and. chan_rf(i) == 1 .and. &
                proj_max(chan_atom(i), chan_l(i)) < SUPPORT_RATIO*supmax) then
              weight(i) = 0.0_R8
              ngated = ngated + 1
           end if
        end do
     end if
     if (is_root .and. ngated > 0) &
          write(unit_out, '(A,I3,A,F4.2,A)') &
          '   P3 support gate: removed', ngated, &
          ' channel(s) with per-state support < ', SUPPORT_RATIO, &
          ' of best'
  end if

  ! Sort channels by weight descending (simple selection sort, nchan is small)
  do i = 1, nchan-1
     j = i
     do w = i+1, nchan
        if (weight(w) > weight(j)) j = w
     end do
     if (j /= i) then
        swap_w = weight(i); weight(i) = weight(j); weight(j) = swap_w
        swap_a = chan_atom(i); chan_atom(i) = chan_atom(j); chan_atom(j) = swap_a
        swap_l = chan_l(i); chan_l(i) = chan_l(j); chan_l(j) = swap_l
        swap_f = mask_flag(i); mask_flag(i) = mask_flag(j); mask_flag(j) = swap_f
        swap_r = chan_rf(i); chan_rf(i) = chan_rf(j); chan_rf(j) = swap_r
        swap_w = chan_sup(i); chan_sup(i) = chan_sup(j); chan_sup(j) = swap_w
     end if
  end do

  ! Compute total weight and find max
  total_weight = sum(weight(1:nchan))
  max_weight = weight(1)

  ! Select channels
  if (config%max_proj > 0) then
     ! User specified maximum number of projections
     nproj_total = 0
     nsel = 0
     do i = 1, nchan
        if (nproj_total + (2*chan_l(i)+1) * stru%mult(chan_atom(i)) &
             > config%max_proj) exit
        nsel = nsel + 1
        nproj_total = nproj_total + (2*chan_l(i)+1) * stru%mult(chan_atom(i))
     end do
  else
     ! Auto-select.  The threshold acts as a NOISE FLOOR: channels with
     ! weight(i)/weight(1) below it are never eligible.  With the P2
     ! energy weight active (recommended), physical channels separate
     ! from ghost/NFE channels by a sharp break in the ranked weights,
     ! so the PRIMARY cutoff is the largest relative drop between
     ! consecutive eligible weights (parameter-free).  Without P2 the
     ! weights are NFE-polluted and the drops are not meaningful, so
     ! the floor alone cuts (legacy threshold behaviour).
     if (config%threshold > 0.0_R8) then
        nelig = 0
        do i = 1, nchan
           if (weight(i) <= 0.0_R8) exit
           if (weight(i)/max_weight < config%threshold) exit
           nelig = nelig + 1
        end do
     else
        ! CLIFF MODE (threshold <= 0): every real candidate is eligible
        ! (the gauge wants breadth); only numerical-noise channels are
        ! excluded.  The mask boundary is detected separately (support
        ! Otsu), not by a weight threshold.
        nelig = count(weight(1:nchan) > 1.0e-6_R8 * max_weight)
     end if
     if (nelig == 0) nelig = 1

     ! Selection = all eligible channels (gate + floor).  The P3 support
     ! gate removes ghost/tail channels structurally, so every survivor
     ! above the noise floor owns bands in the window and is kept —
     ! "embrace larger WF counts".  The largest inter-weight drop is
     ! reported as a diagnostic only (it marks the physical/ghost
     ! boundary when no gate information is available).
     nsel = nelig
     best_gap_pos = nelig
     gap_ratio = 1.0_R8
     do i = 1, min(nelig, nchan-1)
        if (weight(i) <= 0.0_R8) exit
        if (weight(i+1) <= 0.0_R8) exit
        ratio = weight(i+1) / weight(i)
        if (ratio < gap_ratio) then
           gap_ratio = ratio
           best_gap_pos = i
        end if
     end do
     if (config%threshold <= 0.0_R8) then
        ! stable partition WITHIN the eligible (positive-weight) channels:
        ! upper-support (mask) channels first, gauge extras after, weight
        ! order preserved on each side -- the mask is then a leading block
        block
          integer :: iw, jw
          integer :: ord(nelig)
          real(R8) :: wtmp(nelig)
          integer :: atmp(nelig), ltmp(nelig)
          logical :: ftmp(nelig)
          jw = 0
          do iw = 1, nelig
             if (mask_flag(iw)) then
                jw = jw + 1; ord(jw) = iw
             end if
          end do
          do iw = 1, nelig
             if (.not. mask_flag(iw)) then
                jw = jw + 1; ord(jw) = iw
             end if
          end do
          block
            integer :: rtmp(nelig)
            real(R8) :: sptmp(nelig)
            wtmp = weight(1:nelig); atmp = chan_atom(1:nelig)
            ltmp = chan_l(1:nelig); ftmp = mask_flag(1:nelig)
            rtmp = chan_rf(1:nelig); sptmp = chan_sup(1:nelig)
            do iw = 1, nelig
               weight(iw) = wtmp(ord(iw)); chan_atom(iw) = atmp(ord(iw))
               chan_l(iw) = ltmp(ord(iw)); mask_flag(iw) = ftmp(ord(iw))
               chan_rf(iw) = rtmp(ord(iw)); chan_sup(iw) = sptmp(ord(iw))
            end do
          end block
        end block
     end if

     ! Feasibility guard: the downstream pipeline needs Nproj <= Nb
     ! (SCDM gauge and Wannier90 disentanglement both require
     ! num_wann <= num_bands).  Stop the selection before it exceeds
     ! the window.  A single over-complete channel (e.g. a full l
     ! shell wider than a narrow window) is allowed through and left
     ! to the legacy over-complete warning below.
     nproj_total = 0
     do i = 1, nsel
        n_add = (2*chan_l(i)+1) * stru%mult(chan_atom(i))
        if (i > 1 .and. nproj_total + n_add > Nb) then
           if (is_root) write(unit_out, '(A,I3,A,I4,A)') &
                '   feasibility: selection stopped after channel', i-1, &
                ' (next would exceed Nb =', Nb, ')'
           nsel = i - 1
           exit
        end if
        nproj_total = nproj_total + n_add
     end do
  end if

  ! cliff mode: record the mask boundary = columns of the leading
  ! mask-flagged channels among those actually SELECTED (post-feasibility)
  if (config%threshold <= 0.0_R8) then
     block
       integer :: iw, mask_cols
       mask_cols = 0
       do iw = 1, nsel
          if (.not. mask_flag(iw)) exit
          mask_cols = mask_cols + (2*chan_l(iw)+1)*stru%mult(chan_atom(iw))
       end do
       if (is_root) write(unit_out, '(A,I3,A,I3,A,I4,A)') &
            '   cliff (support-otsu): ', count(mask_flag(1:nsel)), &
            ' mask channels of ', nsel, ' selected; mask columns = ', &
            mask_cols, ' (recorded in projmap)'
       inwf%auto_config%mask_cols = mask_cols
     end block
  end if

  ! Ensure at least one channel selected
  if (nsel == 0) then
     nsel = 1
     nproj_total = (2*chan_l(1)+1) * stru%mult(chan_atom(1))
  end if

  ! Export the selection payload for in-process consumers (HYBR one-shot,
  ! capacity gate): ranked channel table + mask identity.  File mirrors
  ! (projmap.dat, diagnostics) are unchanged.
  block
    integer :: iw
    associate (ac => inwf%auto_config)
      if (allocated(ac%sel_atom)) deallocate(ac%sel_atom, ac%sel_l, &
           ac%sel_rf, ac%sel_mask, ac%sel_weight, ac%sel_support)
      ac%nsel = nsel
      allocate(ac%sel_atom(nsel), ac%sel_l(nsel), ac%sel_rf(nsel), &
               ac%sel_mask(nsel), ac%sel_weight(nsel), ac%sel_support(nsel))
      do iw = 1, nsel
         ac%sel_atom(iw)    = chan_atom(iw)
         ac%sel_l(iw)       = chan_l(iw)
         ac%sel_rf(iw)      = chan_rf(iw)
         ac%sel_mask(iw)    = mask_flag(iw)
         ac%sel_weight(iw)  = weight(iw)
         ac%sel_support(iw) = chan_sup(iw)
      end do
    end associate
  end block

  ! Sanity check: warn if Nproj != Nb (entangled or over-complete)
  if (is_root) then
     if (nproj_total < Nb) then
        write(unit_out, '(A)') &
             ' WARNING: Nproj < Nb (entangled bands).'
        write(unit_out, '(A,I4,A,I4)') &
             '   Nproj = ', nproj_total, ', Nb = ', Nb
        write(unit_out, '(A)') &
             '   Wannier90 will require disentanglement parameters'
        write(unit_out, '(A)') &
             '   (dis_win_min/max, dis_froz_min/max) in the .win file.'
        write(unit_out, '(A)') &
             '   Consider adjusting the threshold or band window.'
     else if (nproj_total > Nb) then
        write(unit_out, '(A)') &
             ' WARNING: Nproj > Nb (over-complete projection set).'
        write(unit_out, '(A,I4,A,I4)') &
             '   Nproj = ', nproj_total, ', Nb = ', Nb
        write(unit_out, '(A)') &
             '   Wannier90 will require disentanglement parameters'
        write(unit_out, '(A)') &
             '   or a wider band window to accommodate all projections.'
     end if

     ! Check per-band projectability: warn about poorly represented bands
     band_proj_min = huge(1.0_R8)
     do i = 1, Nb
        band_proj = 0.0_R8
        do j = 1, nsel
           band_proj = band_proj + proj_char(i, chan_atom(j), chan_l(j))
        end do
        band_proj_min = min(band_proj_min, band_proj)
     end do
     if (band_proj_min < 0.05_R8) then
        write(unit_out, '(A)') &
             ' WARNING: Some bands have very low projectability'
        write(unit_out, '(A,F8.4)') &
             '   onto the selected channels (min = ', band_proj_min
        write(unit_out, '(A)') &
             '   ). These bands may not be well represented.'
        write(unit_out, '(A)') &
             '   Consider widening the band window or lowering threshold.'
     end if
  end if

  ! Print diagnostic info
  if (is_root) then
     write(unit_out, *)
     write(unit_out, '(A)') ' AUTO-PROJECTION: Projectability analysis'
     write(unit_out, '(A,F10.3)') '   Threshold (w/w_max): ', config%threshold
     ! (the old negative-max_proj smuggle wrote through the inwf dummy,
     ! so this read of the config copy never saw it: the gap ratio always
     ! printed when the USER max_proj was 0, cliff mode included)
     if (config%max_proj == 0 .and. best_gap_pos > 0) then
        write(unit_out, '(A,F10.4)') '   Largest gap ratio:   ', gap_ratio
        write(unit_out, '(A,I4)')    '   Gap after channel:   ', best_gap_pos
     end if
     write(unit_out, '(A,I4)')   '   Channels found: ', nchan
     write(unit_out, '(A,I4)')   '   Channels selected: ', nsel
     write(unit_out, '(A,I4)')   '   Total projections: ', nproj_total
     write(unit_out, *)
     write(unit_out, '(A)') '   Ranked (atom, l) channels:'
     write(unit_out, '(A)') '   atom   l   weight   selected'
     write(unit_out, '(A)') '   ' // repeat('-', 35)
     do i = 1, min(nchan, nsel + 5)
        if (i <= nsel) then
           write(unit_out, '(I7,I4,F10.4,A,A)') &
                chan_atom(i), chan_l(i), weight(i), '   *', &
                merge('  [LO-shell]', '            ', chan_rf(i) > 1)
        else
           write(unit_out, '(I7,I4,F10.4,A)') &
                chan_atom(i), chan_l(i), weight(i), &
                merge('       [LO-shell]', '                 ', chan_rf(i) > 1)
        end if
     end do
     write(unit_out, *)
  end if

  ! Write diagnostic file if requested
  if (config%write_proj .and. is_root) then
     call write_autoprojection_diagnostics(stru, proj_char, Nb, inwf, &
          nchan, chan_atom, chan_l, weight, nsel, config)
  end if

  ! Generate projection specifications
  inwf%Nproj = nproj_total
  allocate(inwf%projections(nproj_total))

  w = 0
  lfirst = 1
  ! We need lfirst for each atom; recompute from structure
  do i = 1, nsel
     jatom = chan_atom(i)
     l = chan_l(i)

     ! Compute lfirst for this atom
     lfirst = 1
     do j = 1, jatom - 1
        lfirst = lfirst + stru%mult(j)
     end do

     ! Generate one projection per (equivalent atom, m) combination
     do mu = 1, stru%mult(jatom)
        latom = lfirst - 1 + mu
        do m = -l, l
           w = w + 1
           inwf%projections(w)%NY = 1
           allocate(inwf%projections(w)%iat(1))
           allocate(inwf%projections(w)%l(1))
           allocate(inwf%projections(w)%m(1))
           allocate(inwf%projections(w)%coeff(1))
           allocate(inwf%projections(w)%rf(1))
           inwf%projections(w)%iat(1) = latom
           inwf%projections(w)%l(1) = l
           inwf%projections(w)%m(1) = m
           inwf%projections(w)%coeff(1) = cmplx(1.0_DPk, 0.0_DPk, DPk)
           inwf%projections(w)%rf(1) = chan_rf(i)
        end do
     end do
  end do

  deallocate(weight, chan_atom, chan_l, chan_rf, chan_sup)

end subroutine generate_projections


subroutine write_autoprojection_diagnostics(stru, proj_char, Nb, inwf, &
     nchan, chan_atom, chan_l, weight, nsel_auto, config)
  !! Write a diagnostic .autoproj file with:
  !!  - Full channel ranking with gap analysis
  !!  - Per-band projectability breakdown
  !!  - Alternative projection sets at each gap position
  !!  - Ready-to-use .inwf templates for manual fallback
  use param,     only: Lmax2, unit_out
  use const,     only: R8, DPk
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t, auto_config_t

  implicit none

  type(struct_t),      intent(in) :: stru
  real(R8),            intent(in) :: proj_char(:,:,0:)
  integer,             intent(in) :: Nb
  type(inwf_t),        intent(in) :: inwf
  integer,             intent(in) :: nchan
  integer,             intent(in) :: chan_atom(:), chan_l(:)
  real(R8),            intent(in) :: weight(:)
  integer,             intent(in) :: nsel_auto
  type(auto_config_t), intent(in) :: config

  integer, parameter :: UNIT_DIAG = 98
  character(len=6), parameter :: LNAME(0:5) = &
       (/ 's     ', 'p     ', 'd     ', 'f     ', 'g     ', 'h     ' /)
  character(len=80) :: atom_name

  integer  :: i, j, ib, jatom, l, m, mu, lfirst, latom
  integer  :: nproj_at_cut, n_active
  real(R8) :: gap_r, rel_w, cumul, max_weight
  real(R8) :: total_proj_band
  real(R8) :: thresh_needed

  open(UNIT_DIAG, FILE='autoprojection.diag', STATUS='replace', &
       FORM='formatted')

  max_weight = weight(1)

  !-----------------------------------------------------------------
  ! Section 1: Header and full channel ranking
  !-----------------------------------------------------------------
  write(UNIT_DIAG, '(A)') '============================================================'
  write(UNIT_DIAG, '(A)') ' AUTO-PROJECTION DIAGNOSTIC REPORT'
  write(UNIT_DIAG, '(A)') '============================================================'
  write(UNIT_DIAG, *)
  write(UNIT_DIAG, '(A,I4,A,I4,A)') ' Band window: [', inwf%bmin, ', ', &
       inwf%bmax, ']'
  write(UNIT_DIAG, '(A,I4)')   ' Number of bands (Nb): ', Nb
  write(UNIT_DIAG, '(A,I4)')   ' Number of atoms:      ', stru%Nneq
  write(UNIT_DIAG, '(A,F8.4)') ' Threshold used:       ', config%threshold
  write(UNIT_DIAG, '(A,I4)')   ' Channels selected:    ', nsel_auto
  write(UNIT_DIAG, *)
  write(UNIT_DIAG, '(A)') ' CHANNEL RANKING (sorted by weight, descending)'
  write(UNIT_DIAG, '(A)') ' -------------------------------------------------------'
  write(UNIT_DIAG, '(A)') '  rank  atom   l  orb   weight   w/w_max  gap_ratio  sel'
  write(UNIT_DIAG, '(A)') ' -------------------------------------------------------'
  do i = 1, nchan
     if (weight(i) <= 0.0_R8) exit
     rel_w = weight(i) / max_weight
     if (i > 1 .and. weight(i-1) > 0.0_R8) then
        gap_r = weight(i) / weight(i-1)
     else
        gap_r = 1.0_R8
     end if

     ! Get atom name from structure
     atom_name = stru%aname(chan_atom(i))

     if (i <= nsel_auto) then
        write(UNIT_DIAG, '(I5,I6,I4,2X,A2,F10.4,F9.4,F11.4,A)') &
             i, chan_atom(i), chan_l(i), LNAME(chan_l(i))(1:2), &
             weight(i), rel_w, gap_r, '   <--'
     else
        write(UNIT_DIAG, '(I5,I6,I4,2X,A2,F10.4,F9.4,F11.4)') &
             i, chan_atom(i), chan_l(i), LNAME(chan_l(i))(1:2), &
             weight(i), rel_w, gap_r
     end if
  end do
  ! Count active channels (weight > 0)
  n_active = 0
  do i = 1, nchan
     if (weight(i) > 0.0_R8) n_active = n_active + 1
  end do
  write(UNIT_DIAG, '(A)') ' -------------------------------------------------------'
  write(UNIT_DIAG, *)

  !-----------------------------------------------------------------
  ! Section 2: Per-band projectability
  !-----------------------------------------------------------------
  write(UNIT_DIAG, '(A)') ' PER-BAND PROJECTABILITY (summed over all atoms and l)'
  write(UNIT_DIAG, '(A)') ' -------------------------------------------------------'
  write(UNIT_DIAG, '(A)', ADVANCE='no') '  band '
  do jatom = 1, stru%Nneq
     do l = 0, Lmax2
        if (sum(proj_char(1:Nb, jatom, l)) > 0.0_R8) then
           write(UNIT_DIAG, '(I2,A1,A1,A1)', ADVANCE='no') &
                jatom, '-', LNAME(l)(1:1), ' '
        end if
     end do
  end do
  write(UNIT_DIAG, '(A)') '   total'
  write(UNIT_DIAG, '(A)') ' -------------------------------------------------------'
  do ib = 1, Nb
     write(UNIT_DIAG, '(I5,A)', ADVANCE='no') ib + inwf%bmin - 1, '  '
     total_proj_band = 0.0_R8
     do jatom = 1, stru%Nneq
        do l = 0, Lmax2
           if (sum(proj_char(1:Nb, jatom, l)) > 0.0_R8) then
              write(UNIT_DIAG, '(F5.3,A)', ADVANCE='no') &
                   proj_char(ib, jatom, l), ' '
              total_proj_band = total_proj_band + proj_char(ib, jatom, l)
           end if
        end do
     end do
     write(UNIT_DIAG, '(F7.3)') total_proj_band
  end do
  write(UNIT_DIAG, *)

  !-----------------------------------------------------------------
  ! Section 3: Alternative projection sets
  !-----------------------------------------------------------------
  write(UNIT_DIAG, '(A)') ' ALTERNATIVE PROJECTION SETS'
  write(UNIT_DIAG, '(A)') ' ============================================================'
  write(UNIT_DIAG, '(A)') ' Each set includes all channels above a given cutoff.'
  write(UNIT_DIAG, '(A)') ' Use the AUTO threshold or max_proj values shown, or copy'
  write(UNIT_DIAG, '(A)') ' the BOTH-mode .inwf template to use manually.'
  write(UNIT_DIAG, *)

  ! Generate alternatives at each gap position
  do i = 1, min(n_active, 8)
     ! Count projections if we cut after channel i
     nproj_at_cut = 0
     do j = 1, i
        nproj_at_cut = nproj_at_cut + (2*chan_l(j)+1) * stru%mult(chan_atom(j))
     end do

     ! Compute the threshold that would select exactly this many channels
     if (i < n_active .and. weight(i+1) > 0.0_R8) then
        ! Threshold must be between weight(i+1)/max and weight(i)/max
        thresh_needed = (weight(i+1) / max_weight + weight(i) / max_weight) / 2
     else
        thresh_needed = weight(i) / (2 * max_weight)
     end if

     write(UNIT_DIAG, '(A,I2,A)') ' --- Option ', i, &
          ' ---------------------------------------------------'
     write(UNIT_DIAG, '(A,I4)') '  Channels: ', i
     write(UNIT_DIAG, '(A,I4)') '  Projections: ', nproj_at_cut
     write(UNIT_DIAG, '(A)', ADVANCE='no') '  Orbitals: '
     do j = 1, i
        atom_name = stru%aname(chan_atom(j))
        write(UNIT_DIAG, '(A,A1,A,A)', ADVANCE='no') &
             trim(adjustl(atom_name)), '(', &
             trim(LNAME(chan_l(j))), ') '
     end do
     write(UNIT_DIAG, *)
     if (i == nsel_auto) then
        write(UNIT_DIAG, '(A)') '  ** CURRENTLY SELECTED **'
     end if
     write(UNIT_DIAG, '(A,F8.4)') '  AUTO threshold for this set: ', &
          thresh_needed
     write(UNIT_DIAG, '(A,I4)')   '  AUTO max_proj for this set:  ', &
          nproj_at_cut
     write(UNIT_DIAG, *)

     ! Write AUTO .inwf template
     write(UNIT_DIAG, '(A)') '  AUTO .inwf:'
     write(UNIT_DIAG, '(A)') '  ----'
     write(UNIT_DIAG, '(A)') '  AUTO'
     write(UNIT_DIAG, '(2X,I3,I4)') inwf%bmin, inwf%bmax
     write(UNIT_DIAG, '(2X,I3)') inwf%LJmax
     write(UNIT_DIAG, '(2X,F6.4,I3,A2)') thresh_needed, 0, ' T'
     write(UNIT_DIAG, '(A)') '  ----'
     write(UNIT_DIAG, *)

     ! Write BOTH-mode .inwf template
     write(UNIT_DIAG, '(A)') '  BOTH-mode .inwf (for manual use):'
     write(UNIT_DIAG, '(A)') '  ----'
     write(UNIT_DIAG, '(A)') '  BOTH'
     write(UNIT_DIAG, '(2X,I3,I4)') inwf%bmin, inwf%bmax
     write(UNIT_DIAG, '(2X,I3,I4)') inwf%LJmax, nproj_at_cut

     do j = 1, i
        jatom = chan_atom(j)
        l = chan_l(j)
        lfirst = 1
        do mu = 1, jatom - 1
           lfirst = lfirst + stru%mult(mu)
        end do
        do mu = 1, stru%mult(jatom)
           latom = lfirst - 1 + mu
           do m = -l, l
              write(UNIT_DIAG, '(A)') '  1'
              write(UNIT_DIAG, '(2X,I2,I3,I3,2F14.8)') &
                   latom, l, m, 1.0_DPk, 0.0_DPk
           end do
        end do
     end do
     write(UNIT_DIAG, '(A)') '  ----'
     write(UNIT_DIAG, *)
  end do

  close(UNIT_DIAG)

  write(unit_out, '(A)') ' Diagnostic file written: autoprojection.diag'

end subroutine write_autoprojection_diagnostics

end module autoprojection_m


module     planew_m; contains
subroutine planew(stru, NB, num_kpts, Nntot, maxx, maxy, maxz)
  use param,     only: unit_out
  use xa3,       only: gx,gy,gz, a, vecsz
  use amn_mmn,   only: overlap
  use pairs,     only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1
  use lolog,     only: nlo,nlov,nlon
  use const,     only: R8, C16, TAU
  use structmod, only: struct_t
  use mpi_w2w

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
  if (is_root) write(unit_out,*)'MAX Q:',maxx,maxy,maxz
  call planew_check_grid(num_kpts, NNTOT, maxx, maxy, maxz)
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
  if (is_root) &
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
  call planew_mask_extend(a_null_fin, num_kpts, NNTOT, maxx, maxy, maxz, m2x, m2y)

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
  if (is_root) write(unit_out,*)'t_init=',t2-t1

  ! Computation loop: distribute k1 across MPI ranks
  do k1=1,num_kpts
#ifdef MPI
     if (k_owner(k1, num_kpts) /= myrank) cycle
#endif
     pair_index = (k1-1)*Nntot + 1
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
     if (is_root .and. k1 == 1) &
          write(unit_out,'(A,I6,A,I6,A,I10,A,L1)') &
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
        pair_index = (k1-1)*Nntot + k2
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
     end do
     call cpu_time(tt3)
     t_b=t_b+tt2-tt1
     t_prod=t_prod+tt3-tt2
     if (is_root .and. k1 == 1) &
          write(unit_out,*)'PLANEW TIME PER k-point:',t_b+t_prod
  end do
  if (is_root) then
     write(unit_out,*)'t_b=',t_b
     write(unit_out,*)'t_prod=',t_prod
  end if

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

!> Completes planew's interstitial mask.  Step 2 of planew reads B_ at every plane wave of
!> the neighbour, 1..vecsz(kk)-(nlo+nlon+nlov), but the mask loop marks the neighbour's
!> positions only up to the head's count vecsz(kkk).  When a neighbour has more plane waves,
!> the positions no other neighbour marked stay unmarked, B_ is 0 there and the pair loses
!> part of its interstitial overlap.  First pass: report, per
!> pair, the positions left unmarked; second pass: mark them.  A separate, non-inlined
!> routine so that planew's own code generation is unchanged.
subroutine planew_mask_extend(a_null_fin, num_kpts, NNTOT, maxx, maxy, maxz, m2x, m2y)
  use param,     only: unit_out
  use xa3,       only: gx,gy,gz, vecsz
  use pairs,     only: kp, kpb, bqx,bqy,bqz
  use lolog,     only: nlo,nlov,nlon
  use mpi_w2w

  implicit none
!GCC$ ATTRIBUTES NOINLINE :: planew_mask_extend

  logical, intent(inout) :: a_null_fin(0:,:)
  integer, intent(in)    :: num_kpts, NNTOT, maxx, maxy, maxz, m2x, m2y

  integer :: k1, k2, kkk, kk, p, I, index, n_add

  do k1=1,num_kpts
     kkk=KP((k1-1)*NNTOT+1)
     do k2=1,NNTOT
        p=(k1-1)*NNTOT+k2
        kk=KPB(p)
        n_add=0
        do I=vecsz(kkk)+1,vecsz(kk)-(nlo+nlon+nlov)
           index=((GZ(I,kk)+maxz-BQZ(p))*m2y+GY(I,kk)+maxy-BQY(p))*m2x+GX(I,kk)+maxx-BQX(p)
           if (.not. a_null_fin(index,kkk)) n_add=n_add+1
        end do
        if (n_add > 0 .and. is_root) write(unit_out,'(A,I6,A,I6,A,3I3,A,I5,A)') &
             ' planew: head k',kkk,' neighbour k',kk,' G',BQX(p),BQY(p),BQZ(p),':',n_add, &
             ' neighbour plane waves added to the interstitial mask'
     end do
     do k2=1,NNTOT
        p=(k1-1)*NNTOT+k2
        kk=KPB(p)
        do I=vecsz(kkk)+1,vecsz(kk)-(nlo+nlon+nlov)
           index=((GZ(I,kk)+maxz-BQZ(p))*m2y+GY(I,kk)+maxy-BQY(p))*m2x+GX(I,kk)+maxx-BQX(p)
           a_null_fin(index,kkk)=.true.
        end do
     end do
  end do
end subroutine planew_mask_extend

!> Fails closed when an index planew forms falls outside its grid.  planew places a plane wave at
!> component index G + maxx - G_b on a grid of 2*maxx+1 points per direction, where maxx is one more
!> than the largest |G| component read_vec reported.  An index outside 0..2*maxx would put the plane
!> wave at another position (Step 2 checks only the linear index) or write outside a_null, a_null_fin
!> and A_.  Checked for every head's plane waves (no shift) and every neighbour's plane waves with
!> that pair's shift: the positions the mask loops, planew_mask_extend and Step 2 form.
subroutine planew_check_grid(num_kpts, NNTOT, maxx, maxy, maxz)
  use param,     only: unit_out
  use xa3,       only: gx,gy,gz, vecsz
  use pairs,     only: kp, kpb, bqx,bqy,bqz
  use clio,      only: croak
  use mpi_w2w

  implicit none
!GCC$ ATTRIBUTES NOINLINE :: planew_check_grid

  integer, intent(in) :: num_kpts, NNTOT, maxx, maxy, maxz

  integer :: p, k, I, ix, iy, iz, nbad, p1, i1
  character(len=400) :: msg

  nbad=0; p1=0; i1=0
  do k=1,num_kpts
     do I=1,vecsz(k)
        ix=GX(I,k)+maxx; iy=GY(I,k)+maxy; iz=GZ(I,k)+maxz
        if (ix<0 .or. ix>2*maxx .or. iy<0 .or. iy>2*maxy .or. iz<0 .or. iz>2*maxz) nbad=nbad+1
     end do
  end do
  do p=1,num_kpts*NNTOT
     k=KPB(p)
     do I=1,vecsz(k)
        ix=GX(I,k)+maxx-BQX(p); iy=GY(I,k)+maxy-BQY(p); iz=GZ(I,k)+maxz-BQZ(p)
        if (ix<0 .or. ix>2*maxx .or. iy<0 .or. iy>2*maxy .or. iz<0 .or. iz>2*maxz) then
           nbad=nbad+1
           if (p1 == 0) then
              p1=p; i1=I
           end if
        end if
     end do
  end do
  if (nbad == 0) return
  if (p1 > 0) then
     write(msg,'(A,I0,A,I0,A,I0,A,3I4,A,3I3,A,3I4,A)') 'planew: ', nbad, &
          ' plane-wave indices fall outside the interstitial grid; first at pair ', p1, &
          ' (neighbour k ', KPB(p1), '): G', GX(i1,KPB(p1)), GY(i1,KPB(p1)), GZ(i1,KPB(p1)), &
          ', shift', BQX(p1), BQY(p1), BQZ(p1), ', grid half-widths', maxx, maxy, maxz, '.'
  else
     write(msg,'(A,I0,A,3I4,A)') 'planew: ', nbad, &
          ' head plane-wave indices fall outside the interstitial grid; grid half-widths', maxx, maxy, maxz, '.'
  end if
  if (is_root) then
     write(unit_out,'(A)') trim(msg)
     write(unit_out,'(A)') ' planew: the grid half-width is the largest |G| component read_vec reported, plus one.'
     write(unit_out,'(A)') ' planew: alternatives: (1) split vector files (w2w case.def N): read_vec must report the'
     write(unit_out,'(A)') ' planew:     maximum over ALL files; (2) join the vector files and run on one file;'
     write(unit_out,'(A)') ' planew:     (3) a b-vector shift larger than the margin needs a wider grid (report the .nnkp).'
  end if
  call croak(trim(msg)//' Refusing: the interstitial Mmn would be wrong. See case.outputwf.')
end subroutine planew_check_grid
end module planew_m


module     read_vec_m; contains
!!!    Read a 'vector' file
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2015 Elias Assmann
subroutine read_vec(NEmin, NEmax, Nneq, kkk, maxx, maxy, maxz, Efermi, &
                    write_eig)
  use param,  only: unit_out, unit_vector
  use w2w,    only: unit_eig
  use const,  only: R8, Ryd_eV
  use xa3,    only: XK,YK,ZK, GX,GY,GZ, vecsz, A, Eb
  use mpi_w2w

  implicit none

  integer, intent(in)    :: NEmin, NEmax, Nneq
  ! maxx, maxy, maxz accumulate the largest |G| component over every call, like kkk accumulates the
  ! k count: main.F reads k-parallel vector files (w2w case.def N) one call per file, and planew sizes
  ! its grid from these.  The caller zeroes them once before the first vector-file loop.
  integer, intent(inout) :: maxx, maxy, maxz
  integer, intent(inout) :: kkk
  real(R8),intent(in)    :: efermi
  ! suppress the .eig side effect (sigma-2 pass of the SOC loop: the
  ! eigenvalues are shared, a second pass would append duplicates)
  logical, intent(in), optional :: write_eig

  integer :: i, j, n, NB, NE, num, E_size
  logical :: weig

  real(r8), allocatable :: E(:)

  character(10)    BNAME

  E_size = 0
  weig = .true.
  if (present(write_eig)) weig = write_eig

  do I = 1, Nneq
     read(unit_vector)
     read(unit_vector)
  enddo

  kpoint: do
     read(unit_vector,end=998) XK(kkk+1),YK(kkk+1),ZK(kkk+1),BNAME,N,NE
     ! Guard: with NE < bmax the loop below leaves E stale and the A
     ! columns for the missing bands uninitialized -> silent garbage in
     ! amn/mmn.  (main.F clamps bmax from the energy files, but the vector
     ! files can still disagree with them, e.g. after a partial rerun.)
     if (NE < NEmax) then
        if (is_root) write(unit_out,*) 'read_vec: only ', NE, &
             ' bands at k-point ', kkk+1, ' but bmax = ', NEmax
        error stop 'read_vec: bmax exceeds bands in vector file (vector/energy mismatch)'
     end if
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
           Eb(NUM-NEMIN+1, kkk) = E(num)
        else
           read(unit_vector)
        endif
     enddo

     if (is_root .and. weig) then
        do NUM=NEMIN,NEMAX
           NB=NUM-NEMIN+1
           write(unit_eig,"(2I12,F22.16)")NB,kkk, Ryd_eV*(E(NUM)-efermi)
        enddo
     end if
  end do kpoint

998 if (allocated(E)) deallocate(E)
  if (is_root) write(unit_out,*)'vector read in',kkk
end subroutine read_vec
end module read_vec_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
