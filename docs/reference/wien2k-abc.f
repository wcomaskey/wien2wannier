!!! WIEN2k Reference Implementation - SRC_lapw1/abc.f
!!! This is the CORRECT implementation that w2w should match
!!! 
!!! Key differences from w2w:
!!! 1. PI12LO, PE12LO, PILOLO are multi-dimensional arrays
!!! 2. All indexed by (J=l, jlo, JATOM=atom)
!!! 3. PILOLO includes self-overlap, not hardcoded to 1.0
!!!

      SUBROUTINE ABC (JATOM,J,jlo,lapw)
!
      use parallel, only: MYID
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, &
                        DPELO, PEILO, PI12LO, PE12LO, PILOLO
      use atspdt, only  : INS, NL, LMMAX, LQIND, LM, LQNS, DP, DPE, E, P, PE, PEI
      use struk, only : POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      logical            lapw
      INTEGER            J,  JATOM, JLO
!
!     ..................................................................
!
!        ABC calculates the cofficients A,B,C of the local orbitals (LO)
!
!        KEY: Uses multi-dimensional arrays for overlap integrals
!             PI12LO(J,jlo,JATOM) - NOT scalar!
!             PE12LO(J,jlo,JATOM) - NOT scalar!
!             PILOLO(J,jlo,jlo,JATOM) - LO self-overlap
!
!     ..................................................................
!
!        Local Parameters
!
      DOUBLE PRECISION   CUTOFF
      PARAMETER          (CUTOFF = 200.0D+0)
!
!        Local Scalars
!
      DOUBLE PRECISION   XAC, XBC, ALONORM
!
!        Intrinsic Functions
!
      INTRINSIC          SQRT
!

      if(lapw) then
         ! Match energy derivative (u̇) and slope
         XAC = PLO(J,jlo,JATOM)*DPE(J+1,JATOM) -  &
              DPLO(J,jlo,JATOM)*PE(J+1,JATOM)
         XAC = XAC*RMT(JATOM)*RMT(JATOM)
         
         ! Match wavefunction (u) and slope
         XBC = PLO(J,jlo,JATOM)*DP(J+1,JATOM) -  &
              DPLO(J,jlo,JATOM)*P(J+1,JATOM)
         XBC = -XBC*RMT(JATOM)*RMT(JATOM)
         
         ! CRITICAL: Normalization includes ALL overlap terms
         ! Note: Uses PI12LO(J,jlo,JATOM) - indexed by THIS jlo!
         !       Uses PE12LO(J,jlo,JATOM) - indexed by THIS jlo!
         !       Uses PILOLO(J,jlo,jlo,JATOM) - LO self-overlap, NOT 1.0!
         CLO(J,jlo,JATOM) = XAC*(XAC + 2.0D+0*PI12LO(J,jlo,JATOM)) + &
              XBC*(XBC*PEI(J+1,JATOM) + 2.0D+0*PE12LO(J,jlo,JATOM)) + &
              PILOLO(J,jlo,jlo,JATOM)  ! ← w2w had hardcoded 1.0D0 here!
              
         CLO(J,jlo,JATOM) = 1.0D0/SQRT(CLO(J,jlo,JATOM))
         CLO(J,jlo,JATOM) = MIN(CLO(J,jlo,JATOM),CUTOFF)
         
         ! Final coefficients
         ALO(J,jlo,JATOM) = CLO(J,jlo,JATOM)*XAC
         BLO(J,jlo,JATOM) = CLO(J,jlo,JATOM)*XBC
      else
!.....APW+lo definitions
         if(jlo.eq.1) then
            ! First LO: orthogonalize to u only (no u̇ for APW)
            alonorm=sqrt(1.d0+(P(J+1,JATOM)/PE(J+1,JATOM))**2*PEI(J+1,JATOM))
            ALO(J,jlo,JATOM) = 1.d0 /alonorm 
            BLO(J,jlo,JATOM) = -P(J+1,JATOM)/PE(J+1,JATOM)/alonorm
            CLO(J,jlo,JATOM) = 0.d0
         else 
            ! Second and higher LOs: orthogonalize to u and first LO
            xbc=-P(J+1,JATOM)/PLO(J,jlo,JATOM)
            
            ! CRITICAL: For multiple LOs, normalization includes:
            ! - PILOLO(J,jlo,jlo,JATOM): self-overlap
            ! - PI12LO(J,jlo,JATOM): overlap with u
            ! w2w had: sqrt(1+xbc**2+2*xbc*PI12LO) - missing PILOLO factor!
            xac=sqrt(1+PILOLO(J,jlo,jlo,JATOM)*xbc**2+2*xbc*PI12LO(J,jlo,JATOM))
            
            ALO(J,jlo,JATOM) = 1.d0/xac
            BLO(J,jlo,JATOM) = 0.d0
            CLO(J,jlo,JATOM) = xbc/xac
        endif
      endif
      
      ! Output for verification
      if(myid.eq.0) WRITE(6,6000)J,ALO(J,jlo,JATOM),BLO(J,jlo,JATOM),CLO(J,jlo,JATOM)
!
      RETURN
!
 6000 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
!
!        End of 'ABC'
!
      END

!!! Notes on array dimensions in WIEN2k:
!!!
!!! From loabc module:
!!!   ALO(0:LOMAX, NLOAT, nat_1)  - LO coefficients A
!!!   BLO(0:LOMAX, NLOAT, nat_1)  - LO coefficients B  
!!!   CLO(0:LOMAX, NLOAT, nat_1)  - LO coefficients C
!!!   PI12LO(0:LOMAX, NLOAT, nat_1)  - <u | LO>
!!!   PE12LO(0:LOMAX, NLOAT, nat_1)  - <u̇ | LO>
!!!   PILOLO(0:LOMAX, NLOAT, NLOAT, nat_1) - <LO_i | LO_j>
!!!
!!! Where:
!!!   LOMAX = maximum l for LOs (typically 3)
!!!   NLOAT = maximum number of LOs per l (typically 3)
!!!   nat_1 = number of atoms
!!!
!!! The key insight: EVERY overlap integral is stored separately
!!! for each (l, jlo, atom) combination!