# Task 1: Fix Multiple Local Orbital Coefficients

## Problem Statement

### Current Behavior
When multiple Local Orbitals (LOs) exist for the same angular momentum l, w2w produces:
- Infinity or NaN values in LO coefficient calculations
- Error in `abc()` subroutine: `sqrt()` of negative number
- Crash or corrupted output files

### Root Cause
Overlap integrals between radial functions are incorrectly implemented:

**Current (WRONG)**:
```fortran
! In atpar subroutine - modules.f around line 200
real(R8) :: pi12lo, pe12lo   ! SCALARS - get overwritten for each jlo!

! Later in loop:
do jlo=1,ilo(l)
  call RINT13(..., pi12lo)   ! Overwrites previous value!
  call abc(..., pi12lo, ...)  ! Uses wrong value for jlo>1
enddo
```

**Correct (WIEN2k)**:
```fortran
! WIEN2k stores separately for each LO:
PI12LO(J,jlo,JATOM)     ! 3D array indexed by (l, jlo, atom)
PE12LO(J,jlo,JATOM)     ! 3D array
PILOLO(J,jlo,jlo,JATOM) ! 4D array for LO self-overlap
```

### Impact
- **Systems affected**: Any material with multiple LOs per l
  - Transition metals with semicore states (3p + 3d LOs)
  - Rare earths with 4f + 5d LOs
  - Accurate treatment of shallow core states
- **Severity**: CRITICAL - blocks Wannierization of many important materials

## WIEN2k Reference Implementation

See `docs/reference/wien2k-abc.f` for the correct implementation.

**Key differences**:

| Aspect | w2w (wrong) | WIEN2k (correct) |
|--------|-------------|------------------|
| `pi12lo` dimension | scalar | `(l, jlo, atom)` |
| `pe12lo` dimension | scalar | `(l, jlo, atom)` |
| `pilolo` exists? | NO | YES - `(l, jlo1, jlo2, atom)` |
| Calculation | Once, overwritten | Separately for each jlo |
| Usage in abc | Wrong value | Correct value for each jlo |

**Critical WIEN2k code sections**:
```fortran
! From WIEN2k's abc.f:
CLO(J,jlo,JATOM) = XAC*(XAC + 2.0D+0*PI12LO(J,jlo,JATOM)) + &
     XBC*(XBC*PEI(J+1,JATOM) + 2.0D+0*PE12LO(J,jlo,JATOM)) + &
     PILOLO(J,jlo,jlo,JATOM)  ! ← Note: uses PILOLO, not 1.0

! For non-LAPW with jlo>1:
xac=sqrt(1+PILOLO(J,jlo,jlo,JATOM)*xbc**2+2*xbc*PI12LO(J,jlo,JATOM))
```

## Implementation Plan

### Step 1: Update `loabc` Module Data Structures
**File**: `modules.f`, lines ~50-80

**Current code**:
```fortran
module loabc
  use param, only: lomax, Nloat, Nrf
  use const, only: R8

  implicit none
  private; save

  ! abc calculates the cofficients a,b,c of the lo
  real(R8), public :: alo(0:LOmax, Nloat, Nrf)
end module loabc
```

**Required changes**:
1. Add allocatable arrays for overlap integrals
2. Add initialization routine
3. Add cleanup routine

**New code**:
```fortran
module loabc
  use param, only: lomax, Nloat, Nrf
  use const, only: R8

  implicit none
  private; save

  ! abc calculates the cofficients a,b,c of the lo
  real(R8), public :: alo(0:LOmax, Nloat, Nrf)
  
  ! ADD: Overlap integrals for multiple LOs
  real(R8), allocatable, public :: pi12lo(:,:)   ! (0:lomax, nloat)
  real(R8), allocatable, public :: pe12lo(:,:)   ! (0:lomax, nloat)
  real(R8), allocatable, public :: pilolo(:,:,:) ! (0:lomax, nloat, nloat)
  
contains
  ! ADD: Initialization
  subroutine init_loabc()
    implicit none
    
    if (.not. allocated(pi12lo)) allocate(pi12lo(0:lomax, nloat))
    if (.not. allocated(pe12lo)) allocate(pe12lo(0:lomax, nloat))
    if (.not. allocated(pilolo)) allocate(pilolo(0:lomax, nloat, nloat))
    
    pi12lo = 0.0_R8
    pe12lo = 0.0_R8
    pilolo = 0.0_R8
  end subroutine init_loabc
  
  ! ADD: Cleanup
  subroutine cleanup_loabc()
    implicit none
    
    if (allocated(pi12lo)) deallocate(pi12lo)
    if (allocated(pe12lo)) deallocate(pe12lo)
    if (allocated(pilolo)) deallocate(pilolo)
  end subroutine cleanup_loabc
end module loabc
```

**Testing**: Compile and verify no errors. Module should build cleanly.

---

### Step 2: Update `abc` Subroutine Signature
**File**: `modules.f`, lines ~100-150 (abc_m module)

**Current signature**:
```fortran
subroutine abc(l, rmt, pei, pi12lo, pe12lo, jlo, lapw)
  real(R8), intent(in) :: pi12lo, pe12lo  ! WRONG: scalars
```

**Required changes**:
1. Add `pilolo` parameter
2. Use it instead of hardcoded `1.0`

**New signature**:
```fortran
subroutine abc(l, rmt, pei, pi12lo_val, pe12lo_val, pilolo_val, jlo, lapw)
  use param,  only: unit_out, Nrf
  use loabc,  only: alo
  use atspdt, only: P, DP
  use const,  only: R8

  implicit none

  integer,  intent(in) :: l, jlo
  real(R8), intent(in) :: rmt, pei, pi12lo_val, pe12lo_val, pilolo_val
  logical,  intent(in) :: lapw
```

**Update calculations**:
```fortran
! CHANGE 1: LAPW case - use pilolo_val instead of 1.0
if (lapw) then
   irf=2+jlo
   xac=p(l,irf)*dp(l,2)-dp(l,irf)*p(l,2)
   xac= xac * rmt**2
   xbc=p(l,irf)*dp(l,1)-dp(l,irf)*p(l,1)
   xbc= -xbc * rmt**2
   
   ! CHANGED: Use pilolo_val instead of hardcoded 1.0D0
   clo=xac*(xac+2.0D0*pi12lo_val)+xbc* &
        (xbc*pei+2.0D0*pe12lo_val)+pilolo_val  ! ← WAS: +1.0D0
   clo=1.0D0/sqrt(clo)
   ! ... rest unchanged
else
   ! CHANGE 2: Non-LAPW case with jlo>1
   if (jlo.eq.1) then
      ! ... unchanged
   else
      xbc=-P(l,1)/P(l,1+jlo)
      ! CHANGED: Use pilolo_val instead of just xbc**2
      xac=sqrt(1.d0 + pilolo_val*xbc**2 + 2.0d0*xbc*pi12lo_val)
      ! ← WAS: sqrt(1+xbc**2+2*xbc*PI12LO)
      alo(l,jlo,1)=1.d0/xac
      alo(l,jlo,1+jlo)=xbc/xac
   end if
end if
```

**Testing**: Compile both real and complex versions. Check for syntax errors.

---

### Step 3: Update `atpar` Subroutine
**File**: `modules.f`, lines ~200-400 (atpar_m module)

**Current code** (WRONG):
```fortran
subroutine atpar(stru, jatom, itape, jtape)
  ! ...
  real(R8) :: pi12lo, pe12lo  ! WRONG: scalars
  
  ! ... later in LO loop:
  do jlo=1,ilo(l)
     if (lapw(l) .or. jlo/=1) then
        irf=irf+1
        ! ...
        call RINT13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
             &      rf1(:,l,irf), rf2(:,l,irf), pi12lo)  ! Overwrites!
        call RINT13(stru, jatom, rf1(1,l,2), rf2(1,l,2), &
             &      rf1(:,l,irf), rf2(:,l,irf), pe12lo)  ! Overwrites!
     end if
     call abc (l, stru%RMT(jatom), pei(l), pi12lo, pe12lo, jlo, lapw(l))
  end do
```

**Required changes**:
1. Remove local scalar declarations
2. Import module arrays
3. Calculate and store for each jlo separately
4. Calculate LO self-overlap (pilolo)
5. Pass correct values to abc

**New code**:
```fortran
subroutine atpar(stru, jatom, itape, jtape)
  use param,      only: unit_out, Nrad, Nloat, lomax, Lmax2
  use structmod,  only: struct_t
  use lolog,      only: nlo,nlov,nlon,loor,ilo,lapw,n_rad
  use atspdt,     only: P, DP
  use const,      only: R8, clight
  use uhelp,      only: A, B
  use radfu,      only: RF1, RF2
  use loabc,      only: init_loabc, pi12lo, pe12lo, pilolo  ! CHANGED: import

  ! ... existing includes ...

  implicit none

  type(struct_t), intent(in) :: stru
  integer,        intent(in) :: jatom, itape, jtape

  ! ... other variables ...
  
  ! REMOVED: real(R8) :: pi12lo, pe12lo  (now in module)

  ! ... existing code up to LO loop ...
  
  ! ADD: Initialize overlap arrays
  call init_loabc()

  ! ... existing code for reading energies, calculating LAPW functions ...

  ! LO loop - MODIFIED:
  loloop: do l=0,lomax
     irf=2
     iloloop: do jlo=1,ilo(l)
        if (lapw(l) .or. jlo/=1) then
           irf=irf+1
           DELE=2.0D-3
           DELEI=0.25D0/DELE
           FL=L
           EI=elo(l,jlo)/2.d0

           ! Calculate function at EI (existing code)
           if(rlo(jlo,l)) then
              ! ... diracout case ...
           else
              call outwin(stru, jatom, Vr, ei, fl, uv, duv, nodes)
           endif

           call RINT13(stru, jatom, A, B, A, B, OVLP)
           TRX=1.0d0/sqrt(OVLP)
           P(l,irf)=TRX*UV
           DP(l,irf)=TRX*DUV
           IMAX=stru%Npt(JATOM)
           n_rad(l)=irf
           do M=1,IMAX
              rf1(M,l,irf)=TRX*A(M)
              rf2(M,l,irf)=TRX*B(M)
           end do

           ! CHANGED: Calculate and STORE overlap integrals for THIS jlo
           call RINT13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
                &      rf1(:,l,irf), rf2(:,l,irf), pi12lo(l,jlo))
           call RINT13(stru, jatom, rf1(:,l,2), rf2(:,l,2), &
                &      rf1(:,l,irf), rf2(:,l,irf), pe12lo(l,jlo))
           
           ! ADDED: Calculate LO self-overlap
           call RINT13(stru, jatom, rf1(:,l,irf), rf2(:,l,irf), &
                &      rf1(:,l,irf), rf2(:,l,irf), pilolo(l,jlo,jlo))
        else
           ! For jlo=1 with non-LAPW, set defaults
           pi12lo(l,jlo) = 0.0d0
           pe12lo(l,jlo) = 0.0d0
           pilolo(l,jlo,jlo) = 1.0d0
        end if
        
        ! CHANGED: Pass the stored values for THIS jlo
        call abc(l, stru%RMT(jatom), pei(l), pi12lo(l,jlo), &
                 pe12lo(l,jlo), pilolo(l,jlo,jlo), jlo, lapw(l))
     end do iloloop
  end do loloop

  ! ... rest of subroutine unchanged ...
end subroutine atpar
```

**Testing**: 
- Compile and verify arrays are allocated correctly
- Check that pi12lo(l,jlo) values are distinct for different jlo
- Add debug print statements if needed

---

### Step 4: Update `l2Amn` Subroutine
**File**: `modules_rc.F`, lines ~500-800 (l2Amn_m module)

**Current code** (WRONG):
```fortran
! Only uses irf=3 (first LO):
projection(ib, ipr, kkk) = projection(ib, ipr, kkk) &
     & + conjg(Alm(idx, ib, mu, 1))                 &
     &   * inwf%projections(ipr)%coeff(iY)          &
     & + conjg(Alm(idx, ib, mu, 3))                 &  ! Only first LO!
     &   * inwf%projections(ipr)%coeff(iY)          &
     &   * pi12lo(l)  ! Also wrong: should be (l,jlo)
```

**Required changes**:
1. Import module arrays
2. Loop over ALL LOs (not just irf=3)
3. Use correct overlap factors

**New code**:
```fortran
subroutine l2Amn(stru, inwf, num_kpts)
  use param,     only: Lmax2, NRF, lomax, unit_vsp, unit_out, unit_vector
  use w2w,       only: Nmat, iblock, unit_amn
  use const,     only: R8, C16, TAU
  use xa,        only: BK, FJ, DFJ, PHS, BKrot, BKRloc
  use xa3,       only: A, XK,YK,ZK, GX,GY,GZ, BKX,BKY,BKZ, vecsz
  use lolog,     only: nlo,nlov,nlon,loor,ilo,lapw, n_rad
  use loabc,     only: alo, pi12lo  ! CHANGED: import from module
  ! ... other uses ...

  ! ... existing variables ...
  
  ! ADDED: Variables for multiple LO handling
  integer  :: irf_lo, jlo_actual
  
  ! ... existing code up to projection calculation ...

  equiv2: do mu=1,stru%mult(jatom)
     latom = lfirst - 1 + mu

     do ib = 1,Nb
        do ipr = 1, inwf%Nproj
           do iY = 1, inwf%projections(ipr)%NY
              if (inwf%projections(ipr)%iat(iY) /= latom) cycle
              l = inwf%projections(ipr)%l(iY)
              m = inwf%projections(ipr)%m(iY)
              idx = l**2 + m + l + 1

              ! LAPW contribution (irf=1)
              projection(ib, ipr, kkk) = projection(ib, ipr, kkk) &
                   & + conjg(Alm(idx, ib, mu, 1))                 &
                   &   * inwf%projections(ipr)%coeff(iY)
              
              ! CHANGED: Loop over ALL LOs for this l
              do jlo = 1, ilo(l)
                 if (lapw(l)) then
                    irf_lo = 2 + jlo  ! For LAPW: u, u̇, LO1, LO2, ...
                 else
                    irf_lo = 1 + jlo  ! For APW: u, LO1, LO2, ...
                 endif
                 
                 if (irf_lo <= n_rad(l)) then
                    projection(ib, ipr, kkk) = projection(ib, ipr, kkk) &
                         & + conjg(Alm(idx, ib, mu, irf_lo))            &
                         &   * inwf%projections(ipr)%coeff(iY)          &
                         &   * pi12lo(l, jlo)  ! Use correct overlap for this jlo
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo equiv2
  
  ! ... rest unchanged ...
end subroutine l2Amn
```

**Testing**: 
- Check that all LOs contribute to projections
- Verify pi12lo indices are correct
- Compare output with single LO case (should match)

---

### Step 5: Add Cleanup to Main Program
**File**: `main.f`, near end of program

**Add before final cleanup**:
```fortran
program wf
  ! ... existing use statements ...
  use loabc, only: cleanup_loabc  ! ADD
  
  ! ... all existing code ...
  
  if (inwf%Amn) then
     call ptime(unit_out)
     call l2amn(stru, inwf, Nk)
     call ptime('l2Amn')
  endif

  ! ADD: Clean up allocated arrays
  call cleanup_loabc()

  call ERRCLR(ERRFN)
  print "('W2W END')"
end program wf
```

---

## Testing Plan

### Test Case 1: Regression Test (Single LO)
**System**: Fe with single 3d LO
**Location**: `tests/Fe_1LO/`
**Expected**: Should work exactly as before (no regression)

**Steps**:
1. Run patched w2w: `./w2wc Fe_1LO.def`
2. Compare LO coefficients with reference
3. Check no NaN/Inf in output
4. Verify *.amn and *.mmn files are identical to reference

**Success criteria**:
- No errors or warnings
- Output matches reference within numerical precision (1e-10)

---

### Test Case 2: Critical Test (Multiple LOs)
**System**: SrVO₃ with 2 LOs for d-channel
**Location**: `tests/SrVO3_2LO/`
**Expected**: Should now work (currently fails)

**Setup**:
```
# case.inwf excerpt:
amn
1 100
5 20     # Multiple projections with 2 d-channel LOs
```

**Steps**:
1. Run patched w2w: `./w2wc SrVO3_2LO.def`
2. Check LO coefficients are printed for BOTH LOs
3. Verify no NaN/Inf anywhere
4. Compare with WIEN2k lapw1 output (A,B,C coefficients should match)

**Success criteria**:
- Calculation completes without crash
- LO coefficients are finite and reasonable (|A|,|B|,|C| < 10)
- Both LOs contribute to projections in *.amn

**Debug checks**:
```fortran
! Add to atpar after overlap calculation:
write(unit_out, '("DEBUG: l=",I2," jlo=",I2," pi12lo=",E15.7)') &
     l, jlo, pi12lo(l,jlo)
```

---

### Test Case 3: Comparison with WIEN2k
**System**: Any test case from above
**Method**: Direct comparison of LO coefficients

**WIEN2k output** (from lapw1.scf):
```
LO COEFFICIENT: l,A,B,C    2     0.12345     0.67890    -0.23456
```

**w2w output** (from w2w.outputwf):
```
LO COEFFICIENT: l,A,B,C    2     0.12345     0.67890    -0.23456
```

**Success criteria**:
- A, B, C match to 5 decimal places
- Both jlo=1 and jlo=2 match (if multiple LOs)

---

## Validation Checklist

Before considering this task complete:

- [ ] Code compiles without warnings (both real and complex)
- [ ] `loabc` module has allocatable arrays
- [ ] `init_loabc()` and `cleanup_loabc()` are called
- [ ] `abc()` uses `pilolo_val` parameter
- [ ] `atpar()` calculates overlaps for each jlo separately
- [ ] `atpar()` stores overlaps in module arrays (l, jlo)
- [ ] `l2Amn()` loops over all LOs, not just irf=3
- [ ] Test case 1 (single LO) passes - no regression
- [ ] Test case 2 (multiple LOs) passes - bug is fixed
- [ ] LO coefficients match WIEN2k output
- [ ] No NaN/Inf in any output
- [ ] Documentation updated (this file, code comments)

---

## Known Issues and Future Work

### Potential Issues
1. **Cross-overlaps**: Currently only diagonal pilolo(l,jlo,jlo) is calculated
   - For strong LO-LO coupling, might need off-diagonal pilolo(l,jlo1,jlo2)
   - WIEN2k has full 4D array for this
   - Monitor: if orthogonalization is poor, add cross-overlaps

2. **APW+lo indexing**: Careful with irf indexing differences
   - LAPW: irf=1(u), 2(u̇), 3(LO1), 4(LO2), ...
   - APW: irf=1(u), 2(LO1), 3(LO2), ...
   - Current fix handles this but needs thorough testing

3. **Memory**: New arrays add ~(lomax × nloat)² × 8 bytes per atom
   - For typical lomax=3, nloat=3: ~576 bytes per atom
   - Negligible for reasonable system sizes

### Future Enhancements (separate tasks)
- Implement full pilolo(l,jlo1,jlo2) cross-overlaps if needed
- Add automatic detection of multiple LO configurations
- Improve error messages if overlap integrals look suspicious
- Add unit tests for overlap calculations

---

## References

### Code Files
- WIEN2k reference: `docs/reference/wien2k-abc.f`
- w2w before fix: Use `git diff` to see changes
- Mathematical background: See LAPW method textbooks

### Papers
- Original w2w paper: Comput. Phys. Commun. 181, 1888 (2010)
- LAPW method: Singh & Nordström, "Planewaves, Pseudopotentials and the LAPW Method"
- LO formalism: Phys. Rev. B 64, 195134 (2001)

### Online Resources
- WIEN2k user guide: Section on Local Orbitals
- Wannier90 user guide: Section on projections

---

## Change Log

| Date | Author | Change |
|------|--------|--------|
| 2024-12-13 | Initial | Created task specification |

---

## Questions/Issues
If problems arise during implementation, document here:

1. **Q**: Why does WIEN2k use 4D pilolo(J,jlo1,jlo2,JATOM)?  
   **A**: For full orthogonalization between different LOs. We start with diagonal only.

2. **Q**: Should we handle spin-polarized case differently?  
   **A**: No - overlaps are spin-independent (radial functions only).

---

**Next Steps**: 
1. Create `docs/reference/wien2k-abc.f` with WIEN2k reference code
2. Begin Step 1: Update loabc module
3. Test compilation after each step
4. Add debug output to verify overlap values