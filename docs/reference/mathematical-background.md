# Mathematical Background: Local Orbitals in LAPW

## Overview

This document explains the mathematical formalism for Local Orbitals (LOs) in the LAPW+lo method, which is essential to understanding the bug fix in wien2wannier.

## LAPW Basis Functions

### Inside Muffin-Tin Spheres

The LAPW basis function for wave vector **k** and reciprocal lattice vector **G** is:

$$
\phi_{\mathbf{k+G}}(\mathbf{r}) = 
\begin{cases}
\sum_{lm} [A_{lm} u_l(r) + B_{lm} \dot{u}_l(r)] Y_{lm}(\hat{r}) & r < R_{MT} \\
\frac{1}{\sqrt{\Omega}} e^{i(\mathbf{k+G})\cdot\mathbf{r}} & r > R_{MT}
\end{cases}
$$

Where:
- $u_l(r)$ = radial function at energy $E_l$
- $\dot{u}_l(r)$ = energy derivative of radial function
- $Y_{lm}(\hat{r})$ = spherical harmonic
- $A_{lm}, B_{lm}$ = coefficients ensuring continuity and differentiability at $R_{MT}$

### The Energy Derivative Problem

The energy derivative $\dot{u}_l(r)$ provides variational flexibility, but:
- Only valid near the linearization energy $E_l$
- Poor description of states far from $E_l$
- Causes problems for semi-core states

## Local Orbitals (LOs)

### Definition

A Local Orbital adds a third radial function at a different energy:

$$
\phi^{LO}_{lm}(\mathbf{r}) = 
\begin{cases}
\sum_{lm} [A^{LO}_{lm} u_l(r) + B^{LO}_{lm} \dot{u}_l(r) + C^{LO}_{lm} u^{LO}_l(r)] Y_{lm}(\hat{r}) & r < R_{MT} \\
0 & r > R_{MT}
\end{cases}
$$

Key properties:
1. **Localized**: Zero outside muffin-tin (no plane wave tail)
2. **Orthogonal**: Must be orthogonal to LAPW basis
3. **Normalized**: Self-overlap integral equals 1

### Why Multiple LOs?

For materials with:
- Semi-core states (e.g., 3p in transition metals)
- Multiple valence bands in same l-channel
- High accuracy requirements

We need **multiple LOs per l**:
- LO₁ at energy $E_{lo,1}$ (e.g., for semi-core 3p)
- LO₂ at energy $E_{lo,2}$ (e.g., for valence 3p)

## Normalization and Orthogonalization

### Single LO Case

For one LO with radial function $u^{LO}$ at energy $E_{lo}$:

**Orthogonality conditions**:
$$\langle \phi^{LO} | \phi_{LAPW} \rangle = 0$$

This gives two conditions (match to u and u̇):
$$A^{LO} \langle u^{LO} | u \rangle + B^{LO} \langle u^{LO} | \dot{u} \rangle + C^{LO} \langle u^{LO} | u^{LO} \rangle = 0 \quad ...(1)$$

**Normalization**:
$$\langle \phi^{LO} | \phi^{LO} \rangle = 1$$

This gives:
$$
(A^{LO})^2 \langle u|u \rangle + (B^{LO})^2 \langle \dot{u}|\dot{u} \rangle + (C^{LO})^2 \langle u^{LO}|u^{LO} \rangle 
$$
$$
+ 2 A^{LO} B^{LO} \langle u|\dot{u} \rangle + 2 A^{LO} C^{LO} \langle u|u^{LO} \rangle + 2 B^{LO} C^{LO} \langle \dot{u}|u^{LO} \rangle = 1 \quad ...(2)
$$

### Overlap Integrals

Define:
$$
\begin{align}
\pi_{12} &= \langle u | \dot{u} \rangle \quad &\text{(LAPW u-ú overlap)} \\
\pi_{e} &= \langle \dot{u} | \dot{u} \rangle \quad &\text{(ú self-overlap)} \\
\pi_{1,lo} &= \langle u | u^{LO} \rangle \quad &\text{(u-LO overlap) } \leftarrow \text{ PI12LO in code} \\
\pi_{e,lo} &= \langle \dot{u} | u^{LO} \rangle \quad &\text{(ú-LO overlap) } \leftarrow \text{ PE12LO in code} \\
\pi_{lo} &= \langle u^{LO} | u^{LO} \rangle \quad &\text{(LO self-overlap) } \leftarrow \text{ PILOLO in code}
\end{align}
$$

### LAPW+lo Coefficients

For LAPW+lo (with energy derivative):

From matching conditions at sphere boundary:
$$
\begin{align}
A &= (P^{LO} \frac{\partial P}{\partial r} - \frac{\partial P^{LO}}{\partial r} P) R^2_{MT} \\
B &= -(P^{LO} \frac{\partial P_e}{\partial r} - \frac{\partial P^{LO}}{\partial r} P_e) R^2_{MT}
\end{align}
$$

From normalization:
$$
C^{LO} = \frac{1}{\sqrt{A^2 + A^2 \cdot 2\pi_{1,lo} + B^2 \cdot \pi_e + 2B^2 \cdot \pi_{e,lo} + \pi_{lo}}}
$$

$$
\begin{align}
A^{LO} &= C^{LO} \cdot A \\
B^{LO} &= C^{LO} \cdot B
\end{align}
$$

Where A, B are the temporary values from matching.

### APW+lo Coefficients (No Energy Derivative)

For APW+lo (no $\dot{u}$):

**First LO (jlo=1)**:
$$
\begin{align}
A^{LO} &= \frac{1}{\sqrt{1 + (P/P_e)^2 \pi_e}} \\
B^{LO} &= -\frac{P}{P_e} \cdot A^{LO} \\
C^{LO} &= 0
\end{align}
$$

**Second LO (jlo=2)**:  
Must orthogonalize to BOTH u and LO₁:
$$
x_{bc} = -\frac{P}{P^{LO_2}}
$$

$$
C^{LO_2} = \frac{x_{bc}}{x_{ac}}, \quad A^{LO_2} = \frac{1}{x_{ac}}
$$

where:
$$x_{ac} = \sqrt{1 + \pi_{lo,22} \cdot x_{bc}^2 + 2 \cdot x_{bc} \cdot \pi_{1,lo,2}}$$

Here:
- $\pi_{lo,22} = \langle u^{LO_2} | u^{LO_2} \rangle$ - second LO self-overlap
- $\pi_{1,lo,2} = \langle u | u^{LO_2} \rangle$ - u-LO₂ overlap

## Multiple LO Case

### General Formulation

For n LOs at different energies $E_{lo,1}, E_{lo,2}, ..., E_{lo,n}$:

**Overlap matrix**:
$$\Pi_{ij} = \langle u^{LO_i} | u^{LO_j} \rangle \quad \text{for } i,j = 1...n$$

This is the **PILOLO** array in WIEN2k:
```fortran
PILOLO(l, i, j, atom) = <u^LO_i | u^LO_j>
```

**Orthogonalization**:  
Each LO must be orthogonal to:
1. The LAPW basis (u, u̇)
2. All previous LOs

### Why w2w Failed

**The bug**: w2w treated overlap integrals as scalars:
```fortran
real(R8) :: pi12lo, pe12lo   ! Scalars - WRONG!
```

When calculating LO₂:
1. Computed $\langle u | u^{LO_2} \rangle$ → overwrote pi12lo (lost LO₁ value!)
2. Called abc() with wrong pi12lo
3. Formula: $x_{ac} = \sqrt{1 + x_{bc}^2 + 2 \cdot x_{bc} \cdot \pi_{1,lo}}$ used wrong value
4. Result: negative under sqrt → NaN

**The fix**: Use arrays indexed by (l, jlo):
```fortran
real(R8), allocatable :: pi12lo(:,:)    ! (0:lomax, nloat)
real(R8), allocatable :: pe12lo(:,:)    ! (0:lomax, nloat)
real(R8), allocatable :: pilolo(:,:,:)  ! (0:lomax, nloat, nloat)
```

Now each LO has its own overlap values!

## Computational Procedure

### Step 1: Solve Radial Equation

For each LO energy $E_{lo,j}$:
$$\left[-\frac{1}{r^2}\frac{d}{dr}r^2\frac{d}{dr} + \frac{l(l+1)}{r^2} + V(r) - E_{lo,j}\right] u^{LO_j}(r) = 0$$

Normalize: $\int_0^{R_{MT}} [u^{LO_j}(r)]^2 r^2 dr = 1$

### Step 2: Compute Overlap Integrals

```fortran
call RINT13(rf1(:,l,1), rf1(:,l,irf_j), pi12lo(l,j))  ! <u | u^LO_j>
call RINT13(rf1(:,l,2), rf1(:,l,irf_j), pe12lo(l,j))  ! <u̇ | u^LO_j>
call RINT13(rf1(:,l,irf_j), rf1(:,l,irf_j), pilolo(l,j,j))  ! <u^LO_j | u^LO_j>
```

For cross-overlaps (optional but recommended):
```fortran
do i = 1, j-1
   call RINT13(rf1(:,l,irf_i), rf1(:,l,irf_j), pilolo(l,i,j))  ! <u^LO_i | u^LO_j>
   pilolo(l,j,i) = pilolo(l,i,j)  ! Symmetric
enddo
```

### Step 3: Calculate Coefficients

Call abc() with **this LO's overlap values**:
```fortran
call abc(l, RMT, pei(l), pi12lo(l,jlo), pe12lo(l,jlo), &
         pilolo(l,jlo,jlo), jlo, lapw(l))
```

### Step 4: Build Projection Matrix

In l2Amn, include **all LOs**:
```fortran
do jlo = 1, ilo(l)
   irf = get_radial_function_index(l, jlo, lapw)
   projection = projection + conj(Alm(irf)) * coeff * pi12lo(l, jlo)
enddo
```

## Example: Iron with 2 LOs for l=2

### Configuration
- Atom: Fe (Z=26)
- l = 2 (d-channel)
- LO₁: $E_{lo,1}$ = -3.0 Ry (semi-core 3d)
- LO₂: $E_{lo,2}$ = 0.3 Ry (valence 3d)

### Radial Functions
- $u_2(r)$ at $E_2$ = 0.2 Ry (LAPW linearization energy)
- $\dot{u}_2(r)$ energy derivative
- $u^{LO_1}(r)$ at $E_{lo,1}$ = -3.0 Ry
- $u^{LO_2}(r)$ at $E_{lo,2}$ = 0.3 Ry

### Overlap Values (typical)

| Overlap | Symbol | Value |
|---------|--------|-------|
| $\langle u_2 \| u_2 \rangle$ | - | 1.0000 (normalized) |
| $\langle \dot{u}_2 \| \dot{u}_2 \rangle$ | $\pi_e$ | 0.8234 |
| $\langle u_2 \| \dot{u}_2 \rangle$ | $\pi_{12}$ | 0.0123 |
| $\langle u_2 \| u^{LO_1} \rangle$ | $\pi_{1,lo,1}$ | -0.0045 |
| $\langle \dot{u}_2 \| u^{LO_1} \rangle$ | $\pi_{e,lo,1}$ | 0.0234 |
| $\langle u^{LO_1} \| u^{LO_1} \rangle$ | $\pi_{lo,11}$ | 1.0000 |
| $\langle u_2 \| u^{LO_2} \rangle$ | $\pi_{1,lo,2}$ | 0.0567 |
| $\langle \dot{u}_2 \| u^{LO_2} \rangle$ | $\pi_{e,lo,2}$ | -0.0123 |
| $\langle u^{LO_2} \| u^{LO_2} \rangle$ | $\pi_{lo,22}$ | 1.0000 |
| $\langle u^{LO_1} \| u^{LO_2} \rangle$ | $\pi_{lo,12}$ | 0.0012 (small!) |

### Coefficient Calculation

**LO₁**:
$$
\begin{align}
A_1 &= C_1 \cdot (P^{LO_1} \frac{\partial P}{\partial r} - \frac{\partial P^{LO_1}}{\partial r} P) R^2_{MT} \\
B_1 &= C_1 \cdot (-(P^{LO_1} \frac{\partial P_e}{\partial r} - \frac{\partial P^{LO_1}}{\partial r} P_e) R^2_{MT}) \\
C_1 &= \frac{1}{\sqrt{A_1^2 + 2A_1^2 \cdot \pi_{1,lo,1} + B_1^2 \cdot \pi_e + 2B_1^2 \cdot \pi_{e,lo,1} + \pi_{lo,11}}}
\end{align}
$$

**LO₂** (using stored $\pi_{1,lo,2}$, not $\pi_{1,lo,1}$):  
Similar calculation but with jlo=2 overlap values

### What Happens Without Fix

**Bug scenario**:
1. Calculate LO₁: `pi12lo = -0.0045` ✓
2. Call abc(jlo=1) with pi12lo = -0.0045 ✓
3. Calculate LO₂: `pi12lo = 0.0567` → **overwrites** old value!
4. Call abc(jlo=2) with pi12lo = 0.0567...

**But wait!** If there's any code that still uses the LO₁ value, it's now wrong!

Worse, if abc(jlo=2) expects pi12lo for jlo=1 but gets jlo=2's value:
$$x_{ac} = \sqrt{1 + x_{bc}^2 + 2 \cdot x_{bc} \cdot 0.0567}$$

Should use jlo=2 value. If the formula expects the wrong index, this can easily go negative → NaN.

**With fix**:
```fortran
pi12lo(2,1) = -0.0045  ! Stored separately
pi12lo(2,2) = 0.0567   ! Doesn't overwrite

call abc(..., pi12lo(2,1), ..., jlo=1)  ! Uses correct value
call abc(..., pi12lo(2,2), ..., jlo=2)  ! Uses correct value
```

## Physical Interpretation

### Why Do We Need This?

**Energy window problem**:
- Want Wannier functions for valence d-states (E ~ 0 to 5 eV)
- But semi-core 3d states exist at E ~ -40 eV
- Single $u_d(r)$ at $E_d$ = 2 eV describes valence well
- But describes semi-core poorly → leakage into Wannier functions

**Solution**: Two LOs
- LO₁ at E ~ -40 eV captures semi-core precisely
- LO₂ at E ~ 3 eV provides extra flexibility for valence
- Together they span the full energy range accurately

### Orthogonality is Critical

If LOs aren't properly orthogonal:
- Semi-core and valence states mix
- Wannier functions become non-localized
- Unphysical results (negative occupations, etc.)

The overlap integrals ensure:
$$\langle LO_{semicore} | LO_{valence} \rangle \approx 0$$

This requires **separate** overlap calculations for each LO!

## Summary

1. **Multiple LOs** = multiple radial functions per l at different energies
2. **Each LO** needs its own set of overlap integrals
3. **Overlap arrays** must be indexed by (l, jlo)
4. **Normalization** requires LO self-overlap $\pi_{lo}$ (not 1.0!)
5. **WIEN2k** implements this correctly with multi-dimensional arrays
6. **w2w bug** was using scalars → overwrites → wrong values → NaN

The fix is simple in principle: use arrays instead of scalars.  
The implementation requires care: update every place these overlaps are used.

## References

1. Singh, D. J., & Nordström, L. (2006). *Planewaves, Pseudopotentials, and the LAPW Method* (2nd ed.). Springer.

2. Sjöstedt, E., Nordström, L., & Singh, D. J. (2000). An alternative way of linearizing the augmented plane-wave method. *Solid State Communications*, 114(1), 15-20.

3. Madsen, G. K., Blaha, P., Schwarz, K., Sjöstedt, E., & Nordström, L. (2001). Efficient linearization of the augmented plane-wave method. *Physical Review B*, 64(19), 195134.

4. Posternak, M., Baldereschi, A., Massidda, S., & Marzari, N. (2002). Maximally localized Wannier functions in antiferromagnetic MnO within the FLAPW formalism. *Physical Review B*, 65(18), 184422.
