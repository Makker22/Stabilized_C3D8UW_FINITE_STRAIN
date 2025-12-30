# UEL_C3D8UW

Abaqus UEL implementing a 3D 8-node total-Lagrangian U–W saturated porous media element with Hencky strain and Biot/F-bar stabilization. Hencky strain is fed through a UMAT, and the consistent tangent includes the numerical \((dLlog/dE)^T Q\) term for better convergence. Trial/committed SVARS tracking and UVARM post-processing are supported.

## Files
- `UEL_C3D8UW.for`: UEL, UMAT placeholder, SDVINI, and helper routines.
- `ABA_PARAM.INC`: expected from Abaqus include search path.

## Key inputs (PROPS)
Reserved for UEL (defaults in code; override via input file):
1. `rho_s0` – solid density  
2. `rho_f0` – fluid density  
3. `n0`     – initial porosity  
4. `alphaB` – Biot coefficient  
5. `Qb`     – drained bulk modulus  
6. `kperm`  – permeability  
7. `gamma_w`– fluid unit weight  
8. `NPROPS` – passed verbatim to UMAT

## State variables sizing
- Per Gauss point base: `13 + NSTATV_UMAT`.
- 8 Gauss points per element.
- If trial storage is desired (recommended), allocate double: `NSVARS = 8 * 2 * (13 + NSTATV_UMAT)`. Without trial: `NSVARS = 8 * (13 + NSTATV_UMAT)`.

## Initialization (SDVINI)
- Mode 1: read node/element files (`nodes_soil.inp`, `elements_uel.inp`) and set initial stresses from depth.
- Mode 2: map `UVARM1-6` from an Abaqus `.dat` at a target time to all Gauss points.

## Build/use
1) Compile the FORTRAN with the Abaqus-provided compiler environment (e.g., `abaqus make library=UEL_C3D8UW.for`).  
2) Reference in the input file with `*USER ELEMENT` and matching `NSTATEV` from the sizing rule above.  
3) Supply material props in `*UEL PROPERTY`, UMAT props in the tail of PROPS.  
4) To output post variables, request `*EL PRINT, ELSET=..., FREQUENCY=...` with `UVARM`.
