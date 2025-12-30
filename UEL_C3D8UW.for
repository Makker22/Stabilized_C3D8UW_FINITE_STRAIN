C======================================================================
C  Abaqus UEL for a 3D 8-node total-Lagrangian U-W saturated porous 
C  media element (C3D8-like) coupling S-F displacement.
C  UEL properties:
C   - PROPS(1:7) reserved for UEL parameters (solid/fluid densities,
C     porosity, Biot coefficient, bulk modulus, permeability,
C     fluid unit weight).
C   - PROPS(8:) forwarded verbatim to the linked UMAT.
C   - NSVARS = 8 * 2 * (13 + NSTATV_UMAT)
C Copyright (c) Ke Ma (makefs@163.com)
C 2025/12/30
C======================================================================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,
     &     NDOFEL,NRHS,NSVARS,
     &     PROPS,NPROPS,COORDS,MCRD,NNODE,
     &     U,DU,V,A,JTYPE,TIME,DTIME,KSTEP,KINC,
     &     JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,
     &     PNEWDT,JPROPS,NJPROP,PERIOD)

      INCLUDE 'ABA_PARAM.INC'

      INTEGER MAXE, MAXG, MAXV, NUMV, EOFF
      PARAMETER (MAXE=100000, MAXG=8, MAXV=20, NUMV=13, EOFF=1000000)
      DOUBLE PRECISION trackedPostVariables(MAXE,MAXG,MAXV)
      COMMON /POSTDATA/ trackedPostVariables
      LOGICAL POSTINIT
      COMMON /POSTFLAG/ POSTINIT
      SAVE /POSTDATA/, /POSTFLAG/

C     ---- Track KINC to promote trial -> committed states ----
      INTEGER LASTKINC(MAXE)
      LOGICAL KINCINIT
      COMMON /KINCFLAG/ LASTKINC, KINCINIT
      SAVE /KINCFLAG/

C     ---- Integer arguments ----
      INTEGER NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE
      INTEGER JTYPE, KSTEP, KINC, JELEM
      INTEGER NDLOAD, NPREDF, MLVARX, MDLOAD, NJPROP
      INTEGER LFLAGS(*), JDLTYP(MDLOAD,*), JPROPS(NJPROP)

C     ---- Scalar doubles ----
      DOUBLE PRECISION DTIME, TIME(2), PERIOD, PNEWDT

C     ---- Main argument arrays ----
      DOUBLE PRECISION RHS(MLVARX,NRHS), AMATRX(NDOFEL,NDOFEL)
      DOUBLE PRECISION SVARS(NSVARS), ENERGY(8)
      DOUBLE PRECISION PROPS(NPROPS), COORDS(MCRD,NNODE)
      DOUBLE PRECISION U(NDOFEL), DU(MLVARX,NRHS)
      DOUBLE PRECISION V(NDOFEL), A(NDOFEL)
      DOUBLE PRECISION PARAMS(*)
      DOUBLE PRECISION ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*)
      DOUBLE PRECISION PREDEF(2,NPREDF,NNODE)

C     ---- Element constants ----
      INTEGER NDIM, NN, NU, NW, NSTR
      PARAMETER (NDIM=3, NN=8)
      PARAMETER (NU=NDIM*NN, NW=NDIM*NN)
      PARAMETER (NSTR=6)

C     ---- Local matrices/vectors (block form) ----
      DOUBLE PRECISION KUU(NU,NU), KUW(NU,NW), KWU(NW,NU), KWW(NW,NW)
      DOUBLE PRECISION MUU(NU,NU), MUW(NU,NW), MWU(NW,NU), MWW(NW,NW)
      DOUBLE PRECISION CWW(NW,NW), CUW(NU,NW), CWU(NW,NU)
      DOUBLE PRECISION FintU(NU), FintW(NW)

C     ---- Assembled element matrices/vectors ----
      DOUBLE PRECISION Kmat(NDOFEL,NDOFEL)
      DOUBLE PRECISION Mmat(NDOFEL,NDOFEL)
      DOUBLE PRECISION Cmat(NDOFEL,NDOFEL)
      DOUBLE PRECISION Fint(NDOFEL)

C     ---- Gauss points for 2x2x2 ----
      DOUBLE PRECISION GP(2), GW(2)
      DATA GP /-0.577350269189626D0, 0.577350269189626D0/
      DATA GW / 1.D0,                 1.D0                /

C     ---- Shape function buffers ----
      DOUBLE PRECISION Nloc(NN), dNdr(NN,3), dNdX(NN,3)
      DOUBLE PRECISION J0(3,3), J0inv(3,3), detJ0, weight

C     ---- Kinematics ----
      DOUBLE PRECISION F(3,3), Finv(3,3), FinvT(3,3), detF
      DOUBLE PRECISION Cten(3,3), Etens(3,3), LogC(3,3)
      DOUBLE PRECISION Evec(NSTR), Svec(NSTR)
      DOUBLE PRECISION Hvec(NSTR), logCvec(NSTR), Qvec(NSTR)
      DOUBLE PRECISION Llog(NSTR,NSTR), Dhencky(NSTR,NSTR)
      DOUBLE PRECISION TMP6(NSTR,NSTR)
      DOUBLE PRECISION S2(3,3)

C     ---- Llog finite-difference buffers for tangent correction ----
      DOUBLE PRECISION CtenP(3,3), LogCP(3,3), LlogP(NSTR,NSTR)
      DOUBLE PRECISION DLterm(NSTR,NSTR), DeffFull(NSTR,NSTR)
      DOUBLE PRECISION epsE

C     ---- UMAT coupling buffers ----
      INTEGER MAXSTATV
      PARAMETER (MAXSTATV=500)
      INTEGER MAXPROP
      PARAMETER (MAXPROP=200)
      CHARACTER*80 CMNAME_UM
      INTEGER NDI_UM, NSHR_UM, NTENS_UM
      INTEGER NPROPS_MAT, svStrideTot, svStrideC, svStrideT
      INTEGER NSTATV_MAT, NSTATV_CALL
      CHARACTER*80 CMNAME_M
      INTEGER NDI_M, NSHR_M, NTENS_M
      DOUBLE PRECISION PNEWDT_M
      DOUBLE PRECISION PROPS_M(MAXPROP)
      DOUBLE PRECISION COORDS_M(3), TIME_M(2)
      DOUBLE PRECISION STRESS_M(NSTR), STRAN_M(NSTR), DSTRAN_M(NSTR)
      DOUBLE PRECISION DDSDDE_M(NSTR,NSTR)
      DOUBLE PRECISION STATEV_M(MAXSTATV)
      DOUBLE PRECISION DDSDDT_M(NSTR), DRPLDE_M(NSTR)
      DOUBLE PRECISION DROT_M(3,3), DFGRD0_M(3,3), DFGRD1_M(3,3)
      DOUBLE PRECISION coordsUM(3), PREDEF_M(1), DPRED_M(1)
      DOUBLE PRECISION SSE_M, SPD_M, SCD_M, RPL_M, DRPLDT_M
      DOUBLE PRECISION TEMP_M, DTEMP_M, CELENT_M
      DOUBLE PRECISION S0vec(NSTR), Deff(NSTR,NSTR)
      DOUBLE PRECISION Cinv(3,3), dJdE(NSTR)

C     ---- TL strain-displacement ----
      DOUBLE PRECISION B0(NSTR,NU)
      DOUBLE PRECISION DF(3,3), DC(3,3), DE(3,3)

C     ---- Biot operators ----
      DOUBLE PRECISION QS(NU), QW(NW)
      DOUBLE PRECISION zeta0, logJ, piBiot

C     ---- F-bar accumulators (pure F-bar alpha_s=alpha_w=0) ----
      DOUBLE PRECISION vol_tot, logJ_sum, zeta0_sum
      DOUBLE PRECISION logJ_avg, zeta0_avg
      DOUBLE PRECISION Sbar(NU), Wbar(NW)
      DOUBLE PRECISION KUU_dq(NU,NU)

C     ---- Temp nodal buffers ----
      DOUBLE PRECISION uN(3), gradN(3), gradN2(3), q_s(3), wN(3)

C     ---- Elastic stiffness ----
      DOUBLE PRECISION Es, nus, lam, mu
      DOUBLE PRECISION Dmat(NSTR,NSTR)

C     ---- Material / physical params ----
      DOUBLE PRECISION rho_s0, rho_f0, n0, alphaB, Qb
      DOUBLE PRECISION kperm, Rdarcy, rho_mix0, gamma_w

C     ---- Misc local ----
      DOUBLE PRECISION betaN, gammaN, dt, internal
      DOUBLE PRECISION facUU, facUW, facWW, facC, sGG
      DOUBLE PRECISION postVars(NUMV)
      DOUBLE PRECISION dq, f1
      DOUBLE PRECISION PNEWDT_MIN
      INTEGER I,J,K,ig,jg,kg,node,comp
      INTEGER baseU, baseW, idxU, idxW
      INTEGER rU, cU, rW, cW
      INTEGER node2, comp2, jU, jW, gU, gWidx, gU2, gW2idx
      INTEGER gpIndex, svBase0, svBaseC0, svBaseT0, svBaseC, svBaseT
      LOGICAL hasTrial

C     ---- Init: minimum recommended time step from this routine ----
      PNEWDT_MIN = 1.D0

C     ---- Read UEL properties from PROPS(1:7); PROPS(8:) reserved for UMAT ----
      rho_s0 = PROPS(1)
      rho_f0 = PROPS(2)
      n0     = PROPS(3)
      alphaB = PROPS(4)
      Qb     = PROPS(5)
      kperm  = PROPS(6)
      gamma_w= PROPS(7)

      rho_mix0 = (1.D0 - n0)*rho_s0 + n0*rho_f0
      Rdarcy   = gamma_w / kperm

C     ---- UMAT properties (passed as PROPS(8:)) ----
      NPROPS_MAT = NPROPS - 7
      IF (NPROPS_MAT .LT. 2) NPROPS_MAT = 2

C======================================================================
C     ---- SVARS layout (Point-1 fix: committed + trial) ----
C======================================================================
      svStrideTot = NSVARS / MAXG
      hasTrial = .FALSE.
      svStrideC = svStrideTot
      svStrideT = 0
      IF (svStrideTot .GE. 2*NUMV) THEN
        svStrideC = svStrideTot / 2
        svStrideT = svStrideC
        IF (2*svStrideC .EQ. svStrideTot) THEN
          hasTrial = .TRUE.
        END IF
      END IF

      NSTATV_MAT = svStrideC - NUMV
      IF (NSTATV_MAT .LT. 0) NSTATV_MAT = 0
      NSTATV_CALL = NSTATV_MAT
      IF (NSTATV_CALL .EQ. 0) NSTATV_CALL = 1

C     ---- UMAT call constants ----
      CMNAME_UM = 'UEL_UMAT'
      NDI_UM    = 3
      NSHR_UM   = 3
      NTENS_UM  = 6

C======================================================================
C     Track KINC on entry and promote trial -> committed increments
C======================================================================
      IF (.NOT.KINCINIT) THEN
        DO I=1,MAXE
          LASTKINC(I) = -999999
        END DO
        KINCINIT = .TRUE.
      END IF

      IF (hasTrial) THEN
        IF (JELEM .GE. 1 .AND. JELEM .LE. MAXE) THEN
          IF (KINC .GT. LASTKINC(JELEM)) THEN
C           Promote only after the first visit (LASTKINC already valid)
              DO gpIndex = 1, MAXG
                svBase0  = (gpIndex-1)*svStrideTot
                svBaseC0 = svBase0
                svBaseT0 = svBase0 + svStrideC
                DO I=1, svStrideC
                  SVARS(svBaseC0+I) = SVARS(svBaseT0+I)
                END DO
              END DO
            END IF
            LASTKINC(JELEM) = KINC
          ELSE IF (KINC .LT. LASTKINC(JELEM)) THEN
C           Rewind/restart: update record without promoting
            LASTKINC(JELEM) = KINC
          END IF
        END IF
      END IF

C     ---- Zero local block matrices and forces ----
      DO I=1,NU
        FintU(I) = 0.D0
        DO J=1,NU
          KUU(I,J) = 0.D0
          MUU(I,J) = 0.D0
          KUU_dq(I,J) = 0.D0
        END DO
        DO J=1,NW
          KUW(I,J) = 0.D0
          MUW(I,J) = 0.D0
          CUW(I,J) = 0.D0
        END DO
        Sbar(I) = 0.D0
      END DO

      DO I=1,NW
        FintW(I) = 0.D0
        DO J=1,NW
          KWW(I,J) = 0.D0
          MWW(I,J) = 0.D0
          CWW(I,J) = 0.D0
        END DO
        DO J=1,NU
          KWU(I,J) = 0.D0
          MWU(I,J) = 0.D0
          CWU(I,J) = 0.D0
        END DO
        Wbar(I) = 0.D0
      END DO

      IF (.NOT.POSTINIT) THEN
        DO I=1,MAXE
          DO J=1,MAXG
            DO K=1,MAXV
              trackedPostVariables(I,J,K) = 0.D0
            END DO
          END DO
        END DO
        POSTINIT = .TRUE.
      END IF

C     ---- Initialize F-bar accumulators ----
      vol_tot   = 0.D0
      logJ_sum  = 0.D0
      zeta0_sum = 0.D0

C     ============================================================
C     2x2x2 Gauss integration in reference configuration
C     ============================================================
      DO kg = 1,2
        DO jg = 1,2
          DO ig = 1,2

C           ---- Shape functions N, derivatives wrt natural coords ----
            CALL SHAPE8(GP(ig),GP(jg),GP(kg),Nloc,dNdr)

C           ---- Jacobian dX/dr and its inverse ----
            DO I=1,3
              DO J=1,3
                J0(I,J) = 0.D0
              END DO
            END DO

            DO node=1,NN
              DO I=1,3
                J0(I,1) = J0(I,1) + dNdr(node,1)*COORDS(I,node)
                J0(I,2) = J0(I,2) + dNdr(node,2)*COORDS(I,node)
                J0(I,3) = J0(I,3) + dNdr(node,3)*COORDS(I,node)
              END DO
            END DO

            CALL INV3(J0,J0inv,detJ0)
            weight = detJ0 * GW(ig)*GW(jg)*GW(kg)

C           ---- Gradients wrt reference X ----
            DO node=1,NN
              DO I=1,3
                dNdX(node,I) =
     &              J0inv(1,I)*dNdr(node,1)
     &            + J0inv(2,I)*dNdr(node,2)
     &            + J0inv(3,I)*dNdr(node,3)
              END DO
            END DO

C           =======================================================
C           1) Deformation gradient F = I + sum(U outer gradN)
C           =======================================================
            DO I=1,3
              DO J=1,3
                F(I,J) = 0.D0
                IF (I .EQ. J) F(I,J) = 1.D0
              END DO
            END DO

            DO node=1,NN
              baseU = (node-1)*6
              uN(1) = U(baseU+1)
              uN(2) = U(baseU+2)
              uN(3) = U(baseU+3)
              DO I=1,3
                DO J=1,3
                  F(I,J) = F(I,J) + uN(I)*dNdX(node,J)
                END DO
              END DO
            END DO

C           ---- detF, logJ, F^{-1}, F^{-T} ----
            CALL INV3(F,Finv,detF)
            logJ = LOG(detF)

            DO I=1,3
              DO J=1,3
                FinvT(I,J) = Finv(J,I)
              END DO
            END DO

C           ---- Accumulate averages for F-bar ----
            vol_tot   = vol_tot   + weight
            logJ_sum  = logJ_sum  + logJ*weight

C           ---- Right Cauchy-Green C = F^T F ----
            DO I=1,3
              DO J=1,3
                Cten(I,J) = 0.D0
                DO K=1,3
                  Cten(I,J) = Cten(I,J) + F(K,I)*F(K,J)
                END DO
              END DO
            END DO

C           ===================================================
C           Constitutive update via UMAT with Hencky strain
C           ===================================================

C           ---- log(C) and its derivative wrt C (Voigt) ----
            CALL LOGM_SYM3_AND_DERIV(Cten,LogC,Llog)
            CALL STRAIN_TENSOR_TO_VOIGT(LogC,logCvec)

C           ---- Hencky strain H = 0.5*log(C) ----
            DO I=1,NSTR
              Hvec(I) = 0.5D0*logCvec(I)
            END DO

C           ---- GP index and SVARS bases ----
            gpIndex = (kg-1)*4 + (jg-1)*2 + ig
            svBase0 = (gpIndex-1)*svStrideTot
            svBaseC = svBase0
            svBaseT = svBase0
            IF (hasTrial) svBaseT = svBase0 + svStrideC

C           ---- Initialize UMAT arrays ----
            DO I=1,NSTR
              STRESS_M(I) = 0.D0
              STRAN_M(I)  = 0.D0
              DSTRAN_M(I) = 0.D0
              DO J=1,NSTR
                DDSDDE_M(I,J) = 0.D0
              END DO
            END DO

C           ---- Load committed UMAT stress(Q) and strain(H) from SVARS ----
            IF (NSVARS .GE. svBaseC+NUMV) THEN
              DO I=1,NSTR
                STRESS_M(I) = SVARS(svBaseC+I)
                STRAN_M(I)  = SVARS(svBaseC+6+I)
              END DO
            END IF

C           ---- Strain increment passed to UMAT: DeltaH = H_trial - H_committed ----
            DO I=1,NSTR
              DSTRAN_M(I) = Hvec(I) - STRAN_M(I)
            END DO

C           ---- Load committed UMAT STATEV from SVARS (if present) ----
            IF (NSTATV_MAT .GT. 0) THEN
              IF (NSVARS .GE. svBaseC+NUMV+NSTATV_MAT) THEN
                DO I=1,NSTATV_MAT
                  STATEV_M(I) = SVARS(svBaseC+NUMV+I)
                END DO
              END IF
            END IF

C           ---- UMAT material properties from PROPS(8:) ----
            NPROPS_MAT = NPROPS - 7
            IF (NPROPS_MAT .LT. 1) NPROPS_MAT = 1
            DO I=1,NPROPS_MAT
              PROPS_M(I) = PROPS(7+I)
            END DO

C           ---- Set UMAT call constants ----
            NDI_M   = 3
            NSHR_M  = 3
            NTENS_M = NTENS_UM
            CMNAME_M = CMNAME_UM

C           ---- Dummy values for unused UMAT arguments ----
            DO I=1,NSTR
              DDSDDT_M(I) = 0.D0
              DRPLDE_M(I) = 0.D0
            END DO
            DRPLDT_M = 0.D0
            SSE_M    = 0.D0
            SPD_M    = 0.D0
            SCD_M    = 0.D0
            RPL_M    = 0.D0
            PNEWDT_M = 1.D0
            CELENT_M = 1.D0
            DO I=1,3
              COORDS_M(I) = 0.D0
              DO J=1,3
                DROT_M(I,J)   = 0.D0
                DFGRD0_M(I,J) = 0.D0
                DFGRD1_M(I,J) = 0.D0
              END DO
            END DO
            DO I=1,3
              DROT_M(I,I)   = 1.D0
              DFGRD0_M(I,I) = 1.D0
              DFGRD1_M(I,I) = 1.D0
            END DO
            PREDEF_M(1) = 0.D0
            DPRED_M(1)  = 0.D0
            TEMP_M      = 0.D0
            DTEMP_M     = 0.D0
            TIME_M   = TIME - DTIME

C           ---- Call UMAT: (H, DeltaH) -> (Q, D = dQ/dH) ----
            NSTATV_CALL = NSTATV_MAT
            IF (NSTATV_CALL .LT. 1) NSTATV_CALL = 1
            CALL UMAT(STRESS_M,STATEV_M,DDSDDE_M,SSE_M,SPD_M,SCD_M,RPL_M,
     &        DDSDDT_M,DRPLDE_M,DRPLDT_M,STRAN_M,DSTRAN_M,
     &        TIME_M,DTIME,TEMP_M,DTEMP_M,PREDEF_M,DPRED_M,CMNAME_M,
     &        NDI_M,NSHR_M,NTENS_M,NSTATV_CALL,PROPS_M,NPROPS_MAT,
     &        COORDS_M,DROT_M,PNEWDT_M,CELENT_M,DFGRD0_M,DFGRD1_M,
     &        JELEM,gpIndex,0,0,KSTEP,KINC)

C           ---- Adopt UMAT time-step recommendation (PNEWDT_M) ----
            IF (PNEWDT_M .LT. PNEWDT_MIN) PNEWDT_MIN = PNEWDT_M

C           ---- UMAT outputs: Qvec and Dmat (=DDSDDE_M) ----
            DO I=1,NSTR
              Qvec(I) = STRESS_M(I)
            END DO

C           ---- 2PK stress (reference measure) ----
            DO I=1,NSTR
              S0vec(I) = 0.D0
              DO J=1,NSTR
                S0vec(I) = S0vec(I) + Llog(J,I)*Qvec(J)
              END DO
            END DO
            DO I=1,NSTR
              Svec(I) = S0vec(I)
            END DO

C           ---- Deff = Llog^T * D * Llog  (no J) ----
            DO I=1,NSTR
              DO J=1,NSTR
                TMP6(I,J) = 0.D0
                Deff(I,J) = 0.D0
              END DO
            END DO
            DO I=1,NSTR
              DO J=1,NSTR
                DO K=1,NSTR
                  TMP6(I,J) = TMP6(I,J) + DDSDDE_M(I,K)*Llog(K,J)
                END DO
              END DO
            END DO
            DO I=1,NSTR
              DO J=1,NSTR
                DO K=1,NSTR
                  Deff(I,J) = Deff(I,J) + Llog(K,I)*TMP6(K,J)
                END DO
              END DO
            END DO

C======================================================================
C           (Consistent tangent) add (dLlog/dE)^T * Q term
C             dS0/dE = Llog^T * D * Llog  +  (dLlog/dE)^T * Q
C           Finite-difference dLlog/dE by perturbing each Voigt E component
C           Use C = I + 2E (engineering shear: C12->E4, C23->E5, C13->E6)
            epsE = 1.D-7

            DO I=1,NSTR
              DO J=1,NSTR
                DLterm(I,J) = 0.D0
              END DO
            END DO

            DO K=1,NSTR

C             Copy C -> CtenP
              DO I=1,3
                DO J=1,3
                  CtenP(I,J) = Cten(I,J)
                END DO
              END DO

C             Apply perturbation consistent with Voigt E (engineering shear)
C             normals: Cii += 2*epsE
C             shears : Cij += epsE, Cji += epsE
              IF (K .EQ. 1) THEN
                CtenP(1,1) = CtenP(1,1) + 2.D0*epsE
              ELSE IF (K .EQ. 2) THEN
                CtenP(2,2) = CtenP(2,2) + 2.D0*epsE
              ELSE IF (K .EQ. 3) THEN
                CtenP(3,3) = CtenP(3,3) + 2.D0*epsE
              ELSE IF (K .EQ. 4) THEN
                CtenP(1,2) = CtenP(1,2) + epsE
                CtenP(2,1) = CtenP(2,1) + epsE
              ELSE IF (K .EQ. 5) THEN
                CtenP(2,3) = CtenP(2,3) + epsE
                CtenP(3,2) = CtenP(3,2) + epsE
              ELSE IF (K .EQ. 6) THEN
                CtenP(1,3) = CtenP(1,3) + epsE
                CtenP(3,1) = CtenP(3,1) + epsE
              END IF

              CALL LOGM_SYM3_AND_DERIV(CtenP,LogCP,LlogP)

C             DLterm(:,K) = (dLlog/dE_K)^T * Q
              DO I=1,NSTR
                DLterm(I,K) = 0.D0
                DO J=1,NSTR
                  DLterm(I,K) = DLterm(I,K)
     &                        + ( (LlogP(J,I) - Llog(J,I)) / epsE )
     &                          * Qvec(J)
                END DO
              END DO

            END DO

C           ---- DeffFull = Deff + DLterm ----
            DO I=1,NSTR
              DO J=1,NSTR
                DeffFull(I,J) = Deff(I,J) + DLterm(I,J)
              END DO
            END DO

C           ---- C^{-1} = F^{-1} F^{-T} ----
            DO I=1,3
              DO J=1,3
                Cinv(I,J) = 0.D0
                DO K=1,3
                  Cinv(I,J) = Cinv(I,J) + Finv(K,I)*Finv(K,J)
                END DO
              END DO
            END DO

C           ---- dJ/dE (Voigt) ----
            dJdE(1) = detF * Cinv(1,1)
            dJdE(2) = detF * Cinv(2,2)
            dJdE(3) = detF * Cinv(3,3)
            dJdE(4) = 0.5D0*detF * Cinv(1,2)
            dJdE(5) = 0.5D0*detF * Cinv(2,3)
            dJdE(6) = 0.5D0*detF * Cinv(1,3)

C           ---- Consistent tangent in E-space (reference measure) ----
C                dS/dE = J * dS0/dE + S0 ??dJ/dE
            DO I=1,NSTR
              DO J=1,NSTR
                Dhencky(I,J) = detF*DeffFull(I,J)
              END DO
            END DO

            S2(1,1) = Svec(1)
            S2(2,2) = Svec(2)
            S2(3,3) = Svec(3)
            S2(1,2) = Svec(4)
            S2(2,1) = S2(1,2)
            S2(2,3) = Svec(5)
            S2(3,2) = S2(2,3)
            S2(1,3) = Svec(6)
            S2(3,1) = S2(1,3)

C           ---- Build TL B0 for solid displacements ----
            CALL BUILD_B0_TL(F,dNdX,NN,B0)

C           =======================================================
C           2) Biot volumetric operators QS, QW
C              delta(lnJ) = QS^T delta(dU),  delta(zeta0) = QW^T delta(dW)
C           =======================================================
            DO I=1,NU
              QS(I) = 0.D0
            END DO
            DO I=1,NW
              QW(I) = 0.D0
            END DO

            DO node=1,NN
              gradN(1) = dNdX(node,1)
              gradN(2) = dNdX(node,2)
              gradN(3) = dNdX(node,3)

C             q_s = F^{-T} gradN
              DO I=1,3
                q_s(I) = 0.D0
                DO J=1,3
                  q_s(I) = q_s(I) + FinvT(I,J)*gradN(J)
                END DO
              END DO

              idxU = (node-1)*3
              idxW = (node-1)*3

              DO comp=1,3
                QS(idxU+comp) = q_s(comp)
                QW(idxW+comp) = gradN(comp)
              END DO
            END DO

C           ---- Accumulate element Biot operators (for F-bar tangent and internal force) ----
            DO I=1,NU
              Sbar(I) = Sbar(I) + QS(I)*weight
            END DO
            DO I=1,NW
              Wbar(I) = Wbar(I) + QW(I)*weight
            END DO

C           ---- Flow volumetric strain zeta0 = sum(gradN dot W) ----
            zeta0 = 0.D0
            DO node=1,NN
              baseW = (node-1)*6 + 3
              wN(1) = U(baseW+1)
              wN(2) = U(baseW+2)
              wN(3) = U(baseW+3)
              DO I=1,3
                zeta0 = zeta0 + dNdX(node,I)*wN(I)
              END DO
            END DO

C           ---- Accumulate element-average zeta0 (F-bar) ----
            zeta0_sum = zeta0_sum + zeta0*weight

C           ---- Post variables (pi placeholder replaced later with F-bar value) ----
            postVars(1)  = Svec(1)
            postVars(2)  = Svec(2)
            postVars(3)  = Svec(3)
            postVars(4)  = Svec(4)
            postVars(5)  = Svec(5)
            postVars(6)  = Svec(6)
            postVars(7)  = Hvec(1)
            postVars(8)  = Hvec(2)
            postVars(9)  = Hvec(3)
            postVars(10) = Hvec(4)
            postVars(11) = Hvec(5)
            postVars(12) = Hvec(6)
            postVars(13) = 0.D0

C           ---- Store SVARS: write trial half (or committed if trial storage is absent) ----
            IF (hasTrial) THEN
              IF (NSVARS .GE. svBaseT+NUMV) THEN
                DO I=1,NSTR
                  SVARS(svBaseT+I)     = Qvec(I)
                  SVARS(svBaseT+6+I)   = Hvec(I)
                END DO
                SVARS(svBaseT+13) = postVars(13)
                IF (NSTATV_MAT .GT. 0) THEN
                  IF (NSVARS .GE. svBaseT+NUMV+NSTATV_MAT) THEN
                    DO I=1,NSTATV_MAT
                      SVARS(svBaseT+NUMV+I) = STATEV_M(I)
                    END DO
                  END IF
                END IF
              END IF
            ELSE
              IF (NSVARS .GE. svBaseC+NUMV) THEN
                DO I=1,NSTR
                  SVARS(svBaseC+I)     = Qvec(I)
                  SVARS(svBaseC+6+I)   = Hvec(I)
                END DO
                SVARS(svBaseC+13) = postVars(13)
                IF (NSTATV_MAT .GT. 0) THEN
                  IF (NSVARS .GE. svBaseC+NUMV+NSTATV_MAT) THEN
                    DO I=1,NSTATV_MAT
                      SVARS(svBaseC+NUMV+I) = STATEV_M(I)
                    END DO
                  END IF
                END IF
              END IF
            END IF

            CALL storePostVariables(JELEM+EOFF,gpIndex,
     &           postVars,NUMV)

C           =======================================================
C           3) Internal forces & stiffness
C           =======================================================
            DO I=1,NU
              DO K=1,NSTR
                FintU(I) = FintU(I) + B0(K,I)*Svec(K)*weight
              END DO
            END DO

            CALL ADD_KUU_MATERIAL(B0,Dhencky,NU,NSTR,KUU,weight)

C           Geometric stiffness (initial stress) from 2PK S2
            DO node=1,NN
              DO node2=1,NN
                sGG = 0.D0
                DO I=1,3
                  DO J=1,3
                    sGG = sGG + dNdX(node,I)*S2(I,J)*dNdX(node2,J)
                  END DO
                END DO
                DO comp=1,3
                  rU = (node-1)*3 + comp
                  cU = (node2-1)*3 + comp
                  KUU(rU,cU) = KUU(rU,cU) + sGG*weight
                END DO
              END DO
            END DO

C           =======================================================
C           3.2 Biot volumetric coupling: accumulate dQS/dU
C           =======================================================
            DO node=1,NN
              gradN(1) = dNdX(node,1)
              gradN(2) = dNdX(node,2)
              gradN(3) = dNdX(node,3)
              DO node2=1,NN
                gradN2(1) = dNdX(node2,1)
                gradN2(2) = dNdX(node2,2)
                gradN2(3) = dNdX(node2,3)
                DO comp=1,3
                  DO comp2=1,3
                    dq = 0.D0
                    DO J=1,3
                      f1 = gradN(J)*Finv(J,comp2)
                      DO K=1,3
                        dq = dq - f1*gradN2(K)*Finv(K,comp)
                      END DO
                    END DO
                    rU = (node-1)*3 + comp
                    cU = (node2-1)*3 + comp2
                    KUU_dq(rU,cU) = KUU_dq(rU,cU)
     &                              + dq*weight
                  END DO
                END DO
              END DO
            END DO

C           =======================================================
C           3.3 Mass matrices
C           =======================================================
            DO node=1,NN
              DO J=1,NN
                facUU = rho_mix0 *Nloc(node)*Nloc(J)*weight
                facUW = rho_f0   *Nloc(node)*Nloc(J)*weight
                facWW = rho_f0/n0*Nloc(node)*Nloc(J)*weight

                DO comp=1,3
                  rU = (node-1)*3 + comp
                  cU = (    J-1)*3 + comp
                  rW = rU
                  cW = cU

                  MUU(rU,cU) = MUU(rU,cU) + facUU
                  MWW(rW,cW) = MWW(rW,cW) + facWW
                  MUW(rU,cW) = MUW(rU,cW) + facUW
                  MWU(rW,cU) = MWU(rW,cU) + facUW
                END DO
              END DO
            END DO

C           =======================================================
C           3.4 Darcy damping
C           =======================================================
            DO node=1,NN
              DO J=1,NN
                facC = Rdarcy*detF*Nloc(node)*Nloc(J)*weight

                DO comp=1,3
                  rU = (node-1)*3 + comp
                  cU = (   J-1)*3 + comp
                  rW = rU
                  cW = cU

                  CWW(rW,cW) = CWW(rW,cW) + facC
                END DO
              END DO
            END DO

C           ---- end integration point ----
          END DO
        END DO
      END DO

C     ============================================================
C     F-bar Biot scalar & volumetric coupling
C     ============================================================
      IF (vol_tot .GT. 0.D0) THEN
        logJ_avg  = logJ_sum  / vol_tot
        zeta0_avg = zeta0_sum / vol_tot
      ELSE
        logJ_avg  = 0.D0
        zeta0_avg = 0.D0
      END IF

      piBiot = alphaB*logJ_avg + zeta0_avg

      DO I=1,NU
        FintU(I) = FintU(I) + Qb*piBiot*alphaB * Sbar(I)
      END DO
      DO I=1,NW
        FintW(I) = FintW(I) + Qb*piBiot * Wbar(I)
      END DO

      IF (vol_tot .GT. 0.D0) THEN
        DO I=1,NU
          DO J=1,NU
            KUU(I,J) = KUU(I,J)
     &                 + Qb*alphaB*alphaB
     &                   * Sbar(I)*Sbar(J) / vol_tot
          END DO
        END DO

        DO I=1,NU
          DO J=1,NW
            KUW(I,J) = KUW(I,J)
     &                 + Qb*alphaB
     &                   * Sbar(I)*Wbar(J) / vol_tot
          END DO
        END DO

        DO I=1,NW
          DO J=1,NU
            KWU(I,J) = KWU(I,J)
     &                 + Qb*alphaB
     &                   * Wbar(I)*Sbar(J) / vol_tot
          END DO
        END DO

        DO I=1,NW
          DO J=1,NW
            KWW(I,J) = KWW(I,J)
     &                 + Qb
     &                   * Wbar(I)*Wbar(J) / vol_tot
          END DO
        END DO
      END IF

      DO I=1,NU
        DO J=1,NU
          KUU(I,J) = KUU(I,J)
     &               + Qb*piBiot*alphaB * KUU_dq(I,J)
        END DO
      END DO

C     ---- Update post variables with Qb*piBiot (component 13) ----
      DO gpIndex = 1,MAXG
        svBase0 = (gpIndex-1)*svStrideTot
        svBaseC = svBase0
        svBaseT = svBase0
        IF (hasTrial) svBaseT = svBase0 + svStrideC

        IF (NSVARS .GE. svBase0+NUMV) THEN
          DO I=1,12
            postVars(I) = trackedPostVariables(JELEM,gpIndex,I)
          END DO
          postVars(13) = - Qb*piBiot
          IF (hasTrial) THEN
            SVARS(svBaseT+13) = postVars(13)
          ELSE
            SVARS(svBaseC+13) = postVars(13)
          END IF
          CALL storePostVariables(JELEM+EOFF,gpIndex,
     &         postVars,NUMV)
        END IF
      END DO

C     ============================================================
C     4. Assemble full element matrices and internal force
C        Each node has 6 DOF = [Ux Uy Uz Wx Wy Wz]
C     ============================================================
      DO I=1,NDOFEL
        Fint(I) = 0.D0
        DO J=1,NDOFEL
          Kmat(I,J) = 0.D0
          Mmat(I,J) = 0.D0
          Cmat(I,J) = 0.D0
        END DO
      END DO

      DO node = 1, NN
        DO comp = 1, 3
          idxU = (node-1)*3 + comp
          idxW = (node-1)*3 + comp

          gU     = (node-1)*6 + comp
          gWidx  = (node-1)*6 + 3 + comp

          Fint(gU)    = Fint(gU)    + FintU(idxU)
          Fint(gWidx) = Fint(gWidx) + FintW(idxW)

          DO node2 = 1, NN
            DO comp2 = 1, 3
              jU    = (node2-1)*3 + comp2
              jW    = (node2-1)*3 + comp2

              gU2    = (node2-1)*6 + comp2
              gW2idx = (node2-1)*6 + 3 + comp2

              Kmat(gU,gU2)       = Kmat(gU,gU2)       + KUU(idxU,jU)
              Kmat(gU,gW2idx)    = Kmat(gU,gW2idx)    + KUW(idxU,jW)
              Kmat(gWidx,gU2)    = Kmat(gWidx,gU2)    + KWU(idxW,jU)
              Kmat(gWidx,gW2idx) = Kmat(gWidx,gW2idx) + KWW(idxW,jW)

              Mmat(gU,gU2)       = Mmat(gU,gU2)       + MUU(idxU,jU)
              Mmat(gU,gW2idx)    = Mmat(gU,gW2idx)    + MUW(idxU,jW)
              Mmat(gWidx,gU2)    = Mmat(gWidx,gU2)    + MWU(idxW,jU)
              Mmat(gWidx,gW2idx) = Mmat(gWidx,gW2idx) + MWW(idxW,jW)

              Cmat(gWidx,gW2idx) = Cmat(gWidx,gW2idx) + CWW(idxW,jW)
            END DO
          END DO

        END DO
      END DO

C     ============================================================
C     5. Build effective stiffness and internal force
C     ============================================================
      dt = DTIME

      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMATRX(I,J) = 0.D0
        END DO
        RHS(I,1) = 0.D0
      END DO

      IF (LFLAGS(1).EQ.11 .OR. LFLAGS(1).EQ.12) THEN
        betaN  = PARAMS(2)
        gammaN = PARAMS(3)

        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = Kmat(I,J)
     &                   + gammaN/(betaN*dt)    * Cmat(I,J)
     &                   + 1.D0/(betaN*dt*dt)   * Mmat(I,J)
          END DO
        END DO

        DO I=1,NDOFEL
          internal = Fint(I)
          DO J=1,NDOFEL
            internal = internal
     &               + Cmat(I,J)*V(J)
     &               + Mmat(I,J)*A(J)
          END DO
          RHS(I,1) = -internal
        END DO

      ELSE

        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = Kmat(I,J)
     &                   + 1.D0/dt * Cmat(I,J)
          END DO
        END DO

        DO I=1,NDOFEL
          internal = Fint(I)
          DO J=1,NDOFEL
            internal = internal
     &               + Cmat(I,J)*V(J)
          END DO
          RHS(I,1) = -internal
        END DO

      END IF

      PNEWDT = PNEWDT_MIN

      RETURN
      END




C================================================================
C  8-node hex trilinear shape functions & derivatives
C================================================================
      SUBROUTINE SHAPE8(r,s,t,N,dNdr)
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION r,s,t
      DOUBLE PRECISION N(8), dNdr(8,3)

      DOUBLE PRECISION rp,rm,sp,sm,tp,tm
      rp = 1.D0 + r
      rm = 1.D0 - r
      sp = 1.D0 + s
      sm = 1.D0 - s
      tp = 1.D0 + t
      tm = 1.D0 - t

      N(1) = 0.125D0*rm*sm*tm
      N(2) = 0.125D0*rp*sm*tm
      N(3) = 0.125D0*rp*sp*tm
      N(4) = 0.125D0*rm*sp*tm
      N(5) = 0.125D0*rm*sm*tp
      N(6) = 0.125D0*rp*sm*tp
      N(7) = 0.125D0*rp*sp*tp
      N(8) = 0.125D0*rm*sp*tp

C     dN/dr
      dNdr(1,1) = -0.125D0*sm*tm
      dNdr(2,1) =  0.125D0*sm*tm
      dNdr(3,1) =  0.125D0*sp*tm
      dNdr(4,1) = -0.125D0*sp*tm
      dNdr(5,1) = -0.125D0*sm*tp
      dNdr(6,1) =  0.125D0*sm*tp
      dNdr(7,1) =  0.125D0*sp*tp
      dNdr(8,1) = -0.125D0*sp*tp

C     dN/ds
      dNdr(1,2) = -0.125D0*rm*tm
      dNdr(2,2) = -0.125D0*rp*tm
      dNdr(3,2) =  0.125D0*rp*tm
      dNdr(4,2) =  0.125D0*rm*tm
      dNdr(5,2) = -0.125D0*rm*tp
      dNdr(6,2) = -0.125D0*rp*tp
      dNdr(7,2) =  0.125D0*rp*tp
      dNdr(8,2) =  0.125D0*rm*tp

C     dN/dt
      dNdr(1,3) = -0.125D0*rm*sm
      dNdr(2,3) = -0.125D0*rp*sm
      dNdr(3,3) = -0.125D0*rp*sp
      dNdr(4,3) = -0.125D0*rm*sp
      dNdr(5,3) =  0.125D0*rm*sm
      dNdr(6,3) =  0.125D0*rp*sm
      dNdr(7,3) =  0.125D0*rp*sp
      dNdr(8,3) =  0.125D0*rm*sp

      RETURN
      END

C================================================================
C  3x3 inverse + determinant
C================================================================
      SUBROUTINE INV3(A,Ainv,detA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION A(3,3),Ainv(3,3),detA
      DOUBLE PRECISION invdet

      detA =  A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     &      -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     &      +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

      invdet = 1.D0/detA

      Ainv(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))*invdet
      Ainv(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))*invdet
      Ainv(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))*invdet

      Ainv(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))*invdet
      Ainv(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))*invdet
      Ainv(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))*invdet

      Ainv(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))*invdet
      Ainv(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))*invdet
      Ainv(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))*invdet

      RETURN
      END

C================================================================
C  Strain tensor -> Voigt (E11,E22,E33,2E12,2E23,2E13)
C================================================================
      SUBROUTINE STRAIN_TENSOR_TO_VOIGT(Eten,Evec)
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION Eten(3,3),Evec(6)

      Evec(1) = Eten(1,1)
      Evec(2) = Eten(2,2)
      Evec(3) = Eten(3,3)
      Evec(4) = 2.D0*Eten(1,2)
      Evec(5) = 2.D0*Eten(2,3)
      Evec(6) = 2.D0*Eten(1,3)

      RETURN
      END

C================================================================
C  TL B0: deltaE = B0 * deltaU  (3D, 8-node, Green E)
C================================================================
      SUBROUTINE BUILD_B0_TL(F,dNdX,NN,B0)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NN
      DOUBLE PRECISION F(3,3),dNdX(NN,3),B0(6,3*NN)
      DOUBLE PRECISION DF(3,3),DC(3,3),DE(3,3)
      INTEGER node,comp,I,J,K,col

      DO I=1,6
        DO J=1,3*NN
          B0(I,J) = 0.D0
        END DO
      END DO

      DO node=1,NN
        DO comp=1,3
          col = (node-1)*3 + comp

C         deltaF(i,J) = partial(delta u_i)/partial X_J
C                  = dNdX(node,J) * delta_{i,comp}
          DO I=1,3
            DO J=1,3
              IF (I .EQ. comp) THEN
                DF(I,J) = dNdX(node,J)
              ELSE
                DF(I,J) = 0.D0
              END IF
            END DO
          END DO

C         deltaC = F^T deltaF + deltaF^T F
          DO I=1,3
            DO J=1,3
              DC(I,J) = 0.D0
              DO K=1,3
                DC(I,J) = DC(I,J) + F(K,I)*DF(K,J)
     &                              + DF(K,I)*F(K,J)
              END DO
            END DO
          END DO

C         deltaE = 1/2 deltaC
          DO I=1,3
            DO J=1,3
              DE(I,J) = 0.5D0*DC(I,J)
            END DO
          END DO

C         Voigt components
          B0(1,col) = DE(1,1)
          B0(2,col) = DE(2,2)
          B0(3,col) = DE(3,3)
          B0(4,col) = 2.D0*DE(1,2)
          B0(5,col) = 2.D0*DE(2,3)
          B0(6,col) = 2.D0*DE(1,3)

        END DO
      END DO

      RETURN
      END

C================================================================
C  Add material part: KUU += integral(B0^T D B0 dOmega0)
C================================================================
      SUBROUTINE ADD_KUU_MATERIAL(B0,Dmat,NU,NSTR,KUU,wt)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NU,NSTR
      DOUBLE PRECISION B0(NSTR,NU), Dmat(NSTR,NSTR)
      DOUBLE PRECISION KUU(NU,NU), wt
      DOUBLE PRECISION DB(6,24)
      INTEGER I,J,K

C     DB = D * B0   (6 x NU)
      DO I=1,NSTR
        DO J=1,NU
          DB(I,J) = 0.D0
          DO K=1,NSTR
            DB(I,J) = DB(I,J) + Dmat(I,K)*B0(K,J)
          END DO
        END DO
      END DO

      DO I=1,NU
        DO J=1,NU
          DO K=1,NSTR
            KUU(I,J) = KUU(I,J) + B0(K,I)*DB(K,J)*wt
          END DO
        END DO
      END DO

      RETURN
      END

C================================================================
C  Numerical derivative L = d(logC)_voigt / dC_voigt  (6x6)
C================================================================
      SUBROUTINE DERIV_LOGC_NUM(C,Llog)
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION C(3,3), Llog(6,6)
      DOUBLE PRECISION Cplus(3,3), Cminus(3,3)
      DOUBLE PRECISION logP(3,3), logM(3,3)
      DOUBLE PRECISION vecP(6), vecM(6)
      DOUBLE PRECISION epsFD, baseC
      INTEGER I,J

      baseC = 0.D0
      DO I=1,3
        DO J=1,3
          IF (DABS(C(I,J)) .GT. baseC) baseC = DABS(C(I,J))
        END DO
      END DO
      IF (baseC .LT. 1.D0) baseC = 1.D0
      epsFD = 1.D-8 * baseC

      DO J=1,6

C       Copy original C
        DO I=1,3
          Cplus(I,1)  = C(I,1)
          Cplus(I,2)  = C(I,2)
          Cplus(I,3)  = C(I,3)
          Cminus(I,1) = C(I,1)
          Cminus(I,2) = C(I,2)
          Cminus(I,3) = C(I,3)
        END DO

        IF (J .EQ. 1) THEN
C         Perturb C11 (Voigt 1); step size = epsFD
          Cplus(1,1)  = Cplus(1,1)  + epsFD
          Cminus(1,1) = Cminus(1,1) - epsFD

        ELSEIF (J .EQ. 2) THEN
C         Perturb C22
          Cplus(2,2)  = Cplus(2,2)  + epsFD
          Cminus(2,2) = Cminus(2,2) - epsFD

        ELSEIF (J .EQ. 3) THEN
C         Perturb C33
          Cplus(3,3)  = Cplus(3,3)  + epsFD
          Cminus(3,3) = Cminus(3,3) - epsFD

        ELSEIF (J .EQ. 4) THEN
C         Perturb Voigt 4 term = 2*C12
C         *** Adjust shear perturbation: use 0.5*epsFD so 2*C12 step is epsFD ***
          Cplus(1,2)  = Cplus(1,2)  + 0.5D0*epsFD
          Cplus(2,1)  = Cplus(2,1)  + 0.5D0*epsFD
          Cminus(1,2) = Cminus(1,2) - 0.5D0*epsFD
          Cminus(2,1) = Cminus(2,1) - 0.5D0*epsFD

        ELSEIF (J .EQ. 5) THEN
C         Perturb Voigt 5 term = 2*C23
C         *** Adjust shear perturbation: use 0.5*epsFD ***
          Cplus(2,3)  = Cplus(2,3)  + 0.5D0*epsFD
          Cplus(3,2)  = Cplus(3,2)  + 0.5D0*epsFD
          Cminus(2,3) = Cminus(2,3) - 0.5D0*epsFD
          Cminus(3,2) = Cminus(3,2) - 0.5D0*epsFD

        ELSEIF (J .EQ. 6) THEN
C         Perturb Voigt 6 term = 2*C13
C         *** Adjust shear perturbation: use 0.5*epsFD ***
          Cplus(1,3)  = Cplus(1,3)  + 0.5D0*epsFD
          Cplus(3,1)  = Cplus(3,1)  + 0.5D0*epsFD
          Cminus(1,3) = Cminus(1,3) - 0.5D0*epsFD
          Cminus(3,1) = Cminus(3,1) - 0.5D0*epsFD
        END IF

        CALL LOGM_SYM3(Cplus,logP)
        CALL LOGM_SYM3(Cminus,logM)

        CALL STRAIN_TENSOR_TO_VOIGT(logP,vecP)
        CALL STRAIN_TENSOR_TO_VOIGT(logM,vecM)

        DO I=1,6
          Llog(I,J) = (vecP(I) - vecM(I)) / (2.D0*epsFD)
        END DO

      END DO

      RETURN
      END

C================================================================
C  Symmetric 3x3 eigen-decomposition (Jacobi) for SPD matrix
C================================================================
      SUBROUTINE EIGEN_SYM3(A,V,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION A(3,3), V(3,3), D(3)
      DOUBLE PRECISION B(3,3)
      DOUBLE PRECISION off, tol
      DOUBLE PRECISION app,aqq,apq,phi,c,s,bpj,bqj
      INTEGER i,j,iter,p,q,maxit

      maxit = 50
      tol   = 1.0D-12

C     copy A to B
      DO i=1,3
        DO j=1,3
          B(i,j) = A(i,j)
        END DO
      END DO

C     initialize V as identity
      DO i=1,3
        DO j=1,3
          V(i,j) = 0.D0
        END DO
        V(i,i) = 1.D0
      END DO

      DO iter = 1, maxit

        off = DABS(B(1,2)) + DABS(B(1,3)) + DABS(B(2,3))
        IF (off .LT. tol) GOTO 100

        DO p = 1,2
          DO q = p+1,3

            apq = B(p,q)
            IF (DABS(apq) .GT. 0.D0) THEN

              app = B(p,p)
              aqq = B(q,q)

              phi = 0.5D0*ATAN2(2.D0*apq, (aqq-app))
              c   = DCOS(phi)
              s   = DSIN(phi)

C             rotate B
              DO j = 1,3
                IF (j .NE. p .AND. j .NE. q) THEN
                  bpj = c*B(p,j) - s*B(q,j)
                  bqj = s*B(p,j) + c*B(q,j)
                  B(p,j) = bpj
                  B(j,p) = bpj
                  B(q,j) = bqj
                  B(j,q) = bqj
                END IF
              END DO

              B(p,p) = c*c*app - 2.D0*s*c*apq + s*s*aqq
              B(q,q) = s*s*app + 2.D0*s*c*apq + c*c*aqq
              B(p,q) = 0.D0
              B(q,p) = 0.D0

C             rotate eigenvectors V
              DO j = 1,3
                bpj = c*V(j,p) - s*V(j,q)
                bqj = s*V(j,p) + c*V(j,q)
                V(j,p) = bpj
                V(j,q) = bqj
              END DO

            END IF

          END DO
        END DO

      END DO

 100  CONTINUE

C     eigenvalues on diagonal of B
      DO i=1,3
        D(i) = B(i,i)
      END DO

      RETURN
      END

C================================================================
C  Matrix logarithm for symmetric positive definite 3x3:
C     LogA = log(A) via spectral decomposition
C================================================================
      SUBROUTINE LOGM_SYM3(A,LogA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION A(3,3), LogA(3,3)
      DOUBLE PRECISION V(3,3), D(3)
      INTEGER I,J,K

      CALL EIGEN_SYM3(A,V,D)

C     log of eigenvalues (clip to avoid log(<=0))
      DO I = 1,3
        IF (D(I) .LE. 0.D0) THEN
          D(I) = 1.D-16
        END IF
        D(I) = LOG(D(I))
      END DO

C     LogA = V * diag(D) * V^T
      DO I = 1,3
        DO J = 1,3
          LogA(I,J) = 0.D0
          DO K = 1,3
            LogA(I,J) = LogA(I,J) + V(I,K)*D(K)*V(J,K)
          END DO
        END DO
      END DO

      RETURN
      END
      SUBROUTINE LOGM_SYM3_AND_DERIV(C,LogC,Llog)
C================================================================
C  Compute LogC = log(C) and its Frechet derivative w.r.t. C in
C  Voigt form Llog (6x6) using ONE eigen-decomposition of C.
C
C  This removes the expensive finite-difference DERIV_LOGC_NUM and
C  also avoids doing EIGEN_SYM3 twice (once in LOGM_SYM3, once in
C  derivative), which typically yields a noticeable extra speedup.
C
C  Voigt convention MUST match STRAIN_TENSOR_TO_VOIGT:
C    vec = [A11, A22, A33, 2*A12, 2*A23, 2*A13]
C================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION C(3,3), LogC(3,3), Llog(6,6)
      DOUBLE PRECISION V(3,3), lam(3), loglam(3)
      DOUBLE PRECISION g(3,3)
      DOUBLE PRECISION E(3,3), Ehat(3,3), dLhat(3,3), dL(3,3)
      DOUBLE PRECISION vec(6)
      DOUBLE PRECISION diff, tol, lm, lmin
      INTEGER I,J,A,B,K,L

C     ---- eigen-decomposition C = V*diag(lam)*V^T ----
      CALL EIGEN_SYM3(C,V,lam)

C     ---- logs & safeguards (match LOGM_SYM3 clipping) ----
      lmin = 1.D-16
      DO I=1,3
        IF (lam(I) .LE. 0.D0) lam(I) = lmin
        loglam(I) = LOG(lam(I))
      END DO

C     ---- LogC = V * diag(loglam) * V^T ----
      DO I=1,3
        DO J=1,3
          LogC(I,J) = 0.D0
          DO K=1,3
            LogC(I,J) = LogC(I,J) + V(I,K)*loglam(K)*V(J,K)
          END DO
        END DO
      END DO

C     ---- build Frechet derivative coefficients g_ab ----
C     g_aa = 1/lam_a
C     g_ab = (log(lam_a)-log(lam_b))/(lam_a-lam_b)  (a!=b)
C     stable limit for nearly repeated eigenvalues
      tol = 1.D-12
      DO A=1,3
        DO B=1,3
          IF (A .EQ. B) THEN
            g(A,B) = 1.D0/lam(A)
          ELSE
            diff = lam(A) - lam(B)
            lm   = DMAX1(lam(A),lam(B))
            IF (DABS(diff) .LE. tol*lm) THEN
              g(A,B) = 2.D0/(lam(A) + lam(B))
            ELSE
              g(A,B) = (loglam(A) - loglam(B)) / diff
            END IF
          END IF
        END DO
      END DO

C     ---- columns of Llog: response of vec(logC) to unit vec(C) ----
      DO J=1,6

C       Build symmetric perturbation tensor E for unit Voigt comp J
        DO I=1,3
          DO K=1,3
            E(I,K) = 0.D0
          END DO
        END DO

        IF (J .EQ. 1) THEN
          E(1,1) = 1.D0
        ELSEIF (J .EQ. 2) THEN
          E(2,2) = 1.D0
        ELSEIF (J .EQ. 3) THEN
          E(3,3) = 1.D0
        ELSEIF (J .EQ. 4) THEN
C         Voigt 4 corresponds to 2*C12  => dC12=dC21=0.5
          E(1,2) = 0.5D0
          E(2,1) = 0.5D0
        ELSEIF (J .EQ. 5) THEN
C         Voigt 5 corresponds to 2*C23  => dC23=dC32=0.5
          E(2,3) = 0.5D0
          E(3,2) = 0.5D0
        ELSEIF (J .EQ. 6) THEN
C         Voigt 6 corresponds to 2*C13  => dC13=dC31=0.5
          E(1,3) = 0.5D0
          E(3,1) = 0.5D0
        END IF

C       Transform to eigen-basis: Ehat = V^T * E * V
        DO A=1,3
          DO B=1,3
            Ehat(A,B) = 0.D0
            DO K=1,3
              DO L=1,3
                Ehat(A,B) = Ehat(A,B) + V(K,A)*E(K,L)*V(L,B)
              END DO
            END DO
          END DO
        END DO

C       dLhat = g ??Ehat  (Hadamard product)
        DO A=1,3
          DO B=1,3
            dLhat(A,B) = g(A,B) * Ehat(A,B)
          END DO
        END DO

C       Back to physical basis: dL = V * dLhat * V^T
        DO K=1,3
          DO L=1,3
            dL(K,L) = 0.D0
            DO A=1,3
              DO B=1,3
                dL(K,L) = dL(K,L) + V(K,A)*dLhat(A,B)*V(L,B)
              END DO
            END DO
          END DO
        END DO

C       Symmetrize (guard against tiny asymmetry from numerics)
        DO I=1,3
          DO K=I+1,3
            dL(I,K) = 0.5D0*(dL(I,K) + dL(K,I))
            dL(K,I) = dL(I,K)
          END DO
        END DO

C       Convert to Voigt and store as column J
        CALL STRAIN_TENSOR_TO_VOIGT(dL,vec)
        DO I=1,6
          Llog(I,J) = vec(I)
        END DO

      END DO

      RETURN
      END


      SUBROUTINE UVARM(UVAR, DIRECT, T, TIME, DTIME, CMNAME, ORNAME,
     &                 NUVARM, NOEL, NPT, LAYER, KSPT, KSTEP, KINC,
     &                 NDI, NSHR, COORD, JMAC, JMAYP, MATLAYO, LACCFLA)
C     Post-process user variables retrieval
      IMPLICIT NONE
      CHARACTER*80 CMNAME, ORNAME
      CHARACTER*3  FLGRAY(15)
      INTEGER NUVARM, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      INTEGER NDI, NSHR, JMAC(*), JMAYP(*), MATLAYO, LACCFLA
      INTEGER baseElement, nstore
      INTEGER MAXE, MAXG, MAXV, NUMV, EOFF
      PARAMETER (MAXE=100000, MAXG=8, MAXV=20, NUMV=13, EOFF=1000000)
      DOUBLE PRECISION UVAR(NUVARM), DIRECT(3,3), T(3,3), TIME(2)
      DOUBLE PRECISION DTIME, COORD(*)
      DOUBLE PRECISION ARRAY(15)
      INTEGER JARRAY(15)
      DOUBLE PRECISION trackedPostVariables(MAXE,MAXG,MAXV)
      COMMON /POSTDATA/ trackedPostVariables
      LOGICAL POSTINIT
      COMMON /POSTFLAG/ POSTINIT
      SAVE /POSTDATA/, /POSTFLAG/
      DATA POSTINIT /.FALSE./

      UVAR(1:NUVARM) = 0.D0
      IF (NOEL.LE.EOFF) RETURN
      baseElement = NOEL - EOFF
      IF (baseElement.LT.1 .OR. baseElement.GT.MAXE) RETURN
      IF (NPT.LT.1 .OR. NPT.GT.MAXG) RETURN

      IF (NUVARM.GT.MAXV) THEN
        WRITE(*,*) "UVARM WARNING: NUVARM exceeds storage, truncating"
      END IF
      nstore = MIN(NUVARM,MAXV)
      UVAR(1:nstore) = trackedPostVariables(baseElement,NPT,1:nstore)
      RETURN
      END

      SUBROUTINE storePostVariables(element,npt,values,nvalues)
      IMPLICIT NONE
      INTEGER element,npt,nvalues
      DOUBLE PRECISION values(nvalues)
      INTEGER base, nstore
      INTEGER MAXE, MAXG, MAXV, NUMV, EOFF
      PARAMETER (MAXE=100000, MAXG=8, MAXV=20, NUMV=13, EOFF=1000000)
      DOUBLE PRECISION trackedPostVariables(MAXE,MAXG,MAXV)
      COMMON /POSTDATA/ trackedPostVariables
      LOGICAL POSTINIT
      COMMON /POSTFLAG/ POSTINIT
      SAVE /POSTDATA/, /POSTFLAG/

      base = element - EOFF
      IF (base.LT.1 .OR. base.GT.MAXE) RETURN
      IF (npt.LT.1 .OR. npt.GT.MAXG) RETURN
      nstore = MIN(nvalues,MAXV)
      trackedPostVariables(base,npt,1:nstore) = values(1:nstore)
      RETURN
      END


  
C======================================================================
C  SDVINI: UEL State Variable Initialization Subroutine
C
C  FUNCTIONALITY:
C    Mode 1: Read external mesh files and calculate initial in-situ
C            stress based on Z-coordinate depth.
C    Mode 2: Read an Abaqus generated .dat file, extract UVARM data
C            at a specific time step as initial fields, handling
C            element ID offsets (ID-1000000).
C======================================================================
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1                  LAYER,KSPT)
      INCLUDE 'ABA_PARAM.INC'

      INTEGER NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT
      DOUBLE PRECISION STATEV(NSTATV), COORDS(*)

C======================================================================
C  USER CONFIGURATION SECTION
C======================================================================
C  INI_MODE: 
C    1 = Calculate stress based on node depth (Reads .inp files)
C    2 = Map results from .dat file (Reads .dat file)
      INTEGER INI_MODE
      PARAMETER (INI_MODE = 1)

C  TARGET_TIME: Target analysis step time for extraction (Mode 2)
      DOUBLE PRECISION TARGET_TIME
      PARAMETER (TARGET_TIME = 1.00D0)

C  FILE PATH CONFIGURATION
      CHARACTER*260 FNODES, FELEMS, FDAT
      
C    -- Mode 1 File Paths --
      DATA FNODES /'nodes.inp'/
      DATA FELEMS /'elements.inp'/
      
C    -- Mode 2 File Path --
      DATA FDAT   /'job_output.dat'/

C  MODE 1 PARAMETERS: Initial Stress Profile
      DOUBLE PRECISION SIGV1, Z1, SIGV2, Z2, K0
      PARAMETER (SIGV1=-100.D0, Z1=0.0D0,
     &           SIGV2=   0.D0, Z2=10.D0,
     &           K0  =   0.428D0)

C======================================================================
C  COMMON BLOCKS & DATA CACHE
C======================================================================
      INTEGER MAXNODEID, MAXELEMID, NUM_DAT_VARS
      PARAMETER (MAXNODEID=500000, MAXELEMID=300000)
      PARAMETER (NUM_DAT_VARS=6)

C  -- Cache Arrays --
C  ZNODE:     Node Z coordinates (Mode 1)
C  UVARM_DAT: Variables read from .dat file (Mode 2)
      DOUBLE PRECISION ZNODE(MAXNODEID)
      DOUBLE PRECISION UVARM_DAT(MAXELEMID, NUM_DAT_VARS)
      
C  -- Status Flags --
C  NHAS:  Node existence flag
C  EHAS:  Element existence flag (Mode 1)
C  UHAS:  Element data existence flag (Mode 2)
C  ECONN: Element connectivity
      INTEGER NHAS(MAXNODEID)
      INTEGER ECONN(MAXELEMID,8)
      INTEGER EHAS(MAXELEMID)
      INTEGER UHAS(MAXELEMID)

      INTEGER MAXNID, MAXEID, NREADN, NREADE, NREADU, INITDONE
      
      COMMON /SDV_MESH_CACHE/ ZNODE, UVARM_DAT, NHAS, ECONN, 
     &                        EHAS, UHAS, MAXNID, MAXEID, 
     &                        NREADN, NREADE, NREADU, INITDONE
      SAVE /SDV_MESH_CACHE/
C!$OMP THREADPRIVATE(/SDV_MESH_CACHE/)

C  -- Local Variables --
      INTEGER I, gp, svStride, off, IERR
      LOGICAL IS_ELEM_BLOCK
      DOUBLE PRECISION r,s,t, zgp, sigv, sigh
      INTEGER nid(8)
      DOUBLE PRECISION NN1,NN2,NN3,NN4,NN5,NN6,NN7,NN8
      
C  -- Temporary variables for Mode 2 --
      DOUBLE PRECISION VALS(NUM_DAT_VARS)

C======================================================================
C  EXECUTION LOGIC
C======================================================================

C  1. Execute only for UEL (NCRDS=0)
      IF (NCRDS .NE. 0) RETURN

C  2. Global Initialization (Executed only once)
      IF (INITDONE .NE. 1) THEN
          
C         -- Reset flags --
          DO I=1, MAXELEMID
              EHAS(I) = 0
              UHAS(I) = 0
          END DO
          
          IF (INI_MODE .EQ. 1) THEN
C             --- Mode 1: Read Mesh ---
              CALL SDV_READ_MESH(FNODES,FELEMS,
     &                           ZNODE,NHAS,MAXNODEID,
     &                           ECONN,EHAS,MAXELEMID,
     &                           MAXNID,MAXEID,NREADN,NREADE,IERR)
              IF (IERR .NE. 0) THEN
                  CALL SDV_ABORT_FILES('Mesh Read Failed',NOEL,IERR,
     &                                 FNODES,FELEMS)
              END IF
          
          ELSE IF (INI_MODE .EQ. 2) THEN
C             --- Mode 2: Read DAT File ---
              CALL SDV_READ_DAT(FDAT, TARGET_TIME,
     &                          UVARM_DAT, UHAS, MAXELEMID,
     &                          NUM_DAT_VARS, NREADU, IERR)
              IF (IERR .NE. 0) THEN
                  CALL SDV_ABORT_FILES('DAT Read Failed',NOEL,IERR,
     &                                 FDAT, ' ')
              END IF
          END IF

          INITDONE = 1
      END IF

C  3. Determine SDV storage stride (Block storage check)
      IS_ELEM_BLOCK = .FALSE.
      IF (NSTATV .GE. 104 .AND. MOD(NSTATV,8) .EQ. 0) IS_ELEM_BLOCK=.TRUE.
      svStride = MERGE(NSTATV/8, NSTATV, IS_ELEM_BLOCK)

C  4. Apply initialization based on mode
      IF (INI_MODE .EQ. 1) THEN
C         ---------------------------------------------------
C         Mode 1: Geometric Depth Calculation
C         ---------------------------------------------------
          IF (NOEL .LT. 1 .OR. NOEL .GT. MAXELEMID) RETURN
          IF (EHAS(NOEL) .NE. 1) RETURN

          DO I=1,8
              nid(I) = ECONN(NOEL,I)
          END DO

          DO gp=1,8
              CALL SDV_GP_RST_2X2X2(gp,r,s,t)
              CALL SDV_SHAPE8_N_ONLY_SCALAR(r,s,t,
     &             NN1,NN2,NN3,NN4,NN5,NN6,NN7,NN8)

C             Interpolate Z coordinate at Gauss Point
              zgp = NN1*ZNODE(nid(1)) + NN2*ZNODE(nid(2)) 
     &            + NN3*ZNODE(nid(3)) + NN4*ZNODE(nid(4)) 
     &            + NN5*ZNODE(nid(5)) + NN6*ZNODE(nid(6))
     &            + NN7*ZNODE(nid(7)) + NN8*ZNODE(nid(8))

              CALL SDV_SIGV_LINEAR_CLAMP(SIGV1,Z1,SIGV2,Z2,zgp,sigv)
              sigh = K0*sigv
              
              off = MERGE((gp-1)*svStride, 0, IS_ELEM_BLOCK)

              STATEV(off+1) = sigh
              STATEV(off+2) = sigh
              STATEV(off+3) = sigv
              DO I=4,svStride
                  STATEV(off+I) = 0.D0
              END DO
          END DO

      ELSE IF (INI_MODE .EQ. 2) THEN
C         ---------------------------------------------------
C         Mode 2: Mapping from DAT Data
C         ---------------------------------------------------
          IF (NOEL .LT. 1 .OR. NOEL .GT. MAXELEMID) RETURN
          
C         Check if data exists for this element
          IF (UHAS(NOEL) .NE. 1) THEN
              ! Maintain 0 or handle default behavior if not found
              RETURN
          END IF

C         Retrieve cached data
          DO I=1, NUM_DAT_VARS
              VALS(I) = UVARM_DAT(NOEL, I)
          END DO

C         Assign to all integration points (Assuming C3D8R centroid values)
          DO gp=1,8
              off = MERGE((gp-1)*svStride, 0, IS_ELEM_BLOCK)
              
C             Map UVARM1-6 to STATEV 1-6
              STATEV(off+1) = VALS(1)
              STATEV(off+2) = VALS(2)
              STATEV(off+3) = VALS(3)
              STATEV(off+4) = VALS(4)
              STATEV(off+5) = VALS(5)
              STATEV(off+6) = VALS(6)
              
C             Zero out remaining variables
              DO I=7,svStride
                  STATEV(off+I) = 0.D0
              END DO
          END DO

      END IF

      RETURN
      END

C======================================================================
C  SUBROUTINE: Read DAT File (Core Logic)
C======================================================================
      SUBROUTINE SDV_READ_DAT(FNAME, TARGET_T, CACHE, HAS_FLAG, 
     &                        MAXE, NVARS, NREAD, IERR)
      IMPLICIT NONE
      CHARACTER*(*) FNAME
      DOUBLE PRECISION TARGET_T
      INTEGER MAXE, NVARS, NREAD, IERR
      DOUBLE PRECISION CACHE(MAXE, NVARS)
      INTEGER HAS_FLAG(MAXE)

      INTEGER IU, IOS, I, OFFSET, PT_ID, ID_RAW, ID_MAP
      CHARACTER*512 LINE
      DOUBLE PRECISION T_READ
      LOGICAL TIME_MATCHED, TABLE_FOUND, IN_DATA_BLOCK

C     Element ID Offset (Read ID - OFFSET = Internal Storage ID)
      PARAMETER (OFFSET = 1000000)

      IERR = 0
      NREAD = 0
      TIME_MATCHED = .FALSE.
      TABLE_FOUND = .FALSE.
      IN_DATA_BLOCK = .FALSE.

C     Initialize flags
      DO I=1, MAXE
          HAS_FLAG(I) = 0
      END DO

      IU = 155
      OPEN(UNIT=IU, FILE=FNAME, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
          IERR = 301 ! Unable to open file
          RETURN
      END IF

C     --- Line-by-line Read Loop ---
 100  READ(IU, '(A)', IOSTAT=IOS) LINE
      IF (IOS .NE. 0) GOTO 900 ! End of file

C     Phase 1: Search for Target Time (TOTAL TIME COMPLETED)
      IF (.NOT. TIME_MATCHED) THEN
          I = INDEX(LINE, 'TOTAL TIME COMPLETED')
          IF (I .GT. 0) THEN
C             Found keyword, read value after it
C             Format: ..., TOTAL TIME COMPLETED         1.00
C             Read starting from offset + 20 chars
              READ(LINE(I+20:), *, IOSTAT=IOS) T_READ
              IF (IOS .EQ. 0) THEN
C                 Compare float with tolerance 1e-4
                  IF (ABS(T_READ - TARGET_T) .LT. 1.0D-4) THEN
                      TIME_MATCHED = .TRUE.
                      WRITE(6,*) 'SDVINI: DAT Target Time Found:', T_READ
                  END IF
              END IF
          END IF
          GOTO 100
      END IF

C     Phase 2: Search for Data Table Header (ELEMENT OUTPUT)
      IF (TIME_MATCHED .AND. .NOT. TABLE_FOUND) THEN
C         Check if a new time increment has started (missed the table)
          IF (INDEX(LINE, 'TIME INCREMENT COMPLETED').GT.0) THEN
              GOTO 900 ! Stop reading if we moved to next step
          END IF

C         Search for specific header
          IF (INDEX(LINE, 'ELEMENT   PT FOOT-').GT.0) THEN
              TABLE_FOUND = .TRUE.
              IN_DATA_BLOCK = .TRUE.
C             Read and skip next line (usually the NOTE line)
              READ(IU, '(A)') LINE 
              GOTO 100
          END IF
          GOTO 100
      END IF

C     Phase 3: Read Data Lines
      IF (IN_DATA_BLOCK) THEN
C         Pre-read first integer (ID) to check if it's a valid data line
C         Empty lines, page breaks, or text will cause IOS error or ID=0
          READ(LINE, *, IOSTAT=IOS) ID_RAW
          
          IF (IOS .NE. 0 .OR. ID_RAW .EQ. 0) THEN
C             Read failed, likely reached end of table
C             Check if we have completely exited the table block
              IF (INDEX(LINE, 'MAXIMUM').GT.0 .OR. 
     &            INDEX(LINE, 'MINIMUM').GT.0) THEN
                  IN_DATA_BLOCK = .FALSE.
C                 Task done, exit or continue to other tables
                  GOTO 900
              END IF
              GOTO 100 ! Skip empty line and continue
          END IF

C         Calculate Mapped ID
          ID_MAP = ID_RAW - OFFSET

          IF (ID_MAP .GT. 0 .AND. ID_MAP .LE. MAXE) THEN
C             Read Data: ID, PT, UVARM1, UVARM2, UVARM3, UVARM4, UVARM5, UVARM6
C             Note: There is a FOOTNOTE column. List-directed input (*)
C             skips spaces. If FOOTNOTE has characters, format adjustment needed.
              READ(LINE, *, IOSTAT=IOS) ID_RAW, PT_ID, 
     &             CACHE(ID_MAP, 1), CACHE(ID_MAP, 2), 
     &             CACHE(ID_MAP, 3), CACHE(ID_MAP, 4), 
     &             CACHE(ID_MAP, 5), CACHE(ID_MAP, 6)
              
              IF (IOS .EQ. 0) THEN
                  HAS_FLAG(ID_MAP) = 1
                  NREAD = NREAD + 1
              END IF
          END IF
      END IF

      GOTO 100

 900  CLOSE(IU)
      IF (NREAD .EQ. 0) IERR = 303 ! No valid data read
      RETURN
      END

C======================================================================
C  SUBROUTINE: Read Mesh File (.inp)
C======================================================================
      SUBROUTINE SDV_READ_MESH(FNODES,FELEMS,ZNODE,NHAS,MAXNODE,
     &                         ECONN,EHAS,MAXELEM,MAXNID,MAXEID,
     &                         NREADN,NREADE,IERR)
      IMPLICIT NONE
      CHARACTER*(*) FNODES, FELEMS
      INTEGER MAXNODE, MAXELEM, IERR, MAXNID, MAXEID, NREADN, NREADE
      DOUBLE PRECISION ZNODE(MAXNODE)
      INTEGER NHAS(MAXNODE), ECONN(MAXELEM,8), EHAS(MAXELEM)
      
      INTEGER i, iu, ios, p, nid, eid, n(8)
      DOUBLE PRECISION x, y, z
      CHARACTER*512 line

      IERR=0; MAXNID=0; MAXEID=0; NREADN=0; NREADE=0
      DO i=1,MAXNODE
        NHAS(i)=0; ZNODE(i)=0.D0
      END DO
      DO i=1,MAXELEM
        EHAS(i)=0
      END DO

C     ----- Read Nodes -----
      iu = 151
      OPEN(UNIT=iu,FILE=FNODES,STATUS='OLD',IOSTAT=ios)
      IF (ios.NE.0) THEN 
        IERR=101
        RETURN
      ENDIF
 100  READ(iu,'(A)',IOSTAT=ios) line
      IF (ios.NE.0) GOTO 110
      CALL SDV_FIND_START(line, p)
      IF (p.EQ.0 .OR. line(p:p).EQ.'*') GOTO 100
      READ(line, *, IOSTAT=ios) nid, x, y, z
      IF (ios.EQ.0 .AND. nid.GT.0 .AND. nid.LE.MAXNODE) THEN
        NHAS(nid)=1; ZNODE(nid)=z; NREADN=NREADN+1
        IF(nid.GT.MAXNID) MAXNID=nid
      ENDIF
      GOTO 100
 110  CLOSE(iu)

C     ----- Read Elements -----
      iu = 152
      OPEN(UNIT=iu,FILE=FELEMS,STATUS='OLD',IOSTAT=ios)
      IF (ios.NE.0) THEN 
        IERR=201
        RETURN
      ENDIF
 200  READ(iu,'(A)',IOSTAT=ios) line
      IF (ios.NE.0) GOTO 210
      CALL SDV_FIND_START(line, p)
      IF (p.EQ.0 .OR. line(p:p).EQ.'*') GOTO 200

      READ(line, *, IOSTAT=ios) eid, (n(i),i=1,8)
      IF (ios.EQ.0 .AND. eid.GT.0 .AND. eid.LE.MAXELEM) THEN
        EHAS(eid)=1
        DO i=1,8
          ECONN(eid,i)=n(i)
        END DO
        NREADE=NREADE+1
        IF(eid.GT.MAXEID) MAXEID=eid
      ENDIF
      GOTO 200
 210  CLOSE(iu)
      IF (NREADN.EQ.0) IERR=103
      IF (NREADE.EQ.0) IERR=203
      RETURN
      END

C======================================================================
C  HELPER SUBROUTINES
C======================================================================
      SUBROUTINE SDV_FIND_START(LINE, POS)
      CHARACTER*(*) LINE
      INTEGER POS, I
      POS = 0
      DO I=1, LEN(LINE)
        IF (LINE(I:I).NE.' '.AND.ICHAR(LINE(I:I)).NE.9) THEN
          POS = I
          RETURN
        ENDIF
      ENDDO
      END

      SUBROUTINE SDV_GP_RST_2X2X2(gp,r,s,t)
      INTEGER gp
      DOUBLE PRECISION r,s,t, g
      PARAMETER (g=0.577350269189626D0)
      r=MERGE(-g,g,gp.EQ.1.OR.gp.EQ.4.OR.gp.EQ.5.OR.gp.EQ.8)
      s=MERGE(-g,g,gp.EQ.1.OR.gp.EQ.2.OR.gp.EQ.5.OR.gp.EQ.6)
      t=MERGE(-g,g,gp.LE.4)
      END

      SUBROUTINE SDV_SHAPE8_N_ONLY_SCALAR(r,s,t,N1,N2,N3,N4,N5,N6,N7,N8)
      DOUBLE PRECISION r,s,t,N1,N2,N3,N4,N5,N6,N7,N8,rp,rm,sp,sm,tp,tm
      rp=1.D0+r; rm=1.D0-r; sp=1.D0+s; sm=1.D0-s; tp=1.D0+t; tm=1.D0-t
      N1=0.125D0*rm*sm*tm; N2=0.125D0*rp*sm*tm; N3=0.125D0*rp*sp*tm
      N4=0.125D0*rm*sp*tm; N5=0.125D0*rm*sm*tp; N6=0.125D0*rp*sm*tp
      N7=0.125D0*rp*sp*tp; N8=0.125D0*rm*sp*tp
      END

      SUBROUTINE SDV_SIGV_LINEAR_CLAMP(SIG1,Z1,SIG2,Z2,z,sigv)
      DOUBLE PRECISION SIG1,Z1,SIG2,Z2,z,sigv,den,tt,zmin,zmax
      den = Z2 - Z1
      IF (ABS(den).LT.1D-14) THEN
        sigv=SIG1; RETURN
      ENDIF
      zmin=MIN(Z1,Z2); zmax=MAX(Z1,Z2)
      tt = (MAX(zmin, MIN(zmax, z)) - Z1) / den
      sigv = SIG1 + tt*(SIG2-SIG1)
      END

      SUBROUTINE SDV_ABORT_FILES(MSG,NOEL,IERR,FN,FE)
      CHARACTER*(*) MSG, FN, FE
      INTEGER NOEL, IERR
      WRITE(6,*) '*** SDVINI FATAL ERROR ***'
      WRITE(6,*) 'Message: ', MSG, ' IERR: ', IERR
      WRITE(6,*) 'NOEL: ', NOEL
      WRITE(6,*) 'File 1: ', FN
      WRITE(6,*) 'File 2: ', FE
      CALL XIT
      END


c ======================================================================
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,
     &  ddsddt,drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,
     &  dpred,cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,
     &  pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c ======================================================================
c UMAT: linear-elastic material
c ----------------------------------------------------------------------
      implicit none
      character(80) cmname
      integer(4) ntens,ndi, nshr, nstatv, nprops, noel, npt,
     & layer, kspt, kstep, kinc
      real(8) stress(ntens), statev(nstatv), ddsdde(ntens,ntens),
     &  ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     &  time(2), predef(1), dpred(1), props(nprops), coords(3),
     &  drot(3,3), dfgrd0(3,3), dfgrd1(3,3), sse, spd, scd, rpl,
     &  drpldt, dtime, temp, dtemp, pnewdt, celent

      integer i, j, k
      real(8) E, possion

      E  = PROPS(1)
      possion = PROPS(2)

      if (ntens.eq.3) then
        ddsdde = 0.0d0
        ddsdde(1,1) = 1.0d0
        ddsdde(2,2) = 1.0d0
        ddsdde(1,2) = possion/(1.0d0-possion)
        ddsdde(2,1) = possion/(1.0d0-possion)
        ddsdde(3,3) = (1.0d0-2.0d0*possion)/(2.0d0*(1.0d0-possion))
        ddsdde = ddsdde*E*(1.0d0-possion)/
     &        ((1.0d0+possion)*(1.0d0-2.0d0*possion)) 
      elseif (ntens.eq.6) then
        ddsdde = 0.0d0
        ddsdde(1,1) = 1.0d0
        ddsdde(2,2) = 1.0d0
        ddsdde(3,3) = 1.0d0
        ddsdde(1,2) = possion/(1-possion)
        ddsdde(2,1) = possion/(1-possion)
        ddsdde(2,3) = possion/(1-possion)
        ddsdde(3,2) = possion/(1-possion)
        ddsdde(1,3) = possion/(1-possion)
        ddsdde(3,1) = possion/(1-possion)
        ddsdde(4,4) = (1-2*possion)/(2*(1-possion))
        ddsdde(5,5) = (1-2*possion)/(2*(1-possion))
        ddsdde(6,6) = (1-2*possion)/(2*(1-possion))
        ddsdde = ddsdde*E*(1-possion)/((1+possion)*(1-2*possion))
      endif

      do i=1,ntens
        do j=1,ntens
          stress(i) = stress(i) + ddsdde(i,j)*dstran(j)
        enddo           
      enddo      

      return
      end subroutine umat
