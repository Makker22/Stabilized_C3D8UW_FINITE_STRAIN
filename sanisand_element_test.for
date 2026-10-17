! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
      Subroutine MZEROR(R,K)
C
C***********************************************************************
C
C     Function: To make a real array R with dimension K to zero
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J) = 0.0D0
      End Do

      Return
      End


      Subroutine MZEROI(I,K)
C
C***********************************************************************
C
C     Function: To make an integre array I with Dimension K to zero
C
C***********************************************************************
C
      Dimension I(*)

      Do J=1,K
        I(J)=0
      End Do

      Return
      End

      Subroutine SETRVAL(R,K,V)
C
C***********************************************************************
C
C     Function: To fill a real array R with Dimension K with value V
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J)=V
      End Do

      Return
      End

      Subroutine SETIVAL(I,K,IV)
C
C***********************************************************************
C
C     Function: To fill an integer array I with Dimension K with value IV
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension I(*)

      Do J=1,K
        I(J)=IV
      End Do

      Return
      End

      Subroutine COPYIVEC(I1,I2,K)
C
C***********************************************************************
C
C     Function: To copy an integer array I1 with Dimension K to I2
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension I1(*),I2(*)

      Do  J=1,K
        I2(J)=I1(J)
      End Do

      Return
      End

      Subroutine COPYRVEC(R1,R2,K)
C
C***********************************************************************
C
C     Function: To copy a Double array R1 with Dimension K to R2
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R1(*),R2(*)

      Do J=1,K
        R2(J)=R1(J)
      End Do

      Return
      End


      Logical Function IS0ARR(A,N)
C
C***********************************************************************
C    Function :  To check whether a real array contains only zero values.
C                When an array contains only zero's is might not need to be
C                written to the XXX file.
C                exit Function when first non-zero value occured or when
C                all elements are checked and are zero.
C
C    Input:  A : array to be checked
C            N : number of elements in array that should be checked
C
C    Output : .TRUE.  when all elements are 0
C             .FALSE. when at least one element is not zero
C
C    Called by :  Subroutine TOBXX
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*)
      Is0Arr=.False.
      Do I=1,N
        If ( A(I) .Ne. 0 ) Return
      End Do
      Is0Arr=.True.
      Return
      End

      Logical Function IS0IARR(IARR,N)
C
C***********************************************************************
C    Function :  To check whether a integer array contains only zero values.
C                Similar to IS0ARR
C
C    Input:  IARR : array to be checked
C            N    : number of elements in array that should be checked
C
C    Output : .TRUE.  when all elements are 0
C             .FALSE. when at least one element is not zero
C
C    Called by :  Subroutine TOBXX
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension IARR(*)

      Is0IArr=.False.
      Do I=1,N
        If ( IARR(I) .Ne. 0 ) Return
      End Do
      Is0IArr=.True.
      Return
      End
C***********************************************************************
      Subroutine MulVec(V,F,K)
C***********************************************************************
C
C     Function: To multiply a real vector V with dimension K by F
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*)

      Do J=1,K
        V(J)=F*V(J)
      End Do

      Return
      End     ! Subroutine Mulvec
C***********************************************************************
      Subroutine MatVec(xMat,IM,Vec,N,VecR)
C***********************************************************************
C
C     Calculate VecR = xMat*Vec
C
C I   xMat  : (Square) Matrix (IM,*)
C I   Vec   : Vector
C I   N     : Number of rows/colums
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine MatVec

C***********************************************************************
      Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
C***********************************************************************
C
C     Calculate VecR() = R1*Vec1()+R2*Vec2()
C
C I   Vec1,
C I   Vec2  : Vectors
C I   R1,R2 : Multipliers
C I   N     : Number of rows
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec
C
C***********************************************************************
      Double Precision Function DInProd(A,B,N)
C***********************************************************************
C
C     Returns the Inproduct of two vectors
C
C I   A,B  : Two vectors
C I   N    : Used length of vectors
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*),B(*)
C***********************************************************************

      X = 0
      Do I=1,N
        X = X + A(I)*B(I)
      End Do
      DInProd = X
      Return
      End     ! Function DInProd
C
C***********************************************************************
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2
C
C I   xMat1 : Matrix (Id1,*)
C I   xMat2 : Matrix (Id2,*)
C I   nR1   : Number of rows in resulting matrix    (= No rows in xMat1)
C I   nC2   : Number of columns in resulting matrix (= No cols in xMat2)
C I   nC1   : Number of columns in matrix xMat1
C             = Number  rows    in matrix xMat2
C O   xMatR : Resulting matrix (IdR,*)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(Id1,*),xMat2(Id2,*),xMatR(IdR,*)
C**********************************************************************

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMat

C***********************************************************************
      Subroutine MatMatSq(n, xMat1, xMat2, xMatR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2 for square matrices, size n
C
C I   n     : Dimension of matrices
C I   xMat1 : Matrix (n,*)
C I   xMat2 : Matrix (n,*)
C O   xMatR : Resulting matrix (n,*)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(n,*),xMat2(n,*),xMatR(n,*)
C**********************************************************************

      Do I=1,n
        Do J=1,n
          X=0
          Do K=1,n
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMatSq

C***********************************************************************
      Subroutine WriVal ( io, C , V )
C***********************************************************************
C
C Write (Double) value to file unit io (when io>0)
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,V
    1 Format( A,3x, 1x,1p,e12.5)
      Return
      End
C***********************************************************************
      Subroutine WriIVl ( io, C , I )
C***********************************************************************
C
C Write (integer) value to file unit io (when io>0)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,I
    1 Format( A,3x, 1x,I6)
      Return
      End
C***********************************************************************
      Subroutine WriIVc ( io, C , iV , n )
C***********************************************************************
C
C Write (integer) vector to file unit io (when io>0)
C
C***********************************************************************
      Character C*(*)
      Dimension iV(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) (iv(i),i=1,n)
    1 Format( ( 2(3x,5i4) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec ( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      If (Len_Trim(C).Le.6) Then
        Write(io,2) C,( V(i),i=1,n)
      Else
        Write(io,*) C
        Write(io,1) ( V(i),i=1,n)
      End If
    1 Format( ( 2(1x, 3(1x,1p,e10.3) ) ) )
    2 Format( A, ( T7, 2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec5( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 5 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) ( V(i),i=1,n)
    1 Format( 5(1x,1p,e12.5) )
      Return
      End
C***********************************************************************
      Subroutine WriMat ( io, C , V , nd, nr, nc )
C***********************************************************************
C
C Write (Double) matrix to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(nd,*)

      If (io.Le.0) Return

      Write(io,*) C
      Do j=1,nr
        Write(io,1) j,( V(j,i),i=1,nc)
      End Do
    1 Format(i4, (  T7,2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************

      Subroutine MatInvPiv(Aorig,B,N)
      Implicit Double Precision (A-H,O-Z)
      Dimension Aorig(n,*), B(n,*),A(:,:)
      Allocatable :: A
      Allocate ( A(n,n) ) ! No error checking !!
      Call CopyRVec(AOrig, A, n*n )
      Call MZeroR(B,n*n)
      Do i=1,n
        B(i,i) = 1d0
      End Do
      Do I=1,n
        T=A(I,I)
        iPiv=i
        Do j=i+1,n
          If ( Abs(A(j,i)) .Gt. Abs(A(iPiv,i))  ) iPiv=j
        End Do
        If (iPiv.Ne.i) Then
          Do j=1,n
            x         = A( i  ,j)
            A( i  ,j) = A(iPiv,j)
            A(iPiv,j) = x
            x         = B( i  ,j)
            B( i  ,j) = B(iPiv,j)
            B(iPiv,j) = x
          End Do
          T=A(I,I)
        End If
        Do J=1,n
          A(I,J)=A(I,J)/T
          B(I,J)=B(I,J)/T
        End Do
        Do K=1,n
          If (K.Ne.I) Then
            T=A(K,I)
            Do J=1,n
              A(K,J)=A(K,J)-T*A(I,J)
              B(K,J)=B(K,J)-T*B(I,J)
            End Do
          End If
        End Do
      End Do
      DeAllocate ( A  )
      Return
      End ! MatinvPiv

C***********************************************************************
      Subroutine PrnSig(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
C***********************************************************************
      Subroutine Eig_3(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),
     *          xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And.
     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
!          If (a(ip,iq) .Ne. 0.0) Then
          If (Abs(a(ip,iq)) .Gt. 1d-50) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/Dsqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            v1p=c*v(1,ip)-s*v(1,iq)
            v2p=c*v(2,ip)-s*v(2,iq)
            v3p=c*v(3,ip)-s*v(3,iq)
            v(1,iq)=s*v(1,ip)+c*v(1,iq)
            v(2,iq)=s*v(2,ip)+c*v(2,iq)
            v(3,iq)=s*v(3,ip)+c*v(3,iq)
            v(1,ip)=v1p
            v(2,ip)=v2p
            v(3,ip)=v3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = DSqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      Do i=1,3
        xN1(i) = v(i,is1) ! first  column
        xN2(i) = v(i,is2) ! second column
        xN3(i) = v(i,is3) ! third  column
      End Do
      Return
      End ! Eig_3

      Subroutine Eig_3a(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      Do While ( it.lt.itmax .And.
     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
!          If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
          If (Abs(a(ip,iq)) .Gt. 1d-50) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/Dsqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = DSqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      Return
      End ! Eig_3a

C
C***********************************************************************
      Logical Function LEqual(A,B,Eps)
C***********************************************************************
C
C     Returns .TRUE.  when two real values are (almost) equal,
C             .FALSE. otherwise
C
C I   A,B  : Two real values to be compared
C I   Eps  : Toleration (Magnitude ~= 1E-5)
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
C***********************************************************************
      LEqual =.True.
      If (A .Eq. B) Return
      If (DAbs(A-B) .LT. 0.5D0*Eps*( DAbs(A) + DAbs(B) + Eps ) )Return
      LEqual =.False.
      Return
      End     ! function LEqual
C
C***********************************************************************
      Subroutine CrossProd(xN1,xN2,xN3)
C***********************************************************************
C
C     Returns cross product of xN1 and xN2
C
C I   xN1,xN2 : Two basic vectors
C O   xN3     : Resulting vector
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*)
C***********************************************************************

      xN3(1) = xN1(2)*xN2(3) - xN1(3)*xN2(2)
      xN3(2) = xN1(3)*xN2(1) - xN1(1)*xN2(3)
      xN3(3) = xN1(1)*xN2(2) - xN1(2)*xN2(1)

      Return
      End     ! Subroutine CrossProd
C
C***********************************************************************
      Double Precision Function ArcSin(X,ie)
C***********************************************************************
C
C     Returns the Arc Sine of X
C
C I   X : Input value
C
C     Note : In stead of using default routine DASIN we use this one
C            because ³X³ can be slightly beyond 1 and this will give
C            a RTE using DASIN(X)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
C***********************************************************************
      Ie=0
      S = (1-X*X)
!      If (S .Lt. -1E-10) Ie=1
!      If (S .Lt. -1E-10) Write(*,1) X,S
!      If (S .Lt. -1E-10) Write(2,1) X,S
    1 Format(' ArcSin(',1x,1p,e13.5e3,') , S =',1x,1p,e13.5e3)
      If (S.LT.0) S = 0
      S = DSQRT(S)
      ArcSin = DATan2(X,S)
      Return
      End     ! function ArcSin
C
C***********************************************************************
      Subroutine CarSig(S1,S2,S3,xN1,xN2,xN3,SNew)
C***********************************************************************
C
C     Returns the Cartesian stresses using the principal stresses S1..S3
C     and the principal directions
C
C I   S1..S3   : Principal stresses
C I   xN1..xN3 : Principal directions (xNi for Si)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*),SNew(*)
      Dimension SM(3,3),T(3,3),TT(3,3),STT(3,3)
C***********************************************************************
C
C**** Fill transformation (rotation) matrix
C
      Do I=1,3
        T(I,1) = xN1(I)
        T(I,2) = xN2(I)
        T(I,3) = xN3(I)
        TT(1,I) = T(I,1)
        TT(2,I) = T(I,2)
        TT(3,I) = T(I,3)
      End Do
!      Call MatTranspose(T,3,TT,3,3,3)

      Call MZeroR(SM,9)
      SM(1,1) = S1
      SM(2,2) = S2
      SM(3,3) = S3
C
C**** SMnew = T*SM*TT
C
      Call MatMat(SM ,3,  TT,3 , 3,3,3 ,STT,3)
      Call MatMat( T ,3, STT,3 , 3,3,3 ,SM ,3)
!     Call MatMatSq(3, SM,  TT, STT )   ! STT = SM*TT
!     Call MatMatSq(3,  T, STT, SM  )   ! SM  =  T*STT
C
C**** Extract cartesian stress vector from stress matrix
C
      Do I=1,3
        SNew(I) = SM(I,I)
      End Do
      SNew(4) = SM(2,1)
      SNew(5) = SM(3,2)
      SNew(6) = SM(3,1)

      Return
      End     ! Subroutine CarSig
C**********************************************************************
      subroutine setveclen(xn,n,xl)
C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN(*)
      x=0
      do i=1,n
        x=x+xn(i)**2
      end do
      if (x.Ne.0) Then
        f=xl/Dsqrt(x)
        do i=1,3
          xn(i)=xn(i)*f
        end do
      end if
      return
      end ! setveclen

C**********************************************************************
C End Of file
C**********************************************************************

!DIR$ FREEFORM

! ======== USRADDDF.f90 ========
! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
!----------------------------------------------------------------------------------------
! Subroutines in this file:
!
!  Subroutine GetModelCount( nMod )
!  Subroutine GetModelName ( iMod , ModelName )
!  Subroutine GetParamCount( iMod , nParam )
!  Subroutine GetParamName ( iMod , iParam, ParamName )
!  Subroutine GetParamUnit ( iMod , iParam, Units )
!  Subroutine GetStateVarCount( iMod , nVar )
!  Subroutine GetStateVarName ( iMod , iVar, Name )
!  Subroutine GetStateVarUnit ( iMod , iVar, Unit )
!
! Local:
!  Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!  Subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
!----------------------------------------------------------------------------------------
module extraInfo
      
double precision :: oneAtm = 100d0
double precision :: onekPa =   1d0
      
end module extraInfo
    
Subroutine setDLLExtraInfo( iInfo, rInfo )
use extraInfo, only: oneAtm, onekPa
implicit none      
integer,          intent(in) :: iInfo(*)
double precision, intent(in) :: rInfo(*)
      
!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: setDLLExtraInfo

! Makes it possible to use extra info from kernel
!iInfo ( 1 ) = iVal !Special option from PLAXIS Input Phase Settings
!rInfo ( 1 ) = OneAtm
!iInfo ( 2 ) = TempUnit

oneAtm = rInfo(1)
onekPa = oneAtm / 100.0d0
      
End Subroutine setDLLExtraInfo


Subroutine GetModelCount(nMod)!renamed
!
! Return the maximum model number (iMod) in this DLL
!
implicit none
Integer (Kind=4), intent(out) :: nMod

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetModelCount

nMod = 1 ! Maximum model number (iMod) in current DLL

Return
End ! GetModelCount

Subroutine GetModelName( iMod , ModelName )
!
! Return the name of the different models
!
implicit none
Integer, intent(in) :: iMod
Character (Len= * ), intent(out) :: ModelName
Character (Len=255) tName

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetModelName
!someone ssigne tName
Select Case (iMod)
Case (1)
    ModelName = 'SANISANDMS'
Case Default
    ModelName = 'not in DLL'
End Select
call Add_Str_Length( ModelName )

Return
End ! GetModelName

Subroutine GetParamCount( iMod , nParam )
!
! Return the number of parameters of the different models
!
implicit none
integer, intent(in) :: iMod
integer, intent(out) :: nParam

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamCount

Select Case (iMod)
Case ( 1 )
    nParam = 25
Case Default
    nParam = 0
End Select
Return
End ! GetParamCount

Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!
! Return the parameters name and units of the different models
!
! Units: use F for force unit
!            L for length unit
!            T for time unit
!
implicit none
integer, intent(in) :: iMod , iParam
Character (Len=255), intent(out) :: ParamName, Units

Select Case (iMod)
Case (1)
    ! ModName = 'SANISANDMS'
    Select Case (iParam)
        Case (1)
            ParamName = 'G_0#'       ; Units     = '-'
        Case (2)
            ParamName = '@n'      ; Units     = '-'
        Case (3)
            ParamName = 'M_c#'   ; Units     = '-'
        Case (4)
            ParamName = 'c'    ; Units     = '-'
        Case (5)
            ParamName = '@l_c#'       ; Units     = '-'
        Case (6)
            ParamName = 'e_0#'       ; Units     = '-'
        Case (7)
            ParamName = '@x'       ; Units     = '-'
        Case (8)
            ParamName = 'm'       ; Units     = '-'
        Case (9)
            ParamName = 'h_0#'       ; Units     = '-'
        Case (10)
            ParamName = 'c_h#'       ; Units     = '-'
        Case (11)
            ParamName = 'n^b#'       ; Units     = '-'
        Case (12)
            ParamName = 'A_0#'       ; Units     = '-'
        Case (13)
            ParamName = 'n^d#'       ; Units     = '-'
        Case (14)
            ParamName = '@m_0#'       ; Units     = '-'
        Case (15)
            ParamName = '@z'       ; Units     = '-'
        Case (16)
            ParamName = '@b'       ; Units     = '-'
        Case (17)
            ParamName = 'TolF'       ; Units     = '-'
        Case (18)
            ParamName = 'TolR'       ; Units     = '-'
        Case (19)
            ParamName = 'e_ini#'       ; Units     = '-'
        Case (20)
            ParamName = 'e_max#'       ; Units     = '-'
        Case (21)
            ParamName = 'e_min#'       ; Units     = '-'  
        Case (22)
            ParamName = 'M_c#^b#max'       ; Units     = '-'
        Case (23)
            ParamName = 'ElScheme#'       ; Units     = '-'
        Case (24)
            ParamName = 'h^M#ExcEqu#'     ; Units     = '-' 
        Case (25)
            ParamName = 'h_Max#'     ; Units     = '-'   
        Case Default
            ParamName = '???'     ; Units     = '???'
    End Select  
Case Default
    ! model not in DLL
    ParamName = ' N/A '     ; Units     = ' N/A '
End Select

Return
    End ! GetParamAndUnit

!
! Interfaces to routines to report parameter names, counts etc.
!
Subroutine Add_Str_Length( aString )
Implicit None
Character*(*) aString
Character *255 tString
Integer Lt              ! length of incoming string
!
! routine should add the length of the string as the first character
!
tString = aString
Lt      = Len_Trim(tString)
aString = Char(Lt) // tString(1:Lt)
!
End Subroutine Add_Str_Length

Subroutine GetParamName( iMod , iParam, ParamName )
!
! Return the parameters name of the different models
!
Implicit None
Integer  iMod, iParam
Character (Len=255) ParamName, Units

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamName

Call GetParamAndUnit(iMod,iParam,ParamName,Units)
Call Add_Str_Length( ParamName )


End ! GetParamName

Subroutine GetParamUnit( iMod , iParam, Units )
!
! Return the units of the different parameters of the different models
!
Implicit None
Integer  iMod, iParam
Character (Len=255) ParamName, Units

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamUnit

Call GetParamAndUnit(iMod,iParam,ParamName,Units)
Call Add_Str_Length( Units )

End ! GetParamUnit

Subroutine GetStateVarCount( iMod , nVar )
!
! Return the number of state variables of the different models
!
!DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarCount
implicit none
integer, intent(in) :: iMod
integer, intent(out) :: nVar

Select Case (iMod)
Case (1)
nVar = 31
Case Default
nVar = 0
End Select

Return
End


Subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
!
! Return the name and unit of the different state variables of the different models
!
implicit none
integer, intent(in) :: iMod , iVar
Character (Len=255), intent(out) :: Name, Unit

Select Case (iMod)
Case (1)
    Select Case (iVar)
        Case (1)
            Name = 'alpha_11#'         ; Unit = '-'
        Case (2)
            Name = 'alpha_22#'         ; Unit = '-'
        Case (3)
            Name = 'alpha_33#'         ; Unit = '-'
        Case (4)
            Name = 'alpha_12#'         ; Unit = '-'
        Case (5)
            Name = 'alpha_23#'         ; Unit = '-'
        Case (6)
            Name = 'alpha_31#'         ; Unit = '-'
        Case (7)
            Name = 'alpha^#M_11#'       ; Unit = '-'
        Case (8)
            Name = 'alpha^#M_22#'         ; Unit = '-'
        Case (9)
            Name = 'alpha^#M_33#'         ; Unit = '-'
        Case (10)
            Name = 'alpha^#M_12#'         ; Unit = '-'
        Case (11)
            Name = 'alpha^#M_23#'         ; Unit = '-'
        Case (12)
            Name = 'alpha^#M_31#'         ; Unit = '-'
        Case (13)
            Name = 'MM_+#'         ; Unit = '-'
        Case (14)
            Name = 'MM_-#'         ; Unit = '-'
        Case (15)
            Name = 'alpha^#in_11#'         ; Unit = '-'
        Case (16)
            Name = 'alpha^#in_22#'         ; Unit = '-'
        Case (17)
            Name = 'alpha^#in_33#'         ; Unit = '-'
        Case (18)
            Name = 'alpha^#in_12#'         ; Unit = '-'
        Case (19)
            Name = 'alpha^#in_23#'         ; Unit = '-'
        Case (20)
            Name = 'alpha^#in_31#'         ; Unit = '-'
        Case (21)
            Name = 'void ratio'            ; Unit = '-'
        Case (22)
            Name = 'reinitialize'          ; Unit = '-'
        Case (23)
            Name = 'unused'           ; Unit = '-'
        Case (24)
            Name = 'unused'              ; Unit = '-'
        Case (25)
            Name = 'M^b'                    ; Unit = '-'
        Case (26)
            Name = 'M^d'                    ; Unit = '-'
        Case (27)
            Name = 'cos3@q'             ; Unit = '-'
        Case (28)
            Name = '@y'                   ; Unit = '-'
        Case (29)
            Name = 'q'                     ; Unit = '-'
        Case (30)
            Name = 'p'                     ; Unit = '-'
        Case (31)
            Name = 'm^M'                    ; Unit = '-'
        Case Default
            Name='N/A'              ; Unit = '?'
    End Select
Case Default
Name='N/A'                  ; Unit = '?'
End Select

Return
End
    
Subroutine GetStateVarName( iMod , iVar, Name )
!
! Return the name of the different state variables
! of the different models
!
Implicit None
Integer  iMod, iVar
Character (Len=255) Name, Unit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetStateVarName

Call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
Call Add_Str_Length( Name )

End ! GetStateVarName

Subroutine GetStateVarUnit( iMod , iVar, Unit )
!
! Return the units of the different state variables of the different models
!
Implicit None
Integer  iMod, iVar
Character (Len=255) Name, Unit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetStateVarUnit

Call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
Call Add_Str_Length( Unit )

End ! GetStateVarUnit
    
!----------------------------------------------------------------------------------    
module SANISANDMS_constants
!
double precision, parameter :: one3 = 1.0d0/3.0d0
double precision, parameter :: one6 = 1.0d0/6.0d0
double precision, parameter :: two3 = 2.0d0/3.0d0
double precision, parameter :: three2 = 3.0D0 / 2.0D0
double precision, parameter :: root23   = DSQRT(2.0D0/3.0D0)
double precision, parameter :: root32   = DSQRT(3.0D0/2.0D0)
double precision, parameter :: zeroP5 = 0.5D0
double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0d0
double precision, parameter :: two = 2.0d0
double precision, parameter :: three = 3.0D0
!double precision, parameter :: fou = 4.0d0
!double precision, parameter :: fiv = 5.0d0
!double precision, parameter :: six = 6.0d0
!double precision, parameter :: nin = 9.0d0
!double precision, parameter :: ts  = 27.0d0
!
!double precision, parameter :: tot = 1.5d0
!
!double precision, parameter :: pi = 3.14159265358979d0
!
!double precision, parameter :: tol_for_thrsh = 5.0E-02
!double precision, parameter :: tol_for_small = 1.0E-08
!
! unused in state dependent subroutine
double precision, parameter :: onePM5= 1.0D-5
double precision, parameter :: twoeight = 28.0D0
double precision, parameter :: twoeight625 = 28.0D0/625.0D0
double precision, parameter :: onetwofive = 125.0D0
double precision, parameter :: onetwofive625 = 125.0D0/625.0D0
double precision, parameter :: four6 = 4.0D0/6.0D0
double precision, parameter :: fivefoursix = 546.0D0
double precision, parameter :: fivefoursix625 = 546.0D0/625.0D0
double precision, parameter :: fivefour = 54.0D0
double precision, parameter :: fivefour625 = 54.0D0/625.0D0
double precision, parameter :: threeseveneight = 378.0D0
double precision, parameter :: threeseveneight625 = 378.0D0/625.0D0
double precision, parameter :: sixtwofive = 625.0D0
double precision, parameter :: seven27 = 7.0D0/27.0D0
double precision, parameter :: ten27 = 10.0D0/27.0D0
double precision, parameter :: one27 = 1.0D0/27.0D0
double precision, parameter :: onefour336 = 14.0D0/336.0D0
double precision, parameter :: threefive336 = 35.0D0/336.0D0
double precision, parameter :: onesixtwo336 = 162.0D0/336.0D0
double precision, parameter :: onetwofive336 = 125.0D0/336.0D0
double precision, parameter :: fourtwo336 = 42.0D0/336.0D0
double precision, parameter :: twotwofour336 = 224.0D0/336.0D0
double precision, parameter :: twoone336 = 21.0D0/336.0D0

end module SANISANDMS_constants 
    
Module Parameters_position_in_params_iMod1
!! Purpose: Associate a position of the parameter in the array params
!!
Implicit None
!
        integer, parameter :: i_PAR_G0        = 1
        integer, parameter :: i_PAR_nu        = 2
        integer, parameter :: i_PAR_Mc        = 3
        integer, parameter :: i_PAR_c         = 4
        integer, parameter :: i_PAR_lambda_c  = 5
        integer, parameter :: i_PAR_e0        = 6
        integer, parameter :: i_PAR_ksi       = 7
        integer, parameter :: i_PAR_m         = 8
        integer, parameter :: i_PAR_h0        = 9
        integer, parameter :: i_PAR_ch        = 10
        integer, parameter :: i_PAR_nb        = 11
        integer, parameter :: i_PAR_A0        = 12
        integer, parameter :: i_PAR_nd        = 13
        integer, parameter :: i_PAR_mu0       = 14
        integer, parameter :: i_PAR_zeta      = 15
        integer, parameter :: i_PAR_beta      = 16
        integer, parameter :: i_TolF          = 17   !Yield surface Tolerance
        integer, parameter :: i_TolR          = 18   !RKF23 Error Tolerance
        integer, parameter :: i_PAR_e_init    = 19
        integer, parameter :: i_emax          = 20
        integer, parameter :: i_emin          = 21
        integer, parameter :: i_Mbc_max       = 22   !limit to the opening of the BS in TXCompression
        integer, parameter :: i_El_Scheme     = 23   !integration scheme for the the elastic trial: [0] = Single step Forward Euler, [1] = Automatic sub-step RKF23 for trial only, [2] = Automatic sub-step RKF23 for trial and intersection stress 
        integer, parameter :: i_hM_exc_Eq     = 24
        integer, parameter :: i_Custom_h_max  = 25
    end module Parameters_position_in_params_iMod1

Module State_Var_position_in_StVar_iMod1
!! Purpose: Associate a position of the state variable in the StVar and StVar0 arrays
!!
implicit none
        integer, parameter :: i_Alpha_xx      = 1
        integer, parameter :: i_Alpha_yy      = 2
        integer, parameter :: i_Alpha_zz      = 3
        integer, parameter :: i_Alpha_xy      = 4
        integer, parameter :: i_Alpha_yz      = 5
        integer, parameter :: i_Alpha_zx      = 6
        integer, parameter :: i_AlphaM_xx     = 7
        integer, parameter :: i_AlphaM_yy     = 8
        integer, parameter :: i_AlphaM_zz     = 9
        integer, parameter :: i_AlphaM_xy     = 10
        integer, parameter :: i_AlphaM_yz     = 11
        integer, parameter :: i_AlphaM_zx     = 12
        integer, parameter :: i_sMM_plus      = 13
        integer, parameter :: i_sMM_minus     = 14
        integer, parameter :: i_alpha_in_xx   = 15
        integer, parameter :: i_alpha_in_yy   = 16
        integer, parameter :: i_alpha_in_zz   = 17
        integer, parameter :: i_alpha_in_xy   = 18
        integer, parameter :: i_alpha_in_yz   = 19
        integer, parameter :: i_alpha_in_zx   = 20
        integer, parameter :: i_Void          = 21
        integer, parameter :: i_resetState    = 22 
        integer, parameter :: i_PLNoTension   = 23! deprecated output. plastic strain integrated without accuracy
        integer, parameter :: i_PsPsmall      = 24! deprecated output. indicator of inaccurate solution for low stress levels
        integer, parameter :: i_Mb            = 25
        integer, parameter :: i_Md            = 26
        integer, parameter :: i_cos3theta     = 27
        integer, parameter :: i_psi           = 28
        integer, parameter :: i_q_Norm_New    = 29
        integer, parameter :: i_p             = 30
        integer, parameter :: i_MM            = 31

!    
end  module State_Var_position_in_StVar_iMod1
!---------------------------------------------------------------------------------- 

! ======== TensorOperation.f90 ========
! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
Subroutine GetTrace(v, res)
!
!***********************************************************************
!
!     Function: To return the trace of a 6 vector
!
!***********************************************************************
!
Implicit none

double Precision, intent(in) :: v(6)
double precision, intent(out) :: res

double precision, parameter :: zero  = 0.0D0

res = zero
res = v(1)+v(2)+v(3)

Return
End

Subroutine GetDevPart(v, devP)
!
!***********************************************************************
!
!     Function: To return the deviatoric part of a vector
!
!***********************************************************************
!
implicit none

double precision, intent(in) :: v(6)
double precision, intent(out) :: devP(6)
double precision, parameter :: one3 = 1.0d0/3.0d0

double precision res

call GetTrace(v, res)
res = res * one3
!     call COPYRVEC(v, devP, 6)
devP(1) = v(1) - res
devP(2) = v(2) - res
devP(3) = v(3) - res
devP(4) = v(4)
devP(5) = v(5)
devP(6) = v(6)

Return
End

Subroutine SingleDot(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes v1.v2, v1 and v2 should be both in their "contravariant" form
!
!***********************************************************************
!
implicit none

double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res(6)

double precision, parameter :: one2  = 0.5D0

res(1) = v1(1)*v2(1) + v1(4)*v2(4) + v1(6)*v2(6)
res(2) = v1(4)*v2(4) + v1(2)*v2(2) + v1(5)*v2(5)
res(3) = v1(6)*v2(6) + v1(5)*v2(5) + v1(3)*v2(3)
res(4) = one2*(v1(1)*v2(4) + v1(4)*v2(1) + v1(4)*v2(2) + v1(2)*v2(4) + v1(6)*v2(5) + v1(5)*v2(6))
res(5) = one2*(v1(4)*v2(6) + v1(6)*v2(4) + v1(2)*v2(5) + v1(5)*v2(2) + v1(5)*v2(3) + v1(3)*v2(5))
res(6) = one2*(v1(1)*v2(6) + v1(6)*v2(1) + v1(4)*v2(5) + v1(5)*v2(4) + v1(6)*v2(3) + v1(3)*v2(6))

Return
End

Subroutine SingleDotCube(v1, res)
!
!***********************************************************************
!
!     Function: computes v1.v2, v1 and v2 should be both in their "contravariant" form
!
!***********************************************************************
!
implicit none

double precision, intent(in) :: v1(6)
double precision, intent(out) :: res(6)

double precision, parameter :: one2  = 0.5D0
double precision, parameter :: two  = 2.0D0

res(1) = v1(1)*v1(1)*v1(1) + two*v1(1)*v1(4)*v1(4) + two*v1(1)*v1(6)*v1(6) + v1(2)*v1(4)*v1(4) + two*v1(5)*v1(4)*v1(6) + v1(3)*v1(6)*v1(6)
res(2) = v1(2)*v1(2)*v1(2) + two*v1(2)*v1(4)*v1(4) + two*v1(2)*v1(5)*v1(5) + v1(1)*v1(4)*v1(4) + two*v1(6)*v1(4)*v1(5) + v1(3)*v1(5)*v1(5)
res(3) = v1(3)*v1(3)*v1(3) + two*v1(3)*v1(5)*v1(5) + two*v1(3)*v1(6)*v1(6) + v1(2)*v1(5)*v1(5) + two*v1(4)*v1(5)*v1(6) + v1(1)*v1(6)*v1(6)
res(4) = v1(2)*(v1(1)*v1(4) + v1(2)*v1(4) + v1(5)*v1(6)) + v1(5)*(v1(1)*v1(6) + v1(3)*v1(6) + v1(4)*v1(5)) + v1(4)*(v1(1)*v1(1) + v1(4)*v1(4) + v1(6)*v1(6))
res(5) = v1(3)*(v1(2)*v1(5) + v1(3)*v1(5) + v1(4)*v1(6)) + v1(6)*(v1(1)*v1(4) + v1(2)*v1(4) + v1(5)*v1(6)) + v1(5)*(v1(2)*v1(2) + v1(4)*v1(4) + v1(5)*v1(5))
res(6) = v1(5)*(v1(1)*v1(4) + v1(2)*v1(4) + v1(5)*v1(6)) + v1(3)*(v1(1)*v1(6) + v1(3)*v1(6) + v1(4)*v1(5)) + v1(6)*(v1(1)*v1(1) + v1(4)*v1(4) + v1(6)*v1(6))

Return
End

Subroutine DoubleDot2_2_Contr(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes doubledot product for vector-vector arguments, both "contravariant"
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res

double precision, parameter :: zero = 0.0D0
double precision, parameter :: two = 2.0D0

integer i

res = zero

do i = 1, 6
    if (i > 3) then
        res = res + two * v1(i) * v2(i)
    else
        res = res + v1(i) * v2(i)
    end if
end do

Return
End

Subroutine DoubleDot2_2_Cov(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes doubledot product for vector-vector arguments, both "covariant"
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res

double precision, parameter :: zero = 0.0D0
double precision, parameter :: one2 = 0.5D0

integer i

res = zero

do i = 1, 6
    if (i > 3) then
        res = res + one2 * v1(i) * v2(i)
    else
        res = res + v1(i) * v2(i)
    end if
end do

Return
End

Subroutine DoubleDot2_2_Mixed(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res

double precision, parameter :: zero = 0.0D0

integer i

res = zero

do i = 1, 6
    res = res + v1(i)*v2(i)
end do

Return
End

Subroutine GetNorm_Contr(v, res)
!
!***********************************************************************
!
!     Function: computes contravariant (stress-like) norm of input 6x1 tensor
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v(6)
double precision, intent(out) :: res

call DoubleDot2_2_Contr(v, v, res)
res = Dsqrt(res)

Return
End

Subroutine GetNorm_Cov(v, res)
!
!***********************************************************************
!
!     Function: computes covariant (strain-like) norm of input 6x1 tensor
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v(6)
double precision, intent(out) :: res

call DoubleDot2_2_Cov(v, v, res)
res = Dsqrt(res)

Return
End

subroutine Dyadic2_2(v1, v2, n, res)
!
!***********************************************************************
!
!     Function: computes dyadic product for two vector-storage arguments
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
integer, intent(in) :: n
double precision, intent(out) :: res(6,6)

integer i, j

Do i=1,6
    Do j=1,6
        res(i,j) = v1(i)*v2(j)
    end do
end do

Return
End

subroutine ToContraviant(v1, v2)
!
!***********************************************************************
!
!     Function: Transform to the 'contravariant' form
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6)
double precision, intent(out) :: v2(6)

double precision, parameter :: one2 = 0.5D0

integer i

Do i = 1,6
    if (i > 3) then
        v2(i) = one2 * v1(i)
    else
        v2(i) = v1(i)
    end if
end do

Return
End

subroutine ToCovariant(v1, v2)
!
!***********************************************************************
!
!     Function: Transform to the 'covariant' form
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6)
double precision, intent(out) :: v2(6)

double precision, parameter :: two = 2.0D0

integer i

Do i = 1,6
    if (i > 3) then
        v2(i) = two * v1(i)
    else
        v2(i) = v1(i)
    end if
end do

Return
End

! ======== StatedepSub.f90 ========
! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
SUBROUTINE SignFunction(v1, res)

implicit none
double precision, intent(in) :: v1
double precision, intent(out) :: res

double precision, parameter :: tol = 1.0D-15
double precision, parameter :: zero_plus_tol = 1.0D-30
double precision, parameter :: one = 1.0D0
double precision, parameter :: mone = -1.0D0

res = zero_plus_tol
IF (v1 > zero_plus_tol) THEN
    res = one
ELSEIF (v1 < -zero_plus_tol) THEN
    res = mone
ELSE
    res = zero_plus_tol
END IF

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE Macauley(v1, res)

implicit none
double precision, intent(in) :: v1
double precision, intent(out) :: res

double precision, parameter :: tol = 1.0D-15
double precision, parameter :: zero = 0.0D0

res = zero
IF (v1 > zero) THEN
    res = v1
ELSE
    res = zero
END IF

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE MacauleyIndex(v1, res)

implicit none
double precision, intent(in) :: v1
double precision, intent(out) :: res

double precision, parameter :: tol = 1.0D-15
double precision, parameter :: zero = 1.0D-30
double precision, parameter :: one = 1.0D0

res = zero
IF (v1 > zero) THEN
    res = one
ELSE
    res = zero
END IF

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE MacauleyIndex1(v1, res)

implicit none
double precision, intent(in) :: v1
double precision, intent(out) :: res

double precision, parameter :: tol = 1.0D-15
double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0D0

res = zero
IF (v1 > tol) THEN
    res = one
ELSE
    res = zero
END IF

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE g(cos3theta, c, gth, dgdth)
!
! Purpose: compute the Van Ekelen or Argyris function and their derivative 
!
! Arguments:
!                           I/O   Type
!  cos3theta                 I    R    :
!  c                         I    R    : model paramter = Me/Mc 
!  gth                       O    R    : radius of the Function for the current Lode angle
!  dgdth                     O    R    : derivative of the Function with respect to Teta
!
! Argyris function
!  gth = two * c / ((one + c) - (one - c) * cos3theta)
!  dgdth = 
!
! Van Eekelen function
!  gth = alpha / ((one + beta * cos3theta)^-0.25
!  dgth = 1/g dgd(cos3teta) = -0.25 alpha (one + beta * cos3theta)^-1
      
implicit none
double precision, intent(in) :: cos3theta, c
double precision, intent(out) ::  gth, dgdth

double precision tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, alpha, beta, n_VEm1

double precision, parameter :: one = 1.0D0
double precision, parameter :: two = 2.0D0
double precision, parameter :: n_VE = -0.25d0

!
n_VEm1=one/n_VE

tmp1=one/(two**n_VE)
tmp2=c**n_VEm1
tmp3=one+tmp2
tmp4=one-tmp2

! alpha = 2^0.25/(1+c^-4)^-0.25
alpha=tmp1*(tmp3**n_VE)
! beta = (1-c^-4)/(1+c^-4)
beta=tmp4/tmp3

tmp5=(one+beta*cos3theta)**n_VE
tmp6=one+beta*cos3theta

gth = alpha*tmp5
dgdth=n_VE*beta/tmp6

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE GetF(stress, alpha, Params, res)
!     Yield function
use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: stress(6), alpha(6), Params(25)
double precision, intent(out) :: res

double precision s(6), PAR_m, p

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: two3 = 2.0D0/3.0D0
double precision, parameter :: root23 = DSQRT(two3)

PAR_m = Params(i_PAR_m)

CALL GetTrace(stress, p)
p = one3 * p
CALL GetDevPart(stress, s)
s = s - p * alpha

CALL GetNorm_Contr(s, res)
res = res - root23 * PAR_m * p

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE GetPSI(e, p, Params, res)
!
! Purpose: Compute the state parametr
!          Consider emax and emin
! Arguments:
!            I/O  Type
!  e          I   R    : current void ratio
!  p          I   R    :
!  Params     I   R()  :
!  res        O   R    : state parameter

use Parameters_position_in_params_iMod1
use extrainfo, only: onekPa
use extrainfo, only: oneAtm

implicit none
double precision, intent(in) :: e, p, Params(25)
double precision, intent(out) :: res

double precision PAR_lambda_c, PAR_e0, PAR_ksi, pp

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: oneh = 100.0D0
double precision, parameter :: tol = 1.0D-6
double precision PAR_Pmin, PAR_P_atm, PAR_c, PAR_nb, Mbe, Mbc, psi_maximum, PAR_m, PAR_Mc

PAR_Pmin     = onekPa/oneh
PAR_P_atm    = oneAtm
PAR_Mc       = Params(i_PAR_Mc)
PAR_c        = Params(i_PAR_c)
PAR_lambda_c = Params(i_PAR_lambda_c)
PAR_e0       = Params(i_PAR_e0)
PAR_ksi      = Params(i_PAR_ksi)
PAR_m        = Params(i_PAR_m)
PAR_nb       = Params(i_PAR_nb)

! 25-11-2021 change PAR_Pmin to 0
IF (p > 0.0D0) THEN
    res = e - (PAR_e0 - PAR_lambda_c * ((p / PAR_P_atm) ** PAR_ksi))
ELSE
    !res = e - (PAR_e0 - PAR_lambda_c * ((PAR_Pmin / PAR_P_atm) ** PAR_ksi))
    res = e - PAR_e0
END IF

! maximum of bounding opening
Mbe = 1.5D0 - PAR_m
Mbc = Mbe/(PAR_c * PAR_Mc)
psi_maximum = -dlog(Mbc)/PAR_nb*0.90D0

!IF (res < psi_maximum) THEN
!    res = psi_maximum
!END IF

RETURN
END

!------------------------------------------------------
SUBROUTINE GetLodeAngle(sn, Cos3Theta)
!------------------------------------------------------
! Purpose: Subroutine to calculate lode angle
!------------------------------------------------------
implicit none
double precision, intent(in   )   :: sn(6)
double precision, intent(  out)   :: Cos3Theta

double precision, parameter :: tol = 1.0D-15
double precision, parameter :: rootsix = DSQRT(6.0D0)
double precision, parameter :: one = 1.0D0
double precision res(6), res1(6)

CALL SingleDot(sn, sn, res)
CALL SingleDot(sn, res, res1)
CALL GetTrace(res1, Cos3Theta)
Cos3Theta = rootsix * Cos3Theta

IF (Cos3Theta > one + tol) THEN !tolerance
    !write(100,*), 'problematic cos', Cos3Theta
!    write(100,*), 'sn', sn
    Cos3Theta = one
END IF

IF (Cos3Theta < -one - tol) THEN
        !write(100,*), 'problematic cos2', Cos3Theta
!        write(100,*), 'sn', sn
    Cos3Theta = -one
END IF

RETURN
END
    
!------------------------------------------------------
SUBROUTINE GetElasticModuli(stress, en, Params, E_K, E_G)
!------------------------------------------------------
! Purpose: Calculates G, K
!------------------------------------------------------
! Considering abs(pn) if pn<=0
Use Parameters_position_in_params_iMod1
use extrainfo, only: onekPa
use extrainfo, only: oneAtm

implicit none
double precision, intent(in   ) :: stress(6), Params(25), en
double precision, intent(  out) :: E_K, E_G

double precision PAR_G0, PAR_nu, pn
double precision PAR_P_atm
double precision PAR_Pmin

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: twonineseven = 2.97D0
double precision, parameter :: one = 1.0D0
double precision, parameter :: two = 2.0D0
double precision, parameter :: two3 = 2.0D0 / 3.0D0
double precision, parameter :: oneh = 100.0D0

PAR_Pmin = onekPa / oneh
PAR_P_atm = oneAtm

PAR_G0 = Params(i_PAR_G0)
PAR_nu = Params(i_PAR_nu)

CALL GetTrace(stress, pn)
pn = one3 * pn

IF (pn <= 0.0D0) THEN
    pn = DABS(pn)
END IF

E_G = PAR_G0 * PAR_P_atm * ((twonineseven - en)**2) / (one + en) * DSQRT(pn / PAR_P_atm)
E_K = two3 * (one + PAR_nu) / (one - two * PAR_nu) * E_G

RETURN
END

!------------------------------------------------------
SUBROUTINE GetElasticModuli1(stress, en, Params, E_K, E_G)
! Purpose: Calculate G, K
!------------------------------------------------------
! Consider 0 if pn <= PAR_Pmin

use Parameters_position_in_params_iMod1
use extrainfo, only: onekPa
use extrainfo, only: oneAtm

implicit none
double precision, intent(in   ) :: stress(6), Params(25), en
double precision, intent(  out) :: E_K, E_G

double precision PAR_G0, PAR_nu, pn
double precision PAR_P_atm
double precision PAR_Pmin

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: twonineseven = 2.97D0
double precision, parameter :: one = 1.0D0
double precision, parameter :: two = 2.0D0
double precision, parameter :: two3 = 2.0D0 / 3.0D0
double precision, parameter :: oneh = 100.0D0

PAR_Pmin = onekPa / oneh
PAR_P_atm = oneAtm

PAR_G0 = Params(i_PAR_G0)
PAR_nu = Params(i_PAR_nu)

CALL GetTrace(stress, pn)
pn = one3 * pn

IF (pn <= 0.0D0) THEN
    !WRITE(102, *) 'pn <= 0.0D0'
    pn = 0.0D0
END IF

E_G = PAR_G0 * PAR_P_atm * ((twonineseven - en)**2) / (one + en) * DSQRT(DABS(pn) / PAR_P_atm)
E_K = two3 * (one + PAR_nu) / (one - two * PAR_nu) * E_G

RETURN
END

!------------------------------------------------------
SUBROUTINE GetStiffness(E_K, E_G, C)
!------------------------------------------------------
implicit none
double precision, intent(in   ) :: E_K, E_G
double precision, intent(  out) :: C(6,6)

double precision a, b

double precision, parameter :: four3 = 4.0D0 / 3.0D0
double precision, parameter :: two3 = 2.0D0 / 3.0D0
double precision, parameter :: zero = 0.0D0

!     CALL MZEROR(C, 36)
C = zero
a = E_K + four3 * E_G
b = E_K - two3 * E_G

C(1, 1) = a
C(2, 2) = a
C(3, 3) = a
C(4, 4) = E_G
C(5, 5) = E_G
C(6, 6) = E_G
C(1, 2) = b
C(1, 3) = b
C(2, 3) = b
C(2, 1) = b
C(3, 1) = b
C(3, 2) = b

RETURN
END

!------------------------------------------------------
SUBROUTINE GetCompliance(E_K, E_G, D)
!------------------------------------------------------
implicit none
      
double precision, intent(in) :: E_K, E_G
double precision, intent(out) :: D(6,6)

double precision a, b, c

double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0D0
double precision, parameter :: one9 = 1.0D0 / 9.0D0
double precision, parameter :: one3 = 1.0D0 / 3.0D0
double precision, parameter :: one6 = 1.0D0 / 3.0D0

D = zero
a = one9 * E_K + one3 * E_G
b = one9 * E_K - one6 * E_G
c = one / E_G

D(1, 1) = a
D(2, 2) = a
D(3, 3) = a
D(4, 4) = c
D(5, 5) = c
D(6, 6) = c
D(1, 2) = b
D(1, 3) = b
D(2, 3) = b
D(2, 1) = b
D(3, 1) = b
D(3, 2) = b

RETURN
END

!----------------------------------------------------------------------
SUBROUTINE GetElastoPlasticTangent(Stress1, DGamma1, E_G, E_K, B, C,  &
           D, h, sn, dd, bb, aCep, temp4, HAlpha)
!----------------------------------------------------------------------
! returns the stiffness matrix in its contravarinat-contravariant form
! Consider temp is smaller than small and p is smaller than PAR_Pmin

use extrainfo, only: onekPa
use extrainfo, only: oneAtm

implicit none
double precision, intent(in) :: Stress1(6), DGamma1, E_G, E_K, B, C,  &
    D, sn(6), dd(6), bb(6), temp4
double precision, intent(in   ) :: h
double precision, intent(  out) :: aCep(6, 6), HAlpha(6, 6)

double precision aC(6,6), rr(6), temp1(6), temp2(6), R(6), SI1(6),  &
    res(6), res1(6,6), p, PM, temp3, res2, dfdsig(1, 6), res7(1, 6),  &
    bb1(6)

double precision, parameter :: small  = 1.0D-10
double precision, parameter :: zero  = 0.0D0
double precision, parameter :: one  = 1.0D0
double precision, parameter :: one3  = 1.0D0/3.0D0
double precision, parameter :: oneh  = 100.0D0
double precision, parameter :: two3  = 2.0D0/3.0D0

integer i

double precision PAR_P_atm
double precision PAR_Pmin

PAR_Pmin = onekPa / oneh
PAR_P_atm = oneAtm

SI1(1) = one
SI1(2) = one
SI1(3) = one
SI1(4) = zero
SI1(5) = zero
SI1(6) = zero

CALL GetTrace(Stress1, p)
p = one3 * p
!     IF (p <= PAR_Pmin) THEN
!         p = PAR_Pmin
!     END IF
CALL GetDevPart(Stress1, rr)
IF (p > PAR_Pmin) THEN
    rr = rr / p
ELSE
    rr = zero
END IF

CALL DoubleDot2_2_Contr(bb, sn, PM)
PM = two3 * p * h * PM
bb1 = two3 * h * bb

CALL GetStiffness(E_K, E_G, aC)

CALL SingleDot(sn, sn, temp1)
!     dg/dsig
res = (B * sn) - (C * (temp1 - one3 * SI1)) + (one3 * D * SI1)
CALL ToCovariant(res, R)

!     De * (dg/dsig)
temp1 =  MATMUL(aC, R)

!     df/dsig
CALL DoubleDot2_2_Contr(sn, rr, res2)
res = sn - one3 * res2 * SI1
! 11-10-2021
!write(100,*), 'L = df/dsig = ', res
Do i = 1, 6
    dfdsig(1, i) = res(i)
end do

CALL ToCovariant(res, temp2)
!     (df/dsig) * De
temp2 = MATMUL(temp2, aC)

!     (df/dsig) * De * (dg/dsig) + KP
CALL DoubleDot2_2_Contr(temp2, R, temp3)
temp3 = temp3 + PM

!     DP
IF ((temp3 <= small) .OR. (p <= PAR_Pmin)) THEN
    aCep = aC
ELSE
    CALL Dyadic2_2(temp1, temp2, 6, res1)
    CALL MacauleyIndex(DGamma1, res2)
    aCep = aC - (res2 / temp3 * res1)
END IF

!11-10-2021
! dfdsig * De
res7 = matmul(dfdsig, aC)
! dfdsig * De  /temp4
res7 = res7 / temp4
! dfdsig * De  /temp4 * bb1 (become a 6*6 matrix)
!HAlpha = bb1* res7
CALL Dyadic2_2(bb1, res7, 6, HAlpha)
!write(100, *), 'HH(alpha) :'
!write(100, '(6F14.6)'), HAlpha(1:6, 1:6)
RETURN
END
!
!-----------------------------------------------------------------------------------
SUBROUTINE GetNormalToYield(stress, alpha, iStep, sn, iAbort)
!------------------------------------------------------
! ! Purpose: compute the normal to the yield in the deviatoric plane. It aborts if the noorm is too small
implicit none

integer, intent(in   ) ::  iStep
integer, intent(  out) :: iAbort
double precision, intent(in  ) ::  stress(6), alpha(6)
double precision, intent(  out) :: sn(6)

integer i, izero
double precision devStress(6), p, Snorm, mm, stress_use(6), TRACE, Snorm_check, s_palpha(6)

!double precision, parameter :: small = 1.0D-7
double precision, parameter :: small = 1.0D-12
double precision, parameter :: PAR_P_atm = 100.0d0
double precision, parameter :: one3 = 1.0D0/3.0D0
!------------------------------------------------------------remove
double precision, parameter :: PAR_Pmin = 0.01D0
double precision, parameter :: tol = 1.0D-12
double precision, parameter :: one  = 1.0D0
double precision, parameter :: zero = 0.0D0
!------------------------------------------------------------remove

iAbort = 0
stress_use = stress
!------------------------------------------------------------remove
mm = 0.01D0*DSQRT(2.0D0/3.0D0)
!------------------------------------------------------------remove
CALL GetDevPart(stress_use, devStress)
CALL GetTrace(stress_use, p)
p = one3 * p
s_palpha = devStress/p - alpha
CALL GetNorm_Contr(s_palpha, Snorm)
!------------------------------------------------------------remove
Snorm_check = DSQRT( (s_palpha(1)*s_palpha(1)) + (s_palpha(2)*s_palpha(2)) + (s_palpha(3)*s_palpha(3)))
TRACE = s_palpha(1) + s_palpha(2) + s_palpha(3)
!------------------------------------------------------------remove

IF (Snorm <= small) THEN
    write(100, *) 'set iAbort = 1 in GetNormalToYield'
    write(100, *) 'istep', 'Snorm'
    write(100, *) istep, Snorm
    write(100, *) 'stress', 'alpha'
    write(100, '(6F10.6)') stress
    write(100, '(6F10.6)') alpha
    call flush(100)
    iAbort = 1
    !iAbort = 18
    return
    !stop
END IF

!IF ((p <= PAR_Pmin) .OR. (izero == 1)) THEN
!IF (Snorm <= PAR_Pmin) THEN
!    sn = zero
!    WRITE(100,*),'SN SET TO ZERO'
!ELSE
sn = s_palpha/Snorm
!END IF

RETURN
END
!
!---------------------------------------------------------------------------------------------------------
SUBROUTINE GetNormalToYield2(stress, alpha, iStep, sn)
! Purpose: compute the normal to the yield in the deviatoric plane. ! Not a unit vector if snorm is small
!---------------------------------------------------------------------------------------------------------
implicit none
integer, intent(in) ::  iStep
double precision, intent(in) ::  stress(6), alpha(6)
double precision, intent(out) :: sn(6)

double precision devStress(6), p, Snorm, mm, stress_use(6), TRACE, Snorm_check, s_palpha(6)
integer i, izero

double precision, parameter :: small = 1.0D-7
double precision, parameter :: PAR_P_atm = 100.0d0
double precision, parameter :: PAR_Pmin = 0.01D0
double precision, parameter :: tol = 1.0D-12
double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: one  = 1.0D0
double precision, parameter :: zero = 0.0D0

stress_use = stress
mm = 0.01D0*DSQRT(2.0D0/3.0D0)
CALL GetDevPart(stress_use, devStress)
CALL GetTrace(stress_use, p)
p = one3 * p

s_palpha = devStress/p - alpha

CALL GetNorm_Contr(s_palpha, Snorm)
Snorm_check = DSQRT( (s_palpha(1)*s_palpha(1)) + (s_palpha(2)*s_palpha(2)) + (s_palpha(3)*s_palpha(3)))
TRACE = s_palpha(1) + s_palpha(2) + s_palpha(3)

IF (Snorm <= small) THEN
    Snorm = one
END IF

sn = s_palpha/Snorm

RETURN
END

!-----------------------------------------------------------------------------------
SUBROUTINE IntersectionFactor(CurStress, CurVoid, strainInc,    &
           CurAlpha, a00, a11, Params, a, iAbort)   
! Purpose: computing the intersection with the yield surface
!-----------------------------------------------------------------------------------
! Arguments:
!               I/O  Type
!  CurStress    I    R()  : 
!  CurVoid      I    R    :
!  strainInc    I    R()  : 
!  CurAlpha     I    R()  : 
!  a00          I    R    : lower bound of the strain increment interval
!  a11          I    R    : upper bound of the strain increment interval
!  Params       I    R()  : 
!  a            O    R    : Intersection with the yield surface
!  iAbort       O    I    : X-Position of integration point
!
Use Parameters_position_in_params_iMod1
implicit none

integer, intent(out) :: iAbort
double precision, intent(in) :: CurStress(6), CurVoid, strainInc(6),    &
    CurAlpha(6), a00, a11, Params(25)
double precision, intent(out) :: a

integer i
double precision dSigma(6), dSigma0(6), dSigma1(6), C(6,6), a0, a1,  &
    TolF, f0, f1, f, E_K,E_G

double precision, parameter :: small = 1.0D-7
double precision, parameter :: tol = 1.0D-15
double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0D0

iAbort = 0

a0 = a00
a1 = a11
a = a0
TolF = Params(i_TolF)

! 24-11-2021 change
CALL GetElasticModuli(CurStress, CurVoid,Params, E_K, E_G)
CALL GetStiffness(E_K, E_G, C)

dSigma0 = a0 * MATMUL(C, strainInc)
CALL GetF(CurStress + dSigma0, CurAlpha, Params, f0)

dSigma1 = a1 * MATMUL(C, strainInc)
CALL GetF(CurStress + dSigma1, CurAlpha, Params, f1)

! change maximum iteration number from 20 to 10. 
! I encounter the problem 'fail to find a intersection factor' one time. Change to 20 then works.
DO i = 1, 30
    a = a1 - f1 * (a1 - a0) / (f1 - f0)
    dSigma = a * MATMUL(C, strainInc)
    CALL GetF(CurStress + dSigma, CurAlpha, Params, f)
    IF (DABS(f) < TolF) THEN
        EXIT
    END IF
    IF (f * f0 < -tol) THEN
        a1 = a
        f1 = f
    ELSE
        f1 = f1 * f0 / (f0 + f)
        a0 = a
        f0 = f
    END IF
    IF (i == 30) THEN
        WRITE(100,*) 'fail to find a intersection factor'
        call flush(100)
        iAbort = 1
        return
        !EXIT
    END IF
END DO
IF (a > one) THEN 
    a = one
END IF
IF (a < zero) THEN
    a = zero
END IF

RETURN
    END

!------------------------------------------------------------------------------------
SUBROUTINE IntersectionFactor_Unloading(CurStress, CurVoid, &
           strainInc, CurAlpha, Params, a, iAbort)
! Purpose: It provides the intersection with the yield surfaces in case of unloading
!           this subroutine calls IntersectionFactor
!-------------------------------------------------------------------------------------
! Arguments:
!             I/O  Type
!  CurStress   I   R()    : 
!  CurVoid     I   R      : 
!  strainInc   I   R()    :
!  CurAlpha    I   R()    : 
!  Params      I   R()    :
!  a           O   R      : 
!  iAbort      I   I      : 

use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: CurStress(6), CurVoid, strainInc(6),  &
    CurAlpha(6), Params(25)
double precision, intent(out) :: a
integer, intent(out) :: iAbort

double precision dSigma(6), dSigma0(6), dSigma1(6), C(6,6), a0, a1,  &
    TolF, f, E_K, E_G, da

double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0D0
double precision, parameter :: two = 2.0D0

integer i

! 8-12-2021 check this function find right solution or not
iAbort = 0

a = zero
a0 = zero
a1 = one
TolF = Params(i_TolF)

! 24-11-2021
CALL GetElasticModuli(CurStress, CurVoid, Params, E_K, E_G)
CALL GetStiffness(E_K, E_G, C)
dSigma =  matmul(C, strainInc)

DO i = 1, 30
    da = (a1 - a0) / two
    a = a1 - da
    CALL GetF(CurStress + a * dSigma, CurAlpha, Params, f)
    IF (f > TolF) THEN
        a1 = a
    ELSE
        IF (f < -TolF) THEN
            a0 = a
            EXIT
        ELSE
            RETURN
        END IF
        IF (i == 30) THEN
            a = zero
            WRITE(100,*) 'fail to find a intersection factor of plastic unloading 1'
            WRITE(100, '(6f25.15)') CurStress
            WRITE(100, '(6f25.15)') CurAlpha
            WRITE(100, '(6f25.15)') strainInc
            WRITE(100, '(1f25.15)') CurVoid
            WRITE(100, '(1f25.15)') a0
            WRITE(100, '(1f25.15)') a1
            call flush(100)
            iAbort = 1
            RETURN
        END IF
    END IF
END DO

IF (a0 < 1.0D-8 .AND. a1 < 1.0D-8) THEN
    a = zero
    RETURN
END IF

IF ((f > TolF) .or. (f < -TolF)) THEN
    CALL IntersectionFactor(CurStress, CurVoid, strainInc,  &
           CurAlpha, a0, a1, Params, a, iAbort)
    IF (iAbort == 1) THEN
        WRITE(100,*) 'fail to find a intersection factor of plastic unloading 2'
        WRITE(100, '(6f25.15)') CurStress
        WRITE(100, '(6f25.15)') CurAlpha
        WRITE(100, '(6f25.15)') strainInc
        WRITE(100, '(1f25.15)') CurVoid
        WRITE(100, '(1f25.15)') a0
        WRITE(100, '(1f25.15)') a1
        call flush(100)
        iAbort = 1
        RETURN
    END IF
END IF

RETURN
END

    
!---------------------------------------------------------------------------------------------------
SUBROUTINE GetStateDependent(stress, alpha, alphaM, SV_MM_plus, &
           SV_MM_minus, alpha_in, Params, e, sn, dd, bb, bM,   &
           cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, B,  &
           C, D, R, Z, CurSigma, indicator, iStep, iAbort, iTer, iEl, Int, AlphaAlphaInDotN1,   &
           NormAlphaAlphaInDotN, AlphaAlphaInDotNNorm)
! Purpose: Update the stress tensor and the state variables, the state of the material
!       it updadated the stiffness elastoplastic stiffness matrix (remove)
!----------------------------------------------------------------------------------------------------
! Arguments:
!                           I/O   Type
!  iEl                       I    I    :
!  Int                       I    I    :
!  iStep                     I    I    :
!  iTer                      I    I    :
!  indicator                 I    I    :
!  Params                    I    R()  :
!  stress                    I    R()  :  intermediate stress 1 2 3
!  alpha                     I    R()  :
!  alphaM                    I    R()  : 
!  SV_MM_plus                I    R()  : 
!  SV_MM_minus               I    R()  :
!  alpha_in                  I    R()  : NB alfa_initial is not updated  during the substepping process
!  CurSigma                  I    R()  : current accepted Sigma for the current substep sNextStress
!  e                         I    R()  :
!  sn                        O    R()  : 
!  dd                        O    R()  :
!  bb                        O    R()  :
!  bM                        O    R()  :
!  cos3Theta                 O    R    :
!  h                         O    R    :
!  hM                        O    R    :
!  psi                       O    R    :
!  rBtheta                   O    R    : 
!  rDtheta                   O    R    :
!  b0                        O    R    :
!  A                         O    R    : coefficient [called Ad] in the equation of D = Ad (alfa_teta_d - alfa):n D
!  B                         O    R    : coefficient in the equation of R'
!  C                         O    R    : coefficient in the equation of R'
!  D                         O    R    :
!  R                         O    R()  : R' + D I
!  Z                         O    R    : mM * f_shr / PAR_zeta
!  AlphaAlphaInDotN1         O    R    :
!  NormAlphaAlphaInDotN      O    R    :
!  AlphaAlphaInDotNNorm      O    R    :
!  iAbort                    O    I    :
use extrainfo, only: onekPa
use extrainfo, only: oneAtm
use SANISANDMS_constants
use Parameters_position_in_params_iMod1

implicit none
integer, intent(in) :: indicator, iStep, iTer, iEl, Int
integer, intent(out) :: iAbort
double precision, intent(in) :: stress(6), alpha(6), alphaM(6), &
SV_MM_plus, SV_MM_minus, alpha_in(6), Params(25),e,CurSigma(6)
double precision, intent(out) :: sn(6), dd(6), bb(6), bM(6),    &
    cos3Theta, psi, rBtheta, rDtheta, b0, A, B, C, D, R(6),    &
    Z, AlphaAlphaInDotN1, NormAlphaAlphaInDotN, AlphaAlphaInDotNNorm

!REAL*16, intent(out) :: h, hM
double precision, intent(out) :: h, hM

double precision snM(6), SI1(6), rr(6), alphaBtheta(6), &
    alphaBthetaPLUSpi(6), alphaDtheta(6), r_alphaM(6),rM_tilde(6), &
    r_alphaM_tilde(6), rM(6), r_tilde(6), alphad_tilde(6), res(6), &
    stress_use(6), PAR_G0, PAR_Mc, PAR_c, PAR_m, PAR_h0, PAR_ch,   &
    PAR_nb, PAR_A0, PAR_nd, PAR_mu0, PAR_zeta, PAR_beta, D_factor, &
    p, x3, x2, SV_MM, res1, res100, rBthetaPLUSpi, gthetaPlusPi,   &
    f_shr, b_rM_rin, b_d_r, b_dM_tilde, b_bM, bref_D, bref,    &
    bM_distance_tilde, bM_distance, res77(6), res7, rDthetaPLUSpi, alphaDthetaPLUSpi(6), &
    bM_Norm, res55, res555, res8, T_rr, T_alfa, T_bb, T_alphaBtheta, T_sn, dgdth, rBtheta1, rDtheta1, rr_Norm, &
    AlphaAlphaInDotN, VectorAlphaAlphaIn(6), res2, res3, Mb_check, Mbc_max, Custom_h_max, &
    xN1(6),xN2(6),xN3(6),sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress, bM_length, bM_equ(6)

!REAL*16 hMM(6)

integer i, itestnan

double precision, parameter :: zerozo = 1.0D-2  !it is zero + tolerance (must be rewritten)
double precision, parameter :: zero_plus_tol = 1.0D-30
double precision, parameter :: twoPthree = 2.3D0
double precision, parameter :: small = 1.0D-7
double precision, parameter :: slarge = 1.0D7
double precision, parameter :: tol = 1.0D-15
double precision, parameter :: sqrt6 = DSQRT(6.0D0)

double precision PAR_P_atm
double precision PAR_Pmin

logical, external :: is_hLimitedWithDefaultStrategy, is_hM_ExceptionsStrategyActive

1070  FORMAT(A, 6F20.8)
1080  FORMAT(A, 1F25.12)

iAbort = 0

!PAR_Pmin = onekPa/100.0D0
PAR_Pmin = 0.0001D0 * onekPa
PAR_P_atm = oneAtm

PAR_G0   = Params(i_PAR_G0)
PAR_Mc   = Params(i_PAR_Mc)
PAR_c    = Params(i_PAR_c)
PAR_m    = Params(i_PAR_m)
PAR_h0   = Params(i_PAR_h0)
PAR_ch   = Params(i_PAR_ch)
PAR_nb   = Params(i_PAR_nb)
PAR_A0   = Params(i_PAR_A0)
PAR_nd   = Params(i_PAR_nd)
PAR_mu0  = Params(i_PAR_mu0)
PAR_zeta = Params(i_PAR_zeta)
PAR_beta = Params(i_PAR_beta)
Mbc_max  = Params(i_Mbc_max)
Custom_h_max = Params(i_Custom_h_max)

D_factor = one
CALL GetTrace(stress, p)
p = one3 * p

stress_use = stress
! move p
! remove PAR_Pmin
!IF (p <= PAR_Pmin) THEN
!    !write(2,*), "correction of the stress for p too low"
!    stress_use(1) = stress(1) + PAR_Pmin-p
!    stress_use(2) = stress(2) + PAR_Pmin-p
!    stress_use(3) = stress(3) + PAR_Pmin-p
!    stress_use(4) = stress(4)
!    stress_use(5) = stress(5)
!    stress_use(6) = stress(6)
!!    p = p + ABS(PAR_Pmin-p)
!    p = PAR_Pmin
!ELSE
!    stress_use = stress
!END IF

SI1(1) = one
SI1(2) = one
SI1(3) = one
SI1(4) = zero_plus_tol
SI1(5) = zero_plus_tol
SI1(6) = zero_plus_tol


Call GetDevPart(stress, rr)
!IF (p < PAR_Pmin) THEN
!    WRITE(1, *) 'p <= PAR_Pmin'
!    rr = rr/PAR_Pmin
!ELSE
rr = rr / p
!
T_rr = rr(1) + rr(2) +rr(3)
T_alfa = alpha(1) + alpha(2) + alpha(3)
!
Call GetNormalToYield(stress, alpha, iStep, sn, iAbort)
If (iAbort == 1) Then
    RETURN
End If
!
T_sn = sn(1) + sn(2) + sn(3)
!
!c---------------------------------------------------------------------------------------
! 10-01-2022 set new output for norm of |alpha_in-alpha| and s_(a_in-a):sn

! (Alpha - Alpha_in) 
res77 = alpha - alpha_in

! |(Alpha - Alpha_in)|
Call GetNorm_Contr(res77, AlphaAlphaInDotNNorm)

! (Alpha - Alpha_in)/|(Alpha - Alpha_in)|
VectorAlphaAlphaIn = res77 / AlphaAlphaInDotNNorm

! (Alpha - Alpha_in):n    [in the denominator of h and hM ]
CALL DoubleDot2_2_Contr(res77, sn, AlphaAlphaInDotN)

!  (Alpha - Alpha_in)/|(Alpha - Alpha_in)| : n
CALL DoubleDot2_2_Contr(VectorAlphaAlphaIn, sn, NormAlphaAlphaInDotN)

! store the original value of AlphaAlphaInDotN
AlphaAlphaInDotN1 = AlphaAlphaInDotN

CALL GetPSI(e, p, Params, psi)

CALL GetLodeAngle(sn, cos3Theta)

CALL g(cos3Theta, PAR_c, res1, dgdth) 

! calculate Mb and check the limit
Mb_check = PAR_Mc * EXP(-PAR_nb * psi)
if (Mb_check > Mbc_max) then 
    Mb_check = Mbc_max
endif

! calculate coeff for Alpha_teta_B
rBtheta = res1 * Mb_check - PAR_m

rBtheta1 = res1 * Mb_check

CALL g(-cos3Theta, PAR_c, res1, dgdth)

! calculate coeff for Alpha_teta_plus_pi_B
rBthetaPLUSpi = res1* Mb_check - PAR_m

! calculate coeff for Alpha_teta_B
res7 = root23 * rBtheta

! calculate Alpha_teta_B
alphaBtheta = res7 * sn

! calculate coeff for Alpha_teta_plus_pi_B
res7 = -root23 * rBthetaPLUSpi

! Alpha_teta_plus_pi_B
alphaBthetaPLUSpi = res7 * sn

CALL g(cos3Theta, PAR_c, res1, dgdth)

! calculate coeff for Alpha_teta_D
rDtheta = res1 * PAR_Mc * EXP(PAR_nd * psi) - PAR_m

rDtheta1 = res1 * PAR_Mc * EXP(PAR_nd * psi)


CALL g(-cos3Theta, PAR_c, res1, dgdth)

! calculate coeff foe Alpha_teta_plus_pi_D
rDthetaPLUSpi = res1 * PAR_Mc * EXP(PAR_nd * psi) - PAR_m

! calculate coeff for Alpha_teta_D
res7 = root23 * rDtheta

! Alpha_teta_D
alphaDtheta = res7 * sn

! calculate coeff foe Alpha_teta_plus_pi_D
res7 = - root23 * rDthetaPLUSpi

! Alpha_teta_plus_pi_D
alphaDthetaPLUSpi = res7 * sn

! To prevent ivalid sqrt use p_min = p_atm/1000 as limit for p
IF (p >= PAR_Pmin) THEN
    b0 = PAR_G0 * PAR_h0 * DABS(one - PAR_ch * e) / DSQRT(p/PAR_P_atm)
ELSE
    b0 = PAR_G0 * PAR_h0 * DABS(one - PAR_ch * e) / DSQRT(PAR_Pmin/PAR_P_atm)
END IF

dd = alphaDtheta - alpha
bb = alphaBtheta - alpha

! The memory surface and image points on it------------------------------------------------------
SV_MM = MAX(PAR_m, SV_MM_plus + SV_MM_minus)! Radious 

! calculate coeff for first term in r_Alpha_M
res7 = root23 * (SV_MM - PAR_m)

! calculate r_Alpha_M 
r_alphaM = alphaM + res7 * sn

! calculate r_Alpha_M tilde ---remove
r_alphaM_tilde = alphaM - res7 * sn

! The rM  is computed using alphaM (for fshr)
res7 = root23 * SV_MM
rM = alphaM + res7 * sn
rM_tilde = alphaM - res7 * sn

! calculate r_Alpha_M - alfa  = [rM - s23 *m * (rM - alfamM)/|| ] - alfa
res77 = r_alphaM - alpha 

! calculate bM =  (r_Alpha_M - alfa) : sn
CALL DoubleDot2_2_Contr(res77, sn, bM_distance)

IF (bM_distance < 0.0d0) THEN
    bM_distance = 0.0d0
END IF

! The rr_tilde is computed using alpha
res7 = root23 * PAR_m
r_tilde = alpha - res7 * sn

!               NB OpenSees does not use snM but it seems to make a difference
! calculate snM = (rr - rM)/|(rr - rM)|
snM = rM - rr

! Check of the denominator and set snM = sn as for virgin loadings 
CALL GetNorm_Contr(snM, res100)
!IF (abs(res100) < 1.0D-5) THEN!tol
IF (abs(res100) < 0.01D0) THEN!tol
    snM = sn
ELSE
    snM = snM / res100
END IF

res77 = rM - rM_tilde
CALL DoubleDot2_2_Contr(snM, res77, x3)

! x2 here is the numerator of fshr of the paper
res77 = r_tilde - rM_tilde                          
CALL DoubleDot2_2_Contr(snM, res77, x2)

IF ((x2 <= tol).OR.(x3<=tol) .OR. (SV_MM < zerozo)) THEN!tol
    f_shr = zero_plus_tol
ELSE
    f_shr = x2 / x3
END IF

IF (f_shr < -tol) THEN
    f_shr = zero_plus_tol
END IF

IF (f_shr > one) THEN
    f_shr = one
END IF
!c---------------------------------------------------------------------------------------
gthetaPlusPi = two * PAR_c / ((one + PAR_c) + (one - PAR_c) * (cos3Theta))     !Unused!!!!
alphad_tilde = root23 * (gthetaPlusPi * PAR_Mc * EXP(PAR_nd * psi) - PAR_m) * (-one * sn)   !Unused!!!!
!CALL DoubleDot2_2_Contr(alphad_tilde - r_alphaM_tilde, sn, b_dM_tilde)

! calculate alfa_D_tilde - alfa_ini for NEW definition of b_dM_tilde
res77 = alphaDthetaPLUSpi - alpha_in

!res77 = alphaDthetaPLUSpi - r_alphaM_tilde

! calculate b_dM_tilde NEW. From the OpenSees code. It was different in Liu 2019
CALL DoubleDot2_2_Contr(res77, sn, b_dM_tilde)

!CALL DoubleDot2_2_Contr(alphaBtheta - alphaBthetaPLUSpi, sn, res1)

! calculate b_ref = root23(Alpha_teta_B +  Alpha_teta_plus_pi_B) r is actually alfa!
res1 = (rBtheta + rBthetaPLUSpi) * root23

bref = MAX(two * root23 * PAR_m, res1)      !b_ref To be used only in h (New equations of D from OpenSees)

CALL GetNorm_Contr(alpha_in, res1)

bref_D = MAX(two * root23 * PAR_m, res1)    !bref_D To be consistent with OpenSees (see slide), it should be Norm alfa_in.
CALL DoubleDot2_2_Contr(dd, sn, b_d_r)

CALL Macauley(b_dM_tilde, res555)
A = EXP(PAR_beta * res555 / bref_D)
A = PAR_A0 * A
!
D = A * b_d_r

CALL g(cos3Theta, PAR_c, res1, dgdth)
! Argyris function
!B = one + three2 * (one - PAR_c) / PAR_c * res1 * cos3Theta
! Van Eekelen function
B = one + three * cos3Theta * dgdth

CALL g(cos3Theta, PAR_c, res1, dgdth)
! Argyris function
!C = three * root32 * (one - PAR_c) / PAR_c * res1
! Van Eekelen function
C = three * sqrt6 * dgdth

CALL SingleDot(sn, sn, res)
res7 = C * one3
res8 = one3 * D
R = B * sn - C * res + res7 * SI1 + res8 * SI1

Z = SV_MM / PAR_zeta * f_shr                                !to be used in hM

! calculate bM tensor! OUTPUT to be used in dalfaM (the equation in the paper is wrong!)  Also used in hM for b_bM
bM = alphaBtheta - r_alphaM
bM_equ = bM
!----------------------------------------------------------------------------------------------!
!-check if new hM equations for exceptions is active-------------------------------------------!
!----------------------------------------------------------------------------------------------!
if(is_hM_ExceptionsStrategyActive(Params))then
    ! calculate alphaM:n
    CALL DoubleDot2_2_Contr(alphaM, sn, bM_Norm)

    ! calculate (Alpha_teta_B - r_Alpha_M):n
    bM_Norm = root23 * rBtheta - (bM_Norm + root23 * (SV_MM - PAR_m))
    bM_Norm = abs(bM_Norm)

    call GetNorm_Contr(bM, bM_length)

    ! check if (Alpha_teta_B - r_Alpha_M):n close to zero
    IF (bM_Norm < 1.0D-7) THEN
        ! check if |Alpha_teta_B - r_Alpha_M| close to zero
	    IF (bM_length < 1.0D-7) THEN   !tol (check small)
	    !   check its sign
	        CALL DoubleDot2_2_Contr(bM, sn, res55)
	        IF (res55 < -tol) THEN
	            bM = -sn
	        ELSE
	            bM = sn
	        END IF
	        bM_length = 1.0d0
        END IF
    ELSE
        bM = bM/bM_length
        bM = bM_equ
    END IF
endif
!----------------------------------------------------------------------------------------------!
CALL DoubleDot2_2_Contr(bM_equ, sn, b_bM)

! calculate r_alphaM:n = alphaM:n + root23 * (SV_MM - PAR_m) * sn
CALL DoubleDot2_2_Contr(alphaM, sn, b_rM_rin)
b_rM_rin = b_rM_rin + root23 * (SV_MM - PAR_m)

CALL DoubleDot2_2_Contr(alpha_in, sn, res1)

! calculate (r_alphaM - alpha_ini):n
b_rM_rin = b_rM_rin - res1

CALL Macauley(-D, res1)

! turn on absolute value
AlphaAlphaInDotN = DABS(AlphaAlphaInDotN)

! add a small value to (Alpha - alpha_ini):n to prevent 0 
AlphaAlphaInDotN = AlphaAlphaInDotN + 0.001D0

! add a small value to (r_alphaM - alpha_ini):n to prevent 0 
b_rM_rin = DABS(b_rM_rin)+0.001D0

!----------------------------------------------------------------------------------------------!
!-check if new hM equations for exceptions is active-------------------------------------------!
!----------------------------------------------------------------------------------------------!
!
if(bM_Norm < 1.0D-7 .and. is_hM_ExceptionsStrategyActive(Params))then
    !new equation for corner cases
    CALL MacauleyIndex(b_bM, res2)
    CALL SignFunction(b_bM, res3)
    !hM = zeroP5 * b0 / b_rM_rin * bM_Norm + zeroP5 / root23 * Z * res1 / b_bM
	hM = zeroP5 * b0 / b_rM_rin * bM_length + 0.5d0 *Z / root23 * res1 / (res2 + one) * res3
else	
	!old equation
	hM = 0.5d0 * b0 / b_rM_rin + 0.5d0 / root23 * Z * res1 / b_bM
endif

if(.not. is_hM_ExceptionsStrategyActive(Params))then
        If (IsNan(hM)) Then
            hM = 1.0D7
        End If

        hM = MIN(1.0D7, hM)
endif
!----------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------!
!-check if customized limit on h is defined----------------------------------------------------!
!----------------------------------------------------------------------------------------------!
!
IF (p < zero_plus_tol) THEN    
    h = b0 / AlphaAlphaInDotN
ELSE
    res2 = PAR_mu0 * ((p/PAR_P_atm)**zeroP5) * ((bM_distance / bref)**2)
    IF (res2 > 9.2103D0 .and. is_hLimitedWithDefaultStrategy(Params)) THEN
        h = b0 / AlphaAlphaInDotN * EXP(9.2103D0)
    ELSE
        h = b0 / AlphaAlphaInDotN * EXP(res2)
    END IF
END IF

If (IsNaN(h)) Then
    Write(100, *), 'h NAN= ', h
    call flush(100)
    iabort = 1
    return    
End If

IF (is_hLimitedWithDefaultStrategy(Params))then
    h = MIN(1.0D7, h)
ELSE
    h = MIN(Custom_h_max, h)
ENDIF

!----------------------------------------------------------------------------------------------!

!     NaN check and output
!     sn, dd, bb, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, B, C, D, R, Z)

Do i = 1, 6
    if(IsNaN( sn( i))) Then
        Write(100,*), 'sn',sn
        Write(100, 1070), 'stress_use = ', stress_use
        Write(100, 1070), 'alpha = ', alpha
        call flush(100)
        iabort = 1
        return
    End If
End Do
!
Do i=1,6
    If( IsNaN( dd( i))) Then
        Write(100,*), 'iEl', iEl
        Write(100, 1070), 'dd = ', dd
        call flush(100)
        iabort = 1
        return
    End If
End Do

If( IsNaN( cos3Theta)) Then
        Write(100,*), 'iEl', iEl
        Write(100, 1080) 'cos3Theta = ', cos3Theta
        Write(100, 1080) 'psi = ', psi
        call flush(100)
        iabort = 1
        return
End If

Do i=1,6
    If( IsNaN( bM( i))) Then
        Write(100,*), 'iEl', iEl
        Write(100,*), 'bM', bM
        Write(100, 1070), 'r_alphaM = ', r_alphaM
        call flush(100)
        iabort = 1
        return
    End If
End Do

If ( IsNaN( h)) Then
    Write(100,*), 'h NaN',h
    Write(100, 1080), 'AlphaAlphaInDotN = ', AlphaAlphaInDotN
    Write(100, 1080), 'bM_distance = ', bM_distance
    Write(100, 1080), 'bref = ', bref
    Write(100, 1080), 'b0 = ', b0
    call flush(100)
    iabort = 1
    return
End If

If ( IsNaN( hM)) Then
    Write(100,*) 'hM NaN'
    Write(100, 1080), 'b_rM_rin = ', b_rM_rin
    Write(100, 1080), 'Z = ', Z
    Write(100, 1080), 'b_bM = ', b_bM
    Write(100, 1070), 'sn = ', sn
    Write(100, 1080), 'bM_Norm = ', bM_Norm
    Write(100, 1070), 'stress = ', stress
    Write(100, 1070), 'alpha = ', alpha
    call flush(100)
    iabort = 1
    return
End If

If ( IsNaN( A)) Then
    Write(100,*), 'A NaN'
    Write(100,*), 'res555 = ', res555
    Write(100,*), 'bref_D = ', bref_D
    Write(100,*), 'res77 = ',res77
    Write(100,*), 'sn = ', sn
    Write(100,*), 'alphaDthetaPLUSpi = ', alphaDthetaPLUSpi 
    Write(100,*), 'r_alphaM_tilde = ', r_alphaM_tilde
    call flush(100)
    iabort = 1
    return
End If

If ( IsNaN( B)) Then
    Write(100,*), 'B NaN',B
    call flush(100)
    iabort = 1
    return
End If

If (IsNan( C)) Then
    Write(100,*), 'C NaN',C
    call flush(100)
    iabort = 1
    return
End If

If ( IsNaN( D)) Then
    Write(100,*), 'D NaN', D
    call flush(100)
    iabort = 1
    return
End If

Do i=1,6
    If ( IsNaN( R(i))) Then
        Write(100,*), 'R(i)', R
        call flush(100)
        iabort = 1
        return
    End If
End Do

If ( IsNaN(Z) ) Then           
    Write(100,*) 'Z NaN'
    Write(100, 1080) 'SV_MM = ', SV_MM
    Write(100, 1080) 'f_shr = ', f_shr
    Write(100, 1080) 'x3 = ', x3
    Write(100, 1080) 'x2 = ', x2
    Write(100, 1070) 'rM = ', rM
    Write(100, 1070) 'rM_tilde = ', rM_tilde
    Write(100, 1070) 'r_tilde = ', r_tilde
    Write(100, 1070) 'sn = ', sn
    Write(100, 1070) 'alphaM = ', alphaM
    Write(100, 1070) 'stress = ', stress
    Write(100, 1070) 'alpha = ', alpha
    call flush(100)
    iabort = 1
    return
End If

RETURN
END

    
!------------------------------------------------------------------------------------
SUBROUTINE StateVariable(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, Params, &
                        thisVoidRatio, rr, dDevStrain, dVolStrain, &
                        dSigma, dAlpha, dAlphaM,  dMM_plus, dMM_minus, aCep, &
                        iEl, Int, iStep, CurSigma, indicator,iRunge, dT, T, iAbort, iTer, FrozenStress, Lindicator, &
                        dSigmaE, dAlpha0)
!
! Purpose: Update the stress tensor and the state variables, the state of the material
!       it updadated the stiffness elastoplastic stiffness matrix (remove)
!------------------------------------------------------------------------------------
! Arguments:
!                  I/O   Type
!  iEl              I    I    :
!  Int              I    I    :
!  iStep            I    I    :
!  iTer             I    I    :
!  indicator        I    I    :
!  iRunge           I    I    :
!  dT               I    R    :
!  T                I    R    :
!  FrozenStress     I    R()  : stress at the beginning of the step NB: based on the iter, it can be equilibrated or no {yes if Iter =1}
!  Params           I    R()  :
!  thisSigma        I    R()  :
!  thisVoidRatio    I    R()  :
!  thisAlpha        I    R()  :
!  thisAlphaM       I    R()  : 
!  thisMM_plus      I    R()  : 
!  thisMM_minus     I    R()  :
!  Curalpha_in      I    R()  : NB alfa_initial is not updated  during the substepping process
!  CurSigma         I    R()  : current accepted Sigma for the current substep sNextStress
!  rr               I    R()  : computed with the intermediate sigma 1 2 or 3 NB it could be this_rr
!  dDevStrain       I    R()  :
!  dVolStrain       I    R    :
!  Lindicator       O    I    : deprecated after debugging. It is set to 1 when the plastic multiplier is <=0
!  dSigma           O    R()  :
!  dAlpha           O    R()  : 
!  dAlphaM          O    R()  :
!  dMM_plus         O    R    :
!  dMM_minus        O    R    :
!  aCep             O    R    : the elastic stiffness matriax is currently used for IDTask = 6 and 3, therefore it is set to o
!  iAbort           O    i    :
!  dSigmaE          O    R    : elastic portion of the stress increment (used for debugging)
!  dAlpha0          O    R    : output dalpha for Lindicator = 1 (used for debugging)

use extrainfo, only: onekPa
use extrainfo, only: oneAtm
use SANISANDMS_constants

implicit none

integer, intent(in   ) :: iEl, Int, iStep, indicator, iRunge, iTer
integer, intent(  out) :: iAbort, Lindicator
double precision, intent(in  ) :: thisSigma(6), thisAlpha(6), thisAlphaM(6), thisMM_plus, thisMM_minus, Curalpha_in(6),  &
    Params(25), thisVoidRatio, rr(6), dDevStrain(6),dVolStrain, CurSigma(6), dT, T, FrozenStress(6)
double precision, intent(  out) :: dSigma(6), dAlpha(6), dAlphaM(6), dMM_plus, dMM_minus, aCep(6,6), dSigmaE(6), dAlpha0(6)

double precision :: p
double precision sn(6), dd(6),bb(6), bM(6), R(6), sI1(6), dPStrain(6), aD(6, 6), res(6), res1(6), &
    dPStrain1(6),dStrainInc(6), E_K, E_G, Halpha(6,6), PM, res3, A, B, C, D, Z, temp4, sNextDGamma, &
    res2, psi, rBtheta, rDtheta, depsilon_pv, Cos3Theta, b0, dSigma_A(6), &
    dSigma_B(6), dSigma_C(6), res_coef, Cf1_dsigma_C, Cf2_dsigma_C, Cf3_dsigma_C, Cf4_dsigma_C, res222, res8,   &
    E_K1, E_G1, dSigma_A1(6), dSigma_B1(6), Cf1_dsigma_C1, Cf2_dsigma_C1, Cf3_dsigma_C1, Cf4_dsigma_C1, &
    dSigma_C1(6), p_elastic, p_plastic, PM1, temp5, temp4_E, temp4_P, pFrozen, AlphaAlphaInDotN,  &
    NormAlphaAlphaInDotN, AlphaAlphaInDotNNorm


double precision h, hM
integer i, itestnan
double precision, parameter :: zero_with_tol = 1.0D-30 !0.0d0 using a very small value to replace zero
!double precision, parameter :: small  = 1.0D-7
!double precision, parameter :: tol  = 1.0D-18
double precision PAR_P_atm
double precision PAR_Pmin

! initialize Lindicator
Lindicator = 0
dSigma = 0.0D0
aCep = 0.0d0

1170  FORMAT(A, 6F35.23)
1180  FORMAT(A, 1F28.18)
1190  FORMAT(A, 6F30.20)
1200  FORMAT(32F30.15)

iAbort = 0

PAR_Pmin = onekPa/100d0
PAR_P_atm = oneAtm

sI1(1) = one
sI1(2) = one
sI1(3) = one
sI1(4) = zero_with_tol
sI1(5) = zero_with_tol
sI1(6) = zero_with_tol

CALL GetTrace(thisSigma, p)
p = one3 * p
CALL GetTrace(FrozenStress, pFrozen)
pFrozen = one3 * pFrozen

! Calculate elastic moduli 
CALL GetElasticModuli1(thisSigma, thisVoidRatio, Params, E_K1, E_G1)

! Check p for the Elastic moduli
If (p <= 0.0D0) Then
    E_K1 = zero_with_tol
    E_G1 = zero_with_tol
End If

CALL GetStateDependent(thisSigma,thisAlpha,thisAlphaM,thisMM_plus, thisMM_minus, Curalpha_in, Params, thisVoidRatio, sn, dd, &
         bb, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, B, C, D, R, Z, CurSigma, indicator, &
         iStep, iAbort, iTer, iEl, Int, AlphaAlphaInDotN, NormAlphaAlphaInDotN, AlphaAlphaInDotNNorm)

If (iAbort == 1) Then
    return
END IF

! Check p for the Plastic modulus
CALL DoubleDot2_2_Contr(bb, sn, PM1)
If (p >= 0.0D0) Then
    PM = two3 * p * h * PM1
Else
    PM = two3 * 0.0D0 * h * PM1
End If

! Set PM = 0 if is negative, as in the PM4Sand report 
If (PM < 1.0D-30) Then
    PM = 1.0D-30
End If

CALL SingleDot(sn, sn, res)
CALL SingleDot(sn, res, res1)
CALL GetTrace(res1, res2)
CALL DoubleDot2_2_Contr(sn, rr, res3)

! Calculate the denom of the plastic multiplier
temp4 = (PM + two * E_G1 * (B - C * res2) - E_K1 * D * res3)

! for debugging only---------------------------------------------------
temp4_E = two * E_G1 * (B - C * res2) - E_K1 * D * res3
temp4_P = PM
!----------------------------------------------------------------------

CALL DoubleDot2_2_Contr(sn, rr, res2)
CALL DoubleDot2_2_Mixed(sn, dDevStrain, res3)

! Calculate the numerator of the plastic multiplier
temp5 = two * E_G1 * res3 - E_K1 * dVolStrain * res2

! Calculate the plastic multiplier L
sNextDGamma = (two * E_G1 * res3 - E_K1 * dVolStrain * res2)/temp4

! Check the for negative L
IF (sNextDGamma < zero_with_tol) THEN
     Lindicator = 1
END IF

CALL ToContraviant(dDevStrain, res)
CALL SingleDot(sn, sn, res1)

! Calculate <L> = sNextDGamma
CALL Macauley(sNextDGamma, res2)

! Calculate isotropic and deviatoric elastic contributes in dSigmaE
dSigma_A1 = (two * E_G1) * res
dSigma_B1 = (E_K1 * dVolStrain) * sI1

! Calculate plastic portion in dSigma
Cf1_dsigma_C1 = res2 * two * E_G1 * B
Cf2_dsigma_C1 = res2 * two * E_G1 * C
Cf3_dsigma_C1 = res2 * two * E_G1 * C * one3
Cf4_dsigma_C1 = res2 * E_K1 * D
dSigma_C1 = (Cf1_dsigma_C1 * sn) - (Cf2_dsigma_C1 * res1) + (Cf3_dsigma_C1 * sI1) + (Cf4_dsigma_C1 * sI1)

CALL GetTrace(dSigma_A1+dSigma_B1, p_elastic)
CALL GetTrace(dSigma_C1, p_plastic)

! Calculate elastic portion of dSigma
dSigmaE = dSigma_A1 + dSigma_B1

! Calculate dSigma
dSigma = dSigma_A1 + dSigma_B1 - dSigma_C1

! Calculate <L> = sNextDGamma
CALL Macauley(sNextDGamma, res222)

! Calculate the hardening of the yield
res8 = res222 * two3 * h
dAlpha = res8 * bb
dAlpha0 = zero_with_tol ! remove it

! Calculate the hardening of the Memory Surface
res8 = res222 * two3 * hM
dAlphaM = res8 * bM             !this bM is computed with the images along n (tilde in the equation is a typo)

CALL ToCovariant(R, dPStrain1)
dPStrain1 = sNextDGamma * dPStrain1
CALL GetTrace(dPStrain1, depsilon_pv)

CALL DoubleDot2_2_Contr(dAlphaM, sn, dMM_plus)
dMM_plus = one / root23 * dMM_plus

CALL Macauley(-depsilon_pv, res2)
dMM_minus = -Z * res2

!Elastic matrix is currently used-----------------------------------------------------------------
!CALL GetElastoPlasticTangent(thisSigma, sNextDGamma, E_G, E_K, B, C, D, h, sn, dd, bb, aCep, temp4, Halpha(6,6))
! 11-10-2021
!write(1, *), 'Dep :'
!write(1, '(6F14.6)'), aCep
!------------------------------------------------------------------------------------------------

! 11-15-2021
! set Halpha = 0, if L < 10^-15
!CALL MacauleyIndex1(sNextDGamma, res222)
!IF (sNextDGamma > 1.0D-11) THEN
!    dAlpha = matmul(Halpha, dStrainInc)
!ELSE
!    dAlpha = 0.0D0
!END IF

DO i = 1, 6
    IF (IsNaN(dSigma(i))) THEN
        Write(100,*), 'dSigma(i)', dSigma(i)
        call flush(100)
        iabort = 1
        !iabort = 100
        return
    END IF
END DO

DO i = 1, 6
    IF (IsNaN(dAlpha(i))) THEN
        WRITE(100, 1180) 'sNextDGamma = ', res222
        WRITE(100, 1180) 'h = ', h
        WRITE(100, 1170) 'bb = ', bb
        WRITE(100, 1170) 'dDevStrain =', dDevStrain
        WRITE(100, 1180) 'dVolStrain =', dVolStrain
        call flush(100)
        iabort = 1
        !iabort = 100
        return
    END IF
END DO

DO i = 1, 6
    IF (IsNaN(dAlphaM(i))) THEN
        WRITE(100, 1180) 'sNextDGamma = ', res222
        WRITE(100, *) 'hM = ', hM
        WRITE(100, *) 'sNextDGamma * hM = ', res222 * hM
        WRITE(100, *) 'h = ', h
        WRITE(100, *) 'hM = ', hM
        WRITE(100, 1170) 'bM = ', bM
        WRITE(100, 1170) 'bb = ', bb
        WRITE(100, 1170) 'sn = ', sn
        WRITE(100, 1170) 'thisSigma = ', thisSigma
        WRITE(100, 1170) 'thisAlpha = ', thisAlpha
        WRITE(100, 1170) 'dDevStrain =', dDevStrain
        WRITE(100, 1180) 'dVolStrain =', dVolStrain
        call flush(100)
        iabort = 1
        !iabort = 100
        return
    END IF
END DO

IF (IsNaN(dMM_plus)) THEN
    WRITE(100, *) 'dMM_plus = ', dMM_plus
    call flush(100)
    iabort = 1
    !iabort = 100
    return
END IF

IF (IsNaN(dMM_minus)) THEN
    WRITE(100, *) 'dMM_minus = ', dMM_minus
    call flush(100)
    iabort = 1
    !iabort = 100
    return
END IF

!IF ((itestnan == 1) .AND. (indicator == 2)) THEN
!    WRITE(1, *) 'Output'
!    WRITE(1,'(6f15.8)') dSigma
!    WRITE(1,'(6f15.8)') dAlpha
!    WRITE(1,'(6f15.8)') dAlphaM
!    WRITE(1,'(1f30.10)') dMM_plus
!    WRITE(1,'(1f30.10)') dMM_minus
!    WRITE(1, *) 'Input'
!    WRITE(1,'(6f15.8)') thisSigma
!    WRITE(1,'(6f15.8)') thisAlpha
!    WRITE(1,'(6f15.8)') thisAlphaM
!    WRITE(1,'(1f15.8)') thisMM_plus
!    WRITE(1,'(1f15.8)') thisMM_minus
!    WRITE(1,'(6f15.8)') Curalpha_in
!END IF

RETURN
END
!-------------------------------------------------------------------------------------   
SUBROUTINE RungeKutta23(CurStress, CurVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, strainInc, Params, &
    sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoidRatio, aCep, &
    IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
!
! Purpose: Update the stress tensor and the state variables, the state of the material
!       
!
! Arguments:
!                  I/O   Type
!  IDTask           I    I()   :
!  iEl              I    I()   :
!  Int              I    I()   :
!  iStep            I    I()   :
!  iTer             I    I()  :
!  indicator        I    I()  :
!  strainInc        I    R() 
!  Params           I    R()  :
!  CurStress        I    R()  : At the beginning of the elastoplastic integration
!  CurVoid          I    R()  :
!  CurAlpha         I    R()  :
!  CurAlphaM        I    R()  : 
!  CurMM_plus       I    R()  : Current size of the Memory that would result only from modifications due to the rotational hardening 
!  CurMM_minus      I    R()  :
!  Curalpha_in      I    R()  : Next alpha_ini is already updated in integrator
!  sNextStress      O    R()  : At the accepted substep
!  sNextVoid        O    R()  :
!  sNextAlpha       O    R()  :
!  sNextAlphaM      O    R()  :
!  sNextMM_plus     O    R()  : 
!  sNextMM_minus    O    R()  :
!  aCep             O    R()  :
!  PLNoTension      O    R    : dismissed parameter used for debugging to track increment of strain computed without respecting the error tolerance
!  iAbort           O    i    :
!  PsPsmall         IO   R    : dismissed parameter used for debugging as warning for a mean effective stress smaller than the minimum
!
use Parameters_position_in_params_iMod1
use extrainfo, only: onekPa
use extrainfo, only: oneAtm
use SANISANDMS_constants

implicit none
integer, intent(in   ) :: IDTask, iEl, Int, iStep, indicator, iTer
integer, intent(  out) :: iAbort
double precision, intent(in  ) :: CurStress(6), CurVoid, CurAlpha(6),  CurAlphaM(6), CurMM_plus, CurMM_minus, Curalpha_in(6), &
    strainInc(6), Params(25)
double precision, intent(  out) :: sNextStress(6), sNextAlpha(6), sNextAlphaM(6), sNextMM_plus, sNextMM_minus, sNextVoidRatio, &
    aCep(6,6), PLNoTension
double precision, intent(inout) :: PsPsmall

double precision sn(6), dd(6), R(6), dDevStrain(6), rr(6), sI1(6),  snStress(6), snAlpha(6), snAlphaM(6), sndPStrain(6),   &
    dSigma1(6), dSigma2(6), dSigma3(6), dSigma(6),  dAlpha1(6), dAlpha2(6), dAlpha3(6), dAlpha(6), &
    dAlphaM1(6), dAlphaM2(6), dAlphaM3(6),  dAlphaM(6),    & 
    dPStrain1(6), dPStrain2(6), dPStrain3(6), dPStrain(6),    &
    aCep1(6, 6), aCep2(6, 6), aCep3(6, 6), aC(6, 6), aD(6, 6), res(6), res1(6), thisSigma(6),     &
    thisAlpha(6), thisAlphaM(6), r_alphaM_trial(6), rM_Alphatilde_trial(6), rM_trial(6), sInterstress(6),  &
    dSigma_H(6), dAlpha_H(6), dAlphaM_H(6), dPStrain_H(6), snStress_H(6), snAlpha_H(6), snAlphaM_H(6), PAR_m, TolE,   &
    TolR, PAR_e_init, T, dT, dT_min, p, thisVoidRatio, thisMM_plus, stressNorm_H, snMM_plus_H, snMM_plus, &
    snMM_minus_H, snMM_minus, res2, dVolStrain, dMM_plus_H, dMM_plus, dMM_minus_H, dMM_minus, dMM3_plus,   &
    dMM3_minus, dMM2_plus, dMM2_minus, dMM1_plus, dMM1_minus,  depsilon_pv, curStepError2, curStepError1, alphaNorm_H,  &
    stressNorm, q, curStepError, thisMM_minus, E_G, E_K, res7(6), strainIncNew(6), sNextStress1(6), sNextAlpha1(6), &
    fnNext, sNewstress(6), sNewAlpha(6), xN1(3),xN2(3),xN3(3),sNextStress1a,sNextStress2a,sNextStress3a,PNextStress, &
    QNextStress,this_trace_Alpha, sn_trace_Alpha, &
    this_trace_AlphaM, sn_trace_AlphaM, curStepError2_alpha, curStepError1_stress, cos3Theta, res24, dgdth, psi, Mb, Md, &
    PAR_Mc, PAR_c, PAR_nb, PAR_nd, Mc, sNextAlpha_Norm, rr_Norm, FrozenStress(6), FrozenAlpha(6), FrozenAlphaM(6), &
    FrozenMM_plus, FrozenMM_minus, q_Norm_New, q_Norm_Old, alphaNorm, FrozenStress_Norm, sNextStress_Norm, emax, emin, pFrozen, &
    Psmall, dSigmaE1(6), dAlpha01(6), dSigmaE2(6), dAlpha02(6), dSigmaE3(6), dAlpha03(6), sNextStressCopy(6), sNextAlphaCopy(6), &
    sNextAlphaMCopy(6), sNextMM_plusCopy, sNextMM_minusCopy, small_NormSig, interstress(6), interalpha(6), Mbc_max, Cos3ThetaN

integer i, izero, SNum_Inter, iNotension, Lindicator1, Lindicator2, Lindicator3

double precision, parameter :: small = 1.0D-7! small Norm alpha for Err alpha
double precision, parameter :: tol = 1.0D-18
double precision, parameter :: small_adim = 1.0D-12
double precision PAR_P_atm
double precision PAR_Pmin

1170  FORMAT(A, 6F28.18)
1180  FORMAT(A, 1F28.18)
1200  FORMAT(2I, 10F15.10)
1201  FORMAT(16F15.10)
1202  FORMAT(17F15.10)
1203  FORMAT(2F30.15)

iAbort = 0

Psmall = 1.0D-2 * onekPa! when dT is < dTmin
small_NormSig = 0.5D0 * onekPa! for Error Sig 
PAR_Pmin = onekPa/100d0
PAR_P_atm = oneAtm

emax        = Params(i_emax)
emin        = Params(i_emin)
PAR_Mc      = Params(i_PAR_Mc)
PAR_c       = Params(i_PAR_c)
PAR_m       = Params(i_PAR_m)
PAR_nb      = Params(i_PAR_nb)
PAR_nd      = Params(i_PAR_nd)
TolE        = Params(i_TolF)
TolR        = Params(i_TolR)
PAR_e_init	= Params(i_PAR_e_init)
Mbc_max     = Params(i_Mbc_max)

PsPsmall = zero !deprecated variable 
PLNoTension = zero !deprecated variable 

T = zero
dT = one
!dT_min = onePM5
! change to 10-7 26-11-2021
dT_min = 1.0D-8

sI1(1) = one
sI1(2) = one
sI1(3) = one
sI1(4) = zero
sI1(5) = zero
sI1(6) = zero

CALL GetElasticModuli1(CurStress, CurVoid, Params, E_K, E_G)
CALL GetStiffness(E_K, E_G, aC)

sNextStress   = CurStress
sNextAlpha    = CurAlpha
sNextAlphaM   = CurAlphaM
sNextMM_plus  = CurMM_plus
sNextMM_minus = CurMM_minus

FrozenStress   = sNextStress
FrozenAlpha    = sNextAlpha
FrozenAlphaM   = sNextAlphaM
FrozenMM_plus  = sNextMM_plus
FrozenMM_minus = sNextMM_minus

CALL GetTrace(sNextStress, p)
p = one3 * p

CALL GetDevPart(sNextStress, rr)
rr = rr / p

CALL GetTrace(FrozenStress, pFrozen)
pFrozen = one3 * pFrozen

SNum_Inter = 0

DO WHILE (T < one)

    iNotension = 0!    currently, it will be switched to 1 if principal stress negative [<-10^-10] for the final solution or if izero = 1 (if iNotension = 0 everything is fine)
    izero = 0!         currently, it will be switched to 1 if a principal stress component is negative (for the partial solutions)
    SNum_Inter = SNum_Inter + 1
    
!   !Calculate strainInc until T (for calculating the void ratio at T)
    res7 = T * strainInc
    
!   !Calculate  !dEpsVol until T
    CALL GetTrace(res7, res2)
    
!   !Calculate  void ratio    
    sNextVoidRatio = CurVoid - (one + CurVoid) * res2
    
!   ! Limit the voidratio to emax and emin
    IF (sNextVoidRatio >= emax) THEN
        sNextVoidRatio = emax
    END IF
    IF (sNextVoidRatio <= emin) THEN
        sNextVoidRatio = emin
    END IF

!   ! Calculate dEpsVol
    CALL GetTrace(strainInc, res2)
    
!   ! Calculate  d_EpsVol_SS
    dVolStrain = dT * res2
    
!   ! Calculate d_strainInc_SS()
    strainIncNew = dT * strainInc
    
!   ! Calculate dev_d_strainInc_SS()
    CALL GetDevPart(strainIncNew, dDevStrain)

    depsilon_pv = zero   !remove it?

!   ! Calculate Delta 1   
    thisSigma = sNextStress
    thisAlpha = sNextAlpha
    thisAlphaM = sNextAlphaM
    thisMM_plus = sNextMM_plus
    thisMM_minus = sNextMM_minus
    thisVoidRatio = sNextVoidRatio
    
!   ! Correction to restore deviatoric property-----------------------------------
    this_trace_Alpha = thisAlpha(1) + thisAlpha(2) + thisAlpha(3)
    call RestoreDviatoricProperty(thisAlpha,this_trace_Alpha)  

    this_trace_AlphaM = thisAlphaM(1) + thisAlphaM(2) +thisAlphaM(3)
    call RestoreDviatoricProperty(thisAlphaM,this_trace_AlphaM)  

!   ! Check of eventually illegal intermediate stress-----------------------------
    CALL GetTrace(thisSigma, p)
    p = one3 * p
   
    call PrnSig(1,thisSigma,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)    
    if (min(sNextStress1a,sNextStress2a,sNextStress3a)<zero) then
        izero = 1
    end if
!--------------------------------------------------------------------------      
    CALL GetDevPart(thisSigma, rr)
    rr = rr / p

    CALL StateVariable(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, Params, &
        thisVoidRatio, rr, dDevStrain, dVolStrain, dSigma1, dAlpha1, dAlphaM1, dMM1_plus, dMM1_minus, aCep1, &
        iEl, Int, iStep, sNextStress, indicator, 1, dT, T, iAbort, iTer, FrozenStress, Lindicator1, dSigmaE1, dAlpha01)
    IF (iAbort == 1) THEN
        WRITE(100, *) 'RKF1', T, dT
        call flush(100)
        RETURN
    END IF

!   ! Calculate Delta 2
    ! checking Lindicator1 before RKF2
    ! you can comment the following four line, if wanna deactivate elastic solution.
    !IF (Lindicator1 == 1) THEN
    !    dSigma1 = dSigmaE1
    !    dAlpha1 = dAlpha01
    !END IF

    thisSigma = sNextStress + (zeroP5 * dSigma1)
    thisAlpha = sNextAlpha + (zeroP5 * dAlpha1)
    thisAlphaM = sNextAlphaM + (zeroP5 * dAlphaM1)
    thisMM_plus	= sNextMM_plus + (zeroP5 * dMM1_plus)
    thisMM_minus = sNextMM_minus + (zeroP5 * dMM1_minus)

!   ! Correction to restore deviatoric property---------------------------------------------------------
    this_trace_Alpha = thisAlpha(1) + thisAlpha(2) + thisAlpha(3)
    call RestoreDviatoricProperty(thisAlpha,this_trace_Alpha)
    
    this_trace_AlphaM = thisAlphaM(1) + thisAlphaM(2) +thisAlphaM(3)
    call RestoreDviatoricProperty(thisAlphaM,this_trace_AlphaM)

!   ! Check of eventually illegal intermediate stress-----------------------------
    call PrnSig(1,thisSigma,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)
    
    if (min(sNextStress1a,sNextStress2a,sNextStress3a)<zero) then
        izero = 1
    end if
!---------------------------------------------------------------------------------------------------       
    CALL GetTrace(thisSigma, p)
    p = one3 * p

    CALL GetDevPart(thisSigma, rr)
    rr = rr / p

    CALL StateVariable(thisSigma, thisAlpha, thisAlphaM, thisMM_plus,   &
        thisMM_minus, Curalpha_in, Params, thisVoidRatio,   &
        rr, dDevStrain, dVolStrain, dSigma2, dAlpha2, dAlphaM2,    &
        dMM2_plus, dMM2_minus, aCep2, iEl, Int, iStep,sNextStress, &
        indicator, 2, dT, T, iAbort, iTer, FrozenStress, Lindicator2, dSigmaE2, dAlpha02)
    IF (iAbort == 1) THEN
        WRITE(100, *) 'RKF2', T, dT
        call flush(100)
        RETURN
    END IF

!   ! Calculate Delta 3
    ! checking Lindicator1
    !IF (Lindicator2 == 1) THEN
    !    dSigma2 = dSigmaE2
    !    dAlpha2 = dAlpha02
    !END IF

    thisSigma = sNextStress - dSigma1 + two * dSigma2
    thisAlpha = sNextAlpha -  dAlpha1 + two * dAlpha2
    thisAlphaM = sNextAlphaM -dAlphaM1 + two * dAlphaM2
    thisMM_plus = sNextMM_plus - dMM1_plus + two * dMM2_plus
    thisMM_minus = sNextMM_minus - dMM1_minus + two * dMM2_minus
        
!   ! Correction to restore deviatoric property---------------------------------------------------------
    this_trace_Alpha = thisAlpha(1) + thisAlpha(2) + thisAlpha(3)
    call RestoreDviatoricProperty(thisAlpha, this_trace_Alpha)
    
    this_trace_AlphaM = thisAlphaM(1) + thisAlphaM(2) +thisAlphaM(3)
    call RestoreDviatoricProperty(thisAlphaM,this_trace_AlphaM)

!   ! Check of eventually illegal intermediate stress-----------------------------    
    call PrnSig(1,thisSigma,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)

    if (min(sNextStress1a,sNextStress2a,sNextStress3a)<zero) then
        izero = 1
    end if
!---------------------------------------------------------------------------------------------------         
    CALL GetTrace(thisSigma, p)
    p = one3 * p

    CALL GetDevPart(thisSigma, rr)
    rr = rr / p

    CALL StateVariable(thisSigma, thisAlpha, thisAlphaM,    &
        thisMM_plus, thisMM_minus, Curalpha_in, Params, thisVoidRatio,   &
        rr, dDevStrain, dVolStrain, dSigma3, dAlpha3, dAlphaM3,    &
        dMM3_plus, dMM3_minus, aCep3, iEl, Int, iStep,sNextStress, &
        indicator, 3, dT, T, iAbort, iTer, FrozenStress, Lindicator3, dSigmaE3, dAlpha03)
    IF (iAbort == 1) THEN
        WRITE(100, *) 'RKF3', T, dT
        WRITE(100, 1170) 'sNextStress = ', sNextStress
        WRITE(100, 1170) 'dSigma1 = ', dSigma1
        WRITE(100, 1170) 'dSigma2 = ', dSigma2
        WRITE(100, 1170) 'strainIncNew = ', strainIncNew
        call flush(100)
        RETURN
    END IF

    ! checking Lindicator before building final solution
    ! comment the following Lindicator checking, when wanna deactivate elastic solution
    !IF (Lindicator1 == 1 .OR. Lindicator2 == 1 .OR. Lindicator3 == 1) THEN
    !    dSigma1 = dSigmaE1
    !    dSigma2 = dSigmaE2
    !    dSigma3 = dSigmaE3
    !    dAlpha1 = dAlpha01
    !    dAlpha2 = dAlpha02
    !    dAlpha3 = dAlpha03
    !END IF
!// Higer order increments
    dSigma_H    = (one6 * dSigma1) + (two3 * dSigma2) + (one6 * dSigma3)
    dAlpha_H    = (one6 * dAlpha1) + (two3 *  dAlpha2) + (one6 * dAlpha3)
    dAlphaM_H   = (one6 * dAlphaM1) + (two3 *  dAlphaM2) + (one6 * dAlphaM3)
    dMM_plus_H  = (one6 * dMM1_plus) + (two3 *  dMM2_plus) + (one6 * dMM3_plus)
    dMM_minus_H = one6 * dMM1_minus + two3 *  dMM2_minus + (one6 * dMM3_minus)
    dPStrain_H  = (one6 * dPStrain1) + (two3 * dPStrain2) + (one6 * dPStrain3)

!// Lower order increments    
    dSigma   = (dSigma2) 
    dAlpha  = (dAlpha2) 
    dAlphaM  = (dAlphaM2)
    dMM_plus = ( dMM2_plus)
    dMM_minus =( dMM2_minus) 
    dPStrain = ( dPStrain2)

    snStress   = sNextStress   + dSigma
    snAlpha    = sNextAlpha    + dAlpha
    snAlphaM   = sNextAlphaM   + dAlphaM
    snMM_plus  = sNextMM_plus  + dMM_plus
    snMM_minus = sNextMM_minus + dMM_minus
    
    call PrnSig(1,snStress,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)! TODO delete
    CALL GetTrace(snStress, p)
    p = one3 * p

    snStress_H   = sNextStress   + dSigma_H
    snAlpha_H    = sNextAlpha    + dAlpha_H
    snAlphaM_H   = sNextAlphaM   + dAlphaM_H
    snMM_plus_H  = sNextMM_plus  + dMM_plus_H
    snMM_minus_H = sNextMM_minus + dMM_minus_H
!
!   !Calculate the norm. 
    CALL GetNorm_Contr(snStress_H, stressNorm_H)
    CALL GetNorm_Contr(snAlpha_H, alphaNorm_H)
    
!   !Calculate stress Error
    res7 = snStress_H - snStress
    CALL GetNorm_Contr(res7, curStepError1_stress)    
!
    IF (stressNorm_H < small_NormSig) THEN
        IF ((curStepError1_stress / stressNorm_H) <= 0.01D0 .AND. (curStepError1_stress / stressNorm_H) > TolR) THEN
            !stressNorm_H = 1.0D0
            stressNorm_H = curStepError1_stress / (0.5d0*TolR)
        END IF
    END IF
    curStepError1 = curStepError1_stress / stressNorm_H

!   !Calculate back-stress ratio Error----------------------------------------------------------
    res7 = snAlpha_H - snAlpha
    CALL GetNorm_Contr(res7, curStepError2_alpha)
    IF (alphaNorm_H < small) THEN
        !write(1, *) , 'CORRECT NORMAL alpha', iStep, iEl, T
        alphaNorm_H = 1.0D0
    END IF
    curStepError2 = curStepError2_alpha / alphaNorm_H
!---------------------------------------------------------------------------------------------
    !curStepError = MAX(curStepError1, curStepError2)
    !curStepError = curStepError1 + curStepError2
    curStepError = curStepError1
!---------------------------------------------------------------------------------------------
!   !Check for drfit correction before checking tensile principal stresses
    interstress = snStress_H
    interalpha = snAlpha_H
    
    CALL GetTrace(interstress, p)
    p = one3 * p
    CALL GetF(snStress_H, snAlpha_H, Params, fnNext)

    IF (abs(fnNext)> TolE .AND. p > zero) THEN
        IF ((dT <= dT_min .AND. p < Psmall) .OR. curStepError <= TolR) THEN
            call StressCorrection300(CurStress, CurAlpha, snStress_H, snAlpha_H, snAlphaM_H, snMM_plus_H, &
                snMM_minus_H, curalpha_in, params, thisvoidratio, interstress, interalpha, indicator, iStep, iAbort, iTer, iEl, Int)
            IF (iAbort == 1) THEN
                WRITE(100,*) 'Checking tensile stress'
                call flush(100)
                RETURN
            END IF
        END IF
    END IF  
!--------------------------------------------------------------------
    call PrnSig(1,interstress,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)
    
!   !Check the higher order solution
    if (min(sNextStress1a,sNextStress2a,sNextStress3a) < -1.0D-5) then!  TODO small_negative_stress = -1.0D-5*1kPa
        izero = 0 !it is reset 
        iNotension = 1
    end if
!--------------------------------------------------------------------    

    IF ((curStepError > TolR) .and. iNotension==0) THEN
        
        IF (dT <= dT_min) THEN

            CALL GetTrace(snStress_H, p)
            p = one3 * p
            
            IF (p >= Psmall) THEN
                write(100, *) 'dt=dt_min and p>psmall --> abort -> inaccurate sol for p>psmall stress not tensional'
                write(100, *) T, dT
                write(100, *) iEl, Int, istep, iTer
                WRITE(100, 1170) 'dDevStrain = ', dDevStrain
                write(100, *) 'dVolStrain =', dVolStrain
                write(100, *) 'curStepError1 =', curStepError1
                write(100, *) 'curStepError1_stress =', curStepError1_stress
                write(100, *) 'curStepError2 =', curStepError2
                write(100, *) 'curStepError2_alpha =', curStepError2_alpha
                write(100, 1170) 'sNextStress_H = ', snStress_H
                write(100, 1170) 'sNextAlpha_H = ', snAlpha_H
                write(100, 1170) 'dAlpha1 = ', dAlpha1
                write(100, 1170) 'dAlpha2 = ', dAlpha2
                write(100, 1170) 'dAlpha3 = ', dAlpha3
                write(100, 1170) 'dSigma1 = ', dSigma1
                write(100, 1170) 'dSigma2 = ', dSigma2
                write(100, 1170) 'dSigma3 = ', dSigma3
                call flush(100)
                iAbort = 1
                return
            ELSE              
!               Use previous stress to prevent p become too small and continue iterations.
                sNextStress	  = sNextStress
                sNextAlpha    = sNextAlpha
                sNextAlphaM   = sNextAlphaM
                sNextMM_plus  = sNextMM_plus
                sNextMM_minus = sNextMM_minus
            
                sNextStress1 = sNextStress
                sNextAlpha1  = sNextAlpha

                !Check and perform the drift correction if needed
                CALL GetF(sNextStress, sNextAlpha, Params, fnNext)
                if(abs(fnNext)> TolE)then

                    call StressCorrection300(curstress, curalpha, snextstress, snextalpha, snextalpham, snextmm_plus, &
                    snextmm_minus, curalpha_in, params, thisvoidratio, snewstress, snewalpha, indicator, iStep, iAbort, iTer, iEl, Int)
                    IF (iAbort == 1) THEN
                        RETURN
                    END IF
               
                    sNextStress = sNewstress
                    sNextAlpha = sNewAlpha
                
                endif
            
                CALL GetF(sNextStress, sNextAlpha, Params, fnNext)
            
                T = T + dT

!               Update the void ratio if T=1
                res7 = T * strainInc !strainInc_at_T
                CALL GetTrace(res7, res2) !dEpsVol_at_T
                sNextVoidRatio = CurVoid - (one + CurVoid) * res2 !void ratio at the EpsVol_T

                RETURN
            END IF
        ELSE
            q = MAX(0.9D0 * ((TolR / curStepError)**(one3)), 0.25D0)
            dT = MAX(q * dT, dT_min)
            dT = MIN(dt, one - t)
        END IF
        
    ELSE IF ((curStepError < TolR) .and. iNotension==0) then
        
!       !Store the previus accurate solution with all positive prinicpal stresses and p>psmall before update
        sNextStressCopy = sNextStress
        sNextAlphaCopy  = sNextAlpha
        sNextAlphaMCopy = sNextAlphaM
        sNextMM_plusCopy= sNextMM_plus
        sNextMM_minusCopy = sNextMM_minus
        
!       !Accept the soulution
        sNextStress	  = snStress_H
        sNextAlpha    = snAlpha_H
        sNextAlphaM   = snAlphaM_H
        sNextMM_plus  = snMM_plus_H
        sNextMM_minus = snMM_minus_H
!
!-----------------------------------remove        
        sNextStress1 = sNextStress
        sNextAlpha1  = sNextAlpha        
        call PrnSig(1,sNextStress,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)
!-----------------------------------remove 

!       !Check and perform the drift correction if needed
        call GetF(sNextStress, sNextAlpha, Params, fnNext)
        if(abs(fnNext)> TolE)then
            
            call StressCorrection300(curstress, curalpha, snextstress, snextalpha, snextalpham, snextmm_plus, &
            snextmm_minus, curalpha_in, params, thisvoidratio, snewstress, snewalpha, indicator, iStep, iAbort, iTer, iEl, Int)
            IF (iAbort == 1) THEN
                RETURN
            END IF
            
            sNextStress = sNewstress
            sNextAlpha = sNewAlpha
        endif
        
!--------------------------------------------------------------------------------------remove
 !      !Calculate and output Mb and Md
        call GetNormalToYield(sNextStress, sNextAlpha, iStep, sn, iAbort)
        call GetTrace(sNextStress, p)
        p = one3 * p
        call GetDevPart(sNextStress, rr)
        call GetNorm_Contr(rr, q_Norm_New)
        call GetNorm_Contr(sNextStress, sNextStress_Norm)
        q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
        rr = rr / p
        call GetNorm_Contr(rr, rr_Norm)
        call GetNorm_Contr(sNextAlpha, sNextAlpha_Norm)
        rr_Norm = dsqrt(3.0D0/2.0D0) * rr_Norm
        sNextAlpha_Norm = dsqrt(3.0D0/2.0D0) * sNextAlpha_Norm
        call GetLodeAngle(sn, cos3Theta)
        call g(cos3Theta, PAR_c, res24, dgdth)
        call GetPSI(sNextVoidRatio, p, Params, psi)
        Mb = PAR_Mc * EXP(-PAR_nb * psi)
        if (Mb > Mbc_max) then
            Mb = Mbc_max
        end if
        
        Mb = res24 * Mb
        Mc = res24 * PAR_Mc
        Md = res24 * PAR_Mc * EXP(PAR_nd * psi)

        call GetDevPart(FrozenStress, rr)
        call GetNorm_Contr(rr, q_Norm_Old)
        call GetNorm_Contr(FrozenStress, FrozenStress_Norm)
        q_Norm_Old = dsqrt(3.0D0/2.0D0) * q_Norm_Old
!--------------------------------------------------------------------------------------end of remove
!       !Check the mean effective stress of the accepted stress after drift correction
        IF (p <= Psmall) THEN
        !    IF (pFrozen >= 0.01D0) THEN
        !!      sNextStress   = FrozenStress
        !!      sNextAlpha    = FrozenAlpha
        !!      sNextAlphaM   = FrozenAlphaM
        !!      sNextMM_plus  = FrozenMM_plus
        !!      sNextMM_minus = FrozenMM_minus
        !       return
        !    ELSE
        !       sNextStress   = FrozenStress
        !       sNextAlpha    = FrozenAlpha
        !       sNextAlphaM   = FrozenAlphaM
        !       sNextMM_plus  = FrozenMM_plus
        !       sNextMM_minus = FrozenMM_minus
            
!               Overwrite and accept the last updated to allow progress of the iterative process !             
                sNextStress = sNextStressCopy
                sNextAlpha = sNextAlphaCopy
                sNextAlphaM = sNextAlphaMCopy
                sNextMM_plus = sNextMM_plusCopy
                sNextMM_minus = sNextMM_minusCopy
                
!           ! debugging output (switch to true if needed)----------------------------------------
            if (.false.) then
                write(100, *) 'return curStepError < TolR) .and. iNotension==0 and p<=psmall -> accept previous'
                write(100, *) T, dT
                write(100, *) sNextStress1a, sNextStress2a, sNextStress3a
                write(100, *) iEl, Int, istep, iTer
                call GetDevPart(interstress, rr)
                call GetNorm_Contr(rr, q_Norm_New)
                call GetNormalToYield(interstress, interAlpha, iStep, sn, iAbort)
                if (iabort == 1) then
                    return
                endif               
                call GetLodeAngle(sn, cos3Theta)
                call GetLodeAngleForStressRatio(rr, Cos3ThetaN,iabort)
                call g(cos3Theta, PAR_c, res24, dgdth)
                call GetTrace(interstress, p)
                p = one3 * p
                q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New/p
                call GetPSI(sNextVoidRatio, p, Params, psi)
                Mb = res24 * PAR_Mc * EXP(-PAR_nb * psi)
                !    Mc = res24 * PAR_Mc
                !    Md = res24 * PAR_Mc * EXP(PAR_nd * psi)
                write(100, *) "cos3tetaN", cos3ThetaN
                write(100, *) "cos3teta", cos3Theta
                write(100, *) "p", p
                write(100, *) "eta", q_Norm_New
                write(100, *) "Mb", Mb
                write(100, *) "sNextVoidRatio", sNextVoidRatio
                write(100, *) "psi", psi
                write(100, 1170) 'dDevStrain = ', dDevStrain
                write(100, *) 'dVolStrain =', dVolStrain
                write(100, *) 'curStepError1 =', curStepError1
                write(100, *) 'curStepError1_stress =', curStepError1_stress
                write(100, *) 'curStepError2 =', curStepError2
                write(100, *) 'curStepError2_alpha =', curStepError2_alpha
                write(100, 1170) 'sNextStress_H = ', snStress_H
                write(100, 1170) 'sNextAlpha_H = ', snAlpha_H
                write(100, 1170) 'dAlpha1 = ', dAlpha1
                write(100, 1170) 'dAlpha2 = ', dAlpha2
                write(100, 1170) 'dAlpha3 = ', dAlpha3
                write(100, 1170) 'dSigma1 = ', dSigma1
                write(100, 1170) 'dSigma2 = ', dSigma2
                write(100, 1170) 'dSigma3 = ', dSigma3
                call flush(100)
            end if
!           !--------------------------------End of debugging output----------------------------------------
            return
        end if
        
!----------------------------------------------------------------------------------------
        aCep = 0.0d0
!----------------------------------------------------------------------------------------
        q = min(0.9D0 * ((TolR / curStepError)**(one3)), 4.0D0)
        T = T + dT

!       Update of the void ratio if T=1        
        res7 = T * strainInc !strainInc_at_T
        call GetTrace(res7, res2) !dEpsVol_at_T
        sNextVoidRatio = CurVoid - (one + CurVoid) * res2 !void ratio at the EpsVol_T
        
        dT = max(q * dT, dT_min)
        dT = min(dT, one - T)

    else
            q=0.5
            if (dT <= dT_min) THEN
!               !Compute and write debug output
                call GetF(snStress, snAlpha, Params, fnNext)
                call GetDevPart(interstress, rr)
                call GetNorm_Contr(rr, q_Norm_New)
                call GetNormalToYield(interstress, interAlpha, iStep, sn, iAbort)
                if (iabort == 1) then
                    return
                endif               
                call GetLodeAngle(sn, cos3Theta)
                call GetLodeAngleForStressRatio(rr, Cos3ThetaN,iabort)
                call g(cos3Theta, PAR_c, res24, dgdth)
                call GetTrace(interstress, p)
                p = one3 * p
                q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New/p
                call GetPSI(sNextVoidRatio, p, Params, psi)
                Mb = res24 * PAR_Mc * EXP(-PAR_nb * psi)

                write(100, 1170) 'snStress (in T=0)= ', snStress
                write(100, 1170) 'snAlpha (in T=0)= ', snAlpha
                write(100, *) 'f (in T=0)= ', fnNext
                write(100, *) 'Tensional stress, dT <= dT_min'
                write(100, *) sNextStress1a, sNextStress2a, sNextStress3a
                write(100, *) T, dT
                write(100, *) iEl, Int, istep, iTer
                write(100, *) "cos3tetaN", cos3ThetaN
                write(100, *) "cos3teta", cos3Theta
                write(100, *) "p", p
                write(100, *) "eta", q_Norm_New
                write(100, *) "Mb", Mb
                write(100, *) 'curStepError1 =', curStepError1
                write(100, *) 'curStepError1_stress =', curStepError1_stress
                write(100, *) 'curStepError2 =', curStepError2
                write(100, *) 'curStepError2_alpha =', curStepError2_alpha
                write(100, 1170) 'sNextStress_H = ', snStress_H
                write(100, 1170) 'sNextAlpha_H = ', snAlpha_H
                write(100, 1170) 'dAlpha1 = ', dAlpha1
                write(100, 1170) 'dAlpha2 = ', dAlpha2
                write(100, 1170) 'dAlpha3 = ', dAlpha3
                write(100, 1170) 'dSigma1 = ', dSigma1
                write(100, 1170) 'dSigma2 = ', dSigma2
                write(100, 1170) 'dSigma3 = ', dSigma3
                call flush(100)
                iAbort = 1
                !iAbort = 220
                return
                
            else
                dT = max(q * dT, dT_min)
                dT = min(dt, one - t)
            end if
        
    end if
end do

return
end

!---------------------------------------------------------------------------
SUBROUTINE NaNcheck(chcknum, itestnan)
!     NaN check function
implicit none

double precision, intent(in) :: chcknum
integer, intent(out) :: itestnan

double precision, parameter :: zero = 0.0D0

IF (.not.((chcknum >= zero) .OR. (chcknum < zero))) THEN
    itestnan = 1
END IF
IF (chcknum -1 == chcknum) THEN
    itestnan = 1
END IF
!     IF (chcknum > 1.d30)   itestnan = 1
!     IF (chcknum < -1.d30)  itestnan = 1
IF (chcknum /= chcknum) THEN
    itestnan = 1
END IF

RETURN
END
!---------------------------------------------------------------------------
SUBROUTINE Hcheck(chcknum, itestnan)
!     NaN check function for h,hM (real* 16)
implicit none
REAL*16, intent(in) :: chcknum
integer, intent(out) :: itestnan

double precision, parameter :: zero = 0.0D0

IF (.not.((chcknum >= zero) .OR. (chcknum < zero))) THEN
    itestnan = 1
END IF
IF (chcknum -1 == chcknum) THEN
    itestnan = 1
END IF
!     IF (chcknum > 1.d30)   itestnan = 1        
!     IF (chcknum < -1.d30)  itestnan = 1        
IF (chcknum /= chcknum) THEN
    itestnan = 1        
END IF

RETURN
END
!---------------------------------------------------------------------------   
SUBROUTINE tension_cutoff( sig, Alpha )

use extrainfo, only: onekPa

implicit none

integer ::  ndime, i

logical :: ten_cutoff

double precision :: p, q, pThresh
! 
double precision, dimension(3)  :: xN1, xN2, xN3, Prs
double precision, dimension(6)  :: sig
double precision, dimension(6)  :: Alpha
! 
ndime = 6
! 
ten_cutoff = .false.
!
pThresh = 0.01d0*onekPa !change it
! 
! -------------------------------------------------------------------
CALL PrnSig( 1, sig, xN1, xN2, xN3, Prs(1), Prs(2), Prs(3), p, q )
! -------------------------------------------------------------------
do i = 1, 3
! 
   if ( Prs(i) < pThresh ) then
      Prs(i) = pThresh
      ten_cutoff = .true.
   end if
!   
end do
!---------------------------------------------------------
if ( ten_cutoff.eqv..true. ) then
!
   CALL CarSig( Prs(1), Prs(2), Prs(3), xN1,xN2,xN3,sig )
!    
end if
!---------------------------------------------------------
! 
return

END SUBROUTINE tension_cutoff    
!---------------------------------------------------------------------------

SUBROUTINE bisectionOfYSCrossing(CurStress, CurAlpha,  Params, sNewstress,  indicator, iStep, iAbort, iTer, iEl, Int)

use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: CurStress(6), CurAlpha(6), Params(25)
integer, intent(in) :: indicator, iStep,  iTer, iEl, Int
integer, intent(out) :: iabort
double precision, intent(out) :: sNewstress(6)

integer i, it

!double precision:: beta0 = 0.0d0
!double precision:: beta1 = 1.0d0
double precision:: beta, beta0, beta1
double precision:: Fc, Fa
double precision:: p0, s0(6), pa, sa(6)
double precision:: s1(6)
double precision:: pC, sC(6), TolF, stressA(6), stressC(6), Fcs, Fas, F, dummy_stress_in_alfa(6)

TolF = Params(i_TolF)
beta0 = 0.0d0
beta1 = 1.0d0

iabort = 0

CALL GetTrace(CurStress, p0)
p0 = 1.0d0/3.0d0*p0
CALL GetDevPart(CurStress, s0)
!write(1,*),'s0',s0
s1 = CurAlpha * p0
!write(1,*),'s1',s1
!write(1,*),  'beta0', beta0
!write(1,*), 'beta1', beta1

pa = p0
pC = p0
!------------------------------------------
!dummy_stress_in_alfa = s1
!dummy_stress_in_alfa(1) = s1(1) + p0
!dummy_stress_in_alfa(2) = s1(2) + p0
!dummy_stress_in_alfa(3) = s1(3) + p0
!CALL GetF(dummy_stress_in_alfa, CurAlpha, Params, F)
!write(1,*), "dummy F ", F
!---------------------------------------------
it = 0;
CALL GetF(CurStress, CurAlpha, Params, F)
!write(2,*), "BISECTION for ys drift ", F
!write(1,*), CurStress
!write(1,*), CurAlpha
!write(1,*), p0

do while (it < 50) 
    beta = (beta0 + beta1) / 2.0; !new midpoint
!write(1,*),it,  'beta', beta
!write(1,*),it, 'beta0', beta0
!write(1,*),it, beta0, 'beta1', beta1
    sa = (s1 - s0) * beta0 + s0;
    stressA = sa
    !write(1,*), 'stressA', stressA
    stressA(1) = stressA(1) + p0
    stressA(2) = stressA(2) + p0
    stressA(3) = stressA(3) + p0
    !write(1,*), 'stressA', stressA

    sC = (s1 - s0) * beta + s0;
    stressC = sC
    stressC(1) = stressC(1) + p0
    stressC(2) = stressC(2) + p0
    stressC(3) = stressC(3) + p0
    !write(1,*), 'stressC', stressC

    CALL GetF(stressA, CurAlpha, Params, Fa)
    CALL GetF(stressC, CurAlpha, Params, Fc)

    !write(1,*), "BISECTION for ys drift ", it, beta, Fc, pC, sC;

    if (abs(Fc) <= TolF) then
      !write(1,*), 'beta found! ', beta,  Fc, stressC
      sNewstress = stressC
      return
    else
        
      it=it+1
      !write(1,*), 'iter! ', it,  Fc, Fa
      !if (sign(Fc,Fcs) == sign(Fa,Fas)) then
      if ((Fc>0.0d0 .and. Fa>0.0d0) .or. (Fc<0.0d0 .and. Fa<0.0d0)) then
          !write(1,*), 'iter! ', it,  Fcs, Fas
          beta0 = beta
      else
          beta1 = beta
      end if
      
      if (it == 49) then
       write(100,*) , 'bisection failed Beta=', beta, Fc
       WRITE(100, *) iEl, Int, iStep, iTer
       call flush(100)
       !iAbort = 1 
       !return        
      end if
      
    end if
    
    
end do
  
END SUBROUTINE
!---------------------------------------------------------------------------   

SUBROUTINE     StressCorrection300(PreStress, PreAlpha, CurStress, CurAlpha, CuralphaM, Cur_MM_plus, Cur_MM_minus,   &
    Curalpha_in, Params, e, sNewstress, sNewAlpha, indicator, iStep, iAbort, iTer, iEl, Int)
!
!
! Purpose: perform the consistent and normal drift correction based on Sloan et al.2021
!       
! Arguments:
!                  I/O   Type
!  iEl              I    I    :
!  Int              I    I    :
!  iStep            I    I    :
!  iTer             I    I    :
!  indicator        I    I    :
!  e                I    R    :
!  Params           I    R()  :
!  PreStress        I    R()  : At the beginning of the elastoplastic integration
!  PreAlpha         I    R()  : At the beginning of the elastoplastic integration
!  CurStress        I    R()  : Uncorrected stress at the end of the accepted substep
!  CurAlpha         I    R()  : Uncorrected alpha at the end of the accepted substep
!  CurAlphaM        I    R()  : 
!  Cur_MM_plus      I    R()  : Current theoric size of the Memory that would result from modifications only due to the hardening 
!  Cur_MM_minus     I    R()  :
!  Curalpha_in      I    R()  : 
!  iAbort           O    I    : 
!  sNewstress       O    R()  : Corrected stress at the end of the accepted substep
!  sNewAlpha        O    R()  : Corrected alpha at the end of the accepted substep

use SANISANDMS_constants
use Parameters_position_in_params_iMod1
implicit none

integer, intent(in   ) :: indicator,iStep, iTer, iEl, Int
integer, intent(  out) :: iAbort
double precision, intent(in   ) :: &
PreStress(6), PreAlpha(6), CurStress(6), CurAlpha(6), CuralphaM(6), Cur_MM_plus, Cur_MM_minus, Curalpha_in(6), Params(25), e
double precision, intent(  out) :: sNewstress(6), sNewAlpha(6)

double precision sInterstress(6), sInterAlpha(6), SI1(6), rr(6), f0, f1, f2
double precision sn(6), dd(6), bb(6), bM(6), cos3Theta, psi, rBtheta, rDtheta, b0, A, B, C, D, R(6), Z, PM, &
    E_K, E_G, temp1(6), temp2(6), temp3, dlambda, temp4(6), TolF, small, res(6), res2, p, aC(6,6), temp5, h, hM, df_ds(6),  &
    M_alfa, xN1(3), xN2(3), xN3(3), sNextStress1a, sNextStress2a, sNextStress3a, PNextStress, QNextStress,PAR_m,  &
    TRACE_RR, TRACE_ALFA, rr_alfa(6), norm_r_a, ss(6), AlphaAlphaInDotN, NormAlphaAlphaInDotN, AlphaAlphaInDotNNorm, &
    sNewstressBisection(6), p_sNewAlpha

!REAL*16 h, hM
integer i,switch

!double precision, parameter :: two3 = 2.0d0/3.0d0
double precision, parameter :: one_two = 1.0d0/2.0d0
double precision, parameter :: sqrt_two3 = DSQRT(two3)
double precision, parameter :: small_adim = 1.0D-12

iAbort = 0
switch = 0
TolF = Params(i_TolF)
PAR_m = Params(i_PAR_m)

SI1(1) = 1.0d0
SI1(2) = 1.0d0
SI1(3) = 1.0d0
SI1(4) = 0.0d0
SI1(5) = 0.0d0
SI1(6) = 0.0d0
small  = 1.0D-7
i = 0

sInterstress = CurStress
sInterAlpha  = CurAlpha

CALL GetF(PreStress, PreAlpha, Params, f2)

!DO WHILE (i <= 10)
DO WHILE (i <= 50)
    
    if(switch==0) then
        CALL GetTrace(sInterstress, p)
        p = 1.0d0/3.0d0*p

        CALL GetDevPart(sInterstress, rr)
        rr = rr / p
! only for debugging----------------------------------------------------------
        TRACE_RR = rr(1) + rr(2) + rr(3)
        TRACE_ALFA = sInterAlpha(1) + sInterAlpha(2) + sInterAlpha(3)
        !WRITE(1,*), 'TRACE_RR' , TRACE_RR
        !WRITE(1,*), 'TRACE_ALFA', TRACE_ALFA
!-----------------------------------------------------------------------------
        CALL GetF(sInterstress, sInterAlpha, Params, f0)

        CALL GetStateDependent(sInterstress, sInterAlpha, CuralphaM, Cur_MM_plus,  &
         Cur_MM_minus, Curalpha_in, Params, e, sn, dd, &
         bb, bM, cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, &
         B, C, D, R, Z, sInterstress, indicator, iStep, iAbort, iTer, iEl, Int, AlphaAlphaInDotN,   &
         NormAlphaAlphaInDotN, AlphaAlphaInDotNNorm)

        CALL DoubleDot2_2_Contr(bb, sn, PM)
        PM = two3 * p * h * PM

        ! 24-11-2021
        CALL GetElasticModuli(sInterstress, e, Params, E_K, E_G)
        CALL GetStiffness(E_K, E_G, aC)

        CALL SingleDot(sn, sn, temp1)
!     dg/dsig
        res = (B * sn) - (C * (temp1 - 1.0d0/3.0d0 * SI1)) + (1.0d0/3.0d0 * D * SI1)
        CALL ToCovariant(res, R)

!     De * (dg/dsig)
        temp1 =  MATMUL(aC, R)

!     df/dsig
        CALL DoubleDot2_2_Contr(sn, rr, res2)
        res = sn - 1.0d0/3.0d0 * res2 * SI1
!     dg/dSig * df/dsig
        CALL ToCovariant(res, temp2)
!     (df/dsig) * De
        temp2 = MATMUL(temp2, aC)

!     (df/dsig) * De * (dg/dsig) + KP
        CALL DoubleDot2_2_Contr(temp2, R, temp3)
        temp3 = temp3 + PM

        dlambda = f0/temp3

!     De * (dg/dsig)
!       CALL DoubleDot2_2_Contr(aC, R, temp4)
        temp4 = MATMUL(aC, R)
        sNewstress = sInterstress - dlambda * temp4
        sNewAlpha  = sInterAlpha + dlambda * two3 * h * bb
        
        CALL GetF(sNewstress, sNewAlpha, Params, f1)

!---- The procedure succeded: return with the stress and hardening lying on the yield    
        IF (abs(f1) < TolF) THEN
            return
        END IF
!---- Abandon previous correction and change method
        IF (abs(f1) > abs(f2)) THEN
            switch = 1
!---- Continue correcting
        else
            sInterstress = sNewstress
            sInterAlpha = sNewAlpha
        END IF
        
    else
                
!     (df/dsig) * De * (df/dsig)

        CALL DoubleDot2_2_Contr(temp2, res, temp5)
        dlambda = f0 / temp5
        sNewstress = sInterstress - dlambda * res
        sNewAlpha  = sInterAlpha
        CALL GetF(sNewstress, sNewAlpha, Params, f1)
        
        sInterstress = sNewstress
        sInterAlpha = sNewAlpha
        
        CALL GetF(sNewstress, sNewAlpha, Params, f0)

!---- The alternative procedure succeded: return with the stress and hardening lying on the yield 
        IF (abs(f0) < TolF) THEN
            return
        END IF
        
    end if

!- Both procedures failed
    if (i == 50) THEN
        if (f1> 0.0d0) then
            !the bisection cannot work if f1 < 0
            call bisectionOfYSCrossing(sNewstress, sNewAlpha,  Params, sNewstressBisection,  indicator, iStep, iAbort, iTer, iEl, Int)
            !if (iabort == 1)then
            !    return
            !endif
            !return
            CALL GetF(sNewstressBisection, sNewAlpha, Params, f0)
            if (f0 < tolf) then
                sNewstress = sNewstressBisection
                return
            end if
            
        end if
!- Perform drastic correction because f1 negative or bisectionOfYSCrossing failed
!------------------------------------------------------------------------------------------------
        call PrnSig(1,sNewAlpha,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)
        !PNextStress = sNewAlpha(1) + sNewAlpha(2) +sNewAlpha(3)
        !M_alfa= QNextStress
        !CALL GetNormalToYield(sNewstress, sNewAlpha, iStep, sn, iAbort)
!------------------------------------------------------------------------------------------------remove
        IF (iAbort == 1) then
            return
        end if
!------------------------------------------------------------------------------------------------move after GetNormalToYeld        
        CALL GetTrace(sNewstress, p) ! 07-03-2022 replace sInterstress with sNewstress
        CALL GetTrace(sNewAlpha, p_sNewAlpha) ! 07-03-2022 check trace of alpha and correct        
        call RestoreDviatoricProperty(sNewAlpha,p_sNewAlpha)
        
        CALL GetNormalToYield(sNewstress, sNewAlpha, iStep, sn, iAbort)
        
        p = 1.0d0/3.0d0*p
        rr = (sNewAlpha + sqrt_two3 * sn * PAR_m)
        ss = rr * p
        rr_alfa = rr - sNewAlpha

        sNewstress = rr * p 
        sNewstress(1) = sNewstress(1) + p
        sNewstress(2) = sNewstress(2) + p
        sNewstress(3) = sNewstress(3) + p

        CALL GetF(sNewstress, sNewAlpha, Params, f1)
        
        if(f1>TolF)then
            write(100, *) '    extreme correction fails after drift correction!', istep, f1
            write(100, '(4A)') "iEl", "Int", "iStep", "iTer"
            write(100, *) iEl, Int, iStep, iTer
            write(100, '(A, 6F35.10)'), 'sNewstress = ', sNewstress(1), sNewstress(2), sNewstress(3),&
                                                        sNewstress(4), sNewstress(5), sNewstress(6)
            write(100, '(A, 6F35.10)'), 'sNewAlpha = ', sNewAlpha(1), sNewAlpha(2), sNewAlpha(3),&
                                                      sNewAlpha(4), sNewAlpha(5), sNewAlpha(6)
            write(100, '(A, 6F35.10)'), 'rr = ', rr(1), rr(2), rr(3),&
                                            rr(4), rr(5), rr(6)
            write(100, '(A, 1F35.10)'), 'p = ', p
            write(100, '(A, 1F35.10)'), 'p_sNewAlpha = ', p_sNewAlpha
            call flush(100)
            iAbort = 1
            !iAbort = 350
            return
        endif

    end if
    i = i + 1

end do
return
end subroutine

!-------------------------------------------------------------------------------------   
Subroutine RestoreDviatoricProperty(Alpha, trace_Alpha)

!If needed, it corrects the back stress ratio components to have tr = 0

implicit none
double precision, intent(in) :: trace_Alpha
double precision, intent(inout) :: Alpha(6) 
    
double precision, parameter :: small_adim = 1.0D-12

if (DABS(trace_Alpha) > small_adim) then
    Alpha(1) = Alpha(1) - trace_Alpha / 3.0d0
    Alpha(2) = Alpha(2) - trace_Alpha / 3.0d0
    Alpha(3) = Alpha(3) - trace_Alpha / 3.0d0 
end if

End Subroutine 
!-------------------------------------------------------------------------------------   
SUBROUTINE DerivativeDfDSig(stress, alpha, Params, iStep, DfDSig, iAbort)
!
! Purpose: It compute the derivativative of the yield surface respect to sigma
!             DfDSig = n - 1/3 * (alfa:n + sqrt(2/3) m ) I
!             DfDSig = n - 1/3 * (alfa:n + sqrt(2/3) m n:n ) I
!             DfDSig = n - 1/3 * (r:n ) I
! Arguments:
!          I/O  Type
!  stress   I   R()    :
!  alpha    I   R()    :
!  Params   I   r()    : 
!  iStep    I   I      :
!  iAbort   I   I      : 
!  DfDSig   O   R()    :
!
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) ::  stress(6), alpha(6), Params(25)
integer, intent(in) ::  iStep
integer, intent(out) :: iAbort
double precision, intent(out) :: DfDSig(6)

double precision :: sn(6), rr(6), res(6), SI1(6), p, res2, res3, PAR_m

iAbort = 0
PAR_m = Params(i_PAR_m)

SI1(1) = 1.0d0
SI1(2) = 1.0d0
SI1(3) = 1.0d0
SI1(4) = 0.0d0
SI1(5) = 0.0d0
SI1(6) = 0.0d0

CALL GetNormalToYield(stress, alpha, iStep, sn, iAbort)
if (iabort == 1) then
    write(100,*), 'normal to the yield in DfDsig'
    call flush(100)
    return
endif


!CALL GetDevPart(stress, rr)
!call GetTrace(stress, p)
!rr = rr / p
!
!!df/dsig
!CALL DoubleDot2_2_Contr(sn, rr, res2)
!DfDSig = sn - 1.0d0/3.0d0 * res2 * SI1

! considering the order of updating of the quantities
CALL DoubleDot2_2_Contr(alpha, sn, res3)
DfDSig = sn - 1.0d0/3.0d0 * (res3  + dsqrt(2.0d0/3.0d0)* PAR_m)* SI1

RETURN
END
!-------------------------------------------------------------------------------------   
  
Subroutine FindIO(ioMin,io)
implicit none
integer ioMax, ioMin, io
Logical IsOpen
Parameter (ioMax=99)
io=ioMin-1
1 io=io+1
Inquire(io,Opened=IsOpen)
!        Write(2,*)'unit ',io,' open? ',IsOpen
If (io.Lt.ioMax .And. IsOpen) Goto 1 ! try next unit
!      Write(2,*)'Resulting free unit ',io

Return
End
   

!------------------------------------------------------
SUBROUTINE GetLodeAngleForStressRatio(rr, Cos3ThetaN,iabort)
!------------------------------------------------------
! Purpose: Subroutine to calculate lode angle
!------------------------------------------------------
implicit none
integer, intent(  out) :: iAbort
double precision, intent(in   )   :: rr(6)
double precision, intent(  out)   :: Cos3ThetaN


double precision, parameter :: tol = 1.0D-15
double precision, parameter :: rootsix = DSQRT(6.0D0)
double precision, parameter :: one = 1.0D0
double precision res1(6), res(6), J3D, J2D, den

call SingleDot(rr, rr, res)
call SingleDot(rr, res, res1)

call GetTrace(res1, J3D)
J3D = 1.0d0/3.0d0 * J3D

call DoubleDot2_2_Contr(rr, rr, J2D)
J2D = 0.5d0 * J2D
den = DSQRT(J2D*J2D*J2D)

if(abs(den)>1.0d-10)then
    Cos3ThetaN = 3.0D0 * DSQRT(3.0D0) * 0.5d0 * J3D/den
else
 iabort=1
 write(100,*) , "denom too small in GetLodeAngleForStressRatio"
 return
endif

IF (Cos3ThetaN > one + tol) THEN !tolerance
    !write(100,*), 'problematic cos', Cos3Theta
!    write(100,*), 'sn', sn
    Cos3ThetaN = one
END IF

IF (Cos3ThetaN < -one - tol) THEN
        !write(100,*), 'problematic cos2', Cos3Theta
!        write(100,*), 'sn', sn
    Cos3ThetaN = -one
END IF

RETURN
END

!-------------------------------------------------------------------------------------   
SUBROUTINE EL_RungeKutta23(CurStress, CurVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, strainInc, Params, &
                            sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoidRatio, aCep, &
                            IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
!
! Purpose: Update the stress tensor 
!       
!
! Arguments:
!                  I/O   Type
!  IDTask           I    I()   :
!  iEl              I    I()   :
!  Int              I    I()   :
!  iStep            I    I()   :
!  iTer             I    I()  :
!  indicator        I    I()  :
!  strainInc        I    R() 
!  Params           I    R()  :
!  CurStress        I    R()  : At the beginning of the elastoplastic integration
!  CurVoid          I    R()  :
!  CurAlpha         I    R()  :
!  CurAlphaM        I    R()  : 
!  CurMM_plus       I    R()  : Current theoric size of the Memory that would result from modifications only due to the hardening 
!  CurMM_minus      I    R()  :
!  Curalpha_in      I    R()  : NB the next alpha_ini has been already updated in integrator
!  sNextStress      O    R()  : At the accepted substep
!  sNextVoid        O    R()  :
!  sNextAlpha       I    R()  :
!  sNextAlphaM      I    R()  :
!  sNextMM_plus     I    R()  : 
!  sNextMM_minus    I    R()  :
!  aCep             I    R()  :
!  PLNoTension      I    R    : dismissed parameter used for debugging to track increment of strain computed without respecting the error tolerance
!  iAbort           O    i    :
!  PsPsmall         IO   R    : dismissed parameter used for debugging as warning for a mean effective stress smaller than the minimum
!
use Parameters_position_in_params_iMod1
use extrainfo, only: onekPa
use extrainfo, only: oneAtm
use SANISANDMS_constants

implicit none
integer, intent(in   ) :: IDTask, iEl, Int, iStep, indicator, iTer
integer, intent(  out) :: iAbort
double precision, intent(in  ) :: CurStress(6), CurVoid, CurAlpha(6),  CurAlphaM(6), CurMM_plus, CurMM_minus, Curalpha_in(6), &
    strainInc(6), Params(25), sNextAlpha(6), sNextAlphaM(6), sNextMM_plus, sNextMM_minus, aCep(6,6)
double precision, intent(  out) :: sNextStress(6), sNextVoidRatio, PLNoTension
double precision, intent(inout) :: PsPsmall

double precision sn(6), dDevStrain(6), rr(6), sI1(6), snStress(6), snAlpha(6), snAlphaM(6), sndPStrain(6),   &
    dSigma1(6), dSigma2(6), dSigma3(6), dSigma(6),aC(6, 6), aD(6, 6), res(6), thisSigma(6), thisAlpha(6), thisAlphaM(6),   &
    sInterstress(6), dSigma_H(6), dAlpha_H(6), dAlphaM_H(6), snStress_H(6), snAlpha_H(6), snAlphaM_H(6),   &
    PAR_m, TolE, TolR, PAR_e_init, T, dT, dT_min, p, thisVoidRatio, thisMM_plus, stressNorm_H,  &
    res2, dVolStrain,  depsilon_pv, curStepError1, alphaNorm_H,  &
    stressNorm, q, curStepError, thisMM_minus, E_G, E_K, res7(6), strainIncNew(6), sNextStress1(6),  &
    fnNext, sNewstress(6), sNewAlpha(6), &
    xN1(3),xN2(3),xN3(3),sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress,this_trace_Alpha, sn_trace_Alpha, &
    this_trace_AlphaM, sn_trace_AlphaM, curStepError1_stress, cos3Theta, res24, dgdth, psi, Mb, Md, &
    PAR_Mc, PAR_c, PAR_nb, PAR_nd, Mc, sNextAlpha_Norm, rr_Norm, FrozenStress(6),  &
    q_Norm_New, q_Norm_Old, FrozenStress_Norm, sNextStress_Norm, emax, emin, pFrozen, &
    Psmall, sNextStressCopy(6), small_NormSig, interstress(6), interalpha(6), Mbc_max, Cos3ThetaN

integer i, izero, SNum_Inter, iNotension

double precision, parameter :: small = 1.0D-7! small Norm alpha for Err alpha
double precision, parameter :: tol = 1.0D-18
double precision, parameter :: small_adim = 1.0D-12

double precision PAR_P_atm
double precision PAR_Pmin

1170  FORMAT(A, 6F28.18)
1180  FORMAT(A, 1F28.18)
1200  FORMAT(2I, 10F15.10)
1201  FORMAT(16F15.10)
1202  FORMAT(17F15.10)
1203  FORMAT(2F30.15)

iAbort = 0

Psmall = 1.0D-2 * onekPa! when dT is < dTmin
small_NormSig = 0.5D0 * onekPa! for Error Sig 

PAR_Pmin = onekPa/100d0
PAR_P_atm = oneAtm
PsPsmall = zero !deprecated variable
PLNoTension = zero !deprecated variable

emax        = Params(i_emax)
emin        = Params(i_emin)
PAR_Mc      = Params(i_PAR_Mc)
PAR_c       = Params(i_PAR_c)
PAR_m       = Params(i_PAR_m)
PAR_nb      = Params(i_PAR_nb)
PAR_nd      = Params(i_PAR_nd)
TolE        = Params(i_TolF)
TolR        = Params(i_TolR)
PAR_e_init	= Params(i_PAR_e_init)
Mbc_max     = Params(i_Mbc_max)

T = zero
dT = one
dT_min = 1.0D-8

sI1(1) = one
sI1(2) = one
sI1(3) = one
sI1(4) = zero
sI1(5) = zero
sI1(6) = zero

CALL GetElasticModuli1(CurStress, CurVoid, Params, E_K, E_G)

CALL GetStiffness(E_K, E_G, aC)

sNextStress   = CurStress

FrozenStress   = sNextStress

CALL GetTrace(sNextStress, p)
p = one3 * p

CALL GetDevPart(sNextStress, rr)
rr = rr / p

CALL GetTrace(FrozenStress, pFrozen)
pFrozen = one3 * pFrozen

SNum_Inter = 0


DO WHILE (T < one)

    iNotension = 0!     currently, it will be switched to 1 if principal stress negative [<-10^-10] for the final solution or if izero = 1 (if iNotension = 0 everything is fine)
    izero = 0!         currently, it will be switched to 1 if a principal stress component is negative (for the partial solutions)
    SNum_Inter = SNum_Inter + 1
    
!   ! calculate strainInc until T (for calculating the void ratio at T)
    res7 = T * strainInc
    
!   ! calculate  !dEpsVol until T
    CALL GetTrace(res7, res2)
    
!   ! calculate  void ratio    
    sNextVoidRatio = CurVoid - (one + CurVoid) * res2
    
!   ! limit the voidratio to emax and emin
    IF (sNextVoidRatio >= emax) THEN
        sNextVoidRatio = emax
    END IF
    IF (sNextVoidRatio <= emin) THEN
        sNextVoidRatio = emin
    END IF

!   ! calculate dEpsVol
    CALL GetTrace(strainInc, res2)
    
!   ! calculate  d_EpsVol_SS
    dVolStrain = dT * res2
    
!   ! calculate d_strainInc_SS()
    strainIncNew = dT * strainInc
    
!   ! calculate dev_d_strainInc_SS()
    CALL GetDevPart(strainIncNew, dDevStrain)

    depsilon_pv = zero   !remove it?

! Calc Delta 1    
    thisSigma = sNextStress
    thisAlpha = sNextAlpha
    thisAlphaM = sNextAlphaM
    thisMM_plus = sNextMM_plus
    thisMM_minus = sNextMM_minus
    thisVoidRatio = sNextVoidRatio
    
! Correction to restore deviatoric property-------------------------
    this_trace_Alpha = thisAlpha(1) + thisAlpha(2) + thisAlpha(3)
    call RestoreDviatoricProperty(thisAlpha,this_trace_Alpha)  

    this_trace_AlphaM = thisAlphaM(1) + thisAlphaM(2) +thisAlphaM(3)
    call RestoreDviatoricProperty(thisAlphaM,this_trace_AlphaM)  
! check of eventually illegal intermediate stress-----------------------------
    CALL GetTrace(thisSigma, p)
    p = one3 * p
   
    call PrnSig(1,thisSigma,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)    
    if (min(sNextStress1a,sNextStress2a,sNextStress3a)<zero) then
        izero = 1

    end if
!--------------------------------------------------------------------------      
    CALL GetDevPart(thisSigma, rr)
    rr = rr / p
    
    CALL GetElasticModuli1(thisSigma, thisVoidRatio, Params, E_K, E_G)
    CALL GetStiffness(E_K, E_G, aC)
    
    dSigma1 = MATMUL(aC, strainIncNew)

    IF (iAbort == 1) THEN
        WRITE(100, *) 'RKF1_EL', T, dT
        call flush(100)
        RETURN
    END IF

! Calc Delta 2
    thisSigma = sNextStress + (zeroP5 * dSigma1)
    thisAlpha = sNextAlpha
    thisAlphaM = sNextAlphaM
    thisMM_plus	= sNextMM_plus
    thisMM_minus = sNextMM_minus
    
! correction to restore deviatoric property---------------------------------------------------------
    this_trace_Alpha = thisAlpha(1) + thisAlpha(2) + thisAlpha(3)
    call RestoreDviatoricProperty(thisAlpha,this_trace_Alpha)
    
    this_trace_AlphaM = thisAlphaM(1) + thisAlphaM(2) +thisAlphaM(3)
    call RestoreDviatoricProperty(thisAlphaM,this_trace_AlphaM)
! check of eventually illegal intermediate stress-----------------------------
    call PrnSig(1,thisSigma,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)
    
    if (min(sNextStress1a,sNextStress2a,sNextStress3a)<zero) then
        izero = 1
    end if
!--------------------------------------------------------------------------       
    CALL GetTrace(thisSigma, p)
    p = one3 * p

    CALL GetDevPart(thisSigma, rr)
    rr = rr / p
        
    CALL GetElasticModuli1(thisSigma, thisVoidRatio, Params, E_K, E_G)
    CALL GetStiffness(E_K, E_G, aC)
    
    dSigma2 = MATMUL(aC, strainIncNew)
    IF (iAbort == 1) THEN
        WRITE(100, *) 'RKF2_EL', T, dT
        call flush(100)
        RETURN
    END IF

! Calc Delta 3    
    thisSigma = sNextStress - dSigma1 + two * dSigma2
    thisAlpha = sNextAlpha 
    thisAlphaM = sNextAlphaM 
    thisMM_plus = sNextMM_plus 
    thisMM_minus = sNextMM_minus 
        
! correction to restore deviatoric property
    this_trace_Alpha = thisAlpha(1) + thisAlpha(2) + thisAlpha(3)
    call RestoreDviatoricProperty(thisAlpha, this_trace_Alpha)
    
    this_trace_AlphaM = thisAlphaM(1) + thisAlphaM(2) +thisAlphaM(3)
    call RestoreDviatoricProperty(thisAlphaM,this_trace_AlphaM)
! check of eventually illegal intermediate stress-----------------------------    
    call PrnSig(1,thisSigma,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)

    if (min(sNextStress1a,sNextStress2a,sNextStress3a)<zero) then
        izero = 1
    end if
!--------------------------------------------------------------------------       
    CALL GetTrace(thisSigma, p)
    p = one3 * p

    CALL GetDevPart(thisSigma, rr)
    rr = rr / p

    CALL GetElasticModuli1(thisSigma, thisVoidRatio, Params, E_K, E_G)
    CALL GetStiffness(E_K, E_G, aC)
    dSigma3 = MATMUL(aC, strainIncNew)
    
    IF (iAbort == 1) THEN
        !WRITE(1, *) 'RKF2_EL', T, dT
        !WRITE(1, 1170) 'sNextStress = ', sNextStress
        !WRITE(1, 1170) 'dSigma1 = ', dSigma1
        !WRITE(1, 1170) 'dSigma2 = ', dSigma2
        !WRITE(1, 1170) 'strainIncNew = ', strainIncNew
        RETURN
    END IF

!// Higer order increments
    dSigma_H    = (one6 * dSigma1) + (two3 * dSigma2) + (one6 * dSigma3)

!// Lower order increments    
    dSigma   = (dSigma2) 

    snStress   = sNextStress   + dSigma
    
    call PrnSig(1,snStress,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)

    CALL GetTrace(snStress, p)
    p = one3 * p

    snStress_H   = sNextStress   + dSigma_H
!
! compute the norm. 
    CALL GetNorm_Contr(snStress_H, stressNorm_H)
    
! compute stress Error
    res7 = snStress_H - snStress
    CALL GetNorm_Contr(res7, curStepError1_stress)    
!
    !IF (stressNorm_H < 0.01D0) THEN
    IF (stressNorm_H < small_NormSig) THEN
        !write(1, *) , 'CORRECT NORMAL stress', iStep, iEl, T
        IF ((curStepError1_stress / stressNorm_H) <= 0.01D0 .AND. (curStepError1_stress / stressNorm_H) > TolR) THEN
            !stressNorm_H = 1.0D0
            stressNorm_H = curStepError1_stress / (0.5d0*TolR)
        END IF
    END IF
    curStepError1 = curStepError1_stress / stressNorm_H

!---------------------------------------------------------------------------------------------
    curStepError = curStepError1
!---------------------------------------------------------------------------------------------
!
    interstress = snStress_H
    interalpha = snAlpha_H
    
    CALL GetTrace(interstress, p)
    p = one3 * p
    CALL GetF(snStress_H, snAlpha_H, Params, fnNext)
    
!--------------------------------------------------------------------
    call PrnSig(1,interstress,xN1,xN2,xN3,sNextStress1a,sNextStress2a,sNextStress3a,PNextStress,QNextStress)
    
! check the higher order solution
    if (min(sNextStress1a,sNextStress2a,sNextStress3a) < -1.0D-5) then
        izero = 0 !it is reset 
        iNotension = 1
    end if
!--------------------------------------------------------------------    

    IF ((curStepError > TolR) .and. iNotension==0) THEN
        
        IF (dT <= dT_min) THEN

            CALL GetTrace(snStress_H, p)
            p = one3 * p
            
            IF (p >= Psmall) THEN
                !write(100, *) ' abort -> inaccurate sol not tensional but dt=dt_min and p>psmall'
                !write(100, *) T, dT
                !write(100, *) iEl, Int, istep, iTer
                !WRITE(100, 1170) 'dDevStrain = ', dDevStrain
                !write(100, *) 'dVolStrain =', dVolStrain
                !write(100, *) 'curStepError1 =', curStepError1
                !write(100, *) 'curStepError1_stress =', curStepError1_stress
                !write(100, 1170) 'sNextStress_H = ', snStress_H
                !write(100, 1170) 'sNextAlpha_H = ', snAlpha_H
                !write(100, 1170) 'dSigma1 = ', dSigma1
                !write(100, 1170) 'dSigma2 = ', dSigma2
                !write(100, 1170) 'dSigma3 = ', dSigma3
                !call flush(100)
                iAbort = 1
                return
            ELSE
                !write(100, *) 'return -> inaccurate sol not tensional but dt=dt_min and p<psmall '
                !write(100, *) T, dT
                !write(100, *) iEl, Int, istep, iTer
                !call flush(100)
                sNextStress	  = sNextStress! we accept the last updated just to allow the progress of the iterative process ! 
            
                sNextStress1 = sNextStress
            
                CALL GetF(sNextStress, sNextAlpha, Params, fnNext)
            
                T = T + dT

!               !nedeed update of the void ratio if T=1
                res7 = T * strainInc !strainInc_at_T
                CALL GetTrace(res7, res2) !dEpsVol_at_T
                sNextVoidRatio = CurVoid - (one + CurVoid) * res2 !void ratio at the EpsVol_T

                iabort = 1
                RETURN
            END IF
        ELSE
            !q = MAX(0.8D0 * ((TolR / curStepError)**0.2D0), 0.1D0)
            q = MAX(0.9D0 * ((TolR / curStepError)**(one3)), 0.25D0)
            
            !dT = MAX(q * dT , (one/4.0D0)*dT)
            dT = MAX(q * dT, dT_min)
            dT = MIN(dt, one - t)

        END IF
        
    ELSE IF ((curStepError < TolR) .and. iNotension==0) then
        
!       !copy of the previous accepted substep
        sNextStressCopy = sNextStress
        
!       !accept the soulution
        sNextStress	  = snStress_H
        
!       ! Calculate and output Mb and Md
!--------------------------------------------------------------------------------------remove
        call GetNormalToYield(sNextStress, sNextAlpha, iStep, sn, iAbort)
        call GetTrace(sNextStress, p)
        p = one3 * p
        call GetDevPart(sNextStress, rr)
        call GetNorm_Contr(rr, q_Norm_New)
        call GetNorm_Contr(sNextStress, sNextStress_Norm)
        q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
        rr = rr / p
        call GetNorm_Contr(rr, rr_Norm)
        call GetNorm_Contr(sNextAlpha, sNextAlpha_Norm)
        rr_Norm = dsqrt(3.0D0/2.0D0) * rr_Norm
        sNextAlpha_Norm = dsqrt(3.0D0/2.0D0) * sNextAlpha_Norm
        call GetLodeAngle(sn, cos3Theta)
        call g(cos3Theta, PAR_c, res24, dgdth)
        call GetPSI(sNextVoidRatio, p, Params, psi)
        Mb = PAR_Mc * EXP(-PAR_nb * psi)
        if (Mb > Mbc_max) then
            Mb = Mbc_max
        end if
        
        Mb = res24 * Mb
        Mc = res24 * PAR_Mc
        Md = res24 * PAR_Mc * EXP(PAR_nd * psi)

        call GetDevPart(FrozenStress, rr)
        call GetNorm_Contr(rr, q_Norm_Old)
        call GetNorm_Contr(FrozenStress, FrozenStress_Norm)
        q_Norm_Old = dsqrt(3.0D0/2.0D0) * q_Norm_Old
        ! 29-11-2021 checking q seems not ok, checking norm of stress
        
        !IF (p < 0.01D0 .AND. q_Norm_Old > q_Norm_New) THEN
        !IF (p <= Psmall .AND. FrozenStress_Norm > sNextStress_Norm) THEN
!--------------------------------------------------------------------------------------remove
        IF (p <= Psmall) THEN
!            !--------------------------------New lines for debugging to be removed----------------------------------------
                !write(100, *) 'return curStepError < TolR) .and. iNotension==0 and p<=psmall -> accept previous'
                !write(100, *) T, dT
                !write(100, *) sNextStress1a, sNextStress2a, sNextStress3a
                !write(100, *) iEl, Int, istep, iTer
                !call GetDevPart(interstress, rr)
                !call GetNorm_Contr(rr, q_Norm_New)
                !call GetNormalToYield(interstress, interAlpha, iStep, sn, iAbort)
                !if (iabort == 1) then
                !    return
                !endif               
                !call GetLodeAngle(sn, cos3Theta)
                !call GetLodeAngleForStressRatio(rr, Cos3ThetaN,iabort)
                !call g(cos3Theta, PAR_c, res24, dgdth)
                !call GetTrace(interstress, p)
                !p = one3 * p
                !q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New/p
                !!    p = one3 * p
                !!    ! move p
                !!    !IF (p <= PAR_Pmin) THEN
                !!    !    p = PAR_Pmin
                !!    !END IF
                !call GetPSI(sNextVoidRatio, p, Params, psi)
                !Mb = res24 * PAR_Mc * EXP(-PAR_nb * psi)
                !!    Mc = res24 * PAR_Mc
                !!    Md = res24 * PAR_Mc * EXP(PAR_nd * psi)
                !!    WRITE(2, 1200) iStep, iTer, T, Mb, Mc, Md, sNextStress
                !!    call flush(2)
                !!END IF
                !write(100, *) "cos3tetaN", cos3ThetaN
                !write(100, *) "cos3teta", cos3Theta
                !write(100, *) "p", p
                !write(100, *) "eta", q_Norm_New
                !write(100, *) "Mb", Mb
                !write(100, *) "sNextVoidRatio", sNextVoidRatio
                !write(100, *) "psi", psi
                !write(100, 1170) 'dDevStrain = ', dDevStrain
                !write(100, *) 'dVolStrain =', dVolStrain
                !write(100, *) 'curStepError1 =', curStepError1
                !write(100, *) 'curStepError1_stress =', curStepError1_stress
                !write(100, 1170) 'sNextStress_H = ', snStress_H
                !write(100, 1170) 'sNextAlpha_H = ', snAlpha_H
                !write(100, 1170) 'dSigma1 = ', dSigma1
                !write(100, 1170) 'dSigma2 = ', dSigma2
                !write(100, 1170) 'dSigma3 = ', dSigma3
                !
                !call flush(100)

                sNextStress = sNextStressCopy! we accept the last updated just to allow the progress of the iterative process ! 
                iabort = 1
                return
        !    END IF
        end if

        q = min(0.9D0 * ((TolR / curStepError)**(one3)), 4.0D0)
        T = T + dT

!       !neded update of the void ratio if T=1        
        res7 = T * strainInc !strainInc_at_T
        call GetTrace(res7, res2) !dEpsVol_at_T
        sNextVoidRatio = CurVoid - (one + CurVoid) * res2 !void ratio at the EpsVol_T
        
        dT = max(q * dT, dT_min)
        !dT = MIN(4.0D0*q * dT, dT_min)
        dT = min(dT, one - T)

    else
            q=0.5
            if (dT <= dT_min) THEN

                !write(1, 1170) 'snStress (in T=0)= ', snStress
                !write(1, 1170) 'snAlpha (in T=0)= ', snAlpha
                !call GetF(snStress, snAlpha, Params, fnNext)
                !call GetDevPart(interstress, rr)
                !call GetNorm_Contr(rr, q_Norm_New)
                !call GetNormalToYield(interstress, interAlpha, iStep, sn, iAbort)
                !if (iabort == 1) then
                !    return
                !endif               
                !call GetLodeAngle(sn, cos3Theta)
                !call GetLodeAngleForStressRatio(rr, Cos3ThetaN,iabort)
                !call g(cos3Theta, PAR_c, res24, dgdth)
                !call GetTrace(interstress, p)
                !p = one3 * p
                !q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New/p
                !!    p = one3 * p
                !!    ! move p
                !!    !IF (p <= PAR_Pmin) THEN
                !!    !    p = PAR_Pmin
                !!    !END IF
                !call GetPSI(sNextVoidRatio, p, Params, psi)
                !Mb = res24 * PAR_Mc * EXP(-PAR_nb * psi)
                !!    Mc = res24 * PAR_Mc
                !!    Md = res24 * PAR_Mc * EXP(PAR_nd * psi)
                !!    WRITE(2, 1200) iStep, iTer, T, Mb, Mc, Md, sNextStress
                !!    call flush(2)
                !!END IF
                !write(100, *) 'f (in T=0)= ', fnNext
                !write(100, *) 'Tensional stress, dT <= dT_min'
                !write(100, *) sNextStress1a, sNextStress2a, sNextStress3a
                !write(100, *) T, dT
                !write(100, *) iEl, Int, istep, iTer
                !write(100, *) "cos3tetaN", cos3ThetaN
                !write(100, *) "cos3teta", cos3Theta
                !write(100, *) "p", p
                !write(100, *) "eta", q_Norm_New
                !write(100, *) "Mb", Mb
                !write(100, *) 'curStepError1 =', curStepError1
                !write(100, *) 'curStepError1_stress =', curStepError1_stress
                !write(100, 1170) 'sNextStress_H = ', snStress_H
                !write(100, 1170) 'sNextAlpha_H = ', snAlpha_H
                !write(100, 1170) 'dSigma1 = ', dSigma1
                !write(100, 1170) 'dSigma2 = ', dSigma2
                !write(100, 1170) 'dSigma3 = ', dSigma3
                !call flush(100)

                iAbort = 1
                !iAbort = 220
                return
                
            else

                dT = max(q * dT, dT_min)
                dT = min(dt, one - t)
            end if
        
    end if
end do

return
end

!--------------------------------------------------------------------!
logical function isElIntegrationRKF23(Params) result(res)
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: Params(25)
double precision El_Scheme

El_Scheme = Params(i_El_Scheme)
res = El_Scheme > 0.1d0
end function isElIntegrationRKF23    
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
logical function isElIntegrationOneStepEuler(Params) result(res)
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: Params(25)
double precision El_Scheme

El_Scheme = Params(i_El_Scheme)
res = El_Scheme <= 0.1d0
end function isElIntegrationOneStepEuler 
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
logical function isElIntegrationRKF23alsoForIntersection(Params) result(res)
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: Params(25)
double precision El_Scheme

El_Scheme = Params(i_El_Scheme)
res = El_Scheme > 1.1d0
    end function isElIntegrationRKF23alsoForIntersection 
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
logical function isElIntegrationRKF23NotForIntersection(Params) result(res)
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: Params(25)
double precision El_Scheme

El_Scheme = Params(i_El_Scheme)
res = El_Scheme <+ 1.1d0
    end function isElIntegrationRKF23NotForIntersection 
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
logical function is_hLimitedWithDefaultStrategy(Params) result(res)
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: Params(25)
double precision Custom_h_max

Custom_h_max = Params(i_Custom_h_max)
res = Custom_h_max <+ 0.1d0
end function is_hLimitedWithDefaultStrategy 
!--------------------------------------------------------------------!
!--------------------------------------------------------------------!
logical function is_hM_ExceptionsStrategyActive(Params) result(res)
Use Parameters_position_in_params_iMod1
implicit none
double precision, intent(in) :: Params(25)
double precision hM_exc_Eq

hM_exc_Eq = Params(i_hM_exc_Eq)
res = hM_exc_Eq >+ 0.1d0
end function is_hM_ExceptionsStrategyActive    
!    WRITE(100, 1180) 'sn : dDevStrain = ', res3
!    WRITE(100, 1180) 'sn : rr = ', res2
!    WRITE(100, 1180) 'p = ', p
!    WRITE(100, *) 'h = ', h
!    WRITE(100, *) 'b0 = ', b0
!    WRITE(100, *) 'AlphaAlphaInDotN = ', AlphaAlphaInDotN
!    WRITE(100, *) 'NormAlphaAlphaInDotN = ', NormAlphaAlphaInDotN
!    WRITE(100, *) 'AlphaAlphaInDotNNorm = ', AlphaAlphaInDotNNorm
!    WRITE(100, 1170) 'alpha - alpha_in = ', thisAlpha - Curalpha_in
!    WRITE(100, 1180) 'bb:sn = ', PM1
!    WRITE(1, *) 'denominator = ', temp4
!    !WRITE(102, *) temp4, temp4_E, temp4_P
!    !WRITE(102, *) sNextDGamma, T
!    WRITE(100, 1180) 'numerator = ', temp5
!    WRITE(100, 1180) 'temp4_E = ', temp4_E
!    WRITE(100, *) 'temp4_P = ', temp4_P
!    WRITE(100, 1170) 'dVolStrain = ', dVolStrain
!    WRITE(103, 1200) dT, T, iRunge, thisAlpha, Curalpha_in, sn, bb, AlphaAlphaInDotN, PM, sNextDGamma, h, PM1
!    call flush(100)
!    call flush(102)
!    call flush(103)
!END IF
    
!If (iEl == 1047 .AND. Int == 4 .AND. iStep == 36 .AND. iTer == 5) then
!!If (iEl == 1 .and. Int == 1 .and. indicator == 500) then
!    WRITE(100, 1170) 'thisSigma = ', thisSigma
!    WRITE(100, 1170) 'thisAlpha = ', thisAlpha
!    WRITE(100, 1170) 'dStrainInc = ', dStrainInc
!!    WRITE(100, 1170) 'dSigma_A = ', dSigma_A
!!    WRITE(100, 1170) 'dSigma_B = ', dSigma_B
!    WRITE(100, 1170) 'sn = ', sn
!    WRITE(100, 1170) 'dSigma_A1 = ', dSigma_A1
!    WRITE(100, 1170) 'dp_elastic = ', p_elastic
!    WRITE(100, 1170) 'dSigma_B1 = ', dSigma_B1
!!    WRITE(100, 1170) 'dSigma_C = ', dSigma_C
!    WRITE(100, 1170) 'dSigma_C1 = ', dSigma_C1
!    WRITE(100, 1170) 'dp_plastic = ', p_plastic
!    WRITE(100, 1170) 'bb = ', bb
!    WRITE(100, 1180) 'res222 = ', res222
!    WRITE(100, 1180) 'sNextDGamma = ', sNextDGamma
!    WRITE(100, *) h
!    WRITE(100, 1170) 'dVolStrain = ', dVolStrain
!    call flush(100)
!END IF

! ======== Initialize.f90 ========
! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
!------------------------------------------------------------------------------------------
SUBROUTINE initialize(stress, alpha, alpha_in, alphaM, sMM_plus,    &
           sMM_minus, Void, Mb, Md, cos3theta, psi, q_Norm_New, p, Params)
!
! Purpose: Assigne the initial values to the local variables 
!           representative of the StVar0
! Arguments:
!            I/O  Type
!  stress     I   R()   : 
!  Params     I   R()  :
!  p          O   R    : 
!  q_Norm_New O   R    :
!  alpha      O   R()  :
!  alpha_in   O   R()  :
!  alphaM     O   R()  :
!  sMM_plus   O   R    : current radius resulting from + first  term in the isotropic hardening  [Memory rotation]
!  sMM_minus  O   R    : second term in the isotropic hardening  [volumetric plastic strains]
!  Void       O   R    :
!  psi        O   R    :
!  Mb         O   R    :
!  Md         O   R    :
!  cos3theta  O   R    :

Use Parameters_position_in_params_iMod1    
implicit none

double precision, intent(in) :: stress(6), Params(25)
double precision, intent(out) :: alpha(6), alpha_in(6), alphaM(6),  &
    sMM_plus, sMM_minus, Void, Mb, Md, cos3theta, psi, q_Norm_New, p

double precision rr(6), PAR_m, PAR_nb, PAR_nd, rrNorm, small, snrr(6), PAR_Mc, Mbc_max
    

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: zero = 0.0D0

PAR_Mc = Params(i_PAR_Mc)
PAR_m  = Params(i_PAR_m)
PAR_nb = Params(i_PAR_nb)
PAR_nd = Params(i_PAR_nd)

small = 1.0D-7

CALL GetTrace(stress, p)
p = one3 * p
CALL GetDevPart(stress, rr)
CALL GetNorm_Contr(rr, q_Norm_New)
q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
rr = rr / p

alpha = rr
alpha_in = rr
alphaM = rr
sMM_plus = PAR_m 
sMM_minus = zero

Void = Params(i_PAR_e_init)

CALL GetPSI(Void, p, Params, psi)
Mbc_max = Params(i_Mbc_max)
Mb = PAR_Mc * EXP(-PAR_nb * psi)! this is for tx compression not the current 
if (Mb > Mbc_max) then
    Mb = Mbc_max
end if

Md = PAR_Mc * EXP(PAR_nd * psi)

CALL GetNorm_Contr(rr, rrNorm)
IF (rrNorm < small) THEN
    rrNorm = small
END IF
snrr = rr/rrNorm
CALL GetLodeAngle(snrr, cos3theta)

RETURN
END

! ======== Integration.f90 ========
! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
!------------------------------------------------------------------------------------------
SUBROUTINE Integration(CurStress, CurVoid, CurAlpha, CurAlphaM, &
    CurMM_plus, CurMM_minus, Curalpha_in, sNextStress,  &
    strainInc, sNextVoid, sNextAlpha, sNextAlphaM, sNextMM_plus,   &
    sNextMM_minus, sNextalpha_in, Params, ipl, aCep, IDTask,    &
    iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
!
! Purpose: Update the stress tensor and the state variables
!
! Arguments:
!                  I/O   Type
!  IDTask           I    I()   :
!  iEl              I    I()   :
!  Int              I    I()   :
!  iStep            I    I()   :
!  iTer             I    I()  :
!  indicator        I    I()  :
!  strainInc        I    R() 
!  Params           I    R()  :
!  CurStress        I    R()  :
!  CurVoid          I    R()  :
!  CurAlpha         I    R()  :
!  CurAlphaM        I    R()  : 
!  CurMM_plus       I    R()  : the radius of the memory as affected only by the hardening
!  CurMM_minus      I    R()  : variation of the radius of the memory due to the Dilatancy
!  Curalpha_in      I    R()  :
!  sNextStress      O    R()  :
!  sNextVoid        O    R()  :
!  sNextAlpha       O    R()  :
!  sNextAlphaM      O    R()  :
!  sNextMM_plus     O    R()  : 
!  sNextMM_minus    O    R()  :
!  sNextalpha_in    O    R()  :
!  aCep             O    R()  :
!  PLNoTension      O    R    :
!  iAbort           O    i    :
!  ipl              IO   I    :
!  PsPsmall         IO   R    :
!
use extrainfo, only: onekPa
Use Parameters_position_in_params_iMod1

implicit none
double precision, intent(in) :: CurStress(6), CurAlpha(6), CurAlphaM(6),CurMM_plus, CurMM_minus, Curalpha_in(6),   &
    strainInc(6), Params(25)
integer, intent (in) :: IDTask, iEl, Int, iStep, indicator, iTer
integer, intent (out) :: iAbort
double precision, intent(in) :: CurVoid
double precision, intent(out) :: sNextStress(6), sNextAlpha(6), sNextAlphaM(6),sNextMM_plus, sNextMM_minus, sNextalpha_in(6),   &
    aCep(6,6), PLNoTension, sNextVoid
integer, intent(inout) :: ipl
double precision, intent(inout) :: PsPsmall

double precision dSigma(6), aC(6,6), sn(6), a0, a1, TolF, p, El_Scheme,   &
    sNextDGamma, sNextVoidRatio, res2, res1, pn, f, fn, E_K,E_G,    &
    elasticRatio, strainIncNew(6), InterVoid, res3(6),    &
    xN1(3),xN2(3),xN3(3),sNextStress1,sNextStress2,sNextStress3,PNextStress,QNextStress,sDevNextStress(6),pNext, res4(6)

double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0D0
double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: PAR_P_atm = 100.0D0
double precision, parameter :: small = 1.0D-10
double precision, parameter :: tol = 1.0D-6
double precision PAR_Pmin
!
logical, external :: isElIntegrationRKF23, isElIntegrationOneStepEuler,   &
    isElIntegrationRKF23alsoForIntersection, isElIntegrationRKF23NotForIntersection

iAbort = 0

PAR_Pmin = onekPa/100d0
PLNoTension = zero

a0 = zero
a1 = one

El_Scheme = Params(i_El_Scheme)
TolF = Params(i_TolF)
InterVoid = CurVoid

CALL GetElasticModuli1(CurStress, CurVoid, Params, E_K, E_G)
CALL GetStiffness(E_K, E_G, aC)

! Calculate the elastic trial stress with the integration scheme specified in input via the "El_Scheme" parameter 
if (isElIntegrationRKF23(Params)) then
    call EL_RungeKutta23(CurStress, InterVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, strainInc, Params, &
        sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoid, aCep, &
        IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
endif

if (iabort == 1 .or. isElIntegrationOneStepEuler(Params)) then    
	dSigma = MATMUL(aC, strainInc)
	sNextStress = CurStress + dSigma
else
    dSigma = sNextStress - CurStress
endif
iAbort = 0

pNext = (sNextStress(1) + sNextStress(2) + sNextStress(3))*one3
call PrnSig(1,sNextStress,xN1,xN2,xN3,sNextStress1,sNextStress2,sNextStress3,PNextStress,QNextStress)
!if (min(sNextStress1,sNextStress2,sNextStress3)<small) then
!if (pNext < tol) then
!    !sDevNextStress = PAR_Pmin * CurAlpha
!    !snextstress(1)  = sdevnextstress(1) + par_pmin
!    !snextstress(2)  = sdevnextstress(2) + par_pmin
!    !snextstress(3)  = sdevnextstress(3) + par_pmin
!    !snextstress(4)  = sdevnextstress(4)
!    !snextstress(5)  = sdevnextstress(5)
!    !snextstress(6)  = sdevnextstress(6)
!    
!    !p = PAR_Pmin
!    !CALL GetF(sNextStress, CurAlpha, Params, f)
!end if

CALL GetF(sNextStress, CurAlpha, Params, f)
CALL GetTrace(sNextStress, p)
p = one3 * p

CALL GetF(CurStress, CurAlpha, Params, fn) ! just for a check
CALL GetNorm_Contr(strainInc, res1)

!Check the trial stress
IF (f <= TolF) THEN
    !Pure elastic loading/unloading (trial inside the yield)'
    sNextAlpha    = CurAlpha
    sNextAlphaM   = CurAlphaM
    sNextMM_plus  = CurMM_plus
    sNextMM_minus = CurMM_minus
    sNextalpha_in = Curalpha_in
    sNextDGamma   = zero
    aCep = aC
    CALL GetTrace(strainInc, res2)
    sNextVoid = InterVoid - (1.0d0 + InterVoid) * res2
    ipl = 0
ELSE
    !Chenck the initial stress
    ipl = 1
    CALL GetF(CurStress, CurAlpha, Params, fn)
    CALL GetTrace(CurStress, pn)
    pn = one3 * pn
!   !Check if the initial stress is illegal
    IF (fn > TolF) THEN
        WRITE(100, *) '    trial stress is outside and initial outside the yield surface!', istep, fn
        WRITE(100, *) 'El, Int, iStep, iTer'
        WRITE(100, *) iEl, Int, iStep, iTer
        call flush(100)
        iAbort = 1
        RETURN
        CALL RungeKutta23(CurStress, CurVoid, CurAlpha, CurAlphaM,  &
            CurMM_plus, CurMM_minus, sNextalpha_in,     &
            strainInc, Params, sNextStress, sNextAlpha, &
            sNextAlphaM, sNextMM_plus, sNextMM_minus, &
            sNextVoidRatio, aCep, IDTask, iEl, Int,  &
            iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
            
            sNextVoid = sNextVoidRatio
        
        IF (iAbort == 1) THEN
            RETURN
        END IF
!   !Check if the initial stress is iniside the yield surface
    ELSEIF (fn < -TolF) THEN
!       !Calculate the intersection with the yield
        CALL IntersectionFactor(CurStress, CurVoid, strainInc,  &
                    CurAlpha, a0, a1, Params, elasticRatio, iAbort)
        IF (iAbort == 1) THEN
                RETURN
        END IF
        res3 = elasticRatio * strainInc
        !dSigma = MATMUL(aC, res3)

		! Calculate the stress at the intersection point with the yield
        if(isElIntegrationRKF23alsoForIntersection(Params))then
	        call EL_RungeKutta23(CurStress, InterVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, res3, Params, &
	        sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoid, aCep, &
	        IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
		endif
        if (iabort == 1 .or. isElIntegrationRKF23NotForIntersection(Params)) then
            dSigma = MATMUL(aC, res3)
            sNextStress = CurStress + dSigma
            iAbort = 0
  
        else
            dSigma = sNextStress - CurStress
        endif
        
        CALL GetTrace(res3, res2)
        InterVoid = CurVoid - (one + CurVoid) * res2
        strainIncNew = (one - elasticRatio) * strainInc
        res3 = CurStress + dSigma
        
        !Compute the normal to the yield in the intersection point
        CALL GetNormalToYield(res3, CurAlpha, iStep, sn, iAbort)
        IF (iAbort == 1) THEN
            RETURN
        END IF
        !Check for the updating of alfa_initial---------------------------
        res4 = CurAlpha - Curalpha_in
        CALL DoubleDot2_2_Contr(res4, sn, res1)
        IF (res1 < -tol) THEN
            sNextalpha_in = CurAlpha
        ELSE
            sNextalpha_in = Curalpha_in
        END IF
        !-------------------------------------------------------------------
        CALL RungeKutta23(res3, InterVoid,    &
                    CurAlpha, CurAlphaM, CurMM_plus, &
                    CurMM_minus, sNextalpha_in, strainIncNew,  &
                    Params, sNextStress, sNextAlpha,  &
                    sNextAlphaM, sNextMM_plus, sNextMM_minus, &
                    sNextVoidRatio, aCep, IDTask, iEl, Int, &
                    iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
            
        sNextVoid = sNextVoidRatio
        IF (iAbort == 1) THEN
            RETURN
        END IF
    ELSE
        !Check if unloading occurs
        CALL GetNormalToYield(CurStress, CurAlpha, iStep, sn, iAbort)
        IF (iAbort == 1) THEN
            WRITE(100, *) 'failed to compute sn before chcking for unloading'
            call flush(100)
            RETURN
        END IF
        CALL DoubleDot2_2_Contr(sn, dSigma, res1)
        IF (res1 > -DSQRT(TolF)) THEN
            
            !No unloading (pure plastic loading): Check for the updateing of alfa_initial----------------
            res4 = CurAlpha - Curalpha_in
            CALL DoubleDot2_2_Contr(res4, sn, res1)
            IF (res1 < -tol) THEN
                sNextalpha_in = CurAlpha
            ELSE
                sNextalpha_in = Curalpha_in
            END IF
            !--------------------------------------------------------------------------------------------
            CALL RungeKutta23(CurStress, CurVoid, CurAlpha, &
                CurAlphaM,CurMM_plus, CurMM_minus, sNextalpha_in, &
                strainInc, Params, sNextStress, sNextAlpha,    &
                sNextAlphaM, sNextMM_plus, sNextMM_minus,  &
                sNextVoidRatio, aCep, IDTask, iEl, Int,  &
                iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)               
            sNextVoid = sNextVoidRatio
            IF (iAbort == 1) THEN
                RETURN
            END IF
        ELSE

            !Elastic unloading followed by plastic unloading: search the intersection with the yield surface
            CALL IntersectionFactor_Unloading(CurStress,CurVoid,    &
                            strainInc,CurAlpha, Params, elasticRatio, iAbort)
            IF (iAbort == 1) THEN
                RETURN
            END IF
            strainIncNew = elasticRatio * strainInc
            !dSigma = matmul(aC, strainIncNew)
            
            ! Calculate the stress at the intersection point with the yield
			if(isElIntegrationRKF23alsoForIntersection(Params))then
		        call EL_RungeKutta23(CurStress, InterVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, strainIncNew, Params, &
		        sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoid, aCep, &
		        IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
			endif            
            
	        if (iabort == 1 .or. isElIntegrationRKF23NotForIntersection(Params)) then
	            dSigma = MATMUL(aC, strainIncNew)
	            sNextStress = CurStress + dSigma
	            iAbort = 0
	        else
	            dSigma = sNextStress - CurStress
	        endif            
                       
            CALL GetTrace(strainIncNew, res2)
                InterVoid = CurVoid - (one + CurVoid) * res2
                strainIncNew = (one - elasticRatio) * strainInc
                res3 = CurStress + dSigma
            !Check for the updateing of alfa_initial
            CALL GetNormalToYield(res3, CurAlpha, iStep, sn, iAbort)
            IF (iAbort == 1) THEN
                RETURN
            END IF
            !Check for the updating of alfa_initial------------------------------------------------------
            res4 = CurAlpha - Curalpha_in
            CALL DoubleDot2_2_Contr(res4, sn, res1)
            IF (res1 < -tol) THEN
                sNextalpha_in = CurAlpha
            ELSE
                sNextalpha_in = Curalpha_in
            END IF
            !--------------------------------------------------------------------------------------------
            CALL RungeKutta23(res3, InterVoid,    &
                    CurAlpha, CurAlphaM, CurMM_plus,    &
                    CurMM_minus, sNextalpha_in, strainIncNew, &
                    Params, sNextStress, &
                    sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, &
                    sNextVoidRatio, aCep,    &
                    IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
                
                    sNextVoid = sNextVoidRatio
            IF (iAbort == 1) THEN
                 RETURN
            END IF
        END IF
    END IF
END IF

RETURN
END


! ======== SANISANDMS.f90 ========
! ----------------------------------------------------------------------------------------
!---------------
! 'SANISAND-MS': a user-defined soil model for PLAXIS
! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
Subroutine SANISANDMS ( IDTask, iMod, IsUndr,   &
                        iStep, iTer, iEl, Int,  &
                        X, Y, Z, &
                        Time0, dTime,    &
                        Props, Sig0, Swp0, StVar0,   &
                        dEps, D, BulkW,  &
                        Sig, Swp, StVar, ipl,    &
                        nStat, NonSym, iStrsDep, iTimeDep, iTang,   &
                        iAbort )

! ----------------------------------------------------------------------------------------
! User-defined soil model: SANISAND-MS (iMod=1)
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : RETURN number of state variables
!                       5 : inquire matrix properties
!                           RETURN switch for non-symmetric D-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastic material stiffness matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  Int      I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure at start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulk modulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : =1 for Non-Symmetric D-matrix
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iAbort   O   I    : =1 to force arrest of calculation
!   
! ----------------------------------------------------------------------------------------
! Expected content of Props(1..50)
! 
!  1 : G0           Shear modulus coefficient
!  2 : nu           Poisson's ratio
!  3 : Mc           critical state M value (in triaxial compression) 
!  4 : c            Me/Mc ratio between critical stress ratios in triaxial extention and compression
!  5 : lambda_c     slope of critical state line from void ratio-ln_p'-diagram
!  6 : e0           reference void ratio of the critical state line at pc = 0 kPa
!  7 : ksi          critical state line costant
!  8 : m            radious of the yield surface in the stress ratio space
!  9 : h0           hardening parameter
! 10 : ch           hardening parameter
! 11 : nb           bounding surface parameter
! 12 : A0           dilatancy parameter
! 13 : nd           dilatancy surface parameter
! 14 : mu0          ratcheting parameter
! 15 : zeta         memory surface shrinkage parameter
! 16 : beta         dilatancy memory parameter
! 17 : TolF         yield surface tolerance
! 18 : TolR         error tolerance for the RKF23 adaptive substepping scheme
! 19 : e_init       initial void ratio
! 20 : emax         maximum void ratio
! 21 : emin         minimum void ratio
! 22 : Mbc_max      Maximum value of the Bounding stress ratio (in triaxial compression)  
! 23 : El_Scheme    Toggle to choose the integration scheme of the elastic equations: 
!                     [  0  ] = Single step Forward Euler,
!                     [  1  ] = Automatic sub-step RKF23 for elastic trial only, 
!                     [  2  ] = Automatic sub-step RKF23 for elastic trial and path to the intersection stress
! 24 : hM_exc_Eq    Toggle to activate the new hM equation for tackling corner cases
!                     [  0  ] = inactive. 
!                     [  1  ] = active
! 25 : h_max        Toggle to customize the limit for h:
!                     [  0  ] = Default streategy limit separately the exponential term and h
!                     [value] = h is limited to the input value. It could 1.0d8
!
! Note: Toggles 24 and 25 will exclusively be in place during the testing phase
! ----------------------------------------------------------------------------------------

Use Parameters_position_in_params_iMod1
Use State_Var_position_in_StVar_iMod1
implicit none

integer, intent(in) :: IDTask, iMod, IsUndr, iStep, iTer, iEl, Int

double precision, intent(in) :: X, Y, Z, Time0, dTime, Props(*), Swp0, dEps(*) 

double precision, intent(inout) :: Sig0(*), StVar0(*)
double precision, intent(out) :: D(6,6), BulkW,Sig(*),Swp,StVar(*)
integer, intent(out) :: ipl, nStat, NonSym, iStrsDep, iTimeDep, iTang, iAbort
!
!---  Local variables
!
double precision stress(6), strain(6), strainInc(6), Alpha(6), AlphaM(6), alpha_in(6), sNextStress(6),    &
    sNextAlpha(6), sNextAlphaM(6), sNextalpha_in(6), sn(6), TolF, &
    dd(6), bb(6), bM(6), R(6), rr(6),res4(6),res5(6),dDevStrain(6), A, temp4, sNextDGamma,f,PM, b0, B, ZZ, &
    ddsdde(6,6), xNu, PAR_e_init, E_K,E_G, xNu_U, Void, sNextVoid, sNextMM_plus, &
    sNextMM_minus,sMM_plus,  sMM_minus, res, res3, res2, rBtheta, rDtheta, p, psi,  &
    VoidRatio, Fac, dVolStrain, dSwp, dEpsV, DDD, C, Cos3Theta, strainInc1(6), &
    PLNoTension, alfa_trace, emax, emin, PsPsmall, PAR_c, res24, dgdth, Mb, Md, PAR_Mc, PAR_nb, Mc, PAR_nd, &
    q_Norm_New, q_p, Mbc_max

double precision h, hM, fnNext

integer io, i, j, nStatV, itestnan, indicator

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: two3 = 2.0D0/3.0D0
double precision, parameter :: tol = 1.0D-10
double precision, parameter :: ZERO2D = 1.0D-30


ipl=0

! indicator = 1, open warning output; indicator = 0, close warning output
! 200 presents norm of stress and small stress problem
indicator = 2

io=0
!     IF (iEl==-87 .And. Int==12) io=1 ! possibly write some debug info for a point
!     Number of state variables

nStatV = 31
1060 FORMAT(A, I5, A, I5, A, I5, A, I5)
1070 FORMAT(A, 6F15.10)
1080 FORMAT(A, 1F10.4)
1090 FORMAT(A, 6F35.23)
1110 FORMAT(18F20.15)
1120 FORMAT(42F20.15)

iAbort = 0

xNu         = Props(i_PAR_nu)
PAR_Mc      = Props(i_PAR_Mc)
PAR_c       = Props(i_PAR_c)
PAR_nb      = Props(i_PAR_nb)
PAR_nd      = Props(i_PAR_nd)
TolF        = Props(i_TolF)
PAR_e_init  = Props(i_PAR_e_init)
emax        = props(i_emax)
emin        = props(i_emin)
Mbc_max     = props(i_Mbc_max)

! change the signs (In PLAXIS, compression is negative)	
stress = -Sig0(1:6)
strain = -dEps(7:12)
strainInc = -dEps(1:6)

strainInc1 = strainInc

! ----------------------------------------------------------------------------------------
! get number of state parameters
! ----------------------------------------------------------------------------------------
IF (IDTask == 4) THEN
    nStat = nStatV
    RETURN
End IF

! ----------------------------------------------------------------------------------------
! get matrix attributes
! ----------------------------------------------------------------------------------------
IF (IDTask == 5) THEN ! matrix type
    NonSym   = 0  ! 1 for non-symmetric D-matrix
    iStrsDep = 1  ! 1 for stress dependent D-matrix
    iTang    = 0  ! 1 for tangent D-matrix
    iTimeDep = 0  ! 1 for time dependent D-matrix
End IF

! ----------------------------------------------------------------------------------------
! Initialize state variables
! ----------------------------------------------------------------------------------------
IF (IDTask == 1) THEN
    IF (StVar0(i_resetState) /= 777D0) THEN ! To void reinitialization problem in continuous phases
        CALL initialize(stress, Alpha, alpha_in, AlphaM, sMM_plus,  &
            sMM_minus, Void, Mb, Md, cos3theta, psi, q_Norm_New, p, Props)

        StVar0(i_Alpha_xx)      = Alpha(1)
        StVar0(i_Alpha_yy)      = Alpha(2)
        StVar0(i_Alpha_zz)      = Alpha(3)
        StVar0(i_Alpha_xy)      = Alpha(4)
        StVar0(i_Alpha_yz)      = Alpha(5)
        StVar0(i_Alpha_zx)      = Alpha(6)
        StVar0(i_AlphaM_xx)     = AlphaM(1)
        StVar0(i_AlphaM_yy)     = AlphaM(2)
        StVar0(i_AlphaM_zz)     = AlphaM(3)
        StVar0(i_AlphaM_xy)     = AlphaM(4)
        StVar0(i_AlphaM_yz)     = AlphaM(5)
        StVar0(i_AlphaM_zx)     = AlphaM(6)
        StVar0(i_sMM_plus)      = sMM_plus
        StVar0(i_sMM_minus)     = sMM_minus
        StVar0(i_alpha_in_xx)   = alpha_in(1)
        StVar0(i_alpha_in_yy)   = alpha_in(2)
        StVar0(i_alpha_in_zz)   = alpha_in(3)
        StVar0(i_alpha_in_xy)   = alpha_in(4)
        StVar0(i_alpha_in_yz)   = alpha_in(5)
        StVar0(i_alpha_in_zx)   = alpha_in(6)
        StVar0(i_Void)          = Void
        StVar0(i_resetState)    = 777D0
        StVar0(i_PLNoTension)   = 0.0D0 !deprecated use. Their name is removed from output. They will be repositioned in StVar
        StVar0(i_PsPsmall)      = 0.0D0 !deprecated use. Their name is removed from output. They will be repositioned in StVar
        StVar0(i_Mb)            = Mb
        StVar0(i_Md)            = Md
        StVar0(i_cos3theta)     = cos3theta
        StVar0(i_psi)           = psi
        StVar0(i_q_Norm_New)    = q_Norm_New
        StVar0(i_p)             = p
        StVar0(i_MM)            = sMM_plus
        
        StVar(1:31)             = StVar0(1:31)
    END IF
RETURN
End IF

! ----------------------------------------------------------------------------------------
! calculate stresses
! ----------------------------------------------------------------------------------------
IF (IDTask == 2) THEN
    
    Alpha(1)     = StVar0(i_Alpha_xx)
    Alpha(2)     = StVar0(i_Alpha_yy)
    Alpha(3)     = StVar0(i_Alpha_zz)
    Alpha(4)     = StVar0(i_Alpha_xy)
    Alpha(5)     = StVar0(i_Alpha_yz)
    Alpha(6)     = StVar0(i_Alpha_zx)
    AlphaM(1)    = StVar0(i_AlphaM_xx)
    AlphaM(2)    = StVar0(i_AlphaM_yy)
    AlphaM(3)    = StVar0(i_AlphaM_zz)
    AlphaM(4)    = StVar0(i_AlphaM_xy)
    AlphaM(5)    = StVar0(i_AlphaM_yz)
    AlphaM(6)    = StVar0(i_AlphaM_zx)
    sMM_plus     = StVar0(i_sMM_plus)! current size of the memory as it would be if it was affected only by the expantion due to the hardening
    sMM_minus    = StVar0(i_sMM_minus)! variation of the size of the memory due to the dilatancy
    alpha_in(1)  = StVar0(i_alpha_in_xx)
    alpha_in(2)  = StVar0(i_alpha_in_yy)
    alpha_in(3)  = StVar0(i_alpha_in_zz)
    alpha_in(4)  = StVar0(i_alpha_in_xy)
    alpha_in(5)  = StVar0(i_alpha_in_yz)
    alpha_in(6)  = StVar0(i_alpha_in_zx)
    Void            = StVar0(i_Void)
    PsPsmall        = StVar0(i_PsPsmall)
    Mb              = StVar0(i_Mb)

! 21-12-2021 check input q/p reach Mb or not------------------
    CALL GetTrace(stress, p)
    p = one3 * p
    CALL GetDevPart(stress, rr)
    CALL GetNorm_Contr(rr, q_Norm_New)
    q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
    q_p = q_Norm_New / p
!-------------------------------------------------------------
    !If (iStep == 1 .and. iEl == 1 .and. indicator == 500) then
    !    !WRITE(100, *) 'before integration', iStep, iTer!, iEl, Int
    !    WRITE(100, '(6f25.15)') strainInc
    !    WRITE(100, '(6f25.15)') stress
    !    WRITE(100, '(6f25.15)') Alpha
    !    call flush(100)
    !END IF

    CALL Integration(stress, Void, Alpha, AlphaM, sMM_plus, &
       sMM_minus, alpha_in, sNextStress, strainInc, sNextVoid,   &
       sNextAlpha,sNextAlphaM, sNextMM_plus, sNextMM_minus,    &
       sNextalpha_in, Props, ipl, ddsdde, IDTask, iEl, Int, iStep, &
       indicator, PLNoTension, PsPsmall, iAbort, iTer)

    IF (iAbort == 1) THEN
        RETURN
    END IF
    
    !If (iStep == 1 .and. iEl == 1 .AND. Int = 4.and. indicator == 200) then
        !WRITE(100, *) 'integration succeed', iStep, iTer!, iEl, Int
        !WRITE(100, '(6f25.15)') strainInc
        !WRITE(100, '(6f25.15)') sNextStress
        !WRITE(100, '(6f25.15)') sNextAlpha
        !call flush(100)
    !END IF
    CALL GetNormalToYield2(sNextStress, sNextAlpha, iStep, sn)
    CALL GetTrace(sNextStress, p)
    p = one3 * p
    
    ! Calculate the current stress ratio
    CALL GetDevPart(sNextStress, rr)
    CALL GetNorm_Contr(rr, q_Norm_New)
    q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
    
    ! Calculate the current Bounding, Critical state and Dilatancy stress ratios
    CALL GetLodeAngle(sn, cos3Theta)
    CALL g(cos3Theta, PAR_c, res24, dgdth)
    CALL GetPSI(sNextVoid, p, Props, psi)
    Mb = PAR_Mc * EXP(-PAR_nb * psi)
    if (Mb > Mbc_max) then
        Mb = Mbc_max
    endif
    Mb = res24 * Mb
    Mc = res24 * PAR_Mc
    Md = res24 * PAR_Mc * EXP(PAR_nd * psi)

    Sig(1:6)     = -sNextStress
    StVar(i_Alpha_xx)   = sNextAlpha(1)
    StVar(i_Alpha_yy)   = sNextAlpha(2)
    StVar(i_Alpha_zz)   = sNextAlpha(3)
    StVar(i_Alpha_xy)   = sNextAlpha(4)
    StVar(i_Alpha_yz)   = sNextAlpha(5)
    StVar(i_Alpha_zx)   = sNextAlpha(6)
    StVar(i_AlphaM_xx)  = sNextAlphaM(1)
    StVar(i_AlphaM_yy)  = sNextAlphaM(2)
    StVar(i_AlphaM_zz)  = sNextAlphaM(3)
    StVar(i_AlphaM_xy)  = sNextAlphaM(4)
    StVar(i_AlphaM_yz)  = sNextAlphaM(5)
    StVar(i_AlphaM_zx)  = sNextAlphaM(6)
    StVar(i_sMM_plus)    = sNextMM_plus
    StVar(i_sMM_minus)   = sNextMM_minus
    StVar(i_alpha_in_xx) = sNextalpha_in(1)
    StVar(i_alpha_in_yy) = sNextalpha_in(2)
    StVar(i_alpha_in_zz) = sNextalpha_in(3)
    StVar(i_alpha_in_xy) = sNextalpha_in(4)
    StVar(i_alpha_in_yz) = sNextalpha_in(5)
    StVar(i_alpha_in_zx) = sNextalpha_in(6)
    StVar(i_Void)           = sNextVoid
    StVar(i_resetState)     = 777D0
    StVar(i_PLNoTension)    = PLNoTension
    StVar(i_PsPsmall)       = PsPsmall
    StVar(i_Mb)             = Mb
    StVar(i_Md)             = Md
    StVar(i_cos3theta)      = cos3Theta         !cos(3*Lode angle(sn)) 
    StVar(i_psi)            = psi
    StVar(i_q_Norm_New)     = q_Norm_New
    StVar(i_p)              = p
    StVar(i_MM)            = sNextMM_plus + sNextMM_minus
!------------------------------------------------------    
    CALL GetF(sNextStress, sNextAlpha, Props, fnNext)
    alfa_trace = sNextAlpha(1) + sNextAlpha(2) + sNextAlpha(3)
!------------------------------------------------------   
    IF (IsUndr == 1) THEN
        xNu_U = 0.495d0 ! Undrained Poissons' ratio
        Fac = (1.d0+xNu_U)/(1.0d0-2.0d0*xNu_U)-(1.0d0+xNu)/(1.0d0-2.0d0*xNu)
        ! 24-11-2021
        CALL GetElasticModuli1(sNextStress,sNextVoid,Props, E_K, E_G)
        Fac = 2D0*E_G/3D0*Fac
        BulkW = Fac
        dEpsV = dEps(1) + dEps(2) + dEps(3)
        dSwp  = BulkW * dEpsV
        Swp   = Swp0 + dSwp
    ELSE
        Swp = Swp0
    END IF
!
!     Nan check: Sig, sNextAlpha, sNextAlphaM, sNextalpha_in, sNextMM_plus, sNextMM_minus, sNextVoid
    Do i=1,6
        If ( IsNaN( Sig( i))) Then
            Write(100,*), 'Sig(i)', Sig(i)
            Write(100,*), 'iEl', iEl
            call flush(100)
            iabort = 1
            return
        End If
    End Do

    Do i=1,6
        If (IsNaN( sNextAlpha(i))) Then
                Write(100,*), 'sNextAlpha',sNextAlpha
                Write(100,*), 'iEl',iEl
                call flush(100)
                iabort = 1
                return
        End If
    End Do

! 06-29-2021
    Do i = 1, 6
        If ( IsNaN( sNextAlphaM( i))) Then
        Write(100,*), 'sNextAlphaM', sNextAlphaM
        Write(100,*), 'iEl',iEl
        call flush(100)
        iabort = 1
        return
        endif
    End Do

    Do i = 1,6
        If ( IsNaN( sNextalpha_in(i))) Then
            Write(100,*), 'sNextalpha_in',sNextalpha_in
            Write(100,*), 'iEl',iEl
            call flush(100)
            iabort = 1
            return
        End If
    End Do

    If ( IsNaN( sNextMM_plus)) Then
        Write(100,*), 'sNextMM_plus',sNextMM_plus
        Write(100,*), 'iEl',iEl
        call flush(100)
        iabort = 1
        return
    End If

    If ( IsNaN( sNextMM_minus)) Then
        Write(100,*), 'sNextMM_minus',sNextMM_minus
        Write(100,*), 'iEl',iEl
        call flush(100)
        iabort = 1
        return
    End If

    If ( IsNaN( sNextVoid)) Then 
        iabort = 1
        return
    End If

! NaN output
! 06-29-2021  
!    IF ((itestnan == 1) .AND. indicator == 2) THEN
!        WRITE(100,*) 'IDTask = ', IDTask
!        WRITE(100, 1070) 'Inistress =', stress(1:6)
!        WRITE(100, 1070) 'IniAlpha =', Alpha
!        WRITE(100, 1070) AlphaM
!        WRITE(100, 1080) sMM_plus
!        WRITE(100, 1080) sMM_minus
!        WRITE(100, 1070) 'Inialpha_in =', alpha_in
!        WRITE(100, 1080) 'IniVoid =', Void
!        WRITE(100, 1090) 'strainInc =', strainInc
!        WRITE(100,1060)'iStep =',iStep,' iTer =',iTer,' iEl = ', iEl, ' Int = ', Int
!        WRITE(100, 1070) 'sNextStress = ', sNextStress
!        WRITE(100, 1070) 'sNextAlpha = ', sNextAlpha
!        WRITE(100, 1070) 'sNextAlphaM = ', sNextAlphaM
!        WRITE(100, 1080) 'sNextMM_plus = ', sNextMM_plus
!        WRITE(100, 1080) 'sNextMM_minus = ', sNextMM_minus
!        WRITE(100, 1070) 'sNextalpha_in = ', sNextalpha_in
!        WRITE(100, 1080) 'sNextVoid = ', sNextVoid
!        call flush(100) 
!        iAbort = 1
!         return
!    END IF

RETURN
End IF

! ----------------------------------------------------------------------------------------
! calculate material stiffness matrix (D-matrix)
! ----------------------------------------------------------------------------------------
IF (IDTask == 6 .or. IDTask == 3) THEN ! matrix type
    !write(100,*), IDTask

    Alpha(1)     = StVar0(i_Alpha_xx)
    Alpha(2)     = StVar0(i_Alpha_yy)
    Alpha(3)     = StVar0(i_Alpha_zz)
    Alpha(4)     = StVar0(i_Alpha_xy)
    Alpha(5)     = StVar0(i_Alpha_yz)
    Alpha(6)     = StVar0(i_Alpha_zx)
    AlphaM(1)    = StVar0(i_AlphaM_xx)
    AlphaM(2)    = StVar0(i_AlphaM_yy)
    AlphaM(3)    = StVar0(i_AlphaM_zz)
    AlphaM(4)    = StVar0(i_AlphaM_xy)
    AlphaM(5)    = StVar0(i_AlphaM_yz)
    AlphaM(6)    = StVar0(i_AlphaM_zx)
    sMM_plus     = StVar0(i_sMM_plus)
    sMM_minus    = StVar0(i_sMM_minus)
    alpha_in(1)  = StVar0(i_alpha_in_xx)
    alpha_in(2)  = StVar0(i_alpha_in_yy)
    alpha_in(3)  = StVar0(i_alpha_in_zz)
    alpha_in(4)  = StVar0(i_alpha_in_xy)
    alpha_in(5)  = StVar0(i_alpha_in_yz)
    alpha_in(6)  = StVar0(i_alpha_in_zx)
!------------------------------------------------------------
    Void            = StVar0(i_Void)    
    stress=-Sig0(1:6)
    
    CALL GetElasticModuli1(stress, Void, Props, E_K, E_G)
    CALL GetStiffness(E_K, E_G, D)

! calculate bulk modulus
!        BulkW = 0  
!    IF (IsUndr == 1) THEN
!        xNu_U = 0.495d0
!        Fac = (1+xNu_U)/(1-2*xNu_U) - (1+xNu)/(1-2*xNu)
!        CALL GetTrace(strainInc, res)
!        VoidRatio = Void - (1.0d0 + Void) * res
!        CALL GetElasticModuli(stress, VoidRatio, Props, E_K, E_G)
!        Fac = 2D0*E_G/3D0*Fac
!        BulkW = Fac
!    END IF
    
    !WRITE(100, *) 'IDTask', IDTask, iStep, iTer
    !WRITE(100,'(6F14.6)') D
    !call flush(100)

End IF

Return

End



!=======================================================================
!  Abaqus UMAT interface to the SANISAND-MS PLAXIS-UDSM routine
!  (SANISANDMS in SANISANDMS.f90), called from UEL_C3D8UW.
!
!  Conventions (checked against UEL_C3D8UW and SANISANDMS.f90):
!   - UEL passes Voigt order 11,22,33,12,23,13 with engineering shear
!     strains; this is the PLAXIS order xx,yy,zz,xy,yz,zx -> no reorder.
!   - Tension positive in both codes (SANISANDMS flips signs inside).
!   - Units kPa / m / s / Mg.
!
!  PROPS (= UEL PROPS(8:)):
!    1..25  SANISAND-MS parameters in UDSM order (see USRADDDF.f90)
!    26     p_atm (kPa)
!    27     p_min for the tangent only (kPa); elastic moduli in DDSDDE are
!           evaluated with max(p, p_min) so that the tangent never
!           becomes singular at liquefaction. Stress update is untouched.
!    28     maximum number of sub-steps used if SANISANDMS aborts
!    29     debug flag (>0.5: print TIME/KINC of element 1, GP 1 to .msg)
!
!  STATEV(1..31) = SANISAND-MS StVar (StVar(22)=777 -> initialised).
!  Initialisation is triggered by StVar(22)/=777 and NOT by TIME, so it
!  is independent of the time value passed by the UEL.
!
!  A copy of the last computed (=converged at end of increment) stress
!  and state of every Gauss point is kept in COMMON /SSMSDUMP/ and
!  written to <jobname>_state.txt by UEXTERNALDB at the end of every
!  step (used for the iterative geostatic re-import).
!=======================================================================
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,     &
        drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred, &
        cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
        celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      use extraInfo, only: oneAtm, onekPa
      implicit none
      character(80) cmname
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt, layer, kspt, &
              kstep, kinc
      double precision stress(ntens), statev(nstatv), ddsdde(ntens,ntens), &
        ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),         &
        time(2), predef(1), dpred(1), props(nprops), coords(3), drot(3,3), &
        dfgrd0(3,3), dfgrd1(3,3), sse, spd, scd, rpl, drpldt, dtime, temp, &
        dtemp, pnewdt, celent

      integer, parameter :: NSV = 31, MAXE = 20000, MAXG = 8, NDV = 6+NSV
      double precision stdump(NDV, MAXG, MAXE)
      integer idump(MAXE)
      common /SSMSDUMP/ stdump, idump

      double precision PR(25), sig0(6), sig(6), sv0(NSV), sv(NSV), dEps(12), &
        D(6,6), bulkw, swp, pmin, p, dsub(6), E_K, E_G, sigN(6), svN(NSV)
      integer ipl, nStat, nonsym, istrsdep, itimedep, itang, iabort, i, &
              nsub, isub, maxsub

      PR(1:25) = props(1:25)
      oneAtm  = props(26)
      onekPa  = oneAtm / 100.0d0
      pmin    = props(27)
      maxsub  = max(1, nint(props(28)))

      if (nprops >= 29) then
        if (props(29) > 0.5d0 .and. noel == 1 .and. npt == 1) then
          write(7,'(a,2i6,4es14.6)') 'UMATDBG kstep,kinc,time(1:2),dtime,stran3:', &
            kstep, kinc, time(1), time(2), dtime, stran(3)
        end if
      end if
      sig0 = stress(1:6)
      sv0  = 0.0d0
      sv0(1:min(NSV,nstatv)) = statev(1:min(NSV,nstatv))

!     ---- initialise SANISAND-MS state from the current stress ----
      if (abs(sv0(22) - 777.0d0) > 1.0d-6) then
        call SANISANDMS(1, 1, 0, kstep, 0, noel, npt, 0d0, 0d0, 0d0,    &
          0d0, dtime, PR, sig0, 0d0, sv0, dEps, D, bulkw, sig, swp, sv,   &
          ipl, nStat, nonsym, istrsdep, itimedep, itang, iabort)
      end if

!     ---- stress update, with sub-stepping if the integrator aborts ----
      nsub = 1
      do
        sig = sig0
        sv  = sv0
        iabort = 0
        dsub = dstran(1:6) / dble(nsub)
        do isub = 1, nsub
          dEps(1:6)  = dsub
          dEps(7:12) = stran(1:6) + dble(isub-1)*dsub
          call SANISANDMS(2, 1, 0, kstep, kinc, noel, npt, 0d0, 0d0,    &
            0d0, 0d0, dtime, PR, sig, 0d0, sv, dEps, D, bulkw, sigN, swp, &
            svN, ipl, nStat, nonsym, istrsdep, itimedep, itang, iabort)
          if (iabort /= 0) exit
          sig = sigN
          sv  = svN
        end do
        if (iabort == 0) exit
        nsub = 4*nsub
        if (nsub > maxsub) exit
      end do

      if (iabort /= 0) then
!       give up: keep the old state and ask Abaqus for a smaller step
        sig = sig0
        sv  = sv0
        pnewdt = 0.5d0
        write(7,'(a,i7,i3,a,6es12.4,a,6es12.4)') 'SSMS-ABORT el',noel,npt,' sig0',sig0,' dstran',dstran(1:6)
      end if

!     ---- tangent: elastic moduli at the updated state ----
      p = (sig(1) + sig(2) + sig(3)) / 3.0d0
      call GetElasticModuli1(-sig, sv(21), PR, E_K, E_G)
      if (-p < pmin) then
        call GetElasticModuli1((/pmin,pmin,pmin,0d0,0d0,0d0/), sv(21),  &
                               PR, E_K, E_G)
      end if
      call GetStiffness(E_K, E_G, D)
      ddsdde = D(1:ntens,1:ntens)

      stress(1:6) = sig
      statev(1:min(NSV,nstatv)) = sv(1:min(NSV,nstatv))

      if (noel >= 1 .and. noel <= MAXE .and. npt >= 1 .and. npt <= MAXG) then
        stdump(1:6, npt, noel)     = sig
        stdump(7:NDV, npt, noel)   = sv
        idump(noel) = 1
      end if
      return
      end subroutine umat

!=======================================================================
      subroutine uexternaldb(lop, lrestart, time, dtime, kstep, kinc)
      implicit none
      integer lop, lrestart, kstep, kinc
      double precision time(2), dtime
      integer, parameter :: NSV = 31, MAXE = 20000, MAXG = 8, NDV = 6+NSV
      double precision stdump(NDV, MAXG, MAXE)
      integer idump(MAXE)
      common /SSMSDUMP/ stdump, idump
      character(256) outdir, jobname, fname
      integer lenoutdir, lenjobname, ie, ig, iu

      if (lop == 0) then
        idump = 0
        stdump = 0.0d0
      end if
!     lop=3: end of analysis -> write the last converged state
      if (lop == 3) then
        call getoutdir(outdir, lenoutdir)
        call getjobname(jobname, lenjobname)
        fname = outdir(1:lenoutdir)//'\'//jobname(1:lenjobname)//'_state.txt'
        iu = 157
        open(iu, file=fname, status='replace')
        write(iu,'(a)') '# elem gp  S11 S22 S33 S12 S23 S13  StVar(1..31)  (UEL Voigt order, tension +)'
        do ie = 1, MAXE
          if (idump(ie) /= 1) cycle
          do ig = 1, MAXG
            write(iu,'(i7,1x,i2,37(1x,es17.9e3))') ie, ig, stdump(:, ig, ie)
          end do
        end do
        close(iu)
      end if
      return
      end subroutine uexternaldb

! C-callable shim around the SAME Abaqus UMAT wrapper used in the pile model
subroutine elem_umat(stress, statev, ddsdde, stran, dstran, dtime, props, nprops, noel, npt, kinc, pnewdt) &
           bind(C, name='elem_umat')
  use iso_c_binding
  !DEC$ ATTRIBUTES DLLEXPORT :: elem_umat
  implicit none
  integer(c_int), value :: nprops, noel, npt, kinc
  real(c_double) :: stress(6), statev(31), ddsdde(6,6), stran(6), dstran(6), props(nprops), pnewdt
  real(c_double), value :: dtime
  character(80) :: cmname
  double precision :: sse, spd, scd, rpl, drpldt, temp, dtemp, celent, time(2), predef(1), dpred(1), &
                      ddsddt(6), drplde(6), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
  integer :: kstep, layer, kspt
  cmname = 'SANISANDMS'; time = 0d0; kstep = 1; layer = 1; kspt = 1
  sse=0; spd=0; scd=0; rpl=0; drpldt=0; temp=0; dtemp=0; celent=1
  drot=0; dfgrd0=0; dfgrd1=0; coords=0; ddsddt=0; drplde=0; predef=0; dpred=0
  call umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,drplde,drpldt,stran,dstran,time,dtime,temp, &
            dtemp,predef,dpred,cmname,3,3,6,31,props,nprops,coords,drot,pnewdt,celent,dfgrd0,dfgrd1, &
            noel,npt,layer,kspt,kstep,kinc)
end subroutine
! stubs for Abaqus utilities referenced by UEXTERNALDB (never called here)
subroutine getoutdir(d, l)
  character(*) d; integer l
  d = '.'; l = 1
end subroutine
subroutine getjobname(d, l)
  character(*) d; integer l
  d = 'elem'; l = 4
end subroutine
