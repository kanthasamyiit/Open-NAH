!--------------------------------------------------------------------------------------------------------------------------
!
!	          Author: Dr. Kanthasamy Chelliah,
!
!             Amazon & Illinois Institute of Technology
!
!             Nearfield Acoustic Holography (NAH) Using DFT,ESM,BEM,SONAH
!
!--------------------------------------------------------------------------------------------------------------------------

! Copyright (C) 2018  kanthasamy chelliah

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!                          ESM
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE esm(Me,Nx,Ny,NmX,NmY,NmicNode,NNode,FNAME,HF,Xres,XresMic,ZH,&
pSE,pSg,pSt,pSt2,pSt3,pSt4,pSt5,pSl,pH,vSE,vSg,vSt,vSt2,vSt3,vSt4,vSt5,&
vSl,Sigma,U,VT)
!    implicit none

    double precision, INTENT(IN) :: Xres,XresMic,ZH
    integer, INTENT(IN) :: Nx,Ny,NmX,NmY,NmicNode,NNode,HF,Me
    CHARACTER(LEN=100), INTENT(IN) :: FNAME
    complex*16, DIMENSION(NNode), INTENT(OUT) :: pSE,pSg,pSt,pSt2, &
    pSt3,pSt4,pSt5,pSl
    complex*16, DIMENSION(NNode), INTENT(OUT) :: vSE,vSg,vSt,vSt2, &
    vSt3,vSt4,vSt5,vSl

    double precision, PARAMETER :: ZS0 = 0.0
    double precision :: ZSv
    double precision, PARAMETER :: lambda = 0.0            ! Tikhonov parameter
    double precision, PARAMETER :: lambda_tikh =  sqrt(3.0)/2.0   ! Tikhonov parameter
    double precision, PARAMETER :: lambda_cg =  0.084094337148796

    double precision, PARAMETER :: c = 343.0             ! Speed of sound
    double precision, PARAMETER :: rho = 1.225             ! Density of medium
    double precision, PARAMETER :: pi = 3.141592653589793
    complex*16, PARAMETER :: MINUS_ONE = -1.0

    double precision, DIMENSION(NNode , 3) :: SourLoc  ! Coords of Nodes in the mesh
    double precision, DIMENSION(NNode , 3) :: VirLoc   ! Coords of Nodes in the mesh
    double precision, DIMENSION(NmicNode , 3) :: MicLoc ! Coords of present Mic
    double precision, dimension(NmicNode), INTENT(OUT) :: Sigma

    complex*16, DIMENSION(NmicNode), INTENT(OUT) :: pH    ! Pressure at Hologram
!    complex*16, DIMENSION(NNode) :: pS         ! Pressure at Source
    complex*16, DIMENSION(NNode) :: qv         ! Source strength at virtual surface

    complex*16, DIMENSION(NmicNode, NNode) :: EHS      ! EHS Matrix
    complex*16, DIMENSION(NNode, NNode) :: GHS,GDHS      ! DHS Matrix
    complex*16, DIMENSION(NmicNode, NmicNode), INTENT(OUT) :: U
    complex*16, DIMENSION(NNode, NNode), INTENT(OUT) :: VT

!    complex*16, DIMENSION(NmY,NmX) :: p

    double precision :: k, lambda_min, Yres, YresMic, &
    lambda_tikh2,lambda_tikh3,lambda_tikh4,lambda5
    double precision, DIMENSION(3) :: SL,PL

    CHARACTER(LEN=100) :: FNAME1

!    double precision, DIMENSION(LengthTS) :: TS
!    double precision, PARAMETER :: dT = 1/(real(Fs,kind=8));
!    INTEGER, PARAMETER :: icas = 1;
!    INTEGER, PARAMETER :: nroll = 1;

    INTEGER i,j,M,n
    complex*16 g,dg
    REAL*8, PARAMETER :: Pref = 0.00002;




!    write(*,*) 'Enter Helmholtz Frequency (Integer)'
!    read(*,*) HF
    write(*,*) 'Starting ===>'
!write(*,*) 'Inializing Elements and Nodes'


Yres = Xres
YresMic = XresMic
ZSv = -ZH / 2.0


k = 2.0 * pi * HF / c

if (Me == 1) then
    call readmems(FNAME, HF, pH)
else
    OPEN (66,FILE=FNAME,STATUS='OLD')
    READ(66,*) (pH(i),i=1,NmicNode)
    CLOSE(66)
endif

!! Node Matrix and Mic coordinate Matrix in m

DO j=1,NX
    DO i=1,NY
        VirLoc(NY*(j-1)+i , 1 ) = - Xres * ( (NX-1) / 2. ) + (j-1) * Xres
        VirLoc(NY*(j-1)+i , 2 ) = - Yres * ( (NY-1) / 2. ) + (i-1) * Yres
        VirLoc(NY*(j-1)+i , 3 ) = ZSv
    ENDDO
ENDDO

DO j=1,NX
    DO i=1,NY
        SourLoc(NY*(j-1)+i , 1 ) = - Xres * ( (NX-1) / 2. ) + (j-1) * Xres
        SourLoc(NY*(j-1)+i , 2 ) = - Yres * ( (NY-1) / 2. ) + (i-1) * Yres
        SourLoc(NY*(j-1)+i , 3 ) = ZS0
    ENDDO
ENDDO

DO j=1,NmX
    DO i=1,NmY
        MicLoc(NmY*(j-1)+i , 1 ) = - XresMic * ( (NmX-1) / 2. ) + (j-1) * XresMic
        MicLoc(NmY*(j-1)+i , 2 ) = - YresMic * ( (NmY-1) / 2. ) + (i-1) * YresMic
        MicLoc(NmY*(j-1)+i , 3 ) = ZH
    ENDDO
ENDDO

do m = 1,NmicNode

    PL = MicLoc(m,:)

    do n = 1,NNode

        SL = VirLoc(n,:)
!        r =  sqrt( ( PL(1) - SL(1) )**2 + ( PL(2) - SL(2) )**2 + ( PL(3) - SL(3) )**2 )
        call green(g,PL,SL,k)
        EHS(m,n) = g

    enddo
enddo

write(*,*) 'Calling SVD'

call svd (EHS,U,VT,Sigma,NmicNode,NNode)

U = conjg( transpose(U) )
VT = conjg( transpose(VT) )

!write(*,*) 'Virtual Surface to Source Surface Transformation'

do m = 1,NNode

    PL = VirLoc(m,:)

    do n = 1,NNode

        SL = SourLoc(n,:)
        call green(g,PL,SL,k)
        GHS(m,n) = g

        call dgreen(dg,PL,SL,k)
        GDHS(m,n) = dg / ( sqrt(MINUS_ONE) * rho * c * k * 2.0 )

    enddo
enddo

call svdsolve(U,VT,Sigma,NmicNode,NNode,pH,qv)
pSE = matmul(GHS,qv)
vSE = matmul(GDHS,qv)

 lambda_tikh2 = 0.01 * Sigma(1)
 lambda_tikh3 = 0.001 * Sigma(1)
 lambda_tikh4 = 0.0001 * Sigma(1)
 lambda5 = 0.00001 * Sigma(1)

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh)
pSt = matmul(GHS,qv)
vSt = matmul(GDHS,qv)

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh2)
pSt2 = matmul(GHS,qv)
vSt2 = matmul(GDHS,qv)

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh3)
pSt3 = matmul(GHS,qv)
vSt3 = matmul(GDHS,qv)

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh4)
pSt4 = matmul(GHS,qv)
vSt4 = matmul(GDHS,qv)

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda5)
pSt5 = matmul(GHS,qv)
vSt5 = matmul(GDHS,qv)

call gcv(EHS,U,VT,Sigma,NmicNode,NNode,pH,lambda_min)
write(*,*) 'Lamda_min:',lambda_min
call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_min)
pSg = matmul(GHS,qv)
vSg = matmul(GDHS,qv)

FNAME1='SigmaE.txt'
OPEN (47,FILE=FNAME1,STATUS='REPLACE')
write(47,*) (Sigma(i),i=1,NmicNode)
write(47,*) lambda_min

call lcurve(EHS,U,VT,Sigma,NmicNode,NNode,pH,lambda_min)
write(*,*) 'Lamda_min:',lambda_min
call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_min)
pSl = matmul(GHS,qv)
vSl = matmul(GDHS,qv)

write(47,*) lambda_min
CLOSE(47)
END Subroutine esm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!                          SONAH
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE sonah(Me,Nx,Ny,NmX,NmY,NmicNode,NNode,FNAME,HF,Xres,XresMic,D,&
pS,pSg,pSt,pSt1,pSt2,pSt3,pSt4,pSt5,pSl,vS,vSg,vSt,vSt1,vSt2,vSt3,vSt4,vSt5,vSl)
!    implicit none

    double precision, INTENT(IN) :: Xres,XresMic,D
    integer, INTENT(IN) :: Nx,Ny,NmX,NmY,NmicNode,NNode,HF,Me
    CHARACTER(LEN=100), INTENT(IN) :: FNAME
    complex*16, DIMENSION(NNode), INTENT(OUT) :: pS,pSg,pSt,pSt1,pSt2,&
    pSt3,pSt4,pSt5,pSl
    complex*16, DIMENSION(NNode), INTENT(OUT) :: vS,vSg,vSt,vSt1,vSt2,&
    vSt3,vSt4,vSt5,vSl

    double precision :: ZS,ZV
    double precision, PARAMETER :: ZH = 0.0

    double precision, PARAMETER :: lambda = 0.0 ! Tikhonov parameter
    double precision, PARAMETER :: lambda_tikh = sqrt(3.0)/2.0 ! Tikhonov parameter

    double precision, PARAMETER :: c = 343.0 ! Speed of sound
    double precision, PARAMETER :: rho = 1.225 ! Density of medium
    double precision, PARAMETER :: pi = 3.141592653589793
    complex*16, PARAMETER :: MINUS_ONE = -1.0
    COMPLEX*16, PARAMETER :: sqi = sqrt(MINUS_ONE)
    REAL*8, PARAMETER :: Pref = 0.00002;

    double precision, DIMENSION(NmicNode) :: xij,yij ! Coords of Nodes in the mesh
    double precision, DIMENSION(NNode) :: xsij,ysij ! Coords of Nodes in the mesh


    complex*16, DIMENSION(NmicNode) :: pH ! Pressure at Hologram

    complex*16, DIMENSION(NmicNode) :: q, qt, qt1,qt2,qt3,qt4,qt5, qg, ql !, AHLR, AHBR ! AH * alpha(r)
    complex*16, DIMENSION(NmicNode,NNode) :: AHLR, AHBR ! AH * alpha(r)

    complex*16, DIMENSION(NmicNode, NmicNode) :: AHA !,IAHA      ! AHA Matrix
    double precision, dimension(NmicNode) :: Sigma
    complex*16, DIMENSION(NmicNode, NmicNode) :: U
    complex*16, DIMENSION(NmicNode, NmicNode) :: VT

    double precision :: k, lambda_min, Yres, YresMic, Ri,Rij,theta,divi
    double precision :: lambda_tikh2,lambda_tikh3,lambda_tikh4,lambda5,lambda_tikh1

!    double precision, DIMENSION(LengthTS) :: TS
!    double precision, PARAMETER :: dT = 1/(real(Fs,kind=8));
!    INTEGER, PARAMETER :: icas = 1;
!    INTEGER, PARAMETER :: nroll = 1;

!    Parameters for numerical integration

    INTEGER, PARAMETER :: ngp = 1000;
    double precision, DIMENSION(ngp) :: xx ! gauss quad absi
    double precision, DIMENSION(ngp) :: w ! gauss quad weight

    INTEGER, PARAMETER :: ngp1 = 1000;
    double precision, DIMENSION(ngp) :: xx1 ! gauss quad absi
    double precision, DIMENSION(ngp) :: w1 ! gauss quad weight

    double precision, PARAMETER :: l1 = 0.0
    double precision, PARAMETER :: l2 = pi/2.0 ! limits of integration I
    COMPLEX*16 :: IIp,Ip,Iv,IIv
    CHARACTER(LEN=100) :: FNAME1

    REAL*8, PARAMETER :: F0 = 1.0;


    INTEGER i,j,M,n


    write(*,*) 'Starting ===>'
ZS = -D
ZV = -Xres+ZS

Yres = Xres
YresMic = XresMic

k = 2.0 * pi * HF / c


write(*,*) 'Reading the data'

if (Me == 1) then
    call readmems(FNAME, HF, pH)
else
    OPEN (66,FILE=FNAME,STATUS='OLD')
    READ(66,*) (pH(i),i=1,NmicNode)
    CLOSE(66)
endif


DO j=1,NmX
    DO i=1,NmY
        xij(NmY*(j-1)+i) = - XresMic * ( (NmX-1) / 2. ) + (j-1) * XresMic
        yij(NmY*(j-1)+i) = - YresMic * ( (NmY-1) / 2. ) + (i-1) * YresMic
    ENDDO
ENDDO

DO j=1,NX
    DO i=1,NY
        xsij(NY*(j-1)+i) = - Xres * ( (NX-1) / 2. ) + (j-1) * Xres
        ysij(NY*(j-1)+i) = - Yres * ( (NY-1) / 2. ) + (i-1) * Yres
    ENDDO
ENDDO

!!*********************************************************************************************************************
!!
!!           AH * A calculation
!!
!!********************************************************************************************************************
write(*,*) 'Hologram to virtual surface Transformation'

    divi = k*2.0*(ZH-ZV)

call gauleg(ngp, xx, w) ! Discretizing the integration field
call l_quadrature_rule(ngp1, xx1, w1) ! Discretizing the integration field in an element



!do i = 1,ngp        ! changing the abscissas and weights to up to infinite limit
!
!        xx2(i) = tan(pi/4.0*(1+xx(i)))
!        w2(i) = w(i)*pi / (4.0* ( cos(pi/4.0*(1+xx(i))) )**2 )
!
!enddo

do m = 1,NmicNode

    do n = 1,NmicNode

        Rij = sqrt( ( xij(m) - xij(n) )**2 + ( yij(m) - yij(n) )**2 )

        Ip = 0.0
        IIp = 0.0

        do i = 1,ngp

            theta = xx(i) * (l2-l1)/2.0 +  (l2+l1)/2.0

Ip = Ip + (l2-l1)/2.0 * w(i) * bessel_j0( k*Rij*sin(theta) ) * sin(theta) * cos(theta)

        enddo


        do i = 1,ngp1

IIp = IIp + w1(i) / divi * bessel_j0( k*Rij * sqrt(1.0+(xx1(i)/divi)**2) ) * &
    (xx1(i)/divi) / (1.0+(xx1(i)/divi)**2)

        enddo

        AHA(m,n) = Ip+IIp

    enddo

enddo

write(*,*) 'Calling SVD'

call svd(AHA,U,VT,Sigma,NmicNode,NmicNode)

U = conjg( transpose(U) )
VT = conjg( transpose(VT) )

lambda_tikh1 = Sigma(1) * sqrt(3.0) / 2.0 / 500


lambda_tikh2 = 0.01 * Sigma(1)
lambda_tikh3 = 0.001 * Sigma(1)
lambda_tikh4 = 0.0001 * Sigma(1)
lambda5 = 0.00001 * Sigma(1)

call svdsolve(U,VT,Sigma,NmicNode,NmicNode,pH,q)   ! Regular

call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qt,lambda_tikh)   ! fixed Tikhonov

call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qt1,lambda_tikh1)   ! fixed Tikhonov
call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qt2,lambda_tikh2)   ! fixed Tikhonov
call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qt3,lambda_tikh3)   ! fixed Tikhonov
call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qt4,lambda_tikh4)   ! fixed Tikhonov
call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qt5,lambda5)   ! fixed Tikhonov


call gcv(AHA,U,VT,Sigma,NmicNode,NmicNode,pH,lambda_min)   ! GCV Tikhonov
FNAME1='SigmaS.txt'
OPEN (47,FILE=FNAME1,STATUS='REPLACE')
write(47,*) (Sigma(i),i=1,NmicNode)
write(47,*) lambda_min
call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,qg,lambda_min)


call lcurve(AHA,U,VT,Sigma,NmicNode,NmicNode,pH,lambda_min)   ! GCV Tikhonov
write(47,*) lambda_min
CLOSE(47)
call tikhonov(U,VT,Sigma,NmicNode,NmicNode,pH,ql,lambda_min)

!call gauleg(ngp, xx, w)       ! Discretizing the integration field

    AHLR = 0.0
    AHBR = 0.0

divi = k*(ZH+ZS-2.0*ZV)

do m = 1,NNode

!write(*,*) 'SS for Mic:', m

    do n = 1,NmicNode

        Ri =  sqrt( ( xsij(m) - xij(n) )**2 + ( ysij(m) - yij(n) )**2 )


        Ip = 0.0
        IIp = 0.0
        Iv = 0.0
        IIv = 0.0

        do i = 1,ngp

            theta = xx(i) * (l2-l1)/2.0 +  (l2+l1)/2.0

        Ip = Ip + (l2-l1)/2.0 * w(i) * bessel_j0( k*Ri*sin(theta) ) * exp(sqi*k*(ZH-ZS)*cos(theta)) * &
             sin(theta) * cos(theta)

        Iv = Iv +  (l2-l1)/2.0 * w(i) * bessel_j0( k*Ri*sin(theta) ) * exp(sqi*k*(ZH-ZS)*cos(theta)) * &
             sin(theta) * (cos(theta))**2

        enddo


        do i = 1,ngp1

        IIp = IIp +  w1(i) / divi * bessel_j0( k*Ri * sqrt(1+(xx1(i)/divi)**2) ) *  &
    (xx1(i)/divi) / (1.0+(xx1(i)/divi)**2)

        IIv = IIv +  w1(i) /divi * bessel_j0( k*Ri * sqrt(1+(xx1(i)/divi)**2) ) * &
    (xx1(i)/divi)**2 / (1.0+(xx1(i)/divi)**2)

        enddo

        AHLR(n,m) = Ip+IIp
        AHBR(n,m) = Iv+IIv

    enddo

!    pS(m) = dot_product(q,AHLR)
!    pSt(m) = dot_product(qt,AHLR)
!    pSg(m) = dot_product(qg,AHLR)
!    pSm(m) = dot_product(qm,AHLR)
!
!    vS(m) = dot_product(q,AHBR)
!    vSt(m) = dot_product(qt,AHBR)
!    vSg(m) = dot_product(qg,AHBR)
!    vSm(m) = dot_product(qm,AHBR)

enddo

    pS = matmul(q,AHLR)
    pSt = matmul(qt,AHLR)
    pSt1 = matmul(qt1,AHLR)
    pSg = matmul(qg,AHLR)
    pSl = matmul(ql,AHLR)
!    pSm = matmul(qm,AHLR)

    vS = matmul(q,AHBR)
    vSt = matmul(qt,AHBR)
    vSt1 = matmul(qt1,AHBR)
    vSg = matmul(qg,AHBR)
    vSl = matmul(ql,AHBR)
!    vSm = matmul(qm,AHBR)

    pSt2 = matmul(qt2,AHLR)
    vSt2 = matmul(qt2,AHBR)
    pSt3 = matmul(qt3,AHLR)
    vSt3 = matmul(qt3,AHBR)
    pSt4 = matmul(qt4,AHLR)
    vSt4 = matmul(qt4,AHBR)
    pSt5 = matmul(qt5,AHLR)
    vSt5 = matmul(qt5,AHBR)

    vS = vS/rho/c
    vSg = vSg/rho/c
    vSl = vSl/rho/c
    vSt = vSt/rho/c
    vSt1 = vSt1/rho/c
    vSt2 = vSt2/rho/c
    vSt3 = vSt3/rho/c
    vSt4 = vSt4/rho/c
    vSt5 = vSt5/rho/c

END SUBROUTINE sonah

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!                          BEM
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE bem(Me,Nx,Ny,NmX,NmY,NmicNode,NNode,FNAME,HF,Xres,XresMic,ZH,&
pSL,pSLg,pSLt,pSLt1,pSLt2,pSLt3,pSLt4,pSLt5,pSLl,vSL,vSLg,vSLt,vSLt1,&
vSLt2,vSLt3,vSLt4,vSLt5,vSLl)
!    implicit none

    double precision, INTENT(IN) :: Xres,XresMic,ZH
    integer, INTENT(IN) :: Nx,Ny,NmX,NmY,NmicNode,NNode,HF,Me
    CHARACTER(LEN=100), INTENT(IN) :: FNAME
    complex*16, DIMENSION(NNode), INTENT(OUT) :: pSL,pSLg,pSLt,pSLt2,&
    pSLt3,pSLt4,pSLt5,pSLl,pSLt1
    complex*16, DIMENSION(NNode), INTENT(OUT) :: vSL,vSLg,vSLt,vSLt2,&
    vSLt3,vSLt4,vSLt5,vSLl,vSLt1

    double precision, PARAMETER :: ZS0 = 0.0           ! Source surface Z coordinate

    double precision, PARAMETER :: lambda = 0.0             ! Tikhonov parameter
    double precision, PARAMETER :: lambda_tikh = sqrt(3.0)/2.0             ! Tikhonov parameter
!    double precision, PARAMETER :: lambda_tikh1 = 0.024612673536541976 * sqrt(3.0) / 2.0 / 500            ! Tikhonov parameter
    double precision, PARAMETER :: c = 343.0             ! Speed of sound
    double precision, PARAMETER :: rho = 1.225             ! Density of medium
    double precision, PARAMETER :: pi = 3.141592653589793
    complex*16, PARAMETER :: MINUS_ONE = -1.0

    double precision,dimension(4) :: N
    double precision :: Jac
    complex*16 :: g, dg, gj, dgj


    INTEGER*8, DIMENSION(((NX-1)*(NY-1) ), 4)  :: Elem         ! Element connectivity
    double precision, DIMENSION(NNode , 3) :: Node         ! Coords of Nodes in the mesh
    double precision, DIMENSION(NmicNode , 3) :: NodeMic         ! Coords of present Mic
    double precision, DIMENSION(4,2) :: NL         ! Coords of 4 nodes corresponding the present element
    double precision, dimension(NmicNode) :: Sigma

    complex*16, DIMENSION(NmicNode) :: pH         ! Pressure at Hologram
    complex*16, DIMENSION(NNode) :: qv         ! Pressure at Source

    complex*16, DIMENSION(NmicNode, NNode) :: DHS, SHS, DHS1      ! DHS Matrix
    complex*16, DIMENSION(NNode, NNode) :: DSS, SSS     ! DHS Matrix
    complex*16, DIMENSION(NmicNode, NmicNode) :: U      ! U Matrix
    complex*16, DIMENSION(NNode, NNode) :: VT      ! VT Matrix

    double precision, DIMENSION(NX) :: X0
    double precision, DIMENSION(NY) :: Y0
    double precision, DIMENSION(NmX) :: X01
    double precision, DIMENSION(NmY) :: Y01
    double precision :: xp, yp, zp, zs, k, Yres,YresMic,lambda_min, Xm0, Ym0, &
    lambda_tikh2,lambda_tikh3,lambda_tikh4,lambda5
    double precision :: lambda_tikh1
    CHARACTER(LEN=100) :: FNAME1
    INTEGER i,j,M,e,l, NElem


! Parameters for Gauss Integration
    INTEGER, PARAMETER :: ngp = 6;
    double precision, DIMENSION(ngp) :: xx         ! gauss quad absi
    double precision, DIMENSION(ngp) :: w         ! gauss quad weight

! Parameters for blas ZGEMM
    complex*16, PARAMETER :: alpha = 1.0
    complex*16, PARAMETER :: beta = 0.0
    character,parameter :: TRANS = 'N'
    integer,parameter :: NRHS = 1
!    character :: TRANSA, TRANSB


write(*,*) 'Inializing Elements and Nodes'

Yres = Xres
YresMic = XresMic
NElem = (NX-1)*(NY-1)
Xm0 = ( NX - NmX ) * Xres / 2.0
Ym0 = ( NY - NmY ) * Yres / 2.0

k = 2.0 * pi * HF / c

    DO i = 1,NX
        X0(i) = 0.0 + (i-1) * Xres / 0.0254
    ENDDO
    DO i = 1,NY
        Y0(i) = 0.0 + (i-1) * Yres / 0.0254
    ENDDO

!    X01 = X0
!    Y01 = Y0
    DO i = 1,NmX
        X01(i) = Xm0 + (i-1) * XresMic / 0.0254
    ENDDO
    DO i = 1,NmY
        Y01(i) = Ym0 + (i-1) * YresMic / 0.0254
    ENDDO

!! Node Matrix and Mic coordinate Matrix in m

DO j=1,NX
    DO i=1,NY
        Node(NY*(j-1)+i , 1 ) = - Xres * ( (NX-1) / 2. ) + (j-1) * Xres
        Node(NY*(j-1)+i , 2 ) = - Yres * ( (NY-1) / 2. ) + (i-1) * Yres
        Node(NY*(j-1)+i , 3 ) = ZS0
    ENDDO
ENDDO

DO j=1,NmX
    DO i=1,NmY
        NodeMic(NmY*(j-1)+i , 1 ) = - XresMic * ( (NmX-1) / 2. ) + (j-1) * XresMic
        NodeMic(NmY*(j-1)+i , 2 ) = - YresMic * ( (NmY-1) / 2. ) + (i-1) * YresMic
        NodeMic(NmY*(j-1)+i , 3 ) = ZH
    ENDDO
ENDDO


! Element Connectivity Matrix

DO j=1, (NX-1)
    DO i=1, (NY-1)

        E = (NY-1)  * (j-1) + i

        Elem (E, 1) = Ny * (j-1) + i
        Elem (E, 2) = Ny * (j-1) + i + 1
        Elem (E, 3) = Ny * j     + i + 1
        Elem (E, 4) = Ny * j     + i

    ENDDO
ENDDO

write(*,*) 'Reading the FFT processed data'

if (Me == 1) then
    call readmems(FNAME, HF, pH)
else
    OPEN (66,FILE=FNAME,STATUS='OLD')
    READ(66,*) (pH(i),i=1,NmicNode)
    CLOSE(66)
endif

call gauleg(ngp, xx, w)       ! Discretizing the integration field in an element

write(*,*) 'Building the HS transformation matrices'

DHS = 0.0
SHS = 0.0


do m = 1,NmicNode

!write(*,*) 'HS for Mic:', m

    xp = NodeMic ( m, 1 )
    yp = NodeMic ( m, 2 )
    zp = NodeMic ( m, 3 )

    do e = 1,NElem

        do i = 1,4
            NL(i,1) = Node ( Elem(e,i), 1 )
            NL(i,2) = Node ( Elem(e,i), 2 )
        enddo

        do i = 1,ngp
            do j = 1,ngp

                call shapeGreenL(N, Jac, G, DG, xx(i),xx(j),zs, xp,yp,zp, NL, k)

                GJ = G * cmplx( w(i) * w(j) * Jac / ngp , kind=8 )
                DGJ = DG * cmplx( w(i) * w(j) * Jac / ngp , kind=8 )

                do l = 1,4

                    DHS ( m, Elem(e,l) ) = DHS ( m, Elem(e,l) ) + DGJ * cmplx( N(l) , kind=8 )
                    SHS ( m, Elem(e,l) ) = SHS ( m, Elem(e,l) ) + GJ * cmplx( 2.0 * N(l) , kind=8 )

                enddo

            enddo
        enddo

    enddo
enddo


DSS = 0.0
SSS = 0.0

do m = 1,NNode

!write(*,*) 'SS for Source:', m

    xp = Node ( m, 1 )
    yp = Node ( m, 2 )
    zp = Node ( m, 3 )

    do e = 1,NElem

        do i = 1,4
            NL(i,1) = Node ( Elem(e,i), 1 )
            NL(i,2) = Node ( Elem(e,i), 2 )
        enddo

        do i = 1,ngp
            do j = 1,ngp

                call shapeGreenL(N, Jac, G, DG, xx(i),xx(j),zs, xp,yp,zp, NL, k)

                GJ = G * cmplx( w(i) * w(j) * Jac / ngp , kind=8 )
                DGJ = DG * cmplx( w(i) * w(j) * Jac / ngp , kind=8 )

                do l = 1,4

                    DSS ( m, Elem(e,l) ) = DSS ( m, Elem(e,l) ) + DGJ * cmplx( N(l), kind=8 )
                    SSS ( m, Elem(e,l) ) = SSS ( m, Elem(e,l) ) + GJ * cmplx( N(l), kind=8 )

                enddo

            enddo
        enddo

    enddo
enddo


write(*,*) 'Rearranging the Matrices'

do i = 1,NNode
    DSS(i,i) = DSS(i,i) - 0.5
enddo

do i = 1,NNode
    DSS(i,i) = 1.0/DSS(i,i)
enddo


!write(*,*) DSS(1,1),DSS(15,15)

write(*,*) '1'

!DHS = -2.0 * DHS  ! Multiplying by -2.0 instead of multiplying by DSS (it's a diagonal matrix)

DHS = matmul(DHS,DSS)

!call ZGEMM ( TRANS, TRANS, NmicNode, NNode, NNode, alpha, DHS, NmicNode, DSS, NNode, beta, DHS, NmicNode )
write(*,*) '2'

write(*,*) DHS(1,2),SSS(1,2),SHS(1,2)

DHS = matmul(DHS, SSS)

!write(*,*) DHS(1,2)

!call ZGEMM ( TRANS, TRANS, NmicNode, NNode, NNode, alpha, DHS, NmicNode, SSS, NNode, beta, DHS, NmicNode )
write(*,*) '3'

!write(*,*) DHS

DHS = DHS - SHS
!DHS = DHS
DHS1 = DHS

SSS = matmul(DSS,SSS)


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write(*,*) 'Decomposing the Matrices'

call svd (DHS1,U,VT,Sigma,NmicNode,NNode)

U = conjg( transpose(U) )
VT = conjg( transpose(VT) )


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda)

vSL = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSL = matmul(SSS,qv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 lambda_tikh2 = 0.01 * Sigma(1)
 lambda_tikh3 = 0.001 * Sigma(1)
 lambda_tikh4 = 0.0001 * Sigma(1)
 lambda5 = 0.00001 * Sigma(1)

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh)

vSLt = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLt = matmul(SSS,qv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_tikh1 = Sigma(1) * sqrt(3.0) / 2.0 / 500

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh1)

vSLt1 = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLt1 = matmul(SSS,qv)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh2)

vSLt2 = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLt2 = matmul(SSS,qv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh3)

vSLt3 = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLt3 = matmul(SSS,qv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_tikh4)

vSLt4 = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLt4 = matmul(SSS,qv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda5)

vSLt5 = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLt5 = matmul(SSS,qv)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call gcv(DHS,U,VT,Sigma,NmicNode,NNode,pH,lambda_min)
call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_min)

FNAME1='SigmaB.txt'
OPEN (47,FILE=FNAME1,STATUS='REPLACE')
write(47,*) (Sigma(i),i=1,NmicNode)
write(47,*) lambda_min

vSLg = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLg = matmul(SSS,qv)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call lcurve(DHS,U,VT,Sigma,NmicNode,NNode,pH,lambda_min)
call tikhonov(U,VT,Sigma,NmicNode,NNode,pH,qv,lambda_min)
write(47,*) lambda_min
close(47)

vSLl = qv / ( sqrt(MINUS_ONE) * rho * c * k )
pSLl = matmul(SSS,qv)

END subroutine bem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!                          DFT
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine nah(Me,NmicX,NmicY,NmicNode,FNAME,HF,Xres,D,pDilte,pvz1)

    double precision, INTENT(IN) :: Xres,D
    integer, INTENT(IN) :: NmicX,NmicY,HF,Me,NmicNode
    CHARACTER(LEN=100), INTENT(IN) :: FNAME
    integer :: Nmea,NmeaY

    double precision, PARAMETER :: ZS0 = 0.0

    INTEGER, PARAMETER :: ZPN = 2

    REAL*8, PARAMETER :: alpha = 0.5

    REAL*8, PARAMETER :: c = 343.0
    REAL*8, PARAMETER :: mu = 1.225

    REAL*8, PARAMETER :: pi = 3.141592653589793
    COMPLEX*16, PARAMETER :: MINUS_ONE = -1.0
    COMPLEX*16, PARAMETER :: sqi = sqrt(MINUS_ONE)


!    INTEGER Xmax,Ymax
    INTEGER i,j,l,o,nx,ny,mx,my,NmicXz,NmicYz

    REAL*8 Fsx,Fsy,k2,Yres,kr,Lx,Ly,kc,SNR
!    REAL*8, DIMENSION(Nfft) :: f
    REAL*8 :: k

    REAL*8, DIMENSION(NmicX*ZPN) :: kx
    REAL*8, DIMENSION(NmicY*ZPN) :: ky

    REAL*8, DIMENSION(NmicX*ZPN) :: kx2
    REAL*8, DIMENSION(NmicY*ZPN) :: ky2

    COMPLEX*16, DIMENSION(NmicY*ZPN,NmicX*ZPN) :: kz, Kwindow


    REAL*8, DIMENSION(NmicX)  :: X
    REAL*8, DIMENSION(NmicY)  :: Y
    COMPLEX*16, DIMENSION(NmicY*NmicX) :: pH
    COMPLEX*16, DIMENSION(NmicY,NmicX),INTENT(OUT) :: pDilte
    COMPLEX*16, DIMENSION(NmicY,NmicX),INTENT(OUT) :: pvz1 !,pvx,pvy
    COMPLEX*16, DIMENSION(NmicY*ZPN,NmicX*ZPN) :: pvz,pvzhat
    COMPLEX*16, DIMENSION(NmicY*ZPN,NmicX*ZPN) :: pDzp,phat

    write(*,*) 'Starting ===>'

    Yres = Xres
    Nmea = NmicX
    NmeaY = NmicY
    NmicXz = NmicX * ZPN
    NmicYz = NmicY * ZPN

    DO i = 1,NmicX
        X(i) = 0.0 + (i-1) * Xres / 0.0254;
    ENDDO
    DO i = 1,NmicY
        Y(i) = 0.0 + (i-1) * Yres / 0.0254;
    ENDDO

    SNR = 28

    k = 2.0*pi*real(HF,kind=8)/c
    k2 = k * k

kc = sqrt( k**2 + ( SNR / ( 8.7*D ) )**2 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Lx = (NmicX -1) * Xres;
    Ly = (NmicY -1) * Yres;

    Fsx = 1.0 / Xres
    Fsy = 1.0 / Yres

    DO i = 1,NmicXz
        kx(i) = - pi * Fsx + (i-1) * 2.0 * pi * Fsx / real ( NmicXz-1 ) ;
    ENDDO

    DO i = 1,NmicYz
        ky(i) = - pi * Fsy + (i-1) * 2.0 * pi * Fsy / real ( NmicYz-1 );
    ENDDO

    DO i = 1, NmicXz
        kx2(i) = kx(i) * kx(i)
    ENDDO
    DO i = 1, NmicYz
        ky2(i) = ky(i) * ky(i)
    ENDDO


!  write(*,*) kx

    DO j = 1, NmicXz
        DO i = 1, NmicYz
            kz(i,j) = SQRT( cmplx ( k2 - kx2(j) - ky2(i) , kind = 8) )
        ENDDO
    ENDDO

!write(*,*) kz

Do i=1,NmicYz
    Do j=1,NmicXz
        kr=abs(sqrt(kx2(j)+ky2(i)));
        if(kr <= kc) then
            Kwindow(i,j)=cmplx(1.0-(0.5*exp(((kr/kc)-1.0)/alpha)), kind = 8)
        else
            Kwindow(i,j)=cmplx(0.5*exp((1.0-(kr/kc))/alpha), kind = 8)
        endif
    enddo
enddo

if (Me == 1) then
    call readmems(FNAME, HF, pH)
else
    OPEN (66,FILE=FNAME,STATUS='OLD')
    READ(66,*) (pH(i),i=1,NmicNode)
    CLOSE(66)
endif


pDilte = RESHAPE(pH, (/NmicY, NmicX/))


!! Zero Padding

write(*,*) 'Zero Padding ===>'

pDzp = 0.0

do j = 1,NmicX
    do i = 1,NmicY
        pDzp( NmicY * (ZPN-1) / 2 + i, NmicX * (ZPN-1) / 2 +j) = pDilte(i,j)
    enddo
enddo

!!  Forward FT

write(*,*) 'Forward DFT ===>'

phat = 0.0

do j = 1,NmicXz
    nx = j-(NmicXz/2)

do i = 1,NmicYz
    ny = i-(NmicYz/2)

    do l = 1,NmicXz
        mx = l-(NmicXz/2)
    do o = 1,NmicYz
        my = o-(NmicYz/2)

phat(i,j) = phat(i,j)+Xres*Yres*pDzp(o,l)*exp(2*pi*sqi* (nx*mx/cmplx(NmicXz) + ny*my/cmplx(NmicYz)))

    enddo
enddo
    enddo
enddo

    pvzhat = phat

! Setting up the integrand in the reconstruction equation

    DO j = 1, NmicXz
        DO i = 1, NmicYz
            phat(i,j) = phat(i,j) * EXP ( sqi * kz(i,j) * D )  !* Kwindow(i,j)
pvzhat(i,j) = ( kz(i,j) / ( mu * c * k ) ) * pvzhat(i,j) * EXP ( sqi * kz(i,j) * D ) !* Kwindow(i,j)
        ENDDO
    ENDDO

!! Calling the inverse 2D fft

write(*,*) 'Inverse DFT ===>'

pDzp = 0.0
pvz = 0.0

do j = 1,NmicXz
    mx = j-(NmicXz/2)

do i = 1,NmicYz
    my = i-(NmicYz/2)

    do l = 1,NmicXz
        nx = l-(NmicXz/2)
    do o = 1,NmicYz
        ny = o-(NmicYz/2)

        pDzp(i,j) = pDzp(i,j) + 1.0/(NmicXz*NmicYz*Xres*Yres)*phat(o,l)* &
                exp( -2*pi*sqi* (nx*mx/cmplx(NmicXz) + ny*my/cmplx(NmicYz)) )

        pvz(i,j) = pvz(i,j) + 1.0/(NmicXz*NmicYz*Xres*Yres)*pvzhat(o,l)* &
                exp( -2*pi*sqi* (nx*mx/cmplx(NmicXz) + ny*my/cmplx(NmicYz)) )
    enddo
enddo
    enddo
enddo

do j = 1,NmicX
    do i = 1,NmicY
        pDilte(i,j) = pDzp( NmicY * (ZPN-1) / 2 + i, NmicX * (ZPN-1) / 2 +j)
        pvz1(i,j) = pvz( NmicY * (ZPN-1) / 2 + i, NmicX * (ZPN-1) / 2 +j)
    enddo
enddo


END SUBROUTINE nah

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!                          SUPPORTING SUBROUTINES
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine green(g,PL,SL,k)
! implicit none

    COMPLEX*16, PARAMETER :: MINUS_ONE = -1.0
    COMPLEX*16,INTENT(OUT) :: g
    double precision,dimension(3),INTENT(IN) :: SL, PL
    double precision, PARAMETER :: pi = 3.141592653589793
    double precision :: k, r


!*********************************************************************************************************************
!
!           Freespace Green's function
!
!*********************************************************************************************************************


    r =  sqrt(( SL(1)-PL(1) )**2 + ( SL(2)-PL(2) )**2 + ( SL(3)-PL(3) )**2 )
    g = cmplx( exp( sqrt(MINUS_ONE) * k * r ) / ( 4.0 * pi * r ) )

    RETURN
END subroutine green


Subroutine dgreen(dg,SL,PL,k)
! implicit none

    COMPLEX*16, PARAMETER :: MINUS_ONE = -1.0
    COMPLEX*16,INTENT(OUT) :: dg
    double precision,dimension(3),INTENT(IN) :: SL, PL
    double precision, PARAMETER :: pi = 3.141592653589793
    double precision :: k, r


!*********************************************************************************************************************
!
!           Normal derivative of Freespace Green's function
!
!********************************************************************************************************************

    r =  sqrt(( SL(1)-PL(1) )**2 + ( SL(2)-PL(2) )**2 + ( SL(3)-PL(3) )**2 )
    dg = cmplx( exp( sqrt(MINUS_ONE) * k * r ) * ( sqrt(MINUS_ONE) * k * r - 1.0) * ( SL(3)-PL(3) ) / &
    ( 4.0 * pi * r**3 ) )


    RETURN
END subroutine dgreen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE shapeGreenL(N, Jac, G, DG, x,y,zs, xp,yp,zp, NL, k)

! Author: Kanthasamy Chelliah, Illinois Institute of Technology

! Linear shape functions for quadrilateral elements


double precision,DIMENSION(4), INTENT(OUT) :: N
double precision,DIMENSION(2,4) :: DN
double precision,DIMENSION(4,2), INTENT(IN) :: NL
double precision,DIMENSION(2,2) :: Jmat

COMPLEX*16, PARAMETER :: MINUS_ONE = -1.0
COMPLEX*16, INTENT(OUT) :: g, dg
double precision, PARAMETER :: pi = 3.141592653589793
double precision,INTENT(IN) :: x, y, xp, yp, zp, zs, k
double precision, INTENT(OUT) :: Jac
double precision :: xs, ys, r

!*********************************************************************************************************************
!
!           Shape functions
!
!********************************************************************************************************************

N(1) = 0.25 *  ( 1 - x ) * ( 1 - y )
N(2) = 0.25 *  ( 1 - x ) * ( 1 + y )
N(3) = 0.25 *  ( 1 + x ) * ( 1 + y )
N(4) = 0.25 *  ( 1 + x ) * ( 1 - y )


!*********************************************************************************************************************
!
!           Jacobian
!
!********************************************************************************************************************

DN(1,1) = -0.25 * ( 1 - y )
DN(2,1) = -0.25 * ( 1 - x )

DN(1,2) = -0.25 * ( 1 + y )
DN(2,2) = 0.25 * ( 1 - x )

DN(1,3) = 0.25 * ( 1 + y )
DN(2,3) = 0.25 * ( 1 + x )

DN(1,4) = 0.25 * ( 1 - y )
DN(2,4) = -0.25 * ( 1 + x )

Jmat = matmul(DN,NL)

Jac = Jmat(1,1) * Jmat(2,2) - Jmat(1,2) * Jmat(2,1)

!*********************************************************************************************************************
!
!           Freespace Green's function
!
!********************************************************************************************************************

xs = 0.0
ys = 0.0

do i = 1,4
    xs = xs + N(i) * NL(i,1)
    ys = ys + N(i) * NL(i,2)
enddo

r =  sqrt( ( xp - xs )**2 + ( yp - ys )**2 + ( zp - zs )**2 )

g = cmplx( exp( sqrt(MINUS_ONE) * k * r ) / ( 4.0 * pi * r ) )

!*********************************************************************************************************************
!
!           Normal derivative of Freespace Gree's function
!
!********************************************************************************************************************

dg = cmplx(exp(sqrt(MINUS_ONE)*k*r) * (sqrt(MINUS_ONE)*k*r-1.0) * (zp-zs)/(4.0*pi*r**3))


END SUBROUTINE shapeGreenL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gauleg(ngp, xabsc, weig)

INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!   PRIVATE
!   REAL (dbp) :: newv
REAL(dbp)  :: EPS, M_PI
PARAMETER (EPS=3.0d-15) !EPS is the relative precision
PARAMETER (M_PI=3.141592653589793) ! Pi value      INTEGER  i, j, m
REAL(dbp)  p1, p2, p3, pp, z, z1
INTEGER, INTENT(IN) :: ngp ! # of Gauss Points
REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)


   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

   do i = 1, m ! Loop over the desired roots */

        z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*   we enter the main loop of refinement by NEWTON'S method   */
100     p1 = 1.0d0
        p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        do j = 1, ngp
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z = z1 - p1/pp ! Newton's Method  */

        if (dabs(z-z1) .gt. EPS) GOTO  100

    xabsc(i) =  - z ! Roots will be bewteen -1.0 & 1.0 */
    xabsc(ngp+1-i) =  + z ! and symmetric about the origin  */
    weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
    weig(ngp+1-i) = weig(i) ! symmetric counterpart         */

  end do ! i loop

End subroutine gauleg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine l_quadrature_rule ( n, x, w )

!*****************************************************************************80
!
!! L_QUADRATURE_RULE: Gauss-Laguerre quadrature based on L(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!

  integer,intent(in) :: n

  real*8 bj(n)
  integer i
  real*8, intent(out) :: w(n)
  real*8, intent(out) :: x(n)
  real*8 zemu
!
!  Define the zero-th moment.
!
  zemu = 1.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    x(i) = real ( 2 * i - 1, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end subroutine l_quadrature_rule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!

  integer, intent(in) :: n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real*8, intent(inout) :: d(n)
  real*8, intent(inout) :: e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real*8, intent(inout) :: z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end subroutine imtqlx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine svd (SSS,U,VT,Sigma,N1,N2)
! implicit none

    INTEGER :: N1,N2             ! Number of Nodes (surface S)

    complex*16, PARAMETER :: alpha = 1.0
    complex*16, PARAMETER :: beta = 0.0

    character, parameter :: JOBU = 'A'
    character, parameter :: JOBVT = 'A'
    integer :: LWORK
    integer :: INFO
    complex*16,dimension(:),allocatable :: WORK
    double precision,dimension(5*N2) :: RWORK
    double precision,dimension(N1) :: Sigma
    complex*16, DIMENSION(N1, N2) :: SSS, SSS1
    complex*16, DIMENSION(N1, N1) :: U
    complex*16, DIMENSION(N2, N2) :: VT


SSS1 = SSS

LWORK = 10*max(N1,N2)
allocate(WORK(LWORK))

call ZGESVD(JOBU,JOBVT,N1,N2,SSS1,N1,Sigma,U,N1,VT,N2,WORK,LWORK,RWORK,INFO)

deallocate(WORK)

return
end subroutine svd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine svdsolve (U,VT,Sigma,N1,N2,pH,pS)
! implicit none

    INTEGER :: N1,N2, i             ! Number of Nodes (surface S)

    double precision,dimension(N1) :: Sigma
    complex*16, dimension(N1) :: pH, p1
    complex*16, dimension(N2) :: pS
    complex*16, DIMENSION(N1, N1) :: U
    complex*16, DIMENSION(N2, N2) :: VT

p1 = matmul(U,pH)

pS = 0.0
Do i = 1,N1
    pS(i) = p1(i) /  Sigma(i)
enddo

pS = matmul(VT,pS)

return
end subroutine svdsolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tikhonov (U,VT,Sigma,N1,N2,pH,pS,lambda)
! implicit none

    INTEGER :: N1,N2, i             ! Number of Nodes (surface S)

    double precision,dimension(N1) :: Sigma
    complex*16, dimension(N1) :: pH, p1
    complex*16, dimension(N2) :: pS
    complex*16, DIMENSION(N1, N1) :: U
    complex*16, DIMENSION(N2, N2) :: VT

    double precision :: lambda

p1 = matmul(U,pH)

pS = 0.0
Do i = 1,N1
    pS(i) = p1(i)*cmplx( Sigma(i)/(Sigma(i)**2+lambda**2),kind = 8)
enddo

pS = matmul(VT,pS)

return
end subroutine tikhonov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gcv(A,U,VT,Sigma,N1,N2,pH,lambda_min)
! implicit none

    INTEGER :: N1,N2, i, e, e1             ! Number of Nodes (surface S)
    INTEGER, parameter :: Npoints = 1500

    double precision,dimension(N1) :: Sigma, p2
    double precision,dimension(Npoints) :: reg_param, GF
    complex*16, dimension(N1) :: pH, p1
    complex*16, dimension(N2) :: pS
    complex*16, DIMENSION(N1, N1) :: U
    complex*16, DIMENSION(N2, N2) :: VT
    complex*16, DIMENSION(N1, N2) :: A

    double precision :: lambda_min, fi, smin_ratio, ratio


smin_ratio = epsilon(lambda_min);

!!  % Vector of regularization parameters.

  reg_param(Npoints) = max(Sigma(N1),Sigma(1)*smin_ratio);
  ratio = (Sigma(1)/reg_param(Npoints))**(1/(real(Npoints-1,kind=8)));

Do i = Npoints-1,1,-1
    reg_param(i) = ratio*reg_param(i+1);
enddo

!write(*,*) 'reg_param(Npoints), reg_param(1)',reg_param(Npoints), reg_param(1)

Do e = 1,Npoints

!    write(*,*) 'r1', size(U,1), size(U,2), size(VT,1), size(VT,2)

    call tikhonov(U,VT,Sigma,N1,N2,pH,pS,reg_param(e))

    p1 = matmul(A,pS)

    p2 = abs(p1-pH)
    fi = 0.0
    Do i = 1,N1
        fi = fi + Sigma(i)**2 / ( Sigma(i)**2 + reg_param(e)**2 )
    enddo

GF(e) = norm2(p2)**2 / ( real(N1,kind=8) - fi)**2

enddo

e1 = minloc(GF,1)

lambda_min = reg_param(e1);

!write(*,*) 'reg_param(Npoints), reg_param(1)',GF(Npoints), GF(1)


!call plot(log10(reg_param), log10(GF),xlabel = 'log(Gamma)', ylabel = 'log(Gcv)')



!    FNAME='GCVparam.txt'
!    OPEN (97,FILE=FNAME,STATUS='REPLACE')
!    write(97,*) (GF(i),i=1,Npoints)
!    CLOSE(97)
!
!    FNAME='RegParam.txt'
!    OPEN (98,FILE=FNAME,STATUS='REPLACE')
!    write(98,*) (reg_param(i),i=1,Npoints)
!    CLOSE(98)


return
end subroutine gcv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lcurve(A,U,VT,Sigma,N1,N2,pH,lambda_min)

    INTEGER :: N1,N2, i, e, e1             ! Number of Nodes (surface S)
    INTEGER, parameter :: Npoints = 500

    double precision,dimension(N1) :: Sigma
    double precision,dimension(Npoints) :: reg_param, x, y, kappa
    complex*16, dimension(N1) :: pH, p1
    complex*16, dimension(N2) :: pS
    complex*16, DIMENSION(N1, N1) :: U
    complex*16, DIMENSION(N2, N2) :: VT
    complex*16, DIMENSION(N1, N2) :: A

    double precision :: lambda_min, smin_ratio, ratio, dx,d2x,dy,d2y

smin_ratio = epsilon(lambda_min);

!!  % Vector of regularization parameters.

  reg_param(Npoints) = max(Sigma(N1),Sigma(1)*smin_ratio);
  ratio = (Sigma(1)/reg_param(Npoints))**(1/(real(Npoints-1,kind=8)));

Do i = Npoints-1,1,-1
    reg_param(i) = ratio*reg_param(i+1);
enddo

!write(*,*) 'reg_param(Npoints), reg_param(1)',reg_param(Npoints), reg_param(1)

Do e = 1,Npoints

    call tikhonov(U,VT,Sigma,N1,N2,pH,pS,reg_param(e))

    p1 = matmul(A,pS)

y(e) = norm2(abs(pS))
x(e) = norm2(abs(p1-pH))

enddo

Do e = 2,Npoints-1

dx = ( x(e+1) - x(e-1) ) / ( reg_param(e+1) - reg_param(e-1) )
dy = ( y(e+1) - y(e-1) ) / ( reg_param(e+1) - reg_param(e-1) )

d2x = ( x(e+1) + x(e-1) - 2.0 * x(e) ) / ( ( reg_param(e+1) - &
reg_param(e-1) )**2.0 / 4.0 )
d2y = ( y(e+1) + y(e-1) - 2.0 * y(e) ) / ( ( reg_param(e+1) - &
reg_param(e-1) )**2.0 / 4.0 )

kappa(e) = ( dx*d2y - d2x*dy ) / (dx*dx + dy*dy)**(1.5)

enddo

e1 = maxloc(kappa,1)

lambda_min = reg_param(e1);


!    FNAME='rho.txt'
!    OPEN (97,FILE=FNAME,STATUS='REPLACE')
!    write(97,*) (x(i),i=1,Npoints)
!    CLOSE(97)
!
!    FNAME='eta.txt'
!    OPEN (98,FILE=FNAME,STATUS='REPLACE')
!    write(98,*) (y(i),i=1,Npoints)
!    CLOSE(98)
!
!!    FNAME='d2y.txt'
!!    OPEN (88,FILE=FNAME,STATUS='REPLACE')
!!    write(88,*) (d2ynew(i),i=1,Npoints)
!!    CLOSE(88)
!
!    FNAME='kappa.txt'
!    OPEN (88,FILE=FNAME,STATUS='REPLACE')
!    write(88,*) (kappa(i),i=1,Npoints)
!    CLOSE(88)


return
end subroutine lcurve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE band(s,nd,f1,f2,delt,nroll,icaus)

!c  Butterworth bandpass filter order 2*nroll (nroll<=8) (see Kanasewich,
!c    Time Sequence Analysis in Geophysics, Third Edition,
!c    University of Alberta Press, 1981)
!c  written by W.B. Joyner 01/07/97
!c  causal if icaus.eq.1 - zero phase shift otherwise
!c  s(j) input = the time series to be filtered - output = the
!c    filtered series
!c  dimension of s(j) must be as least as large as the larger of
!c    the following:
!c    nd+3.0*float(nroll)/(f1*delt)
!c    nd+6.0*float(nroll)/((f2-f1)*delt)
!c  nd = the number of points in the time series
!c  f1, f2 = the cutoff frequencies
!c  delt = the timestep

!* Dates: xx/xx/xx - Written by Bill Joyner
!*        09/12/00 - Changed "n" to "nroll" to eliminate confusion with
!*                   Kanesewich, who uses "n" as the order (=2*nroll), and
!*                   cleaned up appearance of code by using spaces, indents, etc.
!*        09/12/00 - double precision statements added so that the output
!*                   series has sufficient precision for double integration.
!*        11/08/00 - Increased dimension of s from 50000 to 100000
!*        02/14/01 - Modified/corrected real numbers function calls to
!*                   double precision - cds
!*        02/22/01 - Removed dimension of s (it is up to the user to specify
!*                   it properly)

      INTEGER, INTENT(in) :: f1, f2, icaus, nroll, nd
      REAL*8, INTENT(INOUT), dimension(nd) :: s
      REAL*8 fact(16),b1(16),b2(16)
      REAL*8 pi,w1,xp,yp,x1,x2,y1,y2,w2
      REAL*8 pre, pim, argre, argim, rho, theta, sjre, sjim
      REAL*8 bj, cj, con
      INTEGER i, j, k, np1, np2
      DOUBLE PRECISION, INTENT(IN):: delt
      double precision npad


      IF(f1.eq.0..or.f1.eq.f2) RETURN

      pi=4.0d0*DATAN(1.0d0)
      w1=2.0d0*pi*f1
      w1=2.0d0*DTAN(w1*delt/2.0d0)/delt
      w2=2.0d0*pi*f2
      w2=2.0d0*DTAN(w2*delt/2.0d0)/delt

      DO k=1,nroll
        pre=-DSIN(pi*DFLOAT(2*k-1)/DFLOAT(4*nroll))
        pim=DCOS(pi*DFLOAT(2*k-1)/DFLOAT(4*nroll))
        argre=(pre**2-pim**2)*(w2-w1)**2/4.0d0-w1*w2
        argim=2.0d0*pre*pim*(w2-w1)**2/4.0d0
        rho=(argre**2+argim**2)**(1.0d0/4.0d0)
        theta=pi+DATAN2(argim,argre)/2.0d0
        DO i=1,2
          sjre=pre*(w2-w1)/2.0d0+(-1)**i*rho*(-DSIN(theta-pi/2.0d0))
          sjim=pim*(w2-w1)/2.0d0+(-1)**i*rho*DCOS(theta-pi/2.0d0)
          bj=-2.0d0*sjre
          cj=sjre**2+sjim**2
          con=1.0d0/(2.0d0/delt+bj+cj*delt/2.0d0)
          fact(2*k+i-2)=(w2-w1)*con
          b1(2*k+i-2)=(cj*delt-4.0d0/delt)*con
          b2(2*k+i-2)=(2.0d0/delt-bj+cj*delt/2.0d0)*con
        END DO
      END DO

      np2=nd

      IF(icaus.ne.1) THEN
        npad=3.0*FLOAT(nroll)/(f1*delt)
        IF( npad .lt. 6.0*FLOAT(nroll)/((f2-f1)*delt) ) THEN
          npad=6.0*FLOAT(nroll)/((f2-f1)*delt)
        END IF
        np1=nd+1
        np2=nd+ NINT(npad)
        DO j=np1,np2
          s(j)=0.0
        END DO
      END IF

      DO k=1,2*nroll
        x1=0.0d0
        x2=0.0d0
        y1=0.0d0
        y2=0.0d0
        DO j=1,np2
          xp=s(j)
          yp=fact(k)*(xp-x2)-b1(k)*y1-b2(k)*y2
          s(j)=yp
          y2=y1
          y1=yp
          x2=x1
          x1=xp
        END DO
      END DO

      IF(icaus.ne.1) THEN
        DO k=1,2*nroll
          x1=0.0d0
          x2=0.0d0
          y1=0.0d0
          y2=0.0d0
          DO j=1,np2
            xp=s(np2-j+1)
            yp=fact(k)*(xp-x2)-b1(k)*y1-b2(k)*y2
            s(np2-j+1)=yp
            y2=y1
            y1=yp
            x2=x1
            x1=xp
          END DO
        END DO
      END IF

      RETURN

END SUBROUTINE band

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE readmems64(FNAME, HF, pDilte1)
! implicit none

    INTEGER, PARAMETER :: Frame = 1024;
    INTEGER, PARAMETER :: Nframe = 96;
    INTEGER, PARAMETER :: Nmea = 64;
    INTEGER, PARAMETER :: LengthTS = Frame*Nframe*Nmea;
    INTEGER, PARAMETER :: LengthTSm = Frame*Nframe;

    INTEGER, PARAMETER :: Transi = 25600;
    INTEGER, PARAMETER :: Fs = 50000;

    REAL*8, PARAMETER :: MicCalib = 1.0;

    REAL*8, PARAMETER :: Pref = 0.00002;
    REAL*8, PARAMETER :: dT = 1.0/real(Fs,kind = 8)

    REAL*8, PARAMETER :: c = 343.0
    REAL*8, PARAMETER :: pi = 3.141592653589793
    COMPLEX*8, PARAMETER :: MINUS_ONE = -1.0

    INTEGER i,j,N,l, fl, fh
    INTEGER , INTENT(IN) :: HF
    double precision MeanTS,w,just
    complex*16 DFT


    double precision, dimension(:), allocatable :: TSm
    integer, dimension(:), allocatable :: TS

    complex*16, DIMENSION(Nmea) :: pDilte
    complex*16, DIMENSION(Nmea), INTENT(OUT) :: pDilte1


    CHARACTER(LEN=100), INTENT(IN) ::   FNAME

    INTEGER, PARAMETER :: icas = 0;
    INTEGER, PARAMETER :: nroll = 2;


allocate(TS(LengthTS))
allocate(TSm(Nframe*Frame))

    fl = HF - 10
    fh = HF + 10

        OPEN (55,FILE=FNAME,STATUS='OLD')
        READ(55,*)  (TS(N),N=1,LengthTS)
        CLOSE(55)

w = 2.0*pi*HF


DO i=1,Nmea

Do j=0, (Nframe-1)
    Do N=1,Frame

        just = real(TS( ( (j*Nmea + i-1) * Frame ) + N ) , kind=8) * 120.0 /( 2.0 ** (24.0) )
        TSm(j*Frame + N) = Pref * 10.0**(just / 20.0);

    end do
END DO


call band(TSm,LengthTSm,fl,fh,dT,nroll,icas)


MeanTS = sum(TSm)/real(LengthTSm, kind = 8 )
DFT = 0.0

    do l = Transi,LengthTSm
            TSm(l) = TSm(l) - MeanTS
    DFT = DFT + cmplx(TSm(l),kind=8) * exp(sqrt(MINUS_ONE) * cmplx(w * (l-Transi) * dT , kind = 8))
    enddo
    pDilte(i) = DFT*cmplx((MicCalib),kind=8)/cmplx(LengthTSm-Transi,kind=8)

End do

deallocate(TS)
deallocate(TSm)


do i=1,Nmea-1,2
pDilte1(i) = pDilte((i+1)/2);
enddo
j=1;
do i=2,Nmea,2
pDilte1(i) = pDilte(32+j);
j=j+1;
enddo

!write(*,*) pDilte1

        OPEN (59,FILE="pH.txt",STATUS='REPLACE')
        WRITE(59,*)  (abs(pDilte1(j)),j=1,Nmea)
        CLOSE(59)

END SUBROUTINE readmems64

SUBROUTINE readmems(FNAME, HF, pDilte1)
! implicit none

!    INTEGER, PARAMETER :: Frame = 1024;
!    INTEGER, PARAMETER :: Nframe = 96;
    INTEGER, PARAMETER :: Nmea = 64;
!    INTEGER, PARAMETER :: LengthTS = Frame*Nframe*Nmea;
    INTEGER, PARAMETER :: LengthTSm = 98304  ! 491520;

    INTEGER, PARAMETER :: Transi = 10000;
    INTEGER, PARAMETER :: Fs = 50000;

    REAL*8, PARAMETER :: MicCalib = 1.0;

    REAL*8, PARAMETER :: Pref = 0.00002;
    REAL*8, PARAMETER :: dT = 1.0/real(Fs,kind = 8)

    REAL*8, PARAMETER :: c = 343.0
    REAL*8, PARAMETER :: pi = 3.141592653589793
    COMPLEX*8, PARAMETER :: MINUS_ONE = -1.0

    INTEGER i,j,l
    INTEGER , INTENT(IN) :: HF
    double precision MeanTS,w,dB
    complex*16 DFT


    double precision, dimension(:), allocatable :: TSm
    integer, dimension(:,:), allocatable :: TS

!    complex*16, DIMENSION(Nmea) :: pDilte
    complex*16, DIMENSION(Nmea), INTENT(OUT) :: pDilte1
    complex*16, DIMENSION(Nmea) :: pDilte


    CHARACTER(LEN=100), INTENT(IN) ::   FNAME

    INTEGER, PARAMETER :: icas = 0;
    INTEGER, PARAMETER :: nroll = 2;


allocate(TS(LengthTSm,Nmea))
allocate(TSm(LengthTSm))

        OPEN (55,FILE=FNAME,STATUS='OLD')
        READ(55,*)  ((TS(i,j),j=1,Nmea),i=1,LengthTSm)
        CLOSE(55)

w = 2.0*pi*HF

write(*,*) TS(1,1)
DO i=1,Nmea


    TSm = real(TS(:,i),kind=8);
    MeanTS = sum(TSm)/real(LengthTSm, kind = 8 )
write(*,*) MeanTS
DFT = 0.0

    do l = 1,LengthTSm
        TSm(l) = TSm(l) - MeanTS
        DFT = DFT + TSm(l) * exp( sqrt(MINUS_ONE) * w * l * dT )
    enddo

    DFT = DFT/real(LengthTSm, kind = 8 )
    dB = 20. * log10( abs(DFT) * 2. / (2.**31.-1.) ) + 120.
    pDilte(i)  = Pref * 10.**(dB/20.) * exp(sqrt(MINUS_ONE) * atan2(dimag(DFT),real(DFT,kind = 8)) );
write(*,*) dB
End do

DO i = 0,31
    pDilte1(2*i+1) = pDilte(i+1);
    pDilte1(2*i+1+1) = pDilte(32+i+1);
ENDDO


pDilte1(17) = ( pDilte1(9) + pDilte1(18) + pDilte1(25) )/3.0

deallocate(TS)
deallocate(TSm)

END SUBROUTINE readmems
