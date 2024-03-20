!C*******************************************************************
!C           LCR THROUGH SIMPLEC IN COLLOCATED GRID                 *
!C   Program to solve numerically, with the Finite volume Method,   *
!C   the two-dimensional lid-driven caity flow with a polynomial-   *
!C   based tangential velocity profile                              *
!C   Viscoelastic fluid: Giesekus mode with LOGARITHM of conforma-  *
!C      tion formulation (mesh B non-staggered)                     *
!c*******************************************************************
C      The user may modify the code and give it to third parties,
C      provided that an acknolwedgement to the source of the original
C      version is retained.
C
C		Disclaimer
C		The code VELUT-LOG is distributed in the hope that it will be useful,
C       but WITHOUT ANY WARRANTY.      
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MODULE BANK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION TITLE
      PARAMETER(NX=210,NY=210,NFMAX=10)
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP
      CHARACTER*16 FNAME
      COMMON F(NX,NY,NFMAX),P(NX,NY),RHO(NX,NY),GAM(NX,NY),CON(NX,NY),
     & AIP(NX,NY),AIM(NX,NY),AJP(NX,NY),AJM(NX,NY),AP(NX,NY),
     & PT(NX),QT(NX) ,GAMP(NX,NY),GAMT(NX,NY),
     & X(NX),XU(NX),XDIF(NX),XCV(NX),XCVS(NX),
     & Y(NY),YV(NY),YDIF(NY),YCV(NY),YCVS(NY),
     & YCVR(NY),YCVRS(NY),ARX(NY),ARXJ(NY),ARXJP(NY),
     & R(NY),RMN(NY),SX(NY),SXMN(NY),XCVI(NY),XCVIP(NY),
     & TPXXE(NX,NY),TPYYN(NX,NY),TPXYN(NX,NY),TPXYE(NX,NY) 
      COMMON DU(NX,NY),DV(NX,NY),FV(NX),FVP(NX),
     & FX(NX),FXM(NX),FY(NY),FYM(NY)
      COMMON /INDX/NF,NP,NRHO,NGAM,L1,L2,L3,M1,M2,M3,
     &  IST,JST,ITER,LAST,TITLE(13),RELAX(13),XL,YL,IPREF,JPREF,
     &  LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),RHOCON,RESISUM,
     & RP,RP0,SMALL1,SMALL,TOL,ALP
      COMMON/CNTL/LSTOP
      COMMON/COD/DGRID,DT,NCODE,NSELECT,NSCHEME,UIN,RESU,RESV,RESTXYMAX
      COMMON/SORC/SMAX,SSUM,NBAK,RESUMAX,RESVMAX,RESTXXMAX,RESTYYMAX
      COMMON/COEF/FLOW,DIFF,ACOF,RESTXX,RESTXY,RESTYY,RES

      DIMENSION TXX(NX,NY),TYY(NX,NY),TXY(NX,NY),AX(NX,NY),AY(NX,NY),
     & P1(NX,NY),P2(NX,NY),CMXX(NX,NY),CMYY(NX,NY),CMXY(NX,NY),
     & CMYX(NX,NY),TB(NX,NY),TBAC(NX,NY),
     & UFX(NX,NY),UFY(NX,NY),VFX(NX,NY),VFY(NX,NY),OMGXY(NX,NY),
     & BXX(NX,NY),BYY(NX,NY),BXY(NX,NY),
     & CXX(NX,NY),CXY(NX,NY),CYY(NX,NY),
     & TPXX(NX,NY),TPXY(NX,NY),TPYY(NX,NY),
     & AX2(NX,NY),AY2(NX,NY),AXY(NX,NY),AXY2(NX,NY),
     & ALB1(NX,NY),ALB2(NX,NY),ALB(NX,NY),
     & WXYXX(NX,NY),WXYXY(NX,NY) ,WXYYY(NX,NY),
     & COLDBXX(NX,NY),COLDBYY(NX,NY),COLDBXY(NX,NY)
      !----------------------------------------------------------------------

      DIMENSION U(NX,NY),V(NX,NY),PC(NX,NY),T(NX,NY),APPU(NX,NY),
     & APPV(NX,NY),pp(NX,NY),FF(NX,NY,NFMAX),UE(NX,NY),VN(NX,NY),
     & UEL(NX,NY),VNL(NX,NY),CONU(NX,NY),CONV(NX,NY),STREAM(NX,NY)
	DIMENSION RESI(NX,NY,10),FLOWX(NX,NY),FLOWY(NX,NY)
	DIMENSION VNN(NX,NY),UEE(NX,NY),ANU1(NY),ANU2(NY)
      EQUIVALENCE (F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
     & ,(F(1,1,4),TXX(1,1)),(F(1,1,5),TYY(1,1)),(F(1,1,6),TXY(1,1))
      END MODULE
C***************************************************************
*     -----------------MAIN PROGRAM------------------------------------
      PROGRAM MAIN
      USE BANK        !Universal variable module
      CALL SETUP0		!Set initial calculation parameters
      CALL GRID		!Divide the grid information 
      CALL SETUP1		!Sets the geometric quantities associated with
                      !the grid system that remain constant during the calculation
      CALL START      !Initial the calculation variables
   10 CALL DENSE		!Density
      CALL BOUND		!Boundary condition
      CALL OUTPUT		!Output the reselut
      IF(.NOT.LSTOP) GO TO 15
      STOP			!STOP
   15 CALL SETUP2		
      CALL SETUP3
      GO TO 10		
      END
C******************************************************************

C******************************************************************
	
      SUBROUTINE USER    
	
*******************************************************************
*        TEST FOR THE SIMPLEC
*******************************************************************  
	USE BANK
      ENTRY GRID          
	NSELECT=1			
      NSCHEME=1			!Deferred correction
	MODE=1				!for mode1 
      RELAX(1)=0.5		!Relaxation factor for u,v,p,TXX,TYY,TXY,p'(11)
	RELAX(2)=0.5
	RELAX(3)=0.5
      RELAX(4)=0.7
      RELAX(5)=0.7
      RELAX(6)=0.7
	RELAX(11)=1.0
      LSOLVE(1)=.TRUE.	!Solve u,v,p
      LSOLVE(2)=.TRUE.
      LSOLVE(3)=.TRUE.
      LSOLVE(4)=.TRUE.	!Solve Txx,Tyy,Txy
      LSOLVE(5)=.TRUE.
      LSOLVE(6)=.TRUE.
	LPRINT(11)=.TRUE.	
      LAST=1000000		!
      XL=1.0   			!the geometric information
      YL=1.0
	D=0.5*XL
      L1=102              !the number of nodes
      M1=102
	XL0=XL/2.0
	YL0=YL/2.0
      CALL UGRID			
      RETURN
C
      ENTRY START		  ! Initialise variables	
      DO 101 J=1,M1  
      DO 101 I=1,L1
             U(I,J)=0.0			
             V(I,J)=0.0			  
             P(I,J)=0.
             TPXX(I,J)=0.0
             TPXY(I,J)=0.0
             TPYY(I,J)=0.0
             TXX(I,J)=0.0
             TXY(I,J)=0.0
             TYY(I,J)=0.0
             CXX(I,J)=0.0
             CXY(I,J)=0.0
             CYY(I,J)=0.0
             AX(I,J)=1.0
             AY(I,J)=0.0
             P1(I,J)=0.0
             P2(I,J)=0.0
  101 CONTINUE
      RETURN
C
      ENTRY DENSE			!Density field
	DO 102 J=1,M1
	DO 102 I=1,L1		
	RHO(I,J)=1.0
  102 CONTINUE
      RETURN
C
      ENTRY BOUND			!Boundary condition
	DO I=1,L1
	U(I,M1)=8*UIN*(1+tanh(8*(TIME-0.5)))*X(I)**2*(1-X(I))**2
	U(I,1)=0.0
	V(I,1)=0.0
	V(I,M1)=0.0 
      END DO

	DO J=1,M1			
	U(1,J)=0.0
	U(L1,J)=0.0
	V(1,J)=0.0
	V(L1,J)=0.0
      END DO


      RETURN
C
      ENTRY OUTPUT		
   
      CALL PRINT
c-------------------------------------------------------
      RETURN
C
      ENTRY GAMSOR
	DO 501 J=1,M1
      DO 501 I=1,L1
      GAM(I,J)=0.5 !1          
      GAMP(I,J)=0.5 
      GAMT(I,J)=1.0
      IF(NF.EQ.4.OR.NF.EQ.5.OR.NF.EQ.6)THEN
          GAM(I,J)=0.0
      ENDIF
  501 CONTINUE

      
	RETURN
	END SUBROUTINE USER
*     *****************************************************************


*     -----------------------------------------------------------------
      SUBROUTINE DIFLOW	
      USE BANK
      ACOF=DIFF
      RETURN
      END SUBROUTINE DIFLOW
*--------------------------------------------------------------------------
      SUBROUTINE SOLVE  
      USE BANK		
c--------------------------------------------------------------------------
      ISTF=IST-1		  
      JSTF=JST-1
      IT1=L2+IST
      IT2=L3+IST
      JT1=M2+JST
      JT2=M3+JST
****************************************************************************
      DO 999 NT=1,NTIMES(NF)	
      DO 999 N=NF,NF
*---------------------------------------------------------------------------
      IF(.NOT. LBLK(NF)) GO TO 10
      PT(ISTF)=0.
      QT(ISTF)=0.
      DO 11 I=IST,L2
      BL=0.
      BLP=0.
      BLM=0.
      BLC=0.
      DO 12 J=JST,M2
      BL=BL+AP(I,J)
      IF(J .NE. M2) BL=BL-AJP(I,J)
      IF(J .NE. JST) BL=BL-AJM(I,J)
      BLP=BLP+AIP(I,J)
      BLM=BLM+AIM(I,J)
      BLC=BLC+CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
     &   +AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)-AP(I,J)*F(I,J,N)
   12 CONTINUE
      DENOM=BL-PT(I-1)*BLM
      DENO=1.E15
      IF(ABS(DENOM/BL) .LT. 1.E-10) DENOM=1.E20*DENO
      PT(I)=BLP/DENOM
      QT(I)=(BLC+BLM*QT(I-1))/DENOM
   11 CONTINUE
      BL=0.
      DO 13 II=IST,L2
      I=IT1-II
      BL=BL*PT(I)+QT(I)
      DO 13 J=JST,M2
   13 F(I,J,N)=F(I,J,N)+BL
*---------------------------------------------------------------------------
      PT(JSTF)=0.
      QT(JSTF)=0.
      DO 21 J=JST,M2
      BL=0.
      BLP=0.
      BLM=0.
      BLC=0.
      DO 22 I=IST,L2
      BL=BL+AP(I,J)
      IF(I .NE. L2) BL=BL-AIP(I,J)
      IF(I .NE. IST) BL=BL-AIM(I,J)
      BLP=BLP+AJP(I,J)
      BLM=BLM+AJM(I,J)
      BLC=BLC+CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
     &   +AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)-AP(I,J)*F(I,J,N)
   22 CONTINUE
      DENOM=BL-PT(J-1)*BLM
      IF (ABS(DENOM/BL) .LT. 1E-10) DENOM=1.E20*DENO
      PT(J)=BLP/DENOM
      QT(J)=(BLC+BLM*QT(J-1))/DENOM
   21 CONTINUE
      BL=0.
      DO 23 JJ=JST,M2
      J=JT1-JJ
      BL=BL*PT(J)+QT(J)
      DO 23 I=IST,L2
   23 F(I,J,N)=F(I,J,N)+BL
   10 CONTINUE
*-----------------------------------------------------------------------
      DO 90 J=JST,M2
      PT(ISTF)=0.
      QT(ISTF)=F(ISTF,J,N)
      DO 70 I=IST,L2
      DENOM=AP(I,J)-PT(I-1)*AIM(I,J)
      PT(I)=AIP(I,J)/DENOM
      TEMP=CON(I,J)+AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)
      QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
   70 CONTINUE
      DO 80 II=IST,L2
      I=IT1-II
   80 F(I,J,N)=F(I+1,J,N)*PT(I)+QT(I)
   90 CONTINUE
*-----------------------------------------------------------------------
      DO 190 JJ=JST,M3
      J=JT2-JJ
      PT(ISTF)=0.
      QT(ISTF)=F(ISTF,J,N)
      DO 170 I=IST,L2
      DENOM=AP(I,J)-PT(I-1)*AIM(I,J)
      PT(I)=AIP(I,J)/DENOM
      TEMP=CON(I,J)+AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)
      QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
  170 CONTINUE
      DO 180 II=IST,L2
      I=IT1-II
  180 F(I,J,N)=F(I+1,J,N)*PT(I)+QT(I)
  190 CONTINUE
*-----------------------------------------------------------------------
      DO 290 I=IST,L2
      PT(JSTF)=0.
      QT(JSTF)=F(I,JSTF,N)
      DO 270 J=JST,M2
      DENOM=AP(I,J)-PT(J-1)*AJM(I,J)
      PT(J)=AJP(I,J)/DENOM
      TEMP=CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
      QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
  270 CONTINUE
      DO 280 JJ=JST,M2
      J=JT1-JJ
  280 F(I,J,N)=F(I,J+1,N)*PT(J)+QT(J)
  290 CONTINUE
*-----------------------------------------------------------------------
      DO 390 II=IST,L3
      I=IT2-II
      PT(JSTF)=0.
      QT(JSTF)=F(I,JSTF,N)
      DO 370 J=JST,M2
      DENOM=AP(I,J)-PT(J-1)*AJM(I,J)
      PT(J)=AJP(I,J)/DENOM
      TEMP=CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
      QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
  370 CONTINUE
      DO 380 JJ=JST,M2
      J=JT1-JJ
  380 F(I,J,N)=F(I,J+1,N)*PT(J)+QT(J)
  390 CONTINUE
************************************************************************
  999 CONTINUE
c-----------------------------------------------------------------------
      IF(N.EQ.3) THEN

	RESISUM0=RESISUM
	RP0=RP
	RESISUM=0.0
      DO J=2,M2
	DO I=2,L2
      RESI(I,J,N)=AIM(I,J)*F(I-1,J,N)*AIP(I,J)*F(I+1,J,N)
     &+AJM(I,J)*F(I,J-1,N)+AJP(I,J)*F(I,J+1,N)+CON(I,J)
     &-AP(I,J)*F(I,J,N)
	RESI(I,J,N)=RESI(I,J,N)*RESI(I,J,N)
	END DO
	END DO

	DO J=2,M2
	DO I=2,L2
	RESISUM=RESISUM+RESI(I,J,N)*XCV(I)*YCV(J)
      END DO
	END DO
	RESISUM=SQRT(RESISUM)
      RP=RESISUM/RESISUM0      
      END IF
c-----------------------------------------------------------------------
      DO 400 J=2,M2
      DO 400 I=2,L2
      CON(I,J)=0.
      AP(I,J)=0.
  400 CONTINUE
      RETURN
      END SUBROUTINE SOLVE
***********************************************************************
      SUBROUTINE SETUP
      USE BANK
************************************************************************
    1 FORMAT(//15X,'COMPUTATION IN CARTISIAN COORDINATES')
    2 FORMAT(//15X,'COMPUTATION FOR AXISYMMETRICAL SITUATION')
    3 FORMAT(//15X,' COMPUTATION IN POLAR COORDINATES  ')
    4 FORMAT(1X,14X,40(1H*),//)
*-----------------------------------------------------------------------
      ENTRY SETUP0
      !NFMAX=10
      NP=11				
      NRHO=12				
      NGAM=13				
      LSTOP=.FALSE.		
      DT=2.0E-4			
      SMALL1=1.D-10       
      SMALL=1.D-30
      ALP=0.5
      UIN=1.0
      TOL=1.D-6
      DO I=1,10		
      LSOLVE(I)=.FALSE.   !when true,we slove for F(I,J,NF)
      LBLK(I)=.TRUE.		!when true,the black correction for F(I,J,NF) is used.
      NTIMES(I)=3	
	END DO				
      DO I=1,13		    
      LPRINT(I)=.FALSE.	!when true,F(I,J,NF) is printed.
      RELAX(I)=1.			
      END DO
	MODE=1				
      LAST=5				
      TIME=0.				
      ITER=0				
      IPREF=1				
      JPREF=1
      RHOCON=1			
      RETURN
*     -----------------------------------------------------------------
      ENTRY SETUP1
      L2=L1-1				
      L3=L2-1
      M2=M1-1
      M3=M2-1

      X(1)=XU(2)
      DO I=2,L2			
      X(I)=0.5*(XU(I+1)+XU(I))
      END DO
	X(L1)=XU(L1)

      Y(1)=YV(2)
      DO J=2,M2
      Y(J)=0.5*(YV(J+1)+YV(J))
      END DO
      Y(M1)=YV(M1)
      !δx
      DO I=2,L1
      XDIF(I)=X(I)-X(I-1)
      END DO
      !Δx
      DO I=2,L2
      XCV(I)=XU(I+1)-XU(I)
	END DO
	!(δx)w+ ,(δx)e-
      DO I=2,L2
      XCVI(I)=0.5*XCV(I)
      XCVIP(I)=XCVI(I)
	END DO
      DO 35 J=2,M1
   35 YDIF(J)=Y(J)-Y(J-1)
      DO 40 J=2,M2
   40 YCV(J)=YV(J+1)-YV(J)
      IF (MODE .NE. 1) GO TO 55
      DO 52 J=1,M1
      RMN(J)=1.
   52 R(J)=1.
      GO TO 56
   55 DO 50 J=2,M1
   50 R(J)=R(J-1)+YDIF(J)
      RMN(2)=R(1)
      DO 60 J=3,M2
   60 RMN(J)=RMN(J-1)+YCV(J-1)
      RMN(M1)=R(M1)
   56 CONTINUE
      DO 57 J=1,M1
      SX(J)=1.
      SXMN(J)=1.
      IF(MODE .NE. 3) GO TO 57
      SX(J)=R(J)
      IF(J .NE. 1) SXMN(J)=RMN(J)
   57 CONTINUE
      DO 62 J=2,M2
      YCVR(J)=R(J)*YCV(J)
      ARX(J)=YCVR(J)
      IF (MODE .NE. 3) GO TO 62
      ARX(J)=YCV(J)
   62 CONTINUE
      IF(MODE .NE. 2) GO TO 67
      DO 65 J=3,M3
      ARXJ(J)=0.25*(1.+RMN(J)/R(J))*ARX(J)
   65 ARXJP(J)=ARX(J)-ARXJ(J)
      GO TO 68
   67 DO 66 J=3,M3
      ARXJ(J)=0.5*ARX(J)
   66 ARXJP(J)=ARXJ(J)
   68 ARXJP(2)=ARX(2)
      ARXJ(M2)=ARX(M2)

      DO 70 J=3,M3
      FV(J)=ARXJP(J)/ARX(J)
   70 FVP(J)=1.-FV(J)

      DO 85 I=3,L2
      FX(I)=0.5*XCV(I-1)/XDIF(I)
   85 FXM(I)=1.-FX(I)
      FX(2)=0.
      FXM(2)=1.
      FX(L1)=1.
      FXM(L1)=0.
      DO 90 J=3,M2
      FY(J)=0.5*YCV(J-1)/YDIF(J)
   90 FYM(J)=1.-FY(J)
      FY(2)=0.
      FYM(2)=1.
      FY(M1)=1.
      FYM(M1)=0.
*     --------CON,AP,U,V,RHO,PC AND P ARRAYS ARE INITIALIZED HERE------
      DO 95 J=1,M1		
      DO 95 I=1,L1
      PC(I,J)=0.
      U(I,J)=0.
      V(I,J)=0.
      CON(I,J)=0.		!Sc
      AP(I,J)=0.		!Sp
      RHO(I,J)=RHOCON

	APPU(I,J)=1.0E+30 !To avoid overflow
	APPV(I,J)=1.0E+30 !To avoid overflow
   95 CONTINUE
      RETURN
*     -----------------------------------------------------------------
      ENTRY SETUP2
*     ---------------COEFFICIENTS FOR ALL EQUATIONS------------------
      RESTXXMAX=0.0
      RESTXYMAX=0.0
      RESTYYMAX=0.0
      RESUMAX=0.0
      RESVMAX=0.0
      RES=0.0
      DO NNN=1,NFMAX
      DO I=1,L1
      DO J=1,M1
      FF(I,J,NNN)=F(I,J,NNN)   
      END DO
      END DO
      END DO

      DO I=1,L1                !the face velocity
      DO J=1,M1
	UEL(I,J)=UE(I,J)
	VNL(I,J)=VN(I,J)
      END DO
      END DO
      DO 898 KIJ=1,5     
      IST=2              
      JST=2              
      DO 600 N=1,NFMAX
      NF=N
      IF(.NOT. LSOLVE(NF)) GO TO 600 
      IF(NF.EQ.3) GO TO 300  
      IF(NF.GT.3) GO TO 898  
      CALL GAMSOR
      REL=1.-RELAX(NF)
      !y=0 As
      DO 602 I=2,L2
      AREA=R(1)*XCV(I)
      FLOW=AREA*V(I,1)*RHO(I,1)
C----------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWY(I,J)=FLOW
	END IF
C------------------------------------------
      DIFF=AREA*GAM(I,1)/YDIF(2)
      CALL DIFLOW
  602 AJM(I,2)=ACOF+DMAX1(0.D0,FLOW)
      !x=0 Aw
      DO 603 J=2,M2
      FLOW=ARX(J)*U(1,J)*RHO(1,J)
C----------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWX(I,J)=FLOW
	END IF
C------------------------------------------
      DIFF=ARX(J)*GAM(1,J)/(XDIF(2)*SX(J))
      CALL DIFLOW
      AIM(2,J)=ACOF+DMAX1(0.D0,FLOW)
      
      DO 603 I=2,L2

      IF(I .EQ. L2) GO TO 604
      FLOW=ARX(J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))*UE(I,J)
c---------------------------------------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWX(I,J)=FLOW
	END IF
C************************************************************************
      DIFF=ARX(J)*2.*GAM(I,J)*GAM(I+1,J)/((XCV(I)*GAM(I+1,J)+
     &   XCV(I+1)*GAM(I,J)+1.0E-30)*SX(J))
      GO TO 605
  604 FLOW=ARX(J)*U(L1,J)*RHO(L1,J)
C---------------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWX(I,J)=FLOW
	END IF
C--------------------------------------
      DIFF=ARX(J)*GAM(L1,J)/(XDIF(L1)*SX(J))
  605 CALL DIFLOW

      AIM(I+1,J)=ACOF+DMAX1(0.D0,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW

      AREA=RMN(J+1)*XCV(I)
      IF(J .EQ. M2) GO TO 606
      FLOW=RMN(J+1)*XCV(I)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
     &	*VN(I,J)
C-----------------------------------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWY(I,J)=FLOW
	END IF
C************************************************************************
      DIFF=AREA*2.*GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+
     &   YCV(J+1)*GAM(I,J)+1.0E-30)
      GO TO 607
  606 FLOW=AREA*V(I,M1)*RHO(I,M1)
c----------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWY(I,J)=FLOW
	END IF
c---------------------------------
      DIFF=AREA*GAM(I,M1)/YDIF(M1)
  607 CALL DIFLOW

      AJM(I,J+1)=ACOF+DMAX1(0.D0,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW

      VOL=YCVR(J)*XCV(I)
      APT=RHO(I,J)/DT

	IF(NF.EQ.1) THEN
	CONU(I,J)=CON(I,J)
	END IF
	IF(NF.EQ.2) THEN
	CONV(I,J)=CON(I,J)
	END IF

      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*FF(I,J,NF)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     &   /RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*F(I,J,NF)


********************COEFFICIENTS FOR VELOCITY U EQUATION***************
      IF(NF.EQ.1)THEN
      DU(I,J)=VOL/(SX(J)*AP(I,J))
      PPX=(P(I+1,J)-P(I-1,J))/((XDIF(I)+XDIF(I+1))*SX(J))
      TTUXX=(TPXX(I+1,J)-TPXX(I-1,J))/((XDIF(I)+XDIF(I+1))*SX(J))
      TTUXY=(TPXY(I,J+1)-TPXY(I,J-1))/(YDIF(J)+YDIF(J+1))
          IF(I.EQ.2) THEN
          PPX=(P(I+1,J)-P(I,J))/(XDIF(I+1)*SX(J))
          TTUXX=(TPXX(I+1,J)-TPXX(I,J))/(XDIF(I+1)*SX(J))
          ENDIF
          IF(I.EQ.L2) THEN
          PPX=(P(I,J)-P(I-1,J))/(XDIF(I)*SX(J))
          TTUXX=(TPXX(I,J)-TPXX(I-1,J))/(XDIF(I)*SX(J))
          END IF
          IF(J.EQ.2) THEN 
          TTUXY=(TPXY(I,J+1)-TPXY(I,J))/YDIF(J+1)
          END IF
          IF(J.EQ.M2) THEN
          TTUXY=(TPXY(I,J)-TPXY(I,J-1))/YDIF(J)  
          END IF
            
      CON(I,J)=CON(I,J)-VOL*PPX
      CON(I,J)=CON(I,J)+VOL*TTUXX+VOL*TTUXY     !The source term for u 
      CONU(I,J)=CONU(I,J)+TTUXX+TTUXY           !
      APPU(I,J)=AP(I,J)

      END IF
********************COEFFICIENTS FOR VELOCITY V EQUATION***************
      IF(NF.EQ.2)THEN
      DV(I,J)=VOL/AP(I,J)
      PPY=(P(I,J+1)-P(I,J-1))/(YDIF(J)+YDIF(J+1))
      TTVXY=(TPXY(I+1,J)-TPXY(I-1,J))/((XDIF(I)+XDIF(I+1))*SX(J))
      TTVYY=(TPYY(I,J+1)-TPYY(I,J-1))/(YDIF(J)+YDIF(J+1))
      IF(J.EQ.2) THEN
          PPY=(P(I,J+1)-P(I,J))/YDIF(J+1)
          TTVYY=(TPYY(I,J+1)-TPYY(I,J))/YDIF(J+1)
      END IF
      IF(J.EQ.M2) THEN
          PPY=(P(I,J)-P(I,J-1))/YDIF(J)
          TTVYY=(TPYY(I,J)-TPYY(I,J-1))/YDIF(J)
      END IF
      IF(I.EQ.2) THEN
          TTVXY=(TPXY(I+1,J)-TPXY(I,J))/(XDIF(I+1)*SX(J))
      END IF
      IF(I.EQ.L2) THEN
          TTVXY=(TPXY(I,J)-TPXY(I-1,J))/(XDIF(I)*SX(J))
      END IF 
      
      CON(I,J)=CON(I,J)-VOL*PPY
      CON(I,J)=CON(I,J)+VOL*TTVYY+VOL*TTVXY    !The source term for v 
      CONV(I,J)=CONV(I,J)+TTVYY+TTVXY          !
      APPV(I,J)=AP(I,J)
      END IF
*---------------------------------------------------------------
  603 CONTINUE
*:------------------------------------------------------------------
	IF(NSCHEME.EQ.1) THEN
      CALL CUBISTA 
	END IF
*********************COME HERE TO SLOVE VELOCITY U EQUATION*************

      IF (NF.EQ.1) THEN
      DO I=2,L3
	DO J=2,M2
	!-----------------------------------------
      VOL1=YCVR(J)*XCV(I)
      VOL2=YCVR(J)*XCV(I+1)
      DUP=VOL1/APPU(I,J)
      DUE=VOL2/APPU(I+1,J)
      IF(NSELECT.EQ.1) THEN
	DUFE=DUE*FX(I+1)+DUP*FXM(I+1)
	ELSE
	DUFE=(VOL2*FX(I+1)+VOL1*FXM(I+1))/
     &  (APPU(I+1,J)*FX(I+1)+APPU(I,J)*FXM(I+1))
      END IF
	!-----------------------------------------
      HUE=AIP(I+1,J)*U(I+2,J)+AIM(I+1,J)*U(I,J)+
     &  AJP(I+1,J)*U(I+1,J+1)+AJM(I+1,J)*U(I+1,J-1)+CONU(I+1,J)*VOL2
      HUP=AIP(I,J)*U(I+1,J)+AIM(I,J)*U(I-1,J)+
     &    AJP(I,J)*U(I,J+1)+AJM(I,J)*U(I,J-1)+CONU(I,J)*VOL1
      IF(NSELECT.EQ.1) THEN
	HUFE=HUE/APPU(I+1,J)*FX(I+1)+HUP/APPU(I,J)*FXM(I+1)
	ELSE
      HUFE=(HUE*FX(I+1)+HUP*FXM(I+1))/
     &	(APPU(I+1,J)*FX(I+1)+APPU(I,J)*FXM(I+1))
	END IF
	!----------------------------------------
	DPFE=(P(I+1,J)-P(I,J))/XDIF(I+1)/SX(J)
	RHOE=FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J)
	UE(I,J)=HUFE+(1-RELAX(1))*UE(I,J)		!Solve the face velocity u
     &  -DUFE*DPFE+RHOE/DT*DUFE*UEL(I,J)
      IF(NSELECT.EQ.3) THEN
      UE(I,J)=U(I+1,J)*FX(I+1)+U(I,J)*FXM(I+1)
	END IF
	!----------------------------------------
	END DO
      END DO
      END IF
************************************************************************
*********************COME HERE TO SLOVE VELOCITY V EQUATION*************
************************************************************************
      IF(NF.EQ.2) THEN
      DO I=2,L2
	DO J=2,M3
	!----------------------------------------------
      VOL1=YCVR(J)*XCV(I)
      VOL2=YCVR(J+1)*XCV(I)
	DVP=VOL1/APPV(I,J)
      DVN=VOL2/APPV(I,J+1)
      IF(NSELECT.EQ.1) THEN
	DVFN=DVN*FY(J+1)+DVP*FYM(J+1)
	ELSE
	DVFN=(VOL2*FY(J+1)+VOL1*FYM(J+1))/
     &  (APPV(I,J+1)*FY(J+1)+APPV(I,J)*FYM(J+1))
	END IF 
      !---------------------------------------------------
      HVN=AIP(I,J+1)*V(I+1,J+1)+AIM(I,J+1)*V(I-1,J+1)+
     &    AJP(I,J+1)*V(I,J+2)+AJM(I,J+1)*V(I,J)+CONV(I,J+1)*VOL2
      HVP=AIP(I,J)*V(I+1,J)+AIM(I,J)*V(I-1,J)+
     &    AJP(I,J)*V(I,J+1)+AJM(I,J)*V(I,J-1)+CONV(I,J)*VOL1
      IF(NSELECT.EQ.1) THEN
      HVFN=HVN/APPV(I,J+1)*FY(J+1)+HVP/APPV(I,J)*FYM(J+1)
	ELSE
      HVFN=(HVN*FY(J+1)+HVP*FYM(J+1))/
     &	(APPV(I,J+1)*FY(J+1)+APPV(I,J)*FYM(J+1))
	END IF
	!---------------------------------------------------
      DPFN=(P(I,J+1)-P(I,J))/YDIF(J+1)
      DTYYFN=(TPYY(I,J+1)-TPYY(I,J))/YDIF(J+1)
      DTXYFN=(TPXY(I+1,J)-TPXY(I,J))/XDIF(I+1)/SX(J)
	RHON=FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J)
	VN(I,J)=HVFN+(1-RELAX(2))*VN(I,J)        !Solve the face velocity v
     &  -DVFN*DPFN+RHON/DT*DVFN*VNL(I,J)
      IF(NSELECT.EQ.3) THEN
      VN(I,J)=V(I,J+1)*FY(J+1)+V(I,J)*FYM(J+1)
	END IF
	!---------------------------------------------------
      END DO
	END DO
      END IF
      
***********************************************************************
      CALL SOLVE   !Solve the momentum equation
      GO TO 600
***********************************************************************
*****************COME HERE TO SLOVE PRESSURE CORRECTION EQUEATION******
***********************************************************************
  300 CALL GAMSOR
      SMAX=0.
      SSUM=0.
      DO 410 J=2,M2
      DO 410 I=2,L2
      VOL=YCVR(J)*XCV(I)
  410 CON(I,J)=CON(I,J)*VOL
      DO 402 I=2,L2
      AR=R(1)*XCV(I)*RHO(I,1)
      CON(I,2)=CON(I,2)+AR*V(I,1)
  402 AJM(I,2)=0.
      DO 403 J=2,M2
      AR=ARX(J)*RHO(1,J)
      CON(2,J)=CON(2,J)+AR*U(1,J)
      AIM(2,J)=0.
      DO 403 I=2,L2
      IF(I.EQ.L2) GO TO 404
      FLOW=ARX(J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))*UE(I,J)

      DUe=ARX(J)*(DU(I+1,J)*RHO(I+1,J)*FX(I+1)+
     &      DU(I,J)*RHO(I,J)*FXM(I+1))
      CON(I,J)=CON(I,J)-FLOW
      CON(I+1,J)=CON(I+1,J)+FLOW
      AIP(I,J)=DUe/XDIF(I+1)
      AIM(I+1,J)=AIP(I,J)
      GO TO 405
  404 CON(I,J)=CON(I,J)-ARX(J)*RHO(I+1,J)*U(I+1,J)
      AIP(I,J)=0.
	!---------------------------------------------------
  405 IF(J.EQ.M2) GO TO 406
      FLOW=RMN(J+1)*XCV(I)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
     &	*VN(I,J)
      DVn=RMN(J+1)*XCV(I)*(DV(I,J)*RHO(I,J)*FYM(J+1)+
     &      DV(I,J+1)*RHO(I,J+1)*FY(J+1))
      CON(I,J)=CON(I,J)-FLOW
      CON(I,J+1)=CON(I,J+1)+FLOW
      AJP(I,J)=DVn/YDIF(J+1)
      AJM(I,J+1)=AJP(I,J)
      GO TO 407
  406 CON(I,J)=CON(I,J)-RMN(M1)*XCV(I)*RHO(I,J+1)*V(I,J+1)
      AJP(I,J)=0.
	!------------------------------------------------
  407 AP(I,J)=AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)
      PC(I,J)=0.
      SMAX=DMAX1(SMAX,DABS(CON(I,J)))
      SSUM=SSUM+CON(I,J)
  403 CONTINUE
      CALL SOLVE
***************************************************************************
****************COME HERE TO CORRECT THE PRESSURE AND VELOCITIES***********
***************************************************************************
	DO I=1,L1
	P(I,1)=P(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -P(I,3)*(YDIF(2))/YDIF(3)
	P(I,M1)=P(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -P(I,M3)*(YDIF(M1))/YDIF(M2)
	END DO
	DO J=1,M1
	P(1,J)=P(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -P(3,J)*(XDIF(2))/XDIF(3)
	P(L1,J)=P(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -P(L3,J)*(XDIF(L1))/XDIF(L2)
	END DO

      DO 501 J=2,M2
      DO 501 I=2,L2
      P(I,J)=P(I,J)+PC(I,J)*RELAX(NP)
      PCX=(PC(I+1,J)-PC(I-1,J))/(XDIF(I)+XDIF(I+1))
      PCY=(PC(I,J+1)-PC(I,J-1))/(YDIF(J)+YDIF(J+1))
	IF(I.EQ.2)  PCX=(PC(I+1,J)-PC(I,J))/XDIF(I+1)
      IF(J.EQ.2)  PCY=(PC(I,J+1)-PC(I,J))/YDIF(J+1)
      IF(I.EQ.L2) PCX=(PC(I,J)-PC(I-1,J))/XDIF(I)
      IF(J.EQ.M2) PCY=(PC(I,J)-PC(I,J-1))/YDIF(J)
      U(I,J)=U(I,J)-DU(I,J)*PCX
      V(I,J)=V(I,J)-DV(I,J)*PCY
  501 CONTINUE

*---------------------------------------------------------------
      DO I=2,L2
          DO J=2,M2
              RESU=U(I,J)-FF(I,J,1)
              RESV=V(I,J)-FF(I,J,2)
          IF(RESU.GT.RESUMAX) THEN 
              RESUMAX=RESU
          END IF
          IF(RESV.GT.RESVMAX) THEN 
              RESVMAX=RESV
          END IF
          END DO
      END DO
      
  600 CONTINUE
  898 CONTINUE
      RETURN
      
      
      ENTRY SETUP3
      DO NNN=1,NFMAX
      DO I=1,L1
      DO J=1,M1
      FF(I,J,NNN)=F(I,J,NNN)   
      END DO
      END DO
      END DO

      
!     Decompose the velocity gradient
          DO I=2,L2
          DO J=2,M2
              UFX(I,J)=(U(I+1,J)-U(I-1,J))/((XDIF(I)+XDIF(I+1))*SX(J))
              VFY(I,J)=(V(I,J+1)-V(I,J-1))/(YDIF(J)+YDIF(J+1))
              UFY(I,J)=(U(I,J+1)-U(I,J-1))/(YDIF(J)+YDIF(J+1))
              VFX(I,J)=(V(I+1,J)-V(I-1,J))/((XDIF(I)+XDIF(I+1))*SX(J))
          IF(I.EQ.2) THEN
          UFX(I,J)=(0.5*(U(I+1,J)+U(I,J))-U(1,J))/XCV(I)
          VFX(I,J)=(0.5*(V(I+1,J)+V(I,J))-V(1,J))/XCV(I)
          ENDIF
          IF(I.EQ.L2) THEN
          UFX(I,J)=(U(L1,J)-0.5*(U(I,J)+U(I-1,J)))/XCV(I)
          VFX(I,J)=(V(L1,J)-0.5*(V(I,J)+V(I-1,J)))/XCV(I)
          END IF  
          IF(J.EQ.2) THEN
          VFY(I,J)=0.5*(V(I,J+1)+V(I,J))/YCVR(J)
          UFY(I,J)=0.5*(U(I,J+1)+U(I,J))/YCVR(J)
          END IF 
          IF(J.EQ.M2) THEN
          VFY(I,J)=(V(I,M1)-0.5*(V(I,J)+V(I,J-1)))/YCVR(J)
          UFY(I,J)=(U(I,M1)-0.5*(U(I,J)+U(I,J-1)))/YCVR(J)
          END IF
          END DO
          END DO  
      
      
          DO J=2,M2
          DO I=2,L2
          AX2(I,J)=AX(I,J)**2
        	AY2(I,J)=AY(I,J)**2
        	AXY(I,J)=AX(I,J)*AY(I,J)
        	AXY2(I,J)=AX2(I,J)-AY2(I,J)
        	ALB1(I,J)=EXP(P1(I,J))
        	ALB2(I,J)=EXP(P2(I,J))
          END DO
          END DO
      DO J=2,M2
	DO I=2,L2
		CMXX(I,J)=AX2(I,J)*UFX(I,J)+AXY(I,J)*VFX(I,J)+
     &              AXY(I,J)*UFY(I,J)+AY2(I,J)*VFY(I,J)
          
      	CMXY(I,J)=AXY(I,J)*VFY(I,J)-AXY(I,J)*UFX(I,J)-
     &              AY2(I,J)*VFX(I,J)+AX2(I,J)*UFY(I,J)
                   
      	CMYX(I,J)=-AXY(I,J)*UFX(I,J)+AX2(I,J)*VFX(I,J)-
     &              AY2(I,J)*UFY(I,J)+AXY(I,J)*VFY(I,J)
          
		CMYY(I,J)=AY2(I,J)*UFX(I,J)-AXY(I,J)*VFX(I,J)-
     &             AXY(I,J)*UFY(I,J)+AX2(I,J)*VFY(I,J)
      
	END DO
      END DO	
      
!     Solve the anti-symmetric matrices and symmetric matrix
      DO J=2,M2
	DO I=2,L2
        IF(ABS(ALB1(I,J)-ALB2(I,J)).LT.SMALL) THEN
        OMGXY(I,J)=0      	
        BXX(I,J)=UFX(I,J)
        BYY(I,J)=VFY(I,J)
        BXY(I,J)=0.5*(UFY(I,J)+VFX(I,J))
        ELSE 
         OMGXY(I,J)=(ALB2(I,J)*CMXY(I,J)+ALB1(I,J)*CMYX(I,J))/
     &               (ALB2(I,J)-ALB1(I,J))

          BXX(I,J)=AX2(I,J)*CMXX(I,J)+AY2(I,J)*CMYY(I,J)
       	BXY(I,J)=AXY(I,J)*(CMXX(I,J)-CMYY(I,J))
      	BYY(I,J)=AY2(I,J)*CMXX(I,J)+AX2(I,J)*CMYY(I,J)   
         END IF
          
      END DO
      END DO	
      


      

      DO 899 KIJ=1,5
      IST=2             
      JST=2             

      DO 700 N=4,NFMAX
      NF=N

      IF(.NOT. LSOLVE(NF)) GO TO 700 
*--------------------------------------------------------------------
***********************************************************************
*****************COME HERE TO SLOVE LOG-CONFORMATION R. EQUEATION******
***********************************************************************

      CALL GAMSOR   
      REL=1.-RELAX(NF)

      DO 702 I=2,L2
      AREA=R(1)*XCV(I)
      FLOW=AREA*V(I,1)
C----------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWY(I,J)=FLOW
	END IF
C------------------------------------------
      DIFF=AREA*GAM(I,1)/YDIF(2)
      CALL DIFLOW
  702 AJM(I,2)=ACOF+DMAX1(0.D0,FLOW)

      DO 703 J=2,M2
      FLOW=ARX(J)*U(1,J)
C----------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWX(I,J)=FLOW
	END IF
C------------------------------------------
      DIFF=ARX(J)*GAM(1,J)/(XDIF(2)*SX(J))
      CALL DIFLOW
      AIM(2,J)=ACOF+DMAX1(0.D0,FLOW)
      
      DO 703 I=2,L2

      IF(I .EQ. L2) GO TO 704
      FLOW=ARX(J)*(FX(I+1)+FXM(I+1))*UE(I,J)
c---------------------------------------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWX(I,J)=FLOW
	END IF
C************************************************************************
      DIFF=ARX(J)*2.*GAM(I,J)*GAM(I+1,J)/((XCV(I)*GAM(I+1,J)+
     &   XCV(I+1)*GAM(I,J)+1.0E-30)*SX(J))
      GO TO 705
  704 FLOW=ARX(J)*U(L1,J)
C---------------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWX(I,J)=FLOW
	END IF
C--------------------------------------
      DIFF=ARX(J)*GAM(L1,J)/(XDIF(L1)*SX(J))
  705 CALL DIFLOW

      AIM(I+1,J)=ACOF+DMAX1(0.D0,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW

      AREA=RMN(J+1)*XCV(I)
      IF(J .EQ. M2) GO TO 706
      FLOW=RMN(J+1)*XCV(I)*(FY(J+1)+FYM(J+1))
     &	*VN(I,J)
C-----------------------------------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWY(I,J)=FLOW
	END IF
C************************************************************************
      DIFF=AREA*2.*GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+
     &   YCV(J+1)*GAM(I,J)+1.0E-30)
      GO TO 707
  706 FLOW=AREA*V(I,M1)
c----------------------------------
	IF(NSCHEME.EQ.1) THEN
	FLOWY(I,J)=FLOW
	END IF
c---------------------------------
      DIFF=AREA*GAM(I,M1)/YDIF(M1)
  707 CALL DIFLOW
      
      AJM(I,J+1)=DMAX1(0.D0,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW

      VOL=YCVR(J)*XCV(I)
      APT=1.0/DT

      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*FF(I,J,NF)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     &   /RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*F(I,J,NF)

      
********************COEFFICIENTS FOR LCR PSIXX EQUATION***************
      IF(NF.EQ.4)THEN
      WXYXX(I,J)=2.*(OMGXY(I,J)*AXY(I,J)*(P1(I,J)-P2(I,J))+BXX(I,J))
      COLDBXX(I,J)=
     & 1./GAMT(I,J)*(AX2(I,J)/ALB1(I,J)+AY2(I,J)/ALB2(I,J)-1.
     & -ALP*(AX2(I,J)/ALB1(I,J)+AY2(I,J)/ALB2(I,J)+
     & AX2(I,J)*ALB1(I,J)+AY2(I,J)*ALB2(I,J)-2.))
      CON(I,J)=CON(I,J)+WXYXX(I,J)*VOL+COLDBXX(I,J)*VOL     
      END IF
********************COEFFICIENTS FOR LCR PSIYY EQUATION***************
      IF(NF.EQ.5)THEN
      WXYYY(I,J)=2.*(OMGXY(I,J)*AXY(I,J)*(P2(I,J)-P1(I,J))+BYY(I,J))
      COLDBYY(I,J)=
     & 1./GAMT(I,J)*(AX2(I,J)/ALB2(I,J)+AY2(I,J)/ALB1(I,J)-1.
     &  -ALP*(AX2(I,J)/ALB2(I,J)+AY2(I,J)/ALB1(I,J)-1.+
     & AX2(I,J)*ALB2(I,J)+AY2(I,J)*ALB1(I,J)-1.))
      CON(I,J)=CON(I,J)+WXYYY(I,J)*VOL+COLDBYY(I,J)*VOL     
      END IF
********************COEFFICIENTS FOR LCR PSIXY EQUATION***************
      IF(NF.EQ.6)THEN
      WXYXY(I,J)=(OMGXY(I,J)*AXY2(I,J)*(P2(I,J)-P1(I,J)))+2*BXY(I,J)
      COLDBXY(I,J)=1./GAMT(I,J)*(AXY(I,J)*(1./ALB1(I,J)-1./ALB2(I,J))
     & -ALP*(AXY(I,J)*(1./ALB1(I,J)-1./ALB2(I,J))+
     & AXY(I,J)*(ALB1(I,J)-ALB2(I,J))))
      CON(I,J)=CON(I,J)+WXYXY(I,J)*VOL+COLDBXY(I,J)*VOL    
      END IF      
*---------------------------------------------------------------
*---------------------------------------------------------------

  703 CONTINUE 
	IF(NSCHEME.EQ.1) THEN
      CALL CUBISTA
	END IF
*:------------------------------------------------------------------  
      CALL SOLVE    

         
      
      
  700 CONTINUE
  899 CONTINUE  
      !Solve the conformation tensor and  inverse transformation to obtain polymeric stress
            
      DO J=2,M2
      DO I=2,L2
      TB(I,J)=TXX(I,J)+TYY(I,J)
      TBAC(I,J)=SQRT((TXX(I,J)-TYY(I,J))**2+4.*TXY(I,J)**2)
      P1(I,J)=0.5*(TB(I,J)+TBAC(I,J))
      P2(I,J)=0.5*(TB(I,J)-TBAC(I,J))
      TBAC(I,J)=SQRT((P1(I,J)-TYY(I,J))**2+TXY(I,J)**2)
      IF(TBAC(I,J).LT.SMALL1) THEN
    	AX(I,J)=1.0
    	AY(I,J)=0.0
      ELSE
    	AX(I,J)=(P1(I,J)-TYY(I,J))/TBAC(I,J)
    	AY(I,J)=(TXY(I,J))/TBAC(I,J)
      END IF 
      END DO
      END DO
      
      DO I=2,L2
      DO J=2,M2
          CXX(I,J)=AX(I,J)**2*EXP(P1(I,J))+AY(I,J)**2*EXP(P2(I,J))
          CYY(I,J)=AY(I,J)**2*EXP(P1(I,J))+AX(I,J)**2*EXP(P2(I,J))
          CXY(I,J)=AX(I,J)*AY(I,J)*(EXP(P1(I,J))-EXP(P2(I,J)))
          TPXX(I,J)=GAMP(I,J)/GAMT(I,J)*(CXX(I,J)-1.)
          TPYY(I,J)=GAMP(I,J)/GAMT(I,J)*(CYY(I,J)-1.)
          TPXY(I,J)=GAMP(I,J)/GAMT(I,J)*CXY(I,J)
      END DO
      END DO


*---------------------插值获得边界上tp--------------------------
	DO I=1,L1
	TPXX(I,1)=TPXX(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -TPXX(I,3)*(YDIF(2))/YDIF(3)
	TPXX(I,M1)=TPXX(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -TPXX(I,M3)*(YDIF(M1))/YDIF(M2)
      
      TPYY(I,1)=TPYY(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -TPYY(I,3)*(YDIF(2))/YDIF(3)
	TPYY(I,M1)=TPYY(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -TPYY(I,M3)*(YDIF(M1))/YDIF(M2)
      
      TPXY(I,1)=TPXY(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -TPXY(I,3)*(YDIF(2))/YDIF(3)
	TPXY(I,M1)=TPXY(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -TXY(I,M3)*(YDIF(M1))/YDIF(M2)
	END DO
	DO J=1,M1
	TPXX(1,J)=TPXX(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -TPXX(3,J)*(XDIF(2))/XDIF(3)
	TPXX(L1,J)=TXX(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -TPXX(L3,J)*(XDIF(L1))/XDIF(L2)
      
      TPYY(1,J)=TPYY(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -TPYY(3,J)*(XDIF(2))/XDIF(3)
	TPYY(L1,J)=TPYY(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -TPYY(L3,J)*(XDIF(L1))/XDIF(L2)
      
      TPXY(1,J)=TPXY(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -TPXY(3,J)*(XDIF(2))/XDIF(3)
	TPXY(L1,J)=TPXY(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -TPXY(L3,J)*(XDIF(L1))/XDIF(L2)
      END DO

            	DO I=1,L1
	TXX(I,1)=TXX(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -TXX(I,3)*(YDIF(2))/YDIF(3)
	TXX(I,M1)=TXX(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -TXX(I,M3)*(YDIF(M1))/YDIF(M2)
      
      TYY(I,1)=TYY(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -TYY(I,3)*(YDIF(2))/YDIF(3)
	TYY(I,M1)=TYY(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -TYY(I,M3)*(YDIF(M1))/YDIF(M2)
      
      TXY(I,1)=TXY(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -TPXY(I,3)*(YDIF(2))/YDIF(3)
	TXY(I,M1)=TXY(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -TXY(I,M3)*(YDIF(M1))/YDIF(M2)
	END DO
	DO J=1,M1
	TXX(1,J)=TXX(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -TXX(3,J)*(XDIF(2))/XDIF(3)
	TXX(L1,J)=TXX(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -TXX(L3,J)*(XDIF(L1))/XDIF(L2)
      
      TYY(1,J)=TYY(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -TYY(3,J)*(XDIF(2))/XDIF(3)
	TYY(L1,J)=TYY(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -TYY(L3,J)*(XDIF(L1))/XDIF(L2)
      
      TXY(1,J)=TXY(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -TXY(3,J)*(XDIF(2))/XDIF(3)
	TXY(L1,J)=TXY(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -TXY(L3,J)*(XDIF(L1))/XDIF(L2)
      END DO
      
	DO I=2,L2
	DO J=3,M2
	TPYYN(I,J)=(TPYY(I,J)+TPYY(I,J-1))/2.0
      TPXYN(I,J)=(TPXY(I,J)+TPXY(I,J-1))/2.0
	END DO
	END DO
     
      DO I=3,L2
	DO J=2,M2
	TPXXE(I-1,J)=(TPXX(I-1,J)+TPXX(I,J))/2.0
      TPXYE(I-1,J)=(TPXY(I-1,J)+TPXY(I,J))/2.0
      END DO
      END DO
      
      DO J=2,L2
      TPXXE(1,J)=TPXX(1,J)
      TPXXE(L2,J)=TPXX(L1,J)
      TPXYE(1,J)=TPXY(1,J)
      TPXYE(L2,J)=TPXY(L1,J)
      END DO
      
      DO I=2,M2
      TPYYN(I,1)=TPYY(I,1)
      TPYYN(I,M2)=TPYY(I,M1)
      TPXYN(I,1)=TPXY(I,1)
      TPXYN(I,M2)=TPXY(I,M1)
      END DO
      DO I=2,L2
          DO J=2,M2
              RESTXX=TXX(I,J)-FF(I,J,4)
              RESTYY=TYY(I,J)-FF(I,J,5)
              RESTXY=TXY(I,J)-FF(I,J,6)
          IF(RESTXX.GT.RESTXXMAX) THEN 
              RESTXXMAX=RESTXX
          END IF
          IF(RESTYY.GT.RESTYYMAX) THEN 
              RESTYYMAX=RESTYY
          END IF
          IF(RESTXY.GT.RESTXYMAX) THEN 
              RESTXYMAX=RESTXY
          END IF
          END DO
      END DO
      
      RES=MAX(RESUMAX,RESVMAX,RESTXXMAX,RESTYYMAX,RESTXYMAX)
      IF(RES.LT.TOL) THEN 
            IF(ITER.GT.1000) STOP
      ENDIF


      ITER=ITER+1		 
	TIME=TIME+DT     !advance in time 
      IF(ITER .GE. LAST) LSTOP=.TRUE. !when true, computations stop
      RETURN
      END SUBROUTINE SETUP
*     -----------------------------------------------------------------
      SUBROUTINE SUPPLY
      USE BANK

*************************************************************************
   10 FORMAT(1X,26(1H*),3X,A10,3X,26(1H*))
   20 FORMAT(1X,4H I =,I6,6I9)
   30 FORMAT(1X,1HJ)
   40 FORMAT(1X,I2,3X,1P7E9.2)
   50 FORMAT(1X,1H )

************************************************************************
      	ENTRY PRINT		
77    FORMAT (I6,4E16.6)  
						
*------------------------CALCULATE THE STREAM FUNCTION----------------------------------------
      IF(MOD(ITER,10).EQ.0) THEN
c-------------------------------------
      KKK=KKK+1

      DO I=1,L1
	VNN(I,1)=0
	VNN(I,M1)=0
	UEE(I,1)=0
	UEE(I,M1)=1.0
	END DO
	DO J=1,M1
	VNN(1,J)=0
	VNN(L1,J)=0
	UEE(1,J)=0
	UEE(L1,J)=0
	END DO
	DO I=2,L2
	DO J=3,M2
	VNN(I,J)=(V(I,J)+V(I,J-1))/2.0
	END DO
	END DO
     
      DO I=3,L2
	DO J=2,M2
	UEE(I-1,J)=(U(I-1,J)+U(I,J))/2.0
      END DO
	END DO
C-------------------------------------
      STREAM(2,2)=0.
      DO 82 I=2,L1
      IF(I .NE. 2) STREAM(I,2)=STREAM(I,2)+RHO(I-1,1)*VNN(I-1,2)  
     &   *R(1)*XCV(I-1)
      DO 82 J=3,M1
      RHOM=FX(I)*RHO(I,J-1)+FXM(I)*RHO(I-1,J-1)
   82 STREAM(I,J)=STREAM(I,J-1)-RHOM*UEE(I,J-1)*ARX(J-1)

	SALAST=STRMAX
      SILAST=STRMIN

	STRMAX=0.0
      STRMIN=0.0 
      END IF

      IF(MOD(ITER,100).EQ.0) THEN	
      WRITE(*,77) ITER,SMAX,SSUM,TIME,RES
	END IF

   !	IF(ITER.GT.8000.AND.SSUM.LT.1.0D-8.AND
   !  &.DIFFMAX.LT.0.1)THEN
   	DO I=1,L1
	P(I,1)=P(I,2)*(YDIF(2)+YDIF(3))/YDIF(3)
     &   -P(I,3)*(YDIF(2))/YDIF(3)
	P(I,M1)=P(I,M2)*(YDIF(M1)+YDIF(M2))/YDIF(M2)
     &   -P(I,M3)*(YDIF(M1))/YDIF(M2)
	END DO
	DO J=1,M1
	P(1,J)=P(2,J)*(XDIF(2)+XDIF(3))/XDIF(3)
     &   -P(3,J)*(XDIF(2))/XDIF(3)
	P(L1,J)=P(L2,J)*(XDIF(L1)+XDIF(L2))/XDIF(L2)
     &   -P(L3,J)*(XDIF(L1))/XDIF(L2)
	END DO
c--------------------------------     
	OPEN(15,FILE='U.DAT')
	DO J=1,M1
	WRITE(15,*) Y(J),U((L1+1)/2,J)
	END DO
	CLOSE(15)

	OPEN(16,FILE='V.DAT')
	DO I=1,L1
	WRITE(16,*) X(I),V(I,(M1+1)/2)
	END DO
	CLOSE(16)
      
      OPEN(16,FILE='TXX.DAT')
	DO J=1,M1
	WRITE(16,*) Y(J),TXX((L1+1)/2,J)
	END DO
	CLOSE(16)
      
      OPEN(16,FILE='TXY.DAT')
	DO J=1,M1
	WRITE(16,*) Y(J),TXY((L1+1)/2,J)
	END DO
	CLOSE(16)
      
      OPEN(16,FILE='TYY.DAT')
	DO J=1,M1
	WRITE(16,*) Y(J),TYY((L1+1)/2,J)
	END DO
	CLOSE(16)
  444 FORMAT(1X,F7.4,F7.4,F10.6,F10.6)   
      OPEN(16,FILE='velocity.DAT')
	DO I=1,L1
          DO J=1,M1
	WRITE(16,444) X(I),Y(J),U(I,J),V(I,J)
          END DO
          ENDDO
	CLOSE(16)
      

C--------------------------------------------------
      OPEN(7,FILE="STREAM.DAT")
      WRITE(7,*) 'TITLE="STREAM"'
 	WRITE(7,*) 'VARIABLE="X","Y","STREAM"'
 	WRITE(7,*) 'ZONE T="STREAM"'," I=",L1," J=",M1, " C=BLACK"
       DO J=1,M1
       DO I=1,L1
       WRITE(7,*) X(I),Y(J),STREAM(I,J)
       END DO
       END DO
 	CLOSE(7)



c-------------------------------------------------------
      RETURN
	
	ENTRY UGRID
      DO 11 I=2,L1
 11   XU(I)=XL*((FLOAT(I)-2.)/(FLOAT(L1)-2.))		
      DO 22 J=2,M1
 22   YV(J)=YL*((FLOAT(J)-2.)/(FLOAT(M1)-2.))		
      RETURN
      END SUBROUTINE SUPPLY

************************************************************************
      SUBROUTINE CUBISTA		
      USE BANK
*-------------------  ----------------------------------------------------
*----- X DIRECTION CUBISTA SCHEME SOURCE ---------------------------------
      DO 888 I=IST,L2
	DO 888 J=JST,M2
      IF (I.EQ.L2) GOTO 128
      IF (FLOWX(I,J).GE.0) THEN  
	UP=F(I-1,J,NF)
 	CENT=F(I,J,NF)
	DOWN=F(I+1,J,NF)
	ELSE
	UP=F(I+2,J,NF)
 	CENT=F(I+1,J,NF)
	DOWN=F(I,J,NF)
	END IF
      C1=(CENT-UP)/(DOWN-UP)
!      IF (FLOWX(I,J).GE.0) THEN
!	X1=-XCV(I)-XCV(I-1)/2.0D+0
!	X2=-XCV(I)/2.0D+0
!	X3=XCV(I+1)/2.0D+0
          IF(C1.LE.0.OR.C1.GE.1) THEN
              ALPHA=1.0
              BATA=0.0
          ELSE IF (C1.GT.0.75) THEN
              ALPHA=0.25
              BATA=0.75
          ELSE IF (C1.GT.0.375) THEN
              ALPHA=0.75
              BATA=0.375
          ELSE
              ALPHA=1.75
              BATA=0
          END IF
!      ELSE
!	X1=-XCV(I+1)-XCV(I+2)/2.0D+0
!	X2=-XCV(I+1)/2.0D+0
!	X3=XCV(I)/2.0D+0
!      END IF

!	A1=X2*X3/((X2-X1)*(X3-X1))
!	A2=X1*X3/((X1-X2)*(X3-X2))
!	A3=X1*X2/((X1-X3)*(X2-X3))
      A1=1.0-ALPHA-BATA
      A2=ALPHA
      A3=BATA
	FAIFACEX=A1*UP+A2*CENT+A3*DOWN
	FACECONX=FLOWX(I,J)*(CENT-FAIFACEX)
	CON(I,J)=CON(I,J)+FACECONX
      CON(I+1,J)=CON(I+1,J)-FACECONX	
128   CONTINUE

*----- Y DIRECTION CUBISTA SCHEME SOURCE ---------------------------------
      IF (J.EQ.M2) GOTO 168
	IF (FLOWY(I,J).GE.0) THEN
	UP=F(I,J-1,NF)
 	CENT=F(I,J,NF)
	DOWN=F(I,J+1,NF)
	ELSE
	UP=F(I,J+2,NF)
 	CENT=F(I,J+1,NF)
	DOWN=F(I,J,NF)
	END IF
      C1=(CENT-UP)/(DOWN-UP)
!	IF (FLOWY(I,J).GE.0) THEN
!	X1=-YCV(J)-YCV(J-1)/2D0
!	X2=-YCV(J)/2D0
!	X3=YCV(J+1)/2D0
!	ELSE
!	X1=-YCV(J+1)-YCV(J+2)/2D0
!	X2=-YCV(J+1)/2D0
!	X3=YCV(J)/2D0
!	END IF
          IF(C1.LE.0.OR.C1.GE.1) THEN
              ALPHA=1.0
              BATA=0.0
          ELSE IF (C1.GT.0.75) THEN
              ALPHA=0.25
              BATA=0.75
          ELSE IF (C1.GT.0.375) THEN
              ALPHA=0.75
              BATA=0.375
          ELSE
              ALPHA=1.75
              BATA=0
      END IF
!	A1=X2*X3/((X2-X1)*(X3-X1))
!	A2=X1*X3/((X1-X2)*(X3-X2))
!	A3=X1*X2/((X1-X3)*(X2-X3))
      A1=1.0-ALPHA-BATA
      A2=ALPHA
      A3=BATA
      FAIFACEY=A1*UP+A2*CENT+A3*DOWN
	FACECONY=FLOWY(I,J)*(CENT-FAIFACEY)
	CON(I,J)=CON(I,J)+FACECONY
      CON(I,J+1)=CON(I,J+1)-FACECONY
168   CONTINUE
888   CONTINUE
      RETURN
	END SUBROUTINE CUBISTA
