C N. R. BADNELL     PROGRAM ADASDR       NRB v2.25              08/03/17
C
C***********************************************************************
C
C          POST-PROCESSOR FOR  ** AUTOSTRUCTURE ** (for ADAS)
C         ***************************************************
C
C CALCULATES CA/LS/IC (+HYBDRID) PARTIAL DR RATE COEFFS IN ADF09 FORMAT
C  *** N.B. PROCESS DN=0 AND DN=1 CORE EXCITATIONS SEPARATELY. ***
C ALSO INCORPORATES THE MOST FREQUENTLY USED MDRCS OPTIONS, E.G. CAN
C GENERATE BINNED CROSS SECTIONS (+/- BYPASSING ADF09) AND CONVOLUTE
C WITH GAUSSIAN, COOLER OR MAXWELLIAN DISTRIBUTIONS. 
C
C***********************************************************************
cparc                                                               !par
cparc                    + Parallel +                               !par
cparc                                                               !par
cparc   Only useful for totals since gives separate adf09 per proc  !par
cparc                                                               !par
cpar!***************************************************************!par
cparc                                                               !par
cpar      module comm_interface                                     !par
cparc                                                               !par
cpar      use mpi                                                   !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      public comm_init          ! Initialize MPI                !par
cpar      public comm_barrier       ! MPI barrier                   !par
cpar      public comm_finalize      ! Terminate MPI                 !par
cpar      integer*4, public  :: iam                                 !par
cpar      integer*4, public  :: nproc                               !par
cparc                                                               !par
cpar      SAVE                                                      !par
cparc                                                               !par
cpar      private                                                   !par
cpar      integer*4 :: mpicom                                       !par
cparc                                                               !par
cpar      CONTAINS                                                  !par
cparc                                                               !par
cpar!---------------------------------------------------------------!par
cpar      subroutine comm_init()                                    !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      integer*4 :: ier                                          !par
cparc                                                               !par
cpar      mpicom = MPI_COMM_WORLD                                   !par
cparc                                                               !par
cpar      call mpi_init(ier)                                        !par
cpar      call mpi_comm_rank(mpicom, iam, ier)                      !par
cpar      call mpi_comm_size(mpicom, nproc, ier)                    !par
cparc                                                               !par
cpar      return                                                    !par
cparc                                                               !par
cpar      end subroutine comm_init                                  !par
cparc                                                               !par
cpar!---------------------------------------------------------------!par
cpar      subroutine comm_barrier()                                 !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      integer*4 :: ier                                          !par
cparc                                                               !par
cpar      call mpi_barrier(mpicom, ier)                             !par
cparc                                                               !par
cpar      return                                                    !par
cparc                                                               !par
cpar      end subroutine comm_barrier                               !par
cpar!---------------------------------------------------------------!par
cparc                                                               !par
cpar      subroutine comm_finalize()                                !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      integer*4 :: ier                                          !par
cparc                                                               !par
cpar      call mpi_finalize(ier)                                    !par
cparc                                                               !par
cpar      return                                                    !par
cparc                                                               !par
cpar      end subroutine comm_finalize                              !par
cpar!---------------------------------------------------------------!par
cparc                                                               !par
cpar      end module comm_interface                                 !par
cparc                                                               !par
cparc***************************************************************!par
C
      PROGRAM MAIN
cparc                                                               !par
cpar      use comm_interface, only : iam,nproc,comm_init,           !par
cpar     A                           comm_barrier,comm_finalize     !par
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*2 NAM0
C
C SUN TIME
      REAL*4 TARRY(2),TIME
C
      NAM0=''
cparc                                                               !par
cpar      call comm_init()                                          !par
cpar      write(0,*)'Starting proc', iam                            !par
cparc                                                               !par
cpar      ic1=iam/10                                                !par
cpar      ic2=iam-10*ic1                                            !par
cpar      ich0=ichar('0')                                           !par
cpar      ic1=ic1+ich0                                              !par
cpar      ic2=ic2+ich0                                              !par
cpar      nam0=char(ic1)//char(ic2)                                 !par
cpar      OPEN(5,FILE='adasin')                                     !par
C
C      OPEN(5,FILE='adasin')                                      !STDIN
      OPEN(6,FILE='adasout'//nam0)                               !STDOUT
C      OPEN(7,FILE='ocs')                         !BINNED CROSS SECTIONS
C      OPEN(9,FILE='CAVES/TERMS/LEVELS')           !OPT TARGET SYMM/ENER
      OPEN(10,FILE='adf09'//nam0)                    !ADAS FORMAT OUTPUT
C      OPEN(14,FILE='XDRTOT')                         !CONVOLUTED TOTALS
C      OPEN(70,FILE='o_str')          !DETAILED TARGET DATA (FOR HYBRID)
C      OPEN(70,FILE='ou_str',FORM='UNFORMATTED')   !DITTO       (UNFORM)
C THEN
C      OPEN(70,FILE='on')                   !AUTOS DATA FILE (FORMATTED)
C OR
C      OPEN(70,FILE='onu',FORM='UNFORMATTED')  !AUTOS DATA FILE (UNFORM)
C
      CALL POSTP
C
C SUN TIME
      DUM=DTIME(TARRY)
      TIME=TARRY(1)
C
C CRAY TIME
CCRAY CALL SECOND(TIME)
C
      TIME=TIME/60.0
      WRITE(6,999)TIME
 999  FORMAT(//1X,'CPU TIME=',F9.3,' MIN')
C
      WRITE(6,*) 'PROGRAM ADASDR: NORMAL END'
C
      CLOSE(6)
C      CLOSE(7)
C      CLOSE(9)
      CLOSE(10)
C      CLOSE(14)
C      CLOSE(70)
cparc                                                               !par
cpar      write(0,*)'Ending proc', iam                              !par
cpar      call comm_barrier()                                       !par
cpar      call comm_finalize()                                      !par
C
      END
C
C***********************************************************************
C
      SUBROUTINE POSTP
cparc                                                               !par
cpar      use comm_interface, only : iam,nproc,comm_init,           !par
cpar     A                           comm_barrier,comm_finalize     !par
cpar      use mpi                                                   !par
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (NDIM1=10001)
      PARAMETER (NDIM2=20)
      parameter (ndim4=100)
      PARAMETER (NDIM5=400)
      PARAMETER (NDIM8=NDIM5)
      PARAMETER (NDIM24=50000)
      PARAMETER (NDIM27=150)
      PARAMETER (NDIM25=NDIM27)
      PARAMETER (NDIM28=19)
      PARAMETER (NDIM31=10)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
C
      CHARACTER NAME*30,DATE*30,DATE8*8
C
      CHARACTER*4 NAMEJ,RAD,COD(20),MERGE
      CHARACTER*5 COREX
      CHARACTER*1 COR(5),LIT(0:9)
c
cpar      real*8 ssend(NDIM1),srecv(NDIM1)                          !par
C
      EQUIVALENCE (COR(1),COREX)
C
      LOGICAL BCA,BLS,BIC,BLSOLD,BLSNEW,BLOG,EX,BHYBRD,bpart
C
      DIMENSION ECORI(NDIM8,NDIM8),IWT(NDIM5),IWS(NDIM5),IWL(NDIM5)
     X         ,EI(NDIM5),IWJ(NDIM5),ILVTM(NDIM5),LCP(NDIM5)
     X         ,ALF(NDIM28,NDIM2),TEMP(NDIM28)
C
      DIMENSION EBIN(NDIM1),SBIN(NDIM1,NDIM2),EET(NDIM2)
C
      COMMON /CORR/ACORN(NDIM25),ACORL(NDIM31),RMIN,NNCOR,NLCOR,NCMN
     X            ,NCMX,LCMN,LCMX,IMATCH,RAD,NFNLMX,FNL(NDIM24)
      COMMON /DITT/A0,B0
      COMMON /ECOR/E1C(NDIM8),E1X(NDIM8),TOLB,TOLB0,TOLBE
      COMMON /JCF/JCFA,JCFR,JCFJ,LSPI,J2PI,JPAR
      COMMON /QDTS/QDTS(0:30),NQDT
      common /part/epart,dee,frake,nparti,w0(ndim4)
C
      DATA LIT /'0','1','2','3','4','5','6','7','8','9'/
C
      NAMELIST/ONE/NTAR1,NTAR2,IPRINT,NCUT,LCUT,JCFR,JCFJ,NRSLMX,RAD
     X     ,LSPI,J2PI,NMIN,LMIN,NMAX,LMAX,IRD,NLMAX,COREX,NCMX,NCMN
     X     ,LCMN,LCMX,IMATCH,UNITS,JPAR,nc,JCFA,ILOG,MR5,iflagw,nxtrp
C
      NAMELIST/TWO/EMIN,EMAX,NR1,NECOR,NNCOR,NLCOR,ACOR,RCOR
     X       ,TOLR,TOLB,TOLI,TOLBE,IOLDE,LVAMX,RMIN,NR2,irwt
     X       ,JTEMP,IRDT,JTHETA,IREL,NRB,IOLDW
     X       ,NBIN,NQDT,NFNLMX,EWIDTH,TPAR,TPER,ESWTCHX,NGAUSS,TC1
     x       ,nbins,epart,dee,frake,nparti,w0
C
C
C READ HEADER TO DETERMINE IF CA, LS OR IC RUN (/CA/, /LS/ OR /IC/), 
C /  / DEFINES THE OLD STYLE (93) LS.
C THE REST OF THE LINE IS FOR COMMENT, IT IS ADDED TO THE END OF ADF09.
C
      READ(5,1000)COD
      WRITE(6,1001)COD
C
      BCA=COD(1).EQ.'/CA/'
      BLSNEW=COD(1).EQ.'/LS/'
      BLSOLD=COD(1).EQ.'/  /'
      BLS=BLSOLD.OR.BLSNEW
      BIC=COD(1).EQ.'/IC/'
      IF(.NOT.BCA.AND..NOT.BLS.AND..NOT.BIC)THEN
        WRITE(6,1002)COD(1)
        STOP 'ERROR: INPUT ERROR ON UNIT5'
      ENDIF
C
C
C**********************NAMELIST-ONE*************************************
C
C           ONLY THE STARRED '***' INPUT IS COMPULSORY
C-----------------------------------------------------------------------
C *** NTAR1 = NO OF INITIAL STATES POPULATED.
C *** NTAR2 = NO OF FINAL (PARENT) STATES 
C           = 0 ALL FOUND IN FIXED FORMAT WEIGHT LIST (NEW DEFAULT)
C           > 0 FIRST NTAR1 OF THEM=INITIAL
C           < 0 THESE ARE CONFIGURATION RESOLVED (REQUIRES TARGET o_str)
C
C     *AFTER* NAMELIST-TWO, READ NTAR2 TARGET/PARENT TERMS/LEVELS INFO.
C     IF   IRD .EQ. 0 (DEFAULT /  /) READ
C   *** IWS(I)=2S+1, IWL(I)=L IF IWS(I).NE.0 ELSE IWL(I)=2J+1 IN IC.
C                       S,L AND J ARE TARGET MOMENTA OF COURSE.
C                       CA CAN USE EITHER TO HOLD W (L=0 IF 2S+1 USED).
C     ELSE IRD .LT. 0 READ
C   *** EI(I),IWS(I),IWL(I) WHERE EI(I) ARE ABSOLUTE TARGET BIN ENERGIES
C       I.E. LT. 0, WHICH SPAN GROUPS OF LEVELS/TERMS.
C        (CAN BE USED TO AVERAGE OVER LEVELS OF A TERM OR TERMS/LEVLS OF
C        A CONFIGURATION.) 
C        THE STAT WEIGHT(S) SHOULD BE THE SUM OF THOSE CONTAINED WITHIN
C        THE BIN(S). 
C        IF EI(I>1) .GT. 0 THEN BIN ENERGY ASSUMED RELATIVE TO EI(1)
C        I.E. EI(1) IS LIKELY THE ACTUAL GROUND STATE ENERGY AND SO
C        SUBSEQUENTLY IT IS LOWERED TO ENSURE IT IS A BIN ENERGY.
C     ELSE IRD .GT. 0   (DEFAULT /CA/, /LS/ OR /IC/) AFTER &TWO READ 
C       NTAR2 CFG/TERM/LEVEL INFO AS OUTPUT BY AUTOS CAVES/TERMS/LEVELS
C       FILES. VIZ.
C   ***            IW  IPAR      CA
C       OR
C   ***      IWS  IWL  IPAR      LS
C       OR
C   ***      IWJ  IPAR IWS  IWL   IC  (HERE IWJ=2J).
C
C   *** COREX='N-M' SELECT CORE SHELL EXCITATION, WHERE N AND M ARE
C                   THE INITIAL AND FINAL PRINCIPAL QUANTUM NUMBERS. 
C                   NOTE, THE CORE L VALUES MUST *NOT* BE SPECIFIED.
C                   HOWEVER:
C       COREX='NL-MK' SELECT CORE SUB-SHELL EXCITATION, WHERE NL AND MK
C                     ARE THE INITIAL AND FINAL PRINCIPAL AND ORBITAL
C                     ANGULAR MOMENTUM QUANTUM NUMBERS.
C                  ***THE INITIAL REFERENCE CONFIGURATION IS CF=IMATCH,
C           IMATCH=1, DEFAULT.
C                     NOTE: THE CORE L VALUES ARE NUMERIC, *NOT* 
C                     SPECTROSCOPIC E.G. '40-52' SELECTS 4s-5d(!)
C                     ALSO, '40-5 ' SELECTS 4s to ALL L, OF N=5.
C
C ***    END OF COMPULSORY INPUT    ***
C
C SINCE, IN GENERAL, THE AUTOSTRUCTURE DATA FILE(S) CONTAIN DATA FOR
C MORE THAN ONE CORE EXCITATION, THESE SHOULD BE PROCESSED SEPARATELY,
C AND IT IS OPTIMUM TO COMPUTE THEM SEPARATELY, HENCE THE ROLE OF COREX.
C       
C    NOTE, 
C      DEFAULT COREX='   ' CAUSES A STOP UNLESS JCFJ, NCMN OR NCMX SET -
C                          SEE BELOW. THIS IS TO STOP NOVICE USERS NOT
C                          SETTING ANYTHING. AN EXPERT USER CAN SET
C                          JCFJ, NCMN/NCMX OR INDEED A DUMMY COREX=' - '
C                          TO OVERRIDE, I.E. CODE ASSUMES YOU KNOW WHAT
C                          YOU ARE DOING THEN.
C                  
C-----------------------------------------------------------------------
C
C     NRSLMX -- IS THE MAX N OF RESOLVED DATA, DEFAULT=8.
C     NLMAX  -- IS THE MAX N OF BUNDLED-NL DATA (IC ONLY), DEFAULT=10.
C
C     JCFJ .GT. 0   NEGLECTS CAPTURE INTO CONFIGS .GT. JCFJ.
C                   CAN BE USED IN PLACE OF COREX IF THE AUTOSTRUCTURE
C                   CONFIGURATION LIST IS ORDERED AMENABLY.
C               *** DEFAULT: INCLUDE ALL. ***
C                   ALSO ATTEMPTS TO INCLUDE ANY CONTRIBUTION FROM
C                   NON-RYDBERG STATES (FROM THE FIRST NL BLOCK) IF
C                   NO EQUIVALENT ELECTRON FILE ALREADY PROCESSED.
C          .LT. 0   SWITCHES OFF THIS ATTEMPT.
C
C     IPRINT=PRINT LEVEL
C           .GE. 0, DETAILED PRINTOUT OF EACH PARTIAL CROSS SECTION
C           .EQ.-1, NL CROSS SECTIONS
C           .EQ.-2,  L CROSS SECTIONS
C           .LE.-3, TOTAL CROSS SECTION ONLY.
C
C     NCUT(OR MAX) .GT. 0 IGNORES CONTRIBUTIONS FROM N .GT. NCUT(OR MAX)
C     LCUT(OR MAX) .GE. 0         "         "        L .GT. LCUT(OR MAX)
C     NMIN .GT. 0, IGNORES CONTRIBUTIONS FROM N .LT. NMIN
C     LMIN .GE. 0,         "         "        L .LT. LMIN
C                            DEFAULT: INCLUDE ALL.
C     NOTE: IF NMIN,NMAX DO NOT MATCH A REPRESENTATIVE-N THEN THE TOTAL
C           IS TRUNCATED AT THE FIRST/LAST REP-N WHICH SATISFIES IT,
C           EVEN IF THERE ARE ADDITIONAL REP-N ON FILE. THUS, IT DIFFERS
C           FROM MDRCS13 IN THIS RESPECT BECAUSE PARTIALS ARE IN CONTROL
C
C     LSPI .GT. 0 INCLUDE PARTIAL WAVE LSPI=10000*(2S+1)+100*L+PI ONLY
C     J2PI .GE. 0 INCLUDE PARTIAL WAVE J2PI=100*(2*J)+PI ONLY
C                            DEFAULT: INCLUDE ALL.
C
C     RAD .EQ. 'YES' INCLUDE RADIATIVE WIDTH (DEFAULT).
C     RAD .EQ. 'BF'  INCLUDE ONLY RADIATIVE TRANSITIONS BETWEEN AUTO-
C                    IONIZING STATES AND TRUE BOUND STATES.
C         .EQ. 'NO'  NEGLECT RADIATIVE WIDTH.
C
C     JCFR .LT. 0 INCLUDES CAPTURE INTO CONFIGURATION NUMBER -JCFR
C               AS OUTPUT BY AUTOSTRUCTURE, FOR L (LV) .GE. 0 ONLY.
C     JCFR .GT. 0 ASSUMES ALL STATES OF FINAL CONFIGURATIONS .LE. JCFR 
C                STABLE AGAINST AUTOIONIZATION, UNLESS
C          .GT. 100 COMPLETELY IGNORES RADIATIVE DECAY INTO AUTOIONIZING
C                   STATES.
C          .GT. 200 NEGLECTS RADIATIVE WIDTH (OR SET RAD='NO')
C          .EQ. 0 DEFAULT, DOES NOTHING.
C
C     JCFA .GT. 0 OMIT CONTINUUM CONFIGURATIONS .GT. JCFA.
C                 THIS CAN BE USED TO OMIT AUTOIONIZATION INTO EXCITED
C                 STATES FOR EXAMPLE, PROVIDED CONFIGS SUITABLY ORDERED.
C                 *** THIS IS A DIFFERENT OPERATION FROM MDRCS13. ***
C          .LT. 0 OMIT CONTINUUM CONFIGURATIONS THAT RESULT FROM
C                 CORE RE-ARRANGEMENT, I.E NL+E-.
C
C     NCMN,NCMX DEFINE THE CORE TRANSITION. THE DEFAULT IS SET
C               BY COREX. BUT IF COREX IS UNDEFINED (DEFAULT)
C               THEN LEAVING ONE AND/OR OTHER UNDEFINED RESULTS IN
C               ALL INITIAL AND/OR FINAL EXCITATIONS BEING INCLUDED,
C               SUBJECT TO ANY RESTRICTIONS BY JCFJ. 
C               UNLESS YOU KNOW BETTER, SET COREX OR JCFJ.
C
C     UNITS ENERGY UNITS USED FOR CALC AND OBS TARGET ENERGY, EMIN,EMAX:
C           13.606 FOR EV,1.0 FOR RYDBERGS (DEFAULT).
C
C     NC .GT. 0 READ OLD ocsp/t FILES.
C
C     NOTE: TESTS INVOLVING LSPI, J2PI AND JCFR MAY HAVE BEEN
C     COMMENTED-OUT BY 'CT' FOR SPEED AS THEY ARE RARELY USED.
C
C
C****************************END-ONE************************************
C
C
      MR5=0
      NTAR1=1
      NTAR2=0
      IPRINT=-1
      NCUT=-66
      LCUT=-77
      JCFA=0
      JCFR=0
      JCFJ=0
      NRSLMX=8
      NLMAX=-9999
      RAD='YES'
      LSPI=0
      J2PI=-1
      NMIN=-10
      LMIN=-10
      NMAX=-10
      LMAX=-10
      IF(BLSOLD)IRD=0
      IF(BCA.OR.BIC.OR.BLSNEW)IRD=1
      NCMX=-1
      NCMN=-1
      LCMX=-1
      LCMN=-1
      IMATCH=0
      COREX='     '
      UNITS=DONE
      JPAR=0
      ILOG=0
      nc=0
      iflagw=1              !set =1 to suppress check on stat weights...
      nxtrp=-1
C
C-----------------------------------------------------------------------
C
      READ(5,ONE)
C
C-----------------------------------------------------------------------
C
      IF(ILOG.GT.0)ILOG=0
C
      IF(COREX.EQ.'     '.AND.NCMN.LT.0.AND.NCMX.LT.0.AND.JCFJ.EQ.0)THEN
        IFLAGX=1                                       !POTENTIAL STOP
      ELSE
        IFLAGX=0
      ENDIF
C
      IF(NCMX.LE.0)THEN
        I3=3
        IF(COR(3).EQ.'-')I3=4
        I2=I3-1
        DO I=1,I3,I2
          DO J=1,9
            IF(LIT(J).EQ.COR(I))THEN
              IF(NCMN.LT.0)THEN
                NCMN=J
              ELSE
                NCMX=J
                GO TO 103
              ENDIF
            ENDIF
          ENDDO
          IF(NCMN.LT.0)NCMN=999
        ENDDO
 103    CONTINUE
        IF(COR(2).EQ.'-')THEN
          IF(NCMX.GT.0)THEN
            WRITE(6,1004)NCMN,NCMX
          ELSE
            WRITE(6,1005)NCMN      !,NCMX
          ENDIF
          IF(COR(1).EQ.' '.OR.COR(1).EQ.'*'.OR.
     X       COR(3).EQ.' '.OR.COR(3).EQ.'*')IMATCH=-IABS(IMATCH)
          IF(IMATCH.GT.0)THEN
            IF(LCMN.LT.0)LCMN=999           !SWITCH TO CMP CF'S
          ENDIF
          IF(IMATCH.NE.0)WRITE(6,1009)IMATCH
        ELSEIF(COR(3).EQ.'-')THEN
          IF(COR(1).EQ.' '.OR.COR(1).EQ.'*'.OR.
     X       COR(4).EQ.' '.OR.COR(4).EQ.'*')THEN
            WRITE(6,1010)
            STOP 'MUST SPECIFY N-VALUES IN COREX FOR NL-SELECTION'
          ENDIF
          IF(IMATCH.LE.0)IMATCH=1           !THEN ASSUME FIRST IS GROUND
          DO I=2,5,3
            DO J=0,9
              IF(LIT(J).EQ.COR(I))THEN
                IF(LCMN.LT.0)THEN
                  LCMN=J
                ELSE
                  LCMX=J
                  GO TO 105
                ENDIF
              ENDIF
            ENDDO
            IF(LCMN.LT.0)LCMN=999
          ENDDO
 105      CONTINUE
          IF(LCMX.GE.0)THEN
            WRITE(6,1006)NCMN,LCMN,NCMX,LCMX
          ELSEIF(COR(5).EQ.' '.OR.COR(5).EQ.'*'.OR.COR(5).EQ.'l')THEN
            IF(LCMN.LT.10)THEN
              WRITE(6,1007)NCMN,LCMN,NCMX
            ELSEIF(COR(2).EQ.' '.OR.COR(2).EQ.'*'.OR.COR(2).EQ.'l')THEN
              WRITE(6,1008)NCMN,NCMX
            ELSE
              WRITE(6,1003)COREX
              STOP 'ERROR: COREX IMPROPERLY DEFINED'
            ENDIF
          ELSE
            WRITE(6,1003)COREX
            STOP 'ERROR: COREX IMPROPERLY DEFINED'
          ENDIF
          WRITE(6,1009)IMATCH
        ELSE
          WRITE(6,1003)COREX
          STOP 'ERROR: COREX IMPROPERLY DEFINED'
        ENDIF
      ENDIF
C
      IF(NTAR1.LT.1) THEN
        WRITE(6,849)NTAR1
        STOP 'ERROR: NTAR1 MUST BE .GT. 0'
      ENDIF
      IF(NTAR1.GT.NDIM2)THEN
        WRITE(6,847)NTAR1
        STOP 'ERROR: INCREASE NDIM2'
      ENDIF
      NBINI=NTAR1+1                                   !HISTORIC INTERNAL
      NBINM=NTAR1
C
      NTAR2O=NTAR2
      IF(NTAR2.LT.0)THEN                              !HYBRID PARENTS
        iflagw=1           !for now, assume incomplete (bundled) o-files
        IF(IRD.LE.0)THEN                              !NEED NEW INPUT
          WRITE(6,*)'HYBRID OPTION NTAR2.LT.0 REQUIRES IRD.GT.0'
          STOP 'ERROR: HYBRID OPTION NTAR2.LT.0 REQUIRES IRD.GT.0'
        ENDIF
        NBINR=NTAR2-1
      ELSE
        IF(NTAR2*NTAR2.LT.NTAR1*NTAR2)NTAR2=NTAR1    !.LT.ALLOWS NTAR2=0
        IF(NTAR2.EQ.0)NTAR2=NDIM5-1
        NBINR=NTAR2+1
      ENDIF
      NBINRM=NTAR2
      IF(IABS(NBINR).GT.NDIM5)THEN
        WRITE(6,848)IABS(NBINR)
 848    FORMAT(/' INCREASE NDIM5 TO AT LEAST',I5)
        STOP 'ERROR: INCREASE NDIM5'
      ENDIF
      BHYBRD=NBINR.LT.0          !FINAL PARENT BY CONFIG
C
      IF(RAD.EQ.'NO')JCFR=222
      IF(JCFA.EQ.0)JCFA=9999
      IF(JCFJ.EQ.0)THEN
        JCFJ=999999
      ELSEIF(JCFJ.LT.0)THEN
        JCFJ=-999999
      ENDIF
      JCFX=0
      NAMEJ='JCF*'
      IF(IABS(JCFJ).NE.999999)THEN
        JCFX=JCFJ
        NAMEJ='JCFJ'
      ENDIF
      IF(JCFR.NE.0)THEN
        JCFX=JCFR
        NAMEJ='JCFR'
      ENDIF
      IF(IABS(JCFA).NE.999999)THEN
        JCFX=JCFA
        NAMEJ='JCFA'
      ENDIF
C
      IF(NMAX.GT.0)NCUT=NMAX
      IF(LMAX.GT.-1)LCUT=LMAX
C
      WRITE(6,11) NTAR1,NTAR2,NMIN,LMIN,NCUT,LCUT,NAMEJ,JCFX
C
      IF(BLS.AND.LSPI.GT.0)WRITE(6,12)LSPI
      IF(BIC.AND.J2PI.GE.0)WRITE(6,6)J2PI
C
      IF(NCUT.LT.1)NCUT=1000
      IF(LCUT.LT.0)LCUT=1000
c
      if(nxtrp.le.0)then
        nxtrp=9999
      else
        write(0,*)'*** Warning, extrapolating from N=',nxtrp
        write(6,*)' '
        write(6,*)'*** Warning, extrapolating from N=',nxtrp
      endif
C
      IF(NLMAX.EQ.-9999)THEN                          !UNSET BY USER, SO
        IF(BCA.OR.NTAR2.LT.0)THEN               !CA OR FINAL CF RESOLVED
          NLMAX=100
        ELSE                                             !LS/IC RESOLVED
          NLMAX=10
        ENDIF
      ENDIF
      NLMAX=MIN(NLMAX,NCUT)
C
      NRSLMX=MIN(NRSLMX,NCUT)
c
      if(nrslmx.gt.nxtrp)then
        write(0,*)'*** Warning, reducing nrslmx=',nrslmx
     x           ,' to nxtrp=',nxtrp
        write(6,*)' '
        write(6,*)'*** Warning, reducing nrslmx=',nrslmx
     x           ,' to nxtrp=',nxtrp
        nrslmx=nxtrp
      endif
C
C*************************NAMELIST-TWO**********************************
C
C  NOTHING IS COMPULSORY BUT NR1 AND NECOR (DN=0 ONLY) ARE RECOMMENDED
C-----------------------------------------------------------------------
C
C    NECOR .GT. 0: READ CALCULATED (E1C) AND OBSERVED (E1X) TARGET/CORE
C                  ENERGIES (IN RYD) IMPORTANT FOR DN=0, NOT SO FOR DN=1
C                  E1C(J) ARE CALCULATED CORE ENERGIES, J=1,NECOR 
C                  E1X(J) ARE OBSERVED CORE ENERGIES, J=1,NECOR.
C                  THESE ARE USED TO CORRECT THE POSITION OF THE DR
C                  RESONANCES - WHICH GOVERNS ITS IMPORTANCE, OR NOT.
C                  E1C(1) (AND E1X(1)) ARE NORMALLY ZERO. HOWEVER, IF
C                  ECOR1,2 WAS SET NON-ZERO IN AUTOSTRUCTURE THEN SET
C                  TC1 OR E1C(1)=-ECOR1,2 TO REMOVE ARTIFICIAL SHIFT.
C                  SINCE V1.17, E1C ARE READ FROM ALONGSIDE TERMS/LEVELS
C                  INFO, BY DEFAULT, AND ONLY OBSERVED ARE READ-IN
C                  AFTERWARDS. IN ADDITION, E1C ARE CHECKED AGAINST
C                  INTERNAL NOW AND IF THEY DIFFER THEN INTERNAL ARE
C                  USED INSTEAD. SEE IOLDW AND TOLBE BELOW TO OVERRIDE.
C
C          .EQ. 0: THEORETICAL POSITIONS ARE USED (OK 4 DN=1)**DEFAULT**
C
C     NR1 (SIMPLE) SET IT EQUAL TO THE LARGEST CORE PRINCIPAL QUANTUM
C                  NUMBER (N) +1. THIS IS OFTEN THE LOWEST NV. 
C                  **DEFAULT** ATTEMPTS TO DETERMINE NR1 INTERNALLY.  
C                  IF IT FAILS THEN NR1 IS SET .EQ. 0; NOT GOOD FOR IC, 
C                  ALSO NO OUTER RADIATION THEN (SEE BELOW).
C
C     NR1 (COMPLEX)
C          .EQ. 0: RECOMBINED PARENTS ARE DETERMINED ON ENERGY GROUNDS
C                  ALONE, BY SUBTRACTING THE ENERGY OF THE RYD. STATE.
C                  THIS IS THE SAME AS THE OLD NR1.EQ.0 OPTION. (IT CAN
C                  ALSO BE FORCED BY SETTING NECOR.LT.0, MAINLY FOR USE
C                  OF ENERGY CORRECTIONS ALONE IN COMPLEX SYSTEMS.)
C                  
C          .NE. 0: THEY ARE DETERMINED BY ENERGY ORDERING AND ANGULAR
C                  MOMENTUM SYMMETRY. THE POSSIBLE RECOMBINED STATES 
C                  ARE BUILT UP FROM THE TARGET (ENERGY ORDERED) SYMMS
C                  BY COUPLING-ON AN EXTRA ELECTRON. FOR A FIXED PARENT,
C                  FIXED RYD. (NV)LV AND FIXED TOTAL SYMMETRY, ONLY 1
C                  TERM OR 2 LEVELS ARE POSSIBLE. WE ASSUME THE ORDERING
C                  OF THE N+1 ELECTRON AUTOSTRUCTURE TERMS/LEVELS TO BE
C                  THE SAME AS THIS SCHEME GENERATES. IN IC, WE ASSUME 
C                  ALSO THAT THE LOWER ENERGY OF A PAIR IS JV=LV-0.5
C                  AND THAT THE HIGHER ONE CORRESPONDS TO JV=LV+0.5.
C           ****** THIS OPERATES FOR ALL RYDBERG STATES.
C                  ALSO, IF NR1 .GT.0 THEN RADIATIVE STABILIZATION
C                  OF THE VALENCE ELECTRON IS GENERATED INTERNALLY. 
C                  NR1 IS THE LOWER BOUND. (THIS IS THE SAME AS THE OLD
C                  NR1.GT.0 OPTION, SO ONLY NR1.LT.0 HAS A NEW MEANING.)
C           ****** NON-ZERO IABS(NR1) SHOULD ALWAYS BE 1 GREATER THAN 
C                  THE LARGEST CORE PRINCIPAL QUANTUM NUMBER.
C
C-----------------------------------------------------------------------
C
C     NBIN .NE. 0 THEN BINNED CROSS SECTIONS WRITTEN, WITH |NBIN| 
C                 ENERGIES BETWEEN EMIN,EMAX (BELOW) I.E. |NBIN|-1 BINS.
C          .GT. 0 NO FINAL-STATE RESOLVED DATA, NO ADF09 WRITTEN (FAST).
C          .EQ. 0 DEFAULT, DOES NOTHING.
C
C     EMIN--MINIMUM ENERGY IN UNITS FOR PROCESSING (DEFAULT ZERO)
C           N.B. ANY CROSS SECTIONS .LT. EMIN ARE OMITTED FROM ocs.
C     EMAX--MAXIMUM ENERGY IN UNITS FOR PROCESSING (DEFAULT HUGE)
C           N.B. ANY CROSS SECTIONS .GT. EMAX ARE OMITTED FROM ocs.
C
C  ***EMIN, EMAX ARE THE SAME FOR *ALL* ELECTRON TARGETS, SINCE THE
C     RELEVANT QUANTITY IS THE ELECTRON ENERGY RELATIVE TO THE INITIAL
C     STATE, *NOT* RELATIVE TO THE GROUND STATE. SO, ENSURE EMIN IS
C     SMALL ENOUGH TO INCLUDE ANY/ALL LOW-LYING RESONANCES ACCESSED 
C     FROM EXCITED ELECTRON TARGETS. I.E. EMIN=ZERO IS THE SAFE CHOICE.
CT
CT    ***IN ADDITION, EMIN,EMAX CAN RESTRICT CROSS SECTIONS TO ADF09***
CT       NOTE: THESE TESTS INVOLVING EMIN AND EMAX MAY HAVE BEEN
CT       COMMENTED-OUT BY 'CT' FOR SPEED AS THEY ARE RARELY USED.
C
C     EWIDTH .GT. 0 THEN CONVOLUTE BINNED CROSS SECTIONS WITH GAUSSIAN
C                   DISTRIBUTION OF FWHM EWIDTH.
C     EWIDTH .EQ. 0 THEN CONVOLUTE BINNED CROSS SECTIONS WITH COOLER
C                   DISTRIBUTION DEFINED BY TPAR, TPER COOLER TEMPS.
C
C     EWIDTH .LT. 0 (AND .GT.-100*UNITS) CONVOLUTE BINNED CROSS SECTIONS
C                   WITH MAXWELLIAN DISTRIBUTION. FOR CASES WHERE 
C                   PARTIAL RATE COEFFICIENTS ARE NOT REQUIRED/TOO 
C                   DEMANDING. USE NBIN0.GT.0 TO BYPASS THEM (NO ADF09).
C                   USES ADAS DEFAULT TEMPS RANGE UNLESS USER OVERRIDES.
C
C                   IF EMIN=0 (DEFAULT) THEN 
C                   ***USES A NON-LINEAR ENERGY MESH *** WITH 
C                      EBIN(1)=0 AND EBIN(2)=-EWDITH
C                      THEREAFTER LOGARITHMIC TO EBIN(NBIN)=EMAX.
C                   EMIN.GT.0 USES LINEAR (AND SO THERE SHOULD BE NO 
C                                            NEAR THRESHOLD RESONANCES.)
c
c     epart .ne. 0 uses Breit-Wigner distribution with width epart to
c                  partition dielectronic capture/autoionization rates
c           .gt. 0 over the usual (user suplied) bin energies
c           .lt. 0 over the autoionizing states present. *** TBD
C
C  ***VARIABLES HERE-ON ARE NORMALLY FOR TESTING ETC. ONLY***
C
C     ILOG .LT. 0 USE LOGARITHMIC BIN ENERGY MESH. IF EMIN=0 THEN
C                 EBIN(1)=0 AND EBIN(2)=10**ILOG.
C
C     NGAUSS .GT. 200, NGAUSS CONVOLUTION POINTS BETWEEN EMIN AND EMAX.
C            .LE. 200 2*NBIN-1 CONVOLUTION POINTS BETWEEN EMIN AND EMAX.
C
C     RCOR .NE. 0.0 CORRECTION FACTOR FOR NON-PP'D RADIATIVE RATES.
C          .LT. ZERO, AND OBSERVED ENERGIES PRESENT, FACTOR-IN CHANGE
C           TO DE**3. SET RCOR=-1.0 TO GET THIS FACTOR ALONE. DEFAULT 1.
C     ACOR .GT. 0.0 CORRECTION FACTOR FOR ALL AUGER RATES.
C     NLCOR .NE.0 CORRECTION FACTORS ACORL(LV+1) FOR LV=0,1,...|NLCOR|-1
C           .GT.0 APPLIED TO EACH AUGER RATE WITH RYDBERG A.M. LV
C           .LT.0 APPLIED TO DR CROSS SECTION WITH RYDBERG A.M. LV.
C     NNCOR .GT. 0 CORRECTION FACTORS; N, ACORN(N) FOR I=1,NNCOR
C                  APPLIED TO ENERGY-AVERAGED CROSS SECTIONS FOR N=NV.
C                  E.G. N-DEPENDENT FIELD-ENHANCEMENT FACTORS FOR DR.
C     ACOR, NLCOR, NNCOR .LE. 0.0 **DEFAULT**, RESET TO 1.0
C
C     RMIN IS THE SMALLEST RADIATIVE RATE RETAINED (DEFAULT -1, ALL)
C
C     NFNLMX.GT.0 MAXIMUM RYDBERG N FOR WHICH DETECTION PROBABILITES
C                 ARE PRESENT IN THE FILE 'FNL' FROM SCHIPPERS MODEL.
C
C     NQDT .GT. 0 READ QUANTUM DEFECTS FOR L=0...NQDT-1, TO BE USED IN
C       DETERMINING ENERGY OF RYDBERG ELECTRON. DEFAULT =0, INTERNAL.
C
C     TOLR  CONTROLS RADIATIVE STABILIZATION, NORMALLY SET INTERNALLY.
C           CASE /  / (BLSOLD) ONLY TO NON-AUTOIONIZING FINAL STATES.
C           CASE /LS/ OR /IC/ TO METASTABLE AUTOIONIZING FINAL STATES,
C                         C-R PARTIAL DATA ONLY, TOTALS STILL CORONAL.
C                         SET .EQ. 0.0 TO FORCE NON-AUTOIONIZING ONLY.
C
C     TOLI CORRECTS POSITION OF BOUND STATES BY COMPARING EXPECTED
C          POSITION OF LOWEST BOUND STATE (EIONMN-(Z/N)**2) WITH
C          ACTUAL. THIS IS TO ALLOW FOR UNBALANCED CI EXPANSION
C          THAT OCCURS WHEN DN.GT.0 PROBLEM IS SPLIT BY PARITY.
C          THE BOUND STATE CAN MOVE ABOVE ITS OWN CONTINUUM, FOR
C          SUFFICIENTLY HIGH-N.
C          ***DEFAULT: TOLI.LT.0, DETERMINED INTERNALLY. CURRENTLY,
C          BOUND STATES ARE NOT SHIFTED UP, ONLY LOWERED.
C          ALSO, NOT APPLIED WHEN FOR DN=0 DR.
C          SET TOLI.EQ.0 IN NAMELIST TWO TO SWITCH-OFF ITS OPERATION.
C          SET POSITIVE, EQUAL TO DOWNWARD SHIFT IN RYD, TO
C          OVERRIDE DEFAULT SHIFT BY USER DETERMINED ONE.
C
C     NR2 IS THE HIGHEST STABLE N FOLLOWING OUTER ELECTRON RADIATION.
C         DEFAULT IS TO DETERMINE IT INTERNALLY. 
C         CAN ONLY SET IT MANUALLY IF TOLR=0.0.
C
C     LVAMX IS MAX LV INCLUDED IN SUM OF REPRESENTATIVE AUGER RATES. 
C                                                              (999)
C
C     TOLB=MAX(1.5D-7,5.0D-9*DZ*NZ),  DEFAULT.
C         SET TOLB COARSER TO HANDLE USER SUPPLIED IMBALANCED CONTINUUM
C         EXPANSIONS, I.E. IF NOT ALL PARTIAL WAVES HAVE SAME TARGET CI.
C
C     JTEMP.EQ.0 USE ADAS DEFAULT TEMPERATURES, IF REQUIRED.
C          .GT.0 READ-IN JTEMP TEMPERATURES IN KELVIN.
C          .LT.0 READ-IN -JTEMP TEMPERATURES IN LOG(K).
C
C     IRDT.EQ.0 HISTORIC TEMP(J=1,JTEMP) READ AFTER ALL TARGET INFO ETC.
C         .NE.0 STRAIGHT AFTER NAMELIST, FOR EASE OF USE WITH SCRIPT.
C
C     JTHETA.GT.0 RESTRICT NUMBER OF ADAS DEFAULT TEMPERATURES TO JTHETA
C
C     IREL.NE.0 THEN APPLY RELATIVISTIC (JUTTNER) CORRECTION
C               TO MAXWELL RATE COEFFICIENTS.
C         .GT.0 DROPS LAST TWO ADAS TEMPS TO SYNC. WITH ADASRR.
C         .LT.0 USES FULL TEMP GRID.
C         .EQ.0 DEFAULT: DOES NOTHING. WE LEAVE JUTTNER TO ADAS TO APPLY
C
C     IOLDW=0 DEFAULT: FIXED FORMAT READ OF STAT. WEIGHTS; AND ANY CALC.
C             ENERGIES ALONGSIDE, AS PRODUCED BY NEW TERMS/LEVELS FILES.
C             ANY OBSERVED ENERGIES THEN FOLLOW DIRECTLY, NO CALCULATED.
C             WILL ALSO READ OLD STAT. WEIGHT (FILES), DETECT ABSENCE
C             OF CALCULATED ENERGIES ALONGSIDE AND THEN ANY READ OF
C             CALCULATED ENERGIES *AFTERWARDS* (THEN OBSERVED) AS OF OLD
C          =1 USE OLD FREE-FORMATTED READ OF STAT. WEIGHTS; ANY
C             CALCULATED ENERGIES ARE THEN READ AFTER. FOR USE WITH OLD
C             FILES THAT ARE MIS-READ BY FIXED FORMAT DUE TO MANUAL
C             ENTRY/SPACING OF STAT. WEIGHTS. *** SEE ALSO IRD ***
C
C     TOLBE=TOLB DEFAULT. IF USER INPUT CALCULATED ENERGIES DIFFER FROM
C                THOSE DETERMINED INTERNALLY BY MORE THAN TOLBE THEN
C                ENERGY CORRECTIONS ARE BASED-ON THE INTERNAL AND THE
C                USER IS WARNED. SET TOLBE LARGE TO FORCE USE OF
C                USER SUPPLIED CALCULATED ENERGIES - THIS WAS HISTORIC
C                OPERATION UNTIL V1.17 STRICTLY SPEAKING, SINCE ONLY THE
C                DIFFERENCE BETWEEN USER INPUT CALCULATED AND "OBSERVED"
C                MATTERS, I.E. THE ENERGY CORRECTION, THE INPUT AND
C                INTERNAL CALCULATED DO NOT HAVE TO MATCH, AS LONG AS 
C                THE "OBSERVED" ARE CONSISTENT WITH THIS. IN PRACTICE,
C                ACTUAL OBSERVED ARE USUALLY INPUT AND THE NEW DEFAULT
C                IS SAFER SINCE IT DOES NOT REQUIRE THE USER TO INPUT
C                THE CORRECT ENERGIES. HOWEVER, THERE ARE SOME HISTORIC
C                OR SPECIAL CASES WHERE ONE MIGHT WANT TO OVERRIDE IT.
C
C     IOLDE: DON'T CHANGE THIS UNLESS YOU KNOW WHAT YOU ARE DOING.
C
C***************************END-TWO*************************************
C
      TOLR=1.D-10
      TOLI=ZERO
      TOLB=-DONE
      TOLBE=-DONE
      LVAMX=999
      RCOR=DONE
      ACOR=-DONE
      RMIN=-DONE
      EMIN=ZERO
      EMAX=99999.0D0
      NLCOR=0
      NNCOR=0
      NR2=-1
      IOLDE=0
      NR1=999                 !NEW DEFAULT
      NECOR=0                 !DEFAULT
      irwt=-1       !set = 1 to recover old unweighted valence radiation
      JTHETA=0
      JTEMP=0
      IRDT=0
      IREL=0
      NRB=0                   !DEBUG SWITCH
      IOLDW=0
      NQDT=-1
      NFNLMX=0
      NBIN=0
      EWIDTH=-999
      TPAR=-DONE
      TPER=-DONE
      ESWTCHX=-DONE
      NGAUSS=0
      TC1=2.D20
      nbins=0
c
      epart=zero       !Breit-Wigner partition width
      dee=done         !Suppression factor on Auger loss
      frake=1.2        !expand target bin
      nparti=0         !number of fully resolved target weights
      w0=0             !config stat weight - should not need to set
C
C-----------------------------------------------------------------------
C
      READ(5,TWO)
C
C-----------------------------------------------------------------------
C
      IF(IFLAGX.NE.0.AND.NBIN.LE.0)THEN
      WRITE(6,*)"*** PLEASE SPECIFY A CORE EXCITATION E.G. COREX='2-3' "
      STOP 'ERROR: *** PLEASE SPECIFY A CORE EXCITATION VIA COREX'
      ENDIF
C
      if(nrb.ne.0)then
      write(0,*)'**** WARNING: OPTION NRB IS FOR EXPERIENCED USERS ONLY'
      write(6,*)'**** WARNING: OPTION NRB IS FOR EXPERIENCED USERS ONLY'
      ENDIF
      IF(NRB.GT.0)OPEN(98,FILE='IMAP',STATUS='OLD')
c
      bpart=epart.ne.zero
      if(bpart)then
        write(6,17)epart,nparti,dee
  17    format(/1x,'epart=',f9.3,3x,'nparti=',i8,3x,'dee=',1pd7.1)
        epart=epart/units
        if(nparti.lt.0)nparti=0
      endif
C
      IF(IRDT*JTEMP.NE.0)THEN            !ALTERNATE TEMP READ FOR SCRIPT
        JJTEMP=ABS(JTEMP)
        IF(JJTEMP.GT.NDIM28)THEN
          WRITE(6,*)'TOO MANY TEMPS; INCREASE NDIM28 TO:',JJTEMP
          STOP 'ERROR: TOO MANY TEMPS; INCREASE NDIM28 TO JTEMP'
        ENDIF
        READ(5,*)(TEMP(K),K=1,JJTEMP)
      ENDIF
C
      IF(NFNLMX.GT.0)THEN         ! nl-specific detection probailities
        IFNLMX=NFNLMX*(NFNLMX-1)/2
        IF(IFNLMX.GT.NDIM24)THEN
          WRITE(6,*)'*** INCREASE NDIM24 TO AT LEAST:', IFNLMX 
          STOP 'ERROR: INCREASE NDIM24'
        ENDIF
        OPEN(80,FILE='fnl')  !n,l SPECIFIC DETECTION EFFICIENCIES (READ)
        DO I=1,3
          READ(80,*)
        ENDDO
        DO IFNL=1,IFNLMX
          READ(80,*,END=802)NFNL,LFNL,FNL(IFNL)
C          WRITE(0,*) NFNL,LFNL,FNL(IFNL)
c          fnl(ifnl)=done-(done-fnl(ifnl))*2.
c          if(nfnl.gt.30)fnl(ifnl)=fnl(ifnl)/2.
c          fnl(ifnl)=max(zero,fnl(ifnl))
c          fnl(ifnl)=min(done,fnl(ifnl))
        ENDDO
        GO TO 803
 802    WRITE(6,*)'END OF *.FNL FILE REACHED TOO SOON'
        STOP 'ERROR: END OF *.FNL FILE REACHED TOO SOON'
 803    CONTINUE
      ENDIF
C
      IF(IRD.LE.0)IOLDW=1
      IF(NTAR2O.EQ.0.AND.IOLDW.EQ.1              !RE-INSTATE OLD DEFAULT
     X                        .OR.NBIN.GT.0)THEN !CASE USER RESETS NTAR2
        NTAR2=NTAR1
        IF(NTAR2O.LT.0)NTAR2=-NTAR2              !SHOULDN'T BE NECESSARY
        NBINRM=NBINM
        NBINR=NBINI
      ENDIF
      TOLB0=TOLB
      IF(BLSOLD.OR.NTAR2.LT.0)TOLR=ZERO
      IF(TOLR.LT.ZERO)TOLR=ZERO
      IF((NCMN.EQ.NCMX.OR.NECOR.NE.0)
     X   .AND.TOLI.LT.ZERO)TOLI=ZERO               !SO, DEFAULT OFF DN=0
C
      IF(MR5.LE.0)THEN
        IF(EMAX.LT.ZERO)THEN
          MR5=9
          EMAX=-EMAX
        ELSE
          MR5=5
        ENDIF
      ELSE
        IF(MR5.NE.5)MR5=9
      ENDIF
C
      WRITE(6,13)EMIN,EMAX,NR1,TOLB,TOLI
C
      EMINC=EMIN/UNITS
      EMAXC=EMAX/UNITS
      if(bpart)then
        eminc=eminc-epart*10
        emaxc=emaxc+epart*10
      endif
      IFLAGE=0
C
C SET-UP BINS: ELECTRON ENERGY RANGE DEFINED BY EMIN, EMAX
C             (RESTRICTS *ALL* OUTPUT)
C
      NBIN0=NBIN
      IF(NBIN0.NE.0)THEN
        IF(NBIN0.GT.0)THEN
          IRWT=1                                  !NO ADF09
          IF(NTAR2.LT.0.AND.NECOR.NE.0)THEN   !ONLY IF USER SETS NTAR2<0
            WRITE(0,*)'ATTENTION: AS RADIATIVE DATA MUST BE AT LEAST '
     X               ,'CONFIGURATION RESOLVED FOR NECOR.NE.0'
            WRITE(6,*)'ATTENTION: AS RADIATIVE DATA MUST BE AT LEAST '
     X               ,'CONFIGURATION RESOLVED FOR NECOR.NE.0'
          ENDIF
          NTAR2=IABS(NTAR2)                     !SWITCH-OFF ANY HYBRID
          NBINR=IABS(NBINR)
          NBINRM=IABS(NBINRM)
        ENDIF
        NBIN=ABS(NBIN0)
        NBIN1=NBIN-1
        IF(NBIN.GT.NDIM1)THEN
          WRITE(6,*)' ***DIMENSION EXCEEDED, INCREASE NDIM1 TO: ',NBIN
          STOP 'ERROR: ***DIMENSION EXCEEDED, INCREASE NDIM1'
        ENDIF
C
        BLOG=EWIDTH.LT.ZERO.AND.EMIN.EQ.ZERO.OR.ILOG.LT.0
C
        IF(.NOT.BLOG)THEN                          !LINEAR
          ERES=(EMAX-EMIN)/NBIN1
          DO N=1,NBIN
            T=N-1
            T=EMIN+T*ERES
            EBIN(N)=T/UNITS
          ENDDO
        ELSE                                       !LOG
          IF(EWIDTH.LT.ZERO)THEN
            EWID=MAX(EWIDTH,-DONE/10**(-ILOG))
          ELSE
            IF(EMIN.LT.ZERO)THEN
              EWID=EMIN
              EMIN=ZERO
              EMINC=ZERO
            ELSE
              EWID=-DONE/10**(-ILOG)
            ENDIF
          ENDIF
          EWID=MAX(EWID,-EMAX/1.E5)
          T0=MAX(-EWID,EMIN)
          T0=LOG10(T0)
          ERES=(LOG10(EMAX)-T0)/NBIN1
          DO N=2,NBIN
            T=N-1
            T=T0+T*ERES
            T=10**T
            EBIN(N)=T/UNITS
          ENDDO
          EBIN(1)=EMIN/UNITS
          ERES=-10**ERES
        ENDIF
C
        WRITE(6,14)NBIN0,ERES
C
        DO L=1,NBINM
          DO N=1,NBIN1
            SBIN(N,L)=ZERO
          ENDDO
        ENDDO
cparc                                                               !par
cpar        if(iam.eq.0)then                                        !par
C
        OPEN(7,FILE='ocs')                         !OPEN FILE
C
        EWIDTH=EWIDTH/UNITS
        IF(ABS(EWIDTH).LT.100)THEN
          IF(ILOG.LT.0)THEN
            IFLAGE=-1
          ELSE
            IFLAGE=1
          ENDIF
          IF(EWIDTH.GT.ZERO)THEN
            WRITE(6,10)EWIDTH*UNITS
          ENDIF
          IF(EWIDTH.EQ.ZERO)THEN
            WRITE(6,9)TPAR,TPER
            IF(TPAR.LE.ZERO.OR.TPER.LE.ZERO)
     X      STOP 'ERROR: ILLEGAL INPUT VALUE FOR TPAR/TPER'
            A0=TPAR/UNITS
            A0=DONE/SQRT(A0)
            B0=TPER/UNITS
            B0=DONE/SQRT(B0)
            ESWTCHX=ESWTCHX/UNITS
            IF(ESWTCHX.LT.ZERO)
     X      ESWTCHX=1.D2*(A0**2-B0**2)/B0**4  !SWITCH TO SCALED GAUSSIAN
          ENDIF
          IF(EWIDTH.LT.ZERO)THEN
            IFLAGE=-1
            WRITE(6,18)
          ENDIF
          OPEN(14,FILE='XDRTOT')
          if(nc.gt.0)then                             !read existing ocs
            lmax=nbinm
            if(ngauss.le.0)ngauss=2*nbin+1
            if(ewidth.ge.zero)go to 77
          endif
        ENDIF
c
cpar        endif                                                   !par
C
      ENDIF
C
COLD      IF(RCOR.LE.ZERO)RCOR=-DONE
      IF(ACOR.LE.ZERO)ACOR=-DONE
      IF(ABS(RCOR*ACOR).NE.DONE)WRITE(6,308)ACOR,RCOR
C
      NECOR0=NECOR
      NECOR=IABS(NECOR)
      IF(NECOR.GT.0.AND.TOLI.NE.ZERO)THEN
        WRITE(6,*)'***ERROR: NECOR AND TOLI BOTH NON ZERO'
        STOP'***ERROR: NECOR AND TOLI BOTH NON ZERO'
      ENDIF
      IF(NECOR.GT.NDIM8)THEN
        WRITE(6,888)NECOR
        STOP 'ERROR: INCREASE NDIM8'
      ENDIF
      IF(NECOR.GT.NDIM5)THEN
        WRITE(6,889)NECOR
        STOP 'ERROR: INCREASE NDIM5'
      ENDIF
C
      IF(NQDT.GT.0)THEN
        READ(5,*)(QDTS(N),N=0,NQDT-1)
        WRITE(6,172)(QDTS(N),N=0,NQDT-1)
      ENDIF
C
C READ ELECTRON TARGET INFO
C
      NBINP=MAX(NECOR,NBINM)
      IF(NTAR2.GT.0)THEN
        NBINRM=MAX(NBINRM,NECOR0)
        NBINR=NBINRM+1
        NBINP=MAX(NBINP,NBINRM)
      ENDIF
C
      LCP(1)=0
      EI(1)=ZERO
      IF(IRD.EQ.0)THEN
        READ(5,*)(IWS(I),IWL(I),I=1,NBINRM)
        WRITE(6,815)(IWS(I),IWL(I),I=1,NBINRM)
      ELSEIF(IRD.LT.0)THEN
        WRITE(6,810)
        READ(5,*)(EI(I),IWS(I),IWL(I),I=1,NBINRM),EI(NBINR)
        WRITE(6,812)(EI(I),IWS(I),IWL(I),I=1,NBINRM),EI(NBINR)
        IF(EI(1).GE.ZERO)THEN
          WRITE(6,*)'*** ERROR: FIRST TARGET BIN ENERGY MUST BE .LT. 0'
          STOP'*** ERROR: FIRST TARGET BIN ENERGY MUST BE .LT. 0'
        ENDIF
        EE1=ZERO
        DO I=2,NBINR
          IF(EI(I).GT.ZERO)THEN
            EI(I)=EI(I)+EI(1)                         !RELATIVE TO FIRST
            EE1=1
          ENDIF
        ENDDO
        EI(1)=EI(1)-EE1                     !ASSUME INPUT AS GROUND THEN
      ENDIF
      IF(IRD.LE.0)THEN
        IF(BIC.AND.IWS(1).NE.0)
     X     STOP 'ERROR:CONFUSION OVER COUPLING SCHEME...'
        IF(BLS.AND.IWS(1).EQ.0)
     X     STOP 'ERROR: CONFUSION OVER COUPLING SCHEME...'
        DO I=1,NBINRM
          E1C(I)=ZERO
          ILVTM(I)=I
          IF(BIC)THEN
            IWJ(I)=IWL(I)-1                        !HISTORIC 2J+1 INPUT
            IWL(I)=0
          ENDIF
        ENDDO
      ELSE
        E1C(2)=ZERO
        IF(BLS)THEN
          IF(IOLDW.EQ.1)WRITE(6,816)
          IF(IOLDW.EQ.0)WRITE(6,814)
          IF(MR5.NE.5)THEN
            INQUIRE(FILE='TERMS',EXIST=EX)
            IF(.NOT.EX)GO TO 822
            OPEN(MR5,FILE='TERMS')
            READ(MR5,*,END=822)
          ENDIF
          DO I=1,NBINP
            LCP(I)=0
            IF(IOLDW.EQ.1)READ(MR5,*)IWS(I),IWL(I),IPAR
            IF(IOLDW.EQ.0)READ(MR5,992)IWS(I),IWL(I),IPAR
     X                                ,LCP(I),NI,E1C(I),MERGE
            IF(IWS(I).EQ.0)THEN
              IF(NBINP+1.LT.NDIM5)THEN
                WRITE(6,*)'*** PREMATURE END OF STAT. WEIGHT INPUT:'
                WRITE(6,*)'*** REQUESTED=',NBINP,' BUT FOUND=',I-1
                STOP 'ERROR: PREMATURE END OF STAT. WEIGHT INPUT'
              ELSEIF(NTAR2.GT.0)THEN      !ASSUME INPUT UNSET (NTAR2=0)
                NBINR=I
                NBINRM=NBINR-1
                NBINP=NBINRM                
                GO TO 994                 !TERMINATOR
              ENDIF
            ENDIF
            IF(IOLDW.EQ.1)WRITE(6,817)IWS(I),IWL(I),IPAR
            IF(IOLDW.EQ.0)WRITE(6,817)IWS(I),IWL(I),IPAR,E1C(I)
          ENDDO
          IF(E1C(2).EQ.ZERO.AND.NBINP.GT.1)IOLDW=1
          IF(IOLDW.EQ.0)THEN
            E1C(NBINP+1)=ZERO
            READ(MR5,992,END=994)ITEST,IDUM,IDUM,IDUM,IDUM,ETARG1
            IF(ITEST.EQ.0)GO TO 994       !TERMINATOR
            E1C(NBINP+1)=ETARG1
            IF(MR5.EQ.5)THEN              !SKIP ANY EXTRA TARGET INFO
              DO I=1,9999
                READ(5,992,END=994)ITEST,IDUM,IDUM,IDUM,IDUM,ETARG1
                IF(ITEST.EQ.0)GO TO 994   !TERMINATOR
              ENDDO
            ENDIF
          ENDIF
        ENDIF
C
        IF(BCA)THEN
          WRITE(6,811)
          IF(MR5.NE.5)THEN
            INQUIRE(FILE='CAVES',EXIST=EX)
            IF(.NOT.EX)GO TO 822
            OPEN(MR5,FILE='CAVES')
            READ(MR5,*,END=822)
          ENDIF
          DO I=1,NBINP
            IWL(I)=0                      !TO USE BLS
            READ(MR5,*)IWS(I),IPAR,LCP(I),E1C(I) !5,* SINCE NO OLD STYLE
            IF(IWS(I).EQ.0)THEN
              IF(NBINP+1.LT.NDIM5)THEN
                WRITE(6,*)'*** PREMATURE END OF STAT. WEIGHT INPUT:'
                WRITE(6,*)'*** REQUESTED=',NBINP,' BUT FOUND=',I-1
                STOP 'ERROR: PREMATURE END OF STAT. WEIGHT INPUT'
              ELSEIF(NTAR2.GT.0)THEN      !ASSUME INPUT UNSET (NTAR2=0)
                NBINR=I
                NBINRM=NBINR-1
                NBINP=NBINRM
                GO TO 994                 !TERMINATOR
              ENDIF
            ENDIF
            WRITE(6,837)IWS(I),IPAR,E1C(I)
          ENDDO
          E1C(NBINP+1)=ZERO
          READ(MR5,*,END=994)ITEST,IDUM,IDUM,ETARG1
          IF(ITEST.EQ.0)GO TO 994         !TERMINATOR
          E1C(NBINP+1)=ETARG1
          IF(MR5.EQ.5)THEN
            DO I=1,9999                   !SKIP ANY EXTRA TARGET INFO
              READ(5,*,END=994)ITEST,IDUM,IDUM,ETARG1
              IF(ITEST.EQ.0)GO TO 994     !TERMINATOR
            ENDDO
          ENDIF
        ENDIF
C
        IF(BIC)THEN
          IF(IOLDW.EQ.1)WRITE(6,818)
          IF(IOLDW.EQ.0)WRITE(6,813)
          IF(MR5.NE.5)THEN
            INQUIRE(FILE='LEVELS',EXIST=EX)
            IF(.NOT.EX)GO TO 822
            OPEN(MR5,FILE='LEVELS')
            READ(MR5,*,END=822)
          ENDIF
          ITM=0
          DO I=1,NBINP
            LCP(I)=0
            IF(IOLDW.EQ.1)READ(MR5,*)IWJ(I),IPAR,IWS(I),IWL(I)
            IF(IOLDW.EQ.0)READ(MR5,993)IWJ(I),IPAR,IWS(I),IWL(I)
     X                                ,LCP(I),NI,E1C(I),MERGE
            IF(IWS(I).EQ.0)THEN
              IF(NBINP+1.LT.NDIM5)THEN
                WRITE(6,*)'*** PREMATURE END OF STAT. WEIGHT INPUT:'
                WRITE(6,*)'*** REQUESTED=',NBINP,' BUT FOUND=',I-1
                STOP 'ERROR: PREMATURE END OF STAT. WEIGHT INPUT'
              ELSEIF(NTAR2.GT.0)THEN      !ASSUME INPUT UNSET (NTAR2=0)
                NBINR=I
                NBINRM=NBINR-1
                NBINP=NBINRM
                GO TO 994                 !TERMINATOR
              ENDIF
            ENDIF
            DO J=1,I-1
              IF(IWS(I).EQ.IWS(J).AND.IWL(I).EQ.IWL(J))GO TO 820
            ENDDO
            ITM=ITM+1                     !NEW TERM
            ILVTM(ITM)=I                  !LOWEST LEVEL OF NEW TERM
 820        CONTINUE
          ENDDO
          DO I=ITM+1,NDIM5
            ILVTM(I)=0
          ENDDO
          DO I=1,NBINP
            IPAR=0
            IF(IWS(I).LT.0)IPAR=1
            IWS(I)=IABS(IWS(I))
            IF(IOLDW.EQ.1)WRITE(6,819)IWJ(I),IPAR,IWS(I),IWL(I)
            IF(IOLDW.EQ.0)WRITE(6,819)IWJ(I),IPAR,IWS(I),IWL(I),E1C(I)
          ENDDO
          IF(E1C(2).EQ.ZERO.AND.NBINP.GT.1)IOLDW=1
          IF(IOLDW.EQ.0)THEN
            E1C(NBINP+1)=ZERO
            READ(MR5,993,END=994)IDUM,IDUM,ITEST,IDUM,IDUM,IDUM,ETARG1
            IF(ITEST.EQ.0)GO TO 994       !TERMINATOR
            E1C(NBINP+1)=ETARG1
            IF(MR5.EQ.5)THEN              !SKIP ANY EXTRA TARGET INFO
              DO I=1,9999
                READ(5,993,END=994)IDUM,IDUM,ITEST,IDUM,IDUM,IDUM,ETARG1
                IF(ITEST.EQ.0)GO TO 994   !TERMINATOR
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C
 994  IF(MR5.NE.5)CLOSE(MR5)                               !=9
c
      if(nbinp.gt.9999)then
        write(6,*)
     x       'Target indexing needs adjusting for ntar.gt.9999'
        stop 'Target indexing needs adjusting for ntar.gt.9999'
      endif
C
C CHECK/SET TOLB
C
      IF(IOLDW.EQ.0)THEN
        TOLB=1.D10
        DO I=2,NBINP
          T=E1C(I)-E1C(I-1)
          IF(T.GT.ZERO.AND.T.LT.TOLB)TOLB=T
        ENDDO
        IF(TOLB.LT.TOLB0)THEN   !WRITE WARNING, BUT ALLOW USER SET VALUE
          WRITE(6,829)TOLB0,TOLB
          TOLB=TOLB0                              !RE-INSTATE USER VALUE
        ELSE
          IF(TOLB0.LE.ZERO.AND.TOLB.GT.ZERO)THEN
            TOLB=TOLB/2
          ELSE                             !ORIGINAL INPUT (MAYBE UNSET)
            TOLB=TOLB0
          ENDIF
        ENDIF 
      ENDIF
C
C CORRECTION FACTORS
C
      
      IF(NLCOR.NE.0)THEN
        NLCOR0=NLCOR
        NLCOR=IABS(NLCOR)
        IF(NLCOR.GT.NDIM31)THEN
          WRITE(6,*)'NLCOR REQUIRES NDIM31 AT LEAST',NLCOR
          STOP 'ERROR: NLCOR REQUIRES NDIM31 INCREASE'
        ENDIF
C
        READ(5,*)(ACORL(I),I=1,NLCOR)
        WRITE(6,177)(ACORL(I),I=1,NLCOR)
C
        DO I=NLCOR,NDIM31
          ACORL(I)=ACORL(NLCOR)
        ENDDO
        NLCOR=NLCOR0
      ELSE
        DO I=1,NDIM31
          ACORL(I)=DONE
        ENDDO
      ENDIF
C
      IF(NNCOR.GT.0)THEN
        IF(NNCOR.GT.NDIM25)THEN
          WRITE(6,*)'NNCOR REQUIRES NDIM25 AT LEAST',NNCOR
          STOP 'ERROR: NNCOR REQUIRES NDIM25 INCREASE'
        ENDIF
        DO I=1,NDIM25
          ACORN(I)=DONE
        ENDDO
        N0=999999
C
        READ(5,*)NNCOR
C
        DO I=1,NNCOR
C
          READ(5,*)N,ACORN(N)
          WRITE(6,181)N,ACORN(N)
C
          IF(N.GT.N0+1)THEN
            T=N-N0
            TT=(ACORN(N)-ACORN(N0))/T
            DO J=N0+1,N-1
              T=J-N0
              ACORN(J)=ACORN(N0)+T*TT
            ENDDO
          ENDIF
          N0=N
        ENDDO
        IF(N.LT.NDIM25)THEN
          DO I=N+1,NDIM25
            ACORN(I)=ACORN(N)
          ENDDO
        ENDIF
      ENDIF
C
C READ NECOR OBSERVED ENERGIES, E1X(I), AND CALCULATED, E1C(I), (IF NOT 
C DONE ALREADY) AND SET ENERGY CORRECTIONS ECORI(I,J) FOR I-J TARGET 
C EXCITATION ENERGY.
C (IF IOLDE=1 THEN READ ECORI CORRECTION ENERGIES DIRECTLY. THEY ARE 
C     THEN APPLIED TO CONTINUUM ATTACHED TO TARGET I - NOT RECOMMENDED.)
C
      IF(IOLDW.EQ.1)E1C(1)=ZERO              !AS UNSET
      TOLE=MAX(ZERO,E1C(1))                  !RYD
C
      IF(IOLDE.NE.1)THEN                     !NEW STYLE
C
        E1C(1)=E1C(1)*UNITS
        IF(NECOR.EQ.0)GO TO 380
        IF(IOLDW.EQ.1)THEN
          READ(5,*)(E1C(J),J=1,NECOR)
        ELSE
          DO J=2,NECOR
            E1C(J)=E1C(J)*UNITS              !AS READ RYD
          ENDDO
        ENDIF
        DO J=1,NECOR
          E1X(J)=-DONE
        ENDDO
C
        READ(5,*,END=379,ERR=379)(E1X(J),J=1,NECOR)
C
        WRITE(6,*)' '
        WRITE(6,377)(E1C(J),J=1,NECOR)
        WRITE(6,378)(E1X(J),J=1,NECOR)
        GO TO 380
 379    DO I=1,NECOR
          IF(E1X(I).LT.ZERO)THEN
            WRITE(0,*)'***MISSING OBSERVED ENERGIES: RESETTING NECOR'
            WRITE(6,*)'***MISSING OBSERVED ENERGIES: RESETTING NECOR='
     X               ,I-1
            NECOR=MAX(0,I-1)
            NECOR0=ISIGN(NECOR,NECOR0)
            GO TO 380
          ENDIF
        ENDDO
C
 380    IF(ABS(TC1).GT.1.D20)TC1=E1C(1)
        TC1=TC1/UNITS
        E1C(1)=ZERO
        TX1=E1X(1)/UNITS
        E1X(1)=ZERO
        DO I=1,NECOR
          DO J=I,NECOR
            TC=E1C(J)-E1C(I)
            TX=E1X(J)-E1X(I)
            ECORI(I,J)=(TC-TX)/UNITS
          ENDDO
        ENDDO
        DO J=2,NECOR
          TOLE=MAX(TOLE,ABS(ECORI(1,J)))
          E1C(J)=E1C(J)/UNITS
          E1X(J)=E1X(J)/UNITS       !NOT USED
        ENDDO
        TOLE=1.1D0*TOLE
        E1C(1)=TC1
        E1X(1)=TX1                  !NOT USED
C
        IF(NTAR2.LT.0.)THEN
          IF(NECOR.GT.NTAR1)THEN
            WRITE(0,*)
     X  ' *** WARNING: BUNDLED AUGER WDITHS CANNOT BE ADJUSTED BY NECOR'
            WRITE(6,*)
     X  ' *** WARNING: BUNDLED AUGER WDITHS CANNOT BE ADJUSTED BY NECOR'
          ENDIF
          IF(TC1.GT.DZERO)THEN
            WRITE(0,*)' *** STRONG WARNING: BUNDLED AUGER WDITHS'
     X               ,' CANNOT BE RAISED BY NECOR'
            WRITE(6,*)' *** STRONG WARNING: BUNDLED AUGER WDITHS'
     X               ,' CANNOT BE RAISED BY NECOR'
          ENDIF
        ENDIF
C
      ELSEIF(NECOR.GT.0)THEN                 !OLD STYLE
C
        READ(5,*)(ECORI(I,1),I=1,NECOR)
        WRITE(6,374)NECOR,(ECORI(I,1),I=1,NECOR)
C
        DO I=1,NECOR
          ECORI(I,1)=ECORI(I,1)/UNITS
        ENDDO
C
      ENDIF
C
C INITIALIZE TARGET WEIGHTS (NOW CHECKED INTERNALLY)
C
      DO I=1,NBINP
C        IF(IRD.LT.0)EI(I)=EI(I)/UNITS       !OLD BIN ENERGY IN UNITS
        IWS(I)=IABS(IWS(I))
        IF(IWS(I).EQ.0.AND.IWL(I).NE.0)THEN
          IWT(I)=IWL(I)                      !2J+1
          IWJ(I)=IWL(I)-1                    !2J
          IWL(I)=0
        ELSE
          IF(BLS.OR.BCA)THEN                 !IWS=W & IWL=0 FOR CA
            IWT(I)=IWS(I)*(2*IWL(I)+1)
          ENDIF
          IF(BIC)THEN
            IWT(I)=IWJ(I)+1                  !2J+1
          ENDIF
        ENDIF
      ENDDO
C
C SET TEMPERATURES (USER INPUT ALWAYS IN KELVIN SINCE ADAS USES KELVIN)
C
      IF(JTEMP.NE.0)THEN
        JJTEMP=ABS(JTEMP)
        IF(JJTEMP.GT.NDIM28)THEN
          WRITE(6,*)'TOO MANY TEMPS; INCREASE NDIM28 TO:',JJTEMP
          STOP 'ERROR: TOO MANY TEMPS; INCREASE NDIM28 TO JTEMP'
        ENDIF
        IF(IRDT.EQ.0)READ(5,*)(TEMP(K),K=1,JJTEMP)  !IF NOT ALREADY READ
        IF(JTEMP.LT.0)THEN
          DO K=1,JJTEMP
            TEMP(K)=10**TEMP(K)
          ENDDO
          JTEMP=-JTEMP
        ENDIF
        if(nc.gt.0)then
          do k=1,jtemp
            temp(k)=temp(k)/1.57895d5
          enddo
          go to 77
        endif
      ELSE
        if(nc.gt.0)then         !note, we must be maxwellian to get here
         write(6,*)'*** error, must specify temperatures for maxwellian'
     x,' convolution when reading previously binned cross sections...'
         stop 'must specify temperatures when nc>0'
        endif
        TEMP(1)=-DONE
        JTEMP=JTHETA                            !MAY RESTRICT ADAS TEMPS
      ENDIF
C
C
      LMAX=NBINM           !.LE.NBINM COULD REDUCE BINNED CROSS SECTIONS
C
C SUM OVER CROSS SECTIONS
C
      CALL CROSSJ(NBINM,NBINR,NMIN,LMIN,NCUT,LCUT,NECOR0,NRSLMX,ECORI
     X     ,EI,IWT,NR1,NR2,IPRINT,TOLR,ACOR,RCOR,EMINC,EMAXC,irwt,TOLI
     X     ,TOLE,LMAX,IWS,IWL,IWJ,LCP,ILVTM,LVAMX,BCA,BLSOLD,BLSNEW,BIC
     X     ,NLMAX,JTEMP,TEMP,IREL,NRB,NBIN0,EBIN,SBIN,EET,IOLDE,iflagw
     x     ,nxtrp)
C
C WRITE-OUT BINNED CROSS SECTIONS
C
      IF(NBIN0.NE.0)THEN
cparc                                                               !par
cpar        call flush(6)                                           !par
cparc                                                               !par
cpar        do l=1,lmax                                             !par
cparc                                                               !par
cpar          do n=1,nbin1                                          !par
cpar            ssend(n)=sbin(n,l)                                  !par
cpar          enddo                                                 !par
cparc                                                               !par
cpar          call comm_barrier()                                   !par
cparc                                                               !par
cpar          call mpi_reduce(ssend,srecv,nbin1,mpi_real8,mpi_sum,  !par
cpar     x                    0,mpi_comm_world,ier)                 !par
cpar          if(ier.ne.0)write(0,*)'mpi_reduce: iam, ier=',iam,ier !par
cparc                                                               !par
cpar          call comm_barrier()                                   !par
cparc                                                               !par
cpar          if(iam.eq.0)then                                      !par
cpar            do n=1,nbin1                                        !par
cpar              sbin(n,l)=srecv(n)                                !par
cpar            enddo                                               !par
cpar          endif                                                 !par
cparc                                                               !par
cpar        enddo                                                   !par
cparc                                                               !par
cpar        if(iam.ne.0)go to 800                                   !par
C
        WRITE(7,16)NBIN
        WRITE(7,15)(EBIN(N),N=1,NBIN)
        DO L=1,LMAX
          WRITE(7,15)(SBIN(N,L),N=1,NBIN1)
        ENDDO
        CLOSE(7)
      ENDIF
c
c read previous ocsp/t files, and convolute.
c
  77  if(nc.gt.0.and.iflage.ne.0)then
        read(7,16,end=78)nbin
        nbin1=nbin-1
        read(7,15)(ebin(n),n=1,nbin)
        do l=1,lmax
          read(7,15)(sbin(n,l),n=1,nbin1)
        enddo
        if(nbins.gt.0)then                !downshift
          nbin1=nbin1-nbins
          do l=1,lmax
            do n=1,nbin1
              t=ebin(n+1+nbins)/ebin(n+1)
              sbin(n,l)=sbin(n+nbins,l)*t
            enddo
          enddo
          nbin=nbin-nbins
        endif
  78    close(7)
      endif
C
C CONVOLUTE
C
      IF(IFLAGE.NE.0)THEN
        N0=1
        IF(IFLAGE.GT.0)THEN
          IF(NGAUSS.GT.200)THEN
            NT=NGAUSS
          ELSE
            NT=2*NBIN1
          ENDIF
          E0=EBIN(1)
          T=NT
          DEG=(EBIN(NBIN)-E0)/T
          NT=NT+1
        ELSE
          IF(EWIDTH.LT.ZERO)THEN
            NT=MAX(NGAUSS,100)
            E0=LOG10(TEMP(1))
            DEG=LOG10(TEMP(JTEMP))
          ELSE
            IF(NGAUSS.GT.200)THEN
              NT=NGAUSS
            ELSE
              NT=2*NBIN1
            ENDIF
            E0=LOG10(EBIN(2))
            DEG=LOG10(EBIN(NBIN))
          ENDIF
          T=NT
          DEG=(DEG-E0)/T
          NT=NT+1
C
          J0=0
          IF(BLOG)THEN     !CHECK SAFE LOW-T, LINEAR YOU'RE ON YOUR OWN!
            IF(EBIN(3)-EBIN(2).GT.EBIN(2))WRITE(6,734)
            IF(EWIDTH.LT.ZERO)THEN
              TMIN=EBIN(2)*100
              DO J=1,JTEMP
                IF(TEMP(J).LT.TMIN)J0=J
              ENDDO
            ELSE
              TMIN=EBIN(2)
            ENDIF
            DO N=1,NT
              E=E0+(N-1)*DEG
              E=10**E
              IF(E.LT.TMIN)N0=N
            ENDDO
          ENDIF
          N0=N0+1
C
          IF(EWIDTH.LT.ZERO)THEN
            J0=J0+1
            IF(J0.GT.1)THEN
              WRITE(6,733)-UNITS*TEMP(1)/100
              WRITE(0,*)'*** LOW TEMPERATURE TABULATION TRUNCATED!'
            ENDIF
            DO L=1,LMAX
              DO J=J0,JTEMP
                E=TEMP(J)
                ALF(J,L)=CONVOLM(E,EBIN,SBIN(1,L),NBIN1)     !ADAS TEMPS
              ENDDO
            ENDDO
            WRITE(6,731)
            DO J=J0,JTEMP
              WRITE(6,732)TEMP(J)*1.5789D5,(1.D-11*ALF(J,L),L=1,LMAX)
            ENDDO 
          ENDIF
C
        ENDIF
C
        DO L=1,LMAX
          WRITE(14,20)EET(L)*UNITS
          DO N=N0,NT
            E=E0+(N-1)*DEG
            IF(EWIDTH.LT.ZERO)THEN
              E=10**E
              SCC=CONVOLM(E,EBIN,SBIN(1,L),NBIN1)            !MAXWELLIAN
            ELSE
              IF(ILOG.LT.0)E=10**E
              IF(EWIDTH.GT.ZERO)THEN
                SCC=CONVOLG(E,EWIDTH,EBIN,SBIN(1,L),NBIN1)     !GAUSSIAN
              ELSE
                IF(E.LT.ESWTCHX)THEN
                  SCC=CONVOLX(E,EBIN,SBIN(1,L),NBIN1)          !COOLER
                ELSE
                  SCC=CONVOLG(E,1.D6,EBIN,SBIN(1,L),NBIN1)     !SCALED G
                  SCC=SCC*21.877D0*SQRT(E)
                ENDIF
              ENDIF
            ENDIF
            IF(SCC.LT.1.D-99)SCC=ZERO
            WRITE(14,19)E*UNITS,SCC
          ENDDO
        ENDDO
        CLOSE(14)
      ENDIF
C
C COMMENTS
C
cpar  800 continue                                                  !par
C
      NAME='NAME: ADASDR'
      CALL DATE_AND_TIME(DATE8)    
      DATE='DATE: '//DATE8(7:7)//DATE8(8:8)//'/'//DATE8(5:5)//
     X               DATE8(6:6)//'/'//DATE8(3:3)//DATE8(4:4)
C
      WRITE(6,770)
      IF(IREL.NE.0)WRITE(6,775)
      WRITE(6,1020)(COD(I),I=2,20)
      WRITE(6,790)NAME,DATE
C
      IF(NBIN0.LE.0)THEN
        WRITE(10,770)
        IF(IREL.NE.0)WRITE(10,775)
        WRITE(10,1020)(COD(I),I=2,20)
        WRITE(10,790)NAME,DATE
      ENDIF
C
      RETURN
C
  822 WRITE(6,735)
      STOP 'ERROR: TARGET SYMM INFO NOT FOUND ON FILE!!!'
C
    6 FORMAT(' J2PI=',I6)
    9 FORMAT(/' AND CONVOLUTED WITH TPAR=',F8.6,3X,'TPER=',F8.6,
     X ' COOLER DISTRIBUTION')
   10 FORMAT(/' AND CONVOLUTED WITH EWIDTH=',F6.2,
     X' FWHM GAUSSIAN DISTRIBUTION')
   11 FORMAT(/' NTAR1=',I3,3X,'NTAR2=',I3,3X,'NMIN=',I4,3X,'LMIN=',I3
     X,3X,'NMAX=',I4,3X,'LMAX=',I3,3X,A4,'=',I5)
   12 FORMAT(' LSPI=',I6)
   13 FORMAT(/1X,'EMIN=',F10.3,3X,'EMAX=',F10.3,3X,'NR1=',I3
     X,3X,'TOLB=',F12.8,3X,'TOLI=',F12.8)
   14 FORMAT(/1X,'BINNED CROSS SECTIONS WRITTEN: NBIN=',I6,
     X' WITH BIN WIDTH=',1PE12.4)
   15 FORMAT(6(1PE12.6))
   16 FORMAT(I5)
   18 FORMAT(/' AND CONVOLUTED WITH A MAXWELLIAN DISTRIBUTION')
   19 FORMAT(1PE12.4,3E14.4)
   20 FORMAT('#',1PE16.6)
  172 FORMAT(/' QUAUNTUM DEFECTS FOR L=0,1,2...:',10F8.3)
  177 FORMAT(/' ACORL',10F10.6)
  181 FORMAT(I5,F10.3)
  308 FORMAT(/1X,'ACOR=',F8.4,3X,'RCOR=',F8.4)
  374 FORMAT(/' NECOR=',I5,5X,'ECORI=',10F12.6)
  377 FORMAT(' E1 THY=',10F12.6)
  378 FORMAT(' E1 EXP=',10F12.6)
  731 FORMAT(//'    T(K) ',4X,'ALFT( 1)',2X,'ALFT( 2)',2X,'ALFT( 3)'
     X,2X,'ALFT( 4)',2X,'ALFT( 5)',2X,'ALFT( 6)'
     X,2X,'ALFT( 7)',2X,'ALFT( 8)',2X,'ALFT( 9)',2X,'ALFT(10)'
     X/4X,'----',3X,10(2X,'--------'))
  732 FORMAT(1PE10.2,1X,(10E10.2))
  733 FORMAT(/'*** LOW TEMPERATURE OUTPUT TRUNCATED: REDUCE EWIDTH TO'
     X      ,1PE9.1,' TO EXTEND TABULATION ***'/
     X'*** MAY NEED TO INCREASE NBIN AS WELL! ***')
  734 FORMAT(/'*** INCREASE NBIN ***',5X,'RESULTS MAY BE INACCURATE'/)
  735 FORMAT(/' *** ERROR: EMAX.LT.0 FLAGS READ OF TARGET SYMMS FROM'
     X,' FILES CAVES/TERMS/LEVELS, BUT NONE FOUND/EMPTY')
  770 FORMAT('C',110('-')/'C')
  775 FORMAT('C     JUTTNER RELATIVISTIC CORRECTION APPLIED TO THE',
     X ' DISTRIBUTION'/'C')
  790 FORMAT('C'/'C',1X,A30/'C',1X,A30/'C',110('-'))
  810 FORMAT(/1X,'TARGET BINS + W(I)'/)
  811 FORMAT(/1X,'    W  P        E1C(RYD)')
  812 FORMAT(8(F15.5,I3,I4))
  813 FORMAT(/1X,'2J P',3X,'(2S+1) L          E1C(RYD)')
  814 FORMAT(/1X,'(2S+1) L  P        E1C(RYD)')
  815 FORMAT(/1X,'(2S+1) L/2J+1'//15(I5,I2))
  816 FORMAT(/1X,'(2S+1) L  P')
  817 FORMAT(I6,2I3,3X,F13.6)
  818 FORMAT(/1X,'2J P',3X,'(2S+1) L')
  819 FORMAT(I3,I2,3X,I5,I3,3X,F15.8)
  829 FORMAT(/' *** WARNING: YOUR INPUT TOLB IS LARGER THAN THE',
     X' MINIMUM TARGET SPLITTING:',1P2E10.3/' *** RECOMMEND',
     X' UNSETTING TOLB AND LET CODE DETERMINE IT!'/)
  837 FORMAT(I6,I3,3X,F13.6)
  847 FORMAT(/' INCREASE NDIM2 TO AT LEAST',I4)
  849 FORMAT(/' NTAR1 MUST BE .GT. 0, BUT INPUT NTAR1=',I3)
  888 FORMAT(/' INCREASE NDIM8 TO AT LEAST',I5)
  889 FORMAT(/' INCREASE NDIM5 TO AT LEAST',I5)
  992 FORMAT(3I2,I5,I5,F18.6,3X,A4)
  993 FORMAT(2I2,1X,I3,I2,2I5,3X,F15.8,3X,A4)
 1000 FORMAT(20A4)
 1001 FORMAT(1X,50('-'),'ADASDR',50('-')//1X,20A4//1X,50('-'),'(V2.25)'
     X,50('-')//)
 1002 FORMAT(' ***INPUT CODE ERROR: ONLY /  /, /CA/, /LS/ OR /IC/ ARE'
     X  ,' ALLOWED, WHILE YOUR INPUT IS "',A4,'"')
 1003 FORMAT(' *** ERROR, COREX IMPROPERLY DEFINED:',A5)
 1004 FORMAT(' *** CORE EXCITATION N=',I2,' TO',I2/)
 1005 FORMAT(' *** CORE EXCITATION N=',I2,' TO',' *'/)
 1006 FORMAT(' *** CORE EXCITATION NL=',2I2,' TO',2I2/)
 1007 FORMAT(' *** CORE EXCITATION NL=',2I2,' TO',I2,' *'/)
 1008 FORMAT(' *** CORE EXCITATION NL=',I2,' *',' TO',I2,' *'/)
 1009 FORMAT(' *** FOR REFERENCE INITIAL CF=',I3/)
 1010 FORMAT(" *** ERROR, MUST AT LEAST SPECIFY THE N-VALUES IN NL-N'L'"
     X      ," SELECTION")
 1020 FORMAT('C',19A4)
C
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION CONVOLG(E,EWIDTH,EBIN,SBIN,NBIN1)
C
C-----------------------------------------------------------------------
C
C  CONVOLUTE BINNED CROSS SECTIONS WITH GAUSSIAN DISTRIBUTION
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DTWO=2.0D0)
C
      DIMENSION EBIN(*),SBIN(*)
      COMMON /DITT/A0,B0
C
      DATA PI/3.14159D0/,E0/9.999D50/
C
      IF(E.LT.E0)THEN
        I0=1
        N1=1
      ENDIF
C
      IF(EWIDTH.LT.1.D3)THEN                                       !FWHM
        A=1.6651092D0/EWIDTH
      ELSE                               !COOLER APPROX: SCALED GAUSSIAN
        A=A0/(SQRT(E)*DTWO)
      ENDIF
C
      EMAX=E+PI/A
      EMIN=E-PI/A
C
      SUM=ZERO
C
      DO I=I0,NBIN1
        IF(SBIN(I).GT.ZERO)THEN
          XI1=EBIN(I+1)
          IF(XI1.GT.EMIN)THEN
            XI=EBIN(I)
            IF(XI.LT.EMAX)THEN
              SUM=SUM+SBIN(I)*(ERF(A*(E-XI))-ERF(A*(E-XI1)))/DTWO
            ELSE
              GO TO 1
            ENDIF
          ELSE
            N1=I
          ENDIF
        ENDIF
      ENDDO
C
   1  I0=N1
      E0=E
C
      CONVOLG=SUM
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION CONVOLM(E,EBIN,SBIN,NBIN1)
C
C-----------------------------------------------------------------------
C
C  CONVOLUTE BINNED CROSS SECTIONS WITH  MAXWELLIAN DISTRIBUTION
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0)
C
      DIMENSION EBIN(*),SBIN(*)
C
      SUM=ZERO
C
      DO I=1,NBIN1
        I1=I+1
        T=-EBIN(I1)/E
        T=EXP(T)
        ERES=EBIN(I1)-EBIN(I)
        TT=EBIN(I1)*ERES
        SUM=SUM+T*SBIN(I)*TT
      ENDDO
C
      T=SQRT(E)
      SUM=24.6854*SUM/(T*E)
C
      CONVOLM=SUM
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION CONVOLX(E,EBIN,SBIN,NBIN1)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C  CONVOLUTE BINNED CROSS SECTIONS WITH COOLER DISTRIBUTION
C
C-----------------------------------------------------------------------
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DTWO=2.0D0)
C
      DIMENSION EBIN(*),SBIN(*)
      COMMON /DITT/A0,B0
C
      DATA PI/3.14159D0/,E0/9.999D50/
C
      IF(E.LT.E0)THEN
        I0=1
        N1=1
      ENDIF
C                                         ESTIMATE WIDTH OF DISTRIBUTION
      T=15.0D0/B0**2
      A=1.0D6
      IF(E.GT.1.0D-10)A=A0/(SQRT(E)*DTWO)
C
      EMIN=MIN(E-PI/A,E-0.8*T)
      EMAX=MAX(E+PI/A,E+T)
C
      SUM=ZERO
C
      DO I=I0,NBIN1
        IF(SBIN(I).GT.ZERO)THEN
          XI1=EBIN(I+1)
          IF(XI1.GT.EMIN)THEN
            XI=EBIN(I)
            IF(XI.LT.EMAX)THEN
              SUM=SUM+SBIN(I)*DITTNER(XI,XI1,E)
            ELSE
              GO TO 1
            ENDIF
          ELSE
            N1=I
          ENDIF
        ENDIF
      ENDDO
C
   1  I0=N1
      E0=E
C
      CONVOLX=SUM*10.938D0
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION DITTNER(EL,EU,E0)
C
C-----------------------------------------------------------------------
C
C COOLER DISTRIBUTION IS CHARACTERIZED BY TWO "TEMPERATURES", THE
C SPREAD OF THE BEAM PARALLEL AND PERPENDICULAR TO THE AXIS.
C A=1.0/SQRT(KTpar(RYD))  B=1.0/SQRT(KTperp(RYD))
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DTWO=2.0D0)
C
      PARAMETER (NMESH=2)                               !TEST 21
C
      DIMENSION E(NMESH),F(NMESH)
C
      COMMON /DITT/A0,B0
C
      NMESH1=NMESH-1
      FNMESH1=NMESH1
      E(1)=EL
      IF(EL.LT.-99.D0)E(1)=ZERO
      EH=(EU-E(1))/FNMESH1
      DO I=2,NMESH
        E(I)=E(1)+(I-1)*EH
      ENDDO
C
      A2=A0**2
      B2=B0**2
      T0=A2-B2
      C1=A0*B2/SQRT(T0)
      T=E0/T0
      C2=A2*B2*T
      Z2=A2*SQRT(T)
C
      IS=1
      IF(EL.LT.-99.D0)THEN
        IS=2
        F(1)=ZERO
      ENDIF
C
      DO I=IS,NMESH
        T2=C2-B2*E(I)
        T1=EXP(T2)
        TE=E(I)
        Z1=SQRT(TE*T0)
        ZP=Z1+Z2
        IF(ZP.GT.0.5D0)THEN
          ZM=Z2-Z1
          T2=ERFC(ZM)
          T3=ERFC(ZP)
          T1=T1*(T2-T3)
          F(I)=SQRT(TE)*T1
        ELSE
          ZM=Z1-Z2
          T2=ERF(ZM)
          T3=ERF(ZP)
          T1=T1*(T2+T3)
          F(I)=SQRT(TE)*T1
        ENDIF
      ENDDO
C
      SUM=EH*(F(1)+F(NMESH))/DTWO
C
      DO I=2,NMESH1
        SUM=SUM+EH*F(I)
      ENDDO
C
      DITTNER=C1*SUM
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE CROSSJ(NBINM,NBINR,NMN,LMN,NCUT,LCUT,NECOR0,NRSLMX
     X     ,ECORI,EI,IWT,NR1,NR2,IPRINT,TOLR,ACOR,RCOR,EMINC,EMAXC,irwt
     X     ,TOLI,TOLE,LMAX,IWS,IWL,IWJ,LCP,ILVTM,LVAMX,BCA,BLSOLD,BLSNEW
     X     ,BIC,NLMAX,JTEMP,TEMP,IREL,NRB,NBIN0,EBIN,SBIN,EET,IOLDE
     X     ,iflagw0,nxtrp)
cparc                                                               !par
cpar      use comm_interface, only : iam                            !par
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NDIM0=1000       !INITIAL METASTABLES * FINAL PARENTS
     X         ,NDIM1=10001      !BIN ENERGIES
     X         ,NDIM2=20         !INITIAL TARGET METASTABLES
C     X         ,NDIM3=100001     !MDRCS ONLY
     x         ,ndim4=100        !parent configs (for epart)
     X         ,NDIM5=400        !FINAL PARENT METASTABLES
     X         ,NDIM6=50         !GROUPS IN A MASTER CONFIG
     X         ,NDIM7=20000000    !RADIATIVE RATES PER NL           ****
     X         ,NDIM8=NDIM5
C     X         ,NDIM9=12         !UNUSED
     X         ,NDIM10=300000    !N+1 TERMS/LEVELS PER NL, INC CORR ****
C     X         ,NDIM11=99999     !MDRCS ONLY
     X         ,NDIM12=9000000   !AUTOIONIZATIONS RATES PER NL      ****
     X         ,NDIM13=200000    !N+1 ENERGIES PER NL, EXC CORR ETC ****
     X         ,NDIM14=16000      !CONFIGS (ALL)
     X         ,NDIM15=101       !HYDROGENIC POST-PROC RAD L-VALUES
C     X         ,NDIM16=50        !MDRCS ONLY
     X         ,NDIM17=30000     !TERM/LEVEL RESOLVED FINAL STATES
     X         ,NDIM18=150000    !DUMMY L-RESOLVED FINAL STATES
     X         ,NDIM19=400       !TERMS/LEVELS IN A (MASTER) GROUP
     X         ,NDIM20=92        !SEQUENCE LABEL
     X         ,NDIM24=50000     !FIELD IONIZATION NL-VALUES
     X         ,NDIM26=150       !NRSLMX NL-VALUES
     X         ,NDIM27=150       !REPRESENTAIVE N-VALUES
     X         ,NDIM28=19        !TEMPERATURES
     X         ,NDIM29=10        !OLD ADF09 STANDARD TEMPERATURES
     X         ,NDIM30=10000       !MASTER CONFIGS
     X         ,NDIM31=10        !MAX L+1 FOR BUNDLED-NL
     X         ,NDIM32=2000      !TERMS/LEVELS IN A MASTER CONFIG
     X         ,NDIM37=19        !NEW ADF09 STANDARD TEMPERATURES
     X         ,NDIM38=750       !TOT NO STATES (CFGS) STRADDLE IP
C     X         ,NDIM66=25        !MDRCS ONLY
     X         )
C
      PARAMETER (NDIM25=NDIM27)
C
      PARAMETER (NSSYM=5)              !NO. OF TOTAL SPINS MAX(INT(S)+1)
      PARAMETER (LMAXZ=31)
      PARAMETER (NLIT=60)
      PARAMETER (MXORB0=50)            !NO. NL READ FROM AS
C
      PARAMETER (TINY=1.D-4)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
      PARAMETER (DFOUR=4.0D0)
      PARAMETER (DTWELV=12.0D0)
      PARAMETER (D1M4=1.D-4)
      PARAMETER (D1M10=1.D-10)
C
      PARAMETER (CON=1.33704D-14)      !2*tau_0*(pi*a_0)**2
      PARAMETER (DKCM=109737.4D0)
C      PARAMETER (PI=3.14159265359D0)
C
      PARAMETER (DFSC=DONE/137.03599976D0)
      PARAMETER (DALF=DFSC*DFSC)
C
      REAL*4 AA,EC,AR,EATOM
      REAL*4 AAN,AANL,AANLJ,ALF,ALFN,AN,BN,BNL
C
      integer*4 MTEST4,MBLNK4              !keep I*4 for backward compat
C
      INTEGER SS,SSR,QS0,QL0,QSB,QLB,QMB,QL,QN,QSP,QLP,QSR,QLR,QNV,QLV
     X,QTI,QTT,QST,QLT,SSZ,QLZ,QSZ,QND,QLD,QTE,QTTE,QTTG,QSH,QLH
C
      CHARACTER 
     X LAB2*2,LSQ*2
     X,O*3
     X,LAB4*4,RAD*4,CMBLNK*4,CMSTAR*4
     X,LAB5*5,O1*5
     X,FILNAM*6,O1U*6
     X,F101*30
      CHARACTER*1 LAB1,CLABL(20),CLIT(0:NLIT)
C
      LOGICAL BPRNT0,BPRNT1,BPRNT2,BCFM,BCFP,BRAD,BINT,BFAST,BRFRST,BBIN
     X,BPASS1,BBNFP,BFORM,BNOT,BCA,BLS,BLSOLD,BLSNEW,BIC,BRADBF,EX,BTEST
     X,BRSLE,BRSLF,BRSLP,BRSLP1,BRSLP2,bfirst,BLOG,BCFA,BHYBRD,BUNIT
     X,BUNA,BUNR1,BUNR2,BCAH,BPRTM1,BPRTM2,BEQN,bpart,bflagp
C
      DIMENSION
     X EBIN(NDIM1),SBIN(NDIM1,NDIM2),EET(NDIM2)
C
      DIMENSION
     X IWT(NDIM5),IWS(NDIM5),IWL(NDIM5),IWJ(NDIM5),LCP(NDIM5)
     X,ILVTM(NDIM5),EI(NDIM5),ECORI(NDIM8,NDIM8)
C
      DIMENSION LIT(0:NLIT),LABL(20),LSQ(NDIM20)
C
      DIMENSION COEF(NDIM28),COFT(NDIM28),FREL(NDIM28)
     X         ,TEMP(NDIM28),THTLS(NDIM29),THTIC(NDIM37)
C
C                                                                   !F95
      ALLOCATABLE                                          !ALL     !F95
     X ITA(:),JTA(:),AA(:),EC(:)                                    !F95
     X,ITR(:),JTR(:),AR(:),EATOM(:)                                 !F95
     X,IK(:),IT(:),SS(:),LL(:),JJ(:)                                !F95
     X,JK(:),LCF(:),ITAG(:),IMAP(:)                                 !F95
     X,JV(:),ITAR(:),IAUTO(:),ILSJ(:)                               !F95
     X,JFIRST(:),JLAST(:),KFIRST(:),KLAST(:)                        !F95
     X,ENERG(:)                                                     !F95
CF77                                                                !F77
CF77      DIMENSION                                        !ALL     !F77
CF77     X ITA(NDIM12),JTA(NDIM12),AA(NDIM12),EC(NDIM12)            !F77
CF77     X,ITR(NDIM7),JTR(NDIM7),AR(NDIM7),EATOM(NDIM7)             !F77
CF77     X,IK(NDIM13),IT(NDIM13),SS(NDIM13),LL(NDIM13),JJ(NDIM13)   !F77
CF77     X,JK(NDIM10),LCF(NDIM13),ITAG(NDIM10),IMAP(NDIM10)         !F77
CF77     X,JV(NDIM13),ITAR(NDIM10),IAUTO(NDIM13),ILSJ(NDIM13)       !F77
CF77     X,JFIRST(NDIM10),JLAST(NDIM10),KFIRST(NDIM10),KLAST(NDIM10)!F77
CF77     X,ENERG(NDIM13)                                            !F77
C                                                                   !F95
      ALLOCATABLE                                          !ALL     !F95
     X TCN(:),TC(:),UB0(:,:),IBN(:)                                 !F95
     X,QN(:),QL(:),QND(:),QLD(:),MXOCC(:),NOCC1(:),NOCC(:)          !F95
CF77                                                                !F77
CF77      DIMENSION                                        !ALL     !F77
CF77     X TCN(NDIM2),TC(NDIM2),UB0(NDIM0,NDIM27),IBN(NDIM27)       !F77
CF77     X,QN(NDIM26),QL(NDIM26),QND(NDIM26),QLD(NDIM26)            !F77
CF77     X,MXOCC(NDIM26),NOCC1(NDIM26),NOCC(NDIM26)                 !F77
C                                                                   !F95
      ALLOCATABLE                                 !BINNED ONLY      !F95
     X UB(:,:,:),TNU(:,:)                                           !F95
CF77                                                                !F77
CF77      DIMENSION                               !BINNED ONLY      !F77
CF77     X UB(NDIM27,NDIM1,NDIM2),TNU(NDIM27,NDIM2)                 !F77
C                                                                   !F95
      ALLOCATABLE                                          !ALL     !F95
     X LMP(:),QSP(:,:),QLP(:,:)                                     !F95
     X,EII(:),WNP(:),ECA(:),SUMAN(:)                                !F95
CF77                                                                !F77
CF77      DIMENSION                                                 !F77
CF77     X LMP(NDIM5),QSP(NDIM5,10),QLP(NDIM5,10)                   !F77
CF77     X,EII(NDIM5),WNP(NDIM5),ECA(NDIM5),SUMAN(NDIM5)            !F77
C                                                                   !F95
      ALLOCATABLE                                          !ALL     !F95
     X ICF(:),LCA(:),NG(:),NII(:),QS0(:)                            !F95
     X,QL0(:),QSB(:,:),QLB(:,:),QMB(:,:)                            !F95
     X,LMX(:),QSH(:,:),QLH(:,:),LMH(:)                              !F95
     X,JKH(:),ITARH(:),KAUTY(:,:),JJH(:)                            !F95
     X,WNH(:)                                                       !F95
CF77                                                                !F77
CF77      DIMENSION                                                 !F77
CF77     X ICF(NDIM14),LCA(NDIM14),NG(NDIM14),NII(NDIM14),QS0(10)   !F77
CF77     X,QL0(10),QSB(NDIM14,10),QLB(NDIM14,10),QMB(NDIM14,10)     !F77
CF77     X,LMX(NDIM14),QSH(NDIM14,10),QLH(NDIM14,10),LMH(NDIM14)    !F77
CF77     X,JKH(NDIM14),ITARH(NDIM14),KAUTY(NDIM14,NSSYM),JJH(NDIM14)!F77
CF77     X,WNH(NDIM14)                                              !F77
C                                                                   !F95
      ALLOCATABLE                                          !ALL     !F95
     X CP(:),CM(:),JDUM(:)                                          !F95
     X,SUMRJ0(:),SUMRN0(:,:)                                        !F95
     X,RSUM(:),RSUMC(:),RWT(:)                                      !F95
CF77                                                                !F77
CF77      DIMENSION                                                 !F77
CF77     X CP(NDIM15),CM(NDIM15),JDUM(NDIM15)                       !F77
CF77     X,SUMRJ0(NDIM25),SUMRN0(NDIM25,NDIM5)                      !F77
CF77     X,RSUM(NDIM27),RSUMC(NDIM27),RWT(NDIM27)                   !F77
C                                                                   !F95
      ALLOCATABLE                                 !ADF09 ONLY       !F95
     X SSR(:),LLR(:),JJR(:)                                         !F95
     X,WNR(:),LMR(:),JVR(:)                                         !F95
     X,QSR(:,:),QLR(:,:)                                            !F95
     X,QNV(:),QLV(:),ITARR(:)                                       !F95
     X,AR0(:),IRSOL0(:)                                             !F95
CF77                                                                !F77
CF77      DIMENSION                               !ADF09 ONLY       !F77
CF77     X SSR(NDIM17),LLR(NDIM17),JJR(NDIM17)                      !F77
CF77     X,WNR(NDIM17),LMR(NDIM17),JVR(0:NDIM17)                    !F77
CF77     X,QSR(NDIM17,10),QLR(NDIM17,10)                            !F77
CF77     X,QNV(NDIM17),QLV(NDIM17),ITARR(NDIM17)                    !F77
CF77     X,AR0(NDIM17),IRSOL0(NDIM17)                               !F77
C                                                                   !F95
      ALLOCATABLE                                 !ADF09 ONLY       !F95
     X ERN(:),ERD(:),JWRN(:),JWRD(:)                                !F95
     X,IAUTY(:)                                                     !F95
     X,QST(:,:),QLT(:,:),LMT(:),IRSOL(:)                            !F95
CF77                                                                !F77
CF77      DIMENSION                               !ADF09 ONLY       !F77
CF77     X ERN(NDIM38),ERD(NDIM38),JWRN(NDIM38),JWRD(NDIM38)        !F77
CF77     X,IAUTY(NDIM38)                                            !F77
CF77     X,QST(NDIM30,10),QLT(NDIM30,10),LMT(NDIM30),IRSOL(NDIM13)  !F77
C                                                                   !F95
      ALLOCATABLE                                 !ADF09 ONLY       !F95
     X RPSL(:),RMSL(:),RPS(:),RMS(:)                                !F95
     X,SUMRNN(:,:),SUMRNL(:,:,:)                                    !F95
     X,SUMRLP(:),SUMRLM(:),SUMRJ(:)                                 !F95
CF77                                                                !F77
CF77      DIMENSION                               !ADF09 ONLY       !F77
CF77     X RPSL(NDIM25),RMSL(NDIM25),RPS(NDIM25),RMS(NDIM25)        !F77
CF77     X,SUMRNN(NDIM25,NDIM5),SUMRNL(NDIM31,NDIM25,NDIM5)         !F77
CF77     X,SUMRLP(NDIM25),SUMRLM(NDIM25),SUMRJ(NDIM27)              !F77
C                                                                   !F95
      ALLOCATABLE                                 !ADF09 ONLY       !F95
     X AAN(:,:),AANL(:,:,:),AANLJ(:,:)                              !F95
     X,ALF(:,:),ALFN(:,:,:)                                         !F95
     X,AN(:,:,:)                                                    !F95
     X,BNL(:,:,:,:,:)                                               !F95
     X,BN(:,:,:,:)                                                  !F95
CF77                                                                !F77
CF77      DIMENSION                               !ADF09 ONLY       !F77
CF77     X AAN(NDIM0,NDIM27),AANL(NDIM0,NDIM31,NDIM27)              !F77
CF77     X,AANLJ(NDIM13,NDIM2)                                      !F77
CF77     X,ALF(NDIM28,NDIM2),ALFN(NDIM28,NDIM27,NDIM2)              !F77
CF77     X,AN(NDIM28,NDIM17,NDIM2)                                  !F77
CF77     X,BNL(NDIM28,NDIM31,NDIM27,2,NDIM0)                        !F77
CF77     X,BN(NDIM28,NDIM27,2,NDIM0)                                !F77
C                                                                   !F95
      ALLOCATABLE              !MASTER SET-UP:     ADF09 ONLY       !F95
     X QTTE(:),QTTG(:,:),ICQTG(:,:,:)                               !F95
     X,QTI(:),QTE(:,:),ICQT(:,:)                                    !F95
     X,NGG(:),QTT(:)                                                !F95
CF77                                                                !F77
CF77      DIMENSION            !MASTER SET-UP:     ADF09 ONLY       !F77
CF77     X QTTE(NDIM13),QTTG(NDIM30,NDIM32)                         !F77
CF77     X,ICQTG(NDIM30,NDIM6,NDIM19)                               !F77
CF77     X,QTI(NDIM30),QTE(NDIM30,0:NDIM6),ICQT(NDIM30,NDIM32)      !F77
CF77     X,NGG(NDIM30),QTT(NDIM13)                                  !F77
C                                                                   !F95
      ALLOCATABLE         !PARENT SET-UP: BRSLP=ADF09 OR NECOR.GT.0 !F95
     X ITARZ(:),SSZ(:),LLZ(:),JJZ(:)                                !F95
     X,LMZ(:),QSZ(:,:),QLZ(:,:)                                     !F95
CF77                                                                !F77
CF77      DIMENSION       !PARENT SET-UP: BRSLP=ADF09 OR NECOR.GT.0 !F77
CF77     X ITARZ(NDIM18),SSZ(NDIM18),LLZ(NDIM18),JJZ(NDIM18)        !F77
CF77     X,LMZ(NDIM18),QSZ(NDIM18,10),QLZ(NDIM18,10)                !F77
c                                                                   !F95
      allocatable                      !partitioned only            !F95
     x eparti(:),iwpart(:,:)                                        !F95
     x,iepart(:,:),ice(:),jepart(:)                                 !F95
     x,icfi(:),sumadi(:),jcai(:)                                    !F95
CF77                                                                !F77
CF77      dimension                    !partitioned only            !F77
CF77     x eparti(0:ndim13),iwpart(0:ndim13,ndim4)                  !F77
CF77     x,iepart(0:ndim13,ndim2),ice(ndim4),jepart(0:ndim1)        !F77
CF77     x,icfi(ndim14),sumadi(0:ndim14),jcai(ndim12)               !F77
c
      COMMON /CORR/ACORN(NDIM25),ACORL(NDIM31),RMIN,NNCOR,NLCOR,NCMN
     X            ,NCMX,LCMN,LCMX,IMATCH,RAD,NFNLMX,FNL(NDIM24)
      COMMON /ECOR/E1C(NDIM8),E1X(NDIM8),TOLB,TOLB0,TOLBE
      COMMON /JCF/JCFA,JCFR,JCFJ,LSPI,J2PI,JPAR
      COMMON /LABEL/IMX(NDIM5,NDIM5),IREV(NDIM5,2)
      common /part/epart,dee,frake,nparti,w0(ndim4)
C
      DATA CLABL /'S','P','D','F','G','H','I','J','K','L','M','N','O'
     X,'P','Q','R','S','T','U','*'/, CMBLNK/'    '/, CMSTAR/'****'/
      DATA CLIT 
     X/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
     X 'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U',
     X 'V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k',
     X 'l','m','n','o','p','q','r','s','t','u','v','w','x','y'/
      DATA LSQ 
     X/'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE','NA','MG','AL'
     X,'SI','P ','S ','CL','AR','K ','CA','SC','TI','V ','CR','MN','FE'
     X,'CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y '
     X,'ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE'
     X,'I ','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB'
     X,'DY','HO','ER','TM','YB','LU','HF','TA',' W','RE','OS','IR','PT'
     X,'AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH','PA'
     X,'U '/       !,'NP','PU','AM','CM','BK','CF','ES','FM','MD','NO'/
      DATA
     XTHTLS/1.0D3,2.0D3,5.0D3,1.0D4,2.0D4,5.0D4,1.0D5,2.0D5,5.0D5,1.0D6/
     XTHTIC/1.0D1,2.0D1,5.0D1,1.0D2,2.0D2,5.0D2,1.0D3,2.0D3,5.0D3,1.0D4
     X     ,2.0D4,5.0D4,1.0D5,2.0D5,5.0D5,1.0D6,2.0D6,5.0D6,1.0D7/
      DATA JTHTLS/10/,JTHTIC/19/
C
      PI=ACOS(-DONE)
C
C FIX FOR FORTRAN 90 COMPILERS THAT DON'T ALLOW ASSIGNMENT OF CHARACTERS
C TO INTEGER VARIABLES, REQUIRED FOR HISTORIC BACKWARDS COMPATIBILITY
C
      OPEN(90,STATUS='SCRATCH',FORM='FORMATTED')
      WRITE(90,1111)CMSTAR,(CLIT(I),I=0,NLIT)
 1111 FORMAT(A4,80A1)
      BACKSPACE(90)
      READ(90,1111)MSTAR,(LIT(I),I=0,NLIT)
      WRITE(90,1111)CMBLNK,(CLABL(I),I=1,20)
      BACKSPACE(90)
      READ(90,1111)MBLNK,(LABL(I),I=1,20)
      BACKSPACE(90)
      READ(90,1111)MBLNK4
      CLOSE(90)
C
C THESE SWITCHES CONTROL RESOLUTION, INC REP AUGERS AND CORRESPONDING
C AUTOIONIZING FINAL STATES
C
      bpart=epart.ne.zero        !partitioned
      bflagp=.false.             !flag if master parent not "found".
C
      BHYBRD=NBINR.LT.0          !FINAL PARENT BY CONFIG
      NBINR=IABS(NBINR)
      NBINRM=NBINR-1
      NSDIM=1
      IF(BHYBRD)THEN
        NRX=0
        NBINM0=0                 !SWITCH-OFF REP AUGERS (FOR NOW)
        NBINRM0=NBINM0           !ONLY TRUE BOUND FINAL
        IF(BLSNEW)NSDIM=NSSYM
      ELSE
        NRX=NBINM
        NBINM0=NBINM             !TO ALLOW ABOVE
        NBINRM0=NBINM            !NBINRM FOR ALL PARENTS TO CONTAIN AUTO
      ENDIF                      !THEN NBINM0=NBINRM0 OR ADD LOSS TERM
      NECOR=IABS(NECOR0)
      IF(IOLDE.NE.1)NECOR=-NECOR
      BRSLE=NECOR0.LT.0               !.TRUE. FORCES PARENTAGE BY ENERGY
      BRSLF=NBIN0.LE.0.AND..NOT.BCA.AND..NOT.BHYBRD
      BRSLP=BRSLF.OR.NECOR0.GT.0.AND.NECOR.LT.0
      BBNFP=.TRUE.          !RESOLVES FINAL PARENTS IN BUNDLED-N PICTURE
      IF(BLSOLD)NLMAX=0
C
      IF(NRSLMX.GT.NDIM25)STOP 'ERROR: INCREASE NDIM25 TO NRSLMX'
      LLMAX=NLMAX-1
      IF(LCUT.LT.NCUT)THEN
        LLMAX=MIN(LCUT+1,LLMAX)                     !+1 AS LCUT IS UPPER
      ELSE
        LLMAX=MIN(9,LLMAX)
      ENDIF
      IF(LLMAX+1.GT.NDIM31)THEN
        WRITE(6,*)'WARNING: REPRESENTATIVE NL REDUCED TO L=',NDIM31-1
     X   ,', INCREASE NDIM31 TO:', LLMAX+1,', TO RETAIN ALL'
        WRITE(0,*)'WARNING: REPRESENTATIVE NL REDUCED TO L=',NDIM31-1
     X   ,', INCREASE NDIM31 TO:', LLMAX+1,', TO RETAIN ALL'
        LLMAX=NDIM31-1
      ENDIF
      MXLL=-1
C
      BLS=BLSOLD.OR.BLSNEW
      NSYM=1
      IF(BLS)THEN
        NSYM=2
      ELSE
        BLS=BCA
      ENDIF
C
C                                                                   !F95
      ALLOCATE(                                            !ALL     !F95
     X ITA(NDIM12),JTA(NDIM12),AA(NDIM12),EC(NDIM12)                !F95
     X,ITR(NDIM7),JTR(NDIM7),AR(NDIM7),EATOM(NDIM7)                 !F95
     X,IK(NDIM13),IT(NDIM13),SS(NDIM13),LL(NDIM13),JJ(NDIM13)       !F95
     X,JV(NDIM13),IAUTO(NDIM13),ILSJ(NDIM13),LCF(NDIM13)            !F95
     X,ENERG(NDIM13) !STRICTLY, ONLY JK NEED BE NDIM10, ELSE NDIM13 !F95
     X,JK(NDIM10),ITAG(NDIM10),IMAP(NDIM10),ITAR(NDIM10)            !F95
     X,JFIRST(NDIM10),JLAST(NDIM10),KFIRST(NDIM10),KLAST(NDIM10)    !F95
     X,STAT=IERR)                                                   !F95
      IF(IERR.NE.0)STOP 'ALLOCATION 1 FAILS'                        !F95
C                                                                   !F95
      ALLOCATE(                                            !ALL     !F95
     X TCN(NDIM2),TC(NDIM2),UB0(NDIM0,NDIM27),IBN(NDIM27)           !F95
     X,QN(NDIM26),QL(NDIM26),MXOCC(NDIM26),QND(NDIM26),QLD(NDIM26)  !F95
     X,NOCC1(NDIM26),NOCC(NDIM26)                                   !F95
     X,STAT=IERR)                                                   !F95
      IF(IERR.NE.0)STOP 'ALLOCATION 2 FAILS'                        !F95
C                                                                   !F95
      ALLOCATE(                                            !ALL     !F95
     X LMP(NDIM5),QSP(NDIM5,10),QLP(NDIM5,10)                       !F95
     X,EII(NDIM5),WNP(NDIM5),ECA(NDIM5),SUMAN(NDIM5)                !F95
     X,STAT=IERR)                                                   !F95
      IF(IERR.NE.0)STOP 'ALLOCATION 3 FAILS'                        !F95
C                                                                   !F95
      ALLOCATE(                                            !ALL     !F95
     X ICF(NDIM14),LCA(NDIM14),NG(NDIM14),NII(NDIM14)               !F95
     X,QS0(10),QL0(10),QSB(NDIM14,10),QLB(NDIM14,10),QMB(NDIM14,10) !F95
     X,LMX(NDIM14),QSH(NDIM14,10),QLH(NDIM14,10),LMH(NDIM14)        !F95
     X,JKH(NDIM14),ITARH(NDIM14),KAUTY(NDIM14,NSDIM),JJH(NDIM14)    !F95
     X,WNH(NDIM14)                                                  !F95
     X,STAT=IERR)                                                   !F95
      IF(IERR.NE.0)STOP 'ALLOCATION 4 FAILS'                        !F95
C                                                                   !F95
      ALLOCATE(                                            !ALL     !F95
     X CP(NDIM15),CM(NDIM15),JDUM(NDIM15)                           !F95
     X,SUMRJ0(NDIM25),SUMRN0(NDIM25,NDIM5)                          !F95
     X,RSUM(NDIM27),RSUMC(NDIM27),RWT(NDIM27)                       !F95
     X,STAT=IERR)                                                   !F95
      IF(IERR.NE.0)STOP 'ALLOCATION 5 FAILS'                        !F95
C                                                                   !F95
      IF(NBIN0.NE.0)THEN                                            !F95
        ALLOCATE(                                 !BINNED ONLY      !F95
     X   UB(NDIM27,NDIM1,NDIM2),TNU(NDIM27,NDIM2)                   !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 6 FAILS'                      !F95
      ENDIF                                                         !F95
C                                                                   !F95
      IF(NBIN0.LE.0)THEN                          !.AND.BRSLF       !F95
        ALLOCATE(                                                   !F95
     X   SSR(NDIM17),LLR(NDIM17),JJR(NDIM17)                        !F95
     X  ,WNR(NDIM17),LMR(NDIM17),JVR(0:NDIM17)                      !F95
     X  ,QSR(NDIM17,10),QLR(NDIM17,10)                              !F95
     X  ,QNV(NDIM17),QLV(NDIM17),ITARR(NDIM17)                      !F95
     X  ,AR0(NDIM17),IRSOL0(NDIM17)                                 !F95
     X  ,AN(NDIM28,NDIM17,NDIM2)                                    !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 7 FAILS'                      !F95
C                                                                   !F95
        ALLOCATE(                                                   !F95
     X   ERN(NDIM38),ERD(NDIM38),JWRN(NDIM38),JWRD(NDIM38)          !F95
     X  ,IAUTY(NDIM38)                                              !F95
     X  ,QST(NDIM30,10),QLT(NDIM30,10),LMT(NDIM30),IRSOL(NDIM13)    !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 8 FAILS'                      !F95
      ENDIF                                                         !F95
C                                                                   !F95
      IF(NBIN0.LE.0)THEN                                            !F95
        ALLOCATE(                                  !ADF09 ONLY      !F95
     X   RPSL(NDIM25),RMSL(NDIM25),RPS(NDIM25),RMS(NDIM25)          !F95
     X  ,SUMRNN(NDIM25,NDIM5),SUMRNL(NDIM31,NDIM25,NDIM5)           !F95
     X  ,SUMRLP(NDIM25),SUMRLM(NDIM25),SUMRJ(NDIM27)                !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 9 FAILS'                      !F95
C                                                                   !F95
        ALLOCATE(                                  !ADF09 ONLY      !F95
     X   AAN(NDIM0,NDIM27),AANL(NDIM0,NDIM31,NDIM27)                !F95
     X  ,AANLJ(NDIM13,NDIM2)                                        !F95
     X  ,ALF(NDIM28,NDIM2),ALFN(NDIM28,NDIM27,NDIM2)                !F95
     X  ,BNL(NDIM28,NDIM31,NDIM27,NSYM,NDIM0)                       !F95
     X  ,BN(NDIM28,NDIM27,NSYM,NDIM0)                               !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 10 FAILS'                     !F95
      ENDIF                                                         !F95
C                                                                   !F95
      IF(BRSLF)THEN                                                 !F95
        ALLOCATE(                                  !ADF09 ONLY      !F95
     X   QTTE(NDIM13),QTTG(NDIM30,NDIM32),ICQTG(NDIM30,NDIM6,NDIM19)!F95
     X  ,QTI(NDIM30),QTE(NDIM30,0:NDIM6),ICQT(NDIM30,NDIM32)        !F95
     X  ,NGG(NDIM30),QTT(NDIM13)                                    !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 11 FAILS'                     !F95
      ENDIF                                                         !F95
C                                                                   !F95
      IF(BRSLP)THEN                                                 !F95
        ALLOCATE(         !PARENT SET-UP: BRSLP=ADF09 OR NECOR.GT.0 !F95
     X   ITARZ(NDIM18),SSZ(NDIM18),LLZ(NDIM18),JJZ(NDIM18)          !F95
     X  ,LMZ(NDIM18),QSZ(NDIM18,10),QLZ(NDIM18,10)                  !F95
     X  ,STAT=IERR)                                                 !F95
        IF(IERR.NE.0)STOP 'ALLOCATION 12 FAILS'                     !F95
      ENDIF                                                         !F95
c                                                                   !F95
      if(bpart)then                                                 !F95
        allocate(                        !partitioned only          !F95
     x   eparti(0:ndim13),iwpart(0:ndim13,ndim4)                    !F95
     x  ,iepart(0:ndim13,ndim2),ice(ndim4),jepart(0:ndim1)          !F95
     x  ,icfi(ndim14),sumadi(0:ndim14),jcai(ndim12)                 !F95
     x  ,stat=ierr)                                                 !F95
        if(ierr.ne.0)stop 'allocation 13 fails'                     !F95
      endif                                                         !F95
c
      iam0=0
cpar      iam0=iam                                                  !par
c
C
C***********
C INITIALIZE
C***********
C
      BUNA=.FALSE.
      BUNR1=.FALSE.
      BUNR2=.FALSE.
      BEQN=.FALSE.
      INR1=IABS(NR1)
      BRAD=NR1.GT.0
      BRADBF=RAD.EQ.'BF'
      RABS=ABS(RCOR)
      NCMX0=-1
      IFLAGB=0
      IFLAGE=0
      IFLAGR=0
      NFLAG2=NDIM5+1
      NBIN=ABS(NBIN0)
      NBIN1=NBIN-1
      NVINT=100
      IF(NCUT.Le.999.AND.NVINT.LT.NCUT)NVINT=NCUT
      NMN0=NMN
      NUMAX=0
      NUMRX=0
      NZOLD=0
      NEOLD=0
C
      IF(BLS)THEN
        IKUN0=2          !HOW FAR TO SEARCH FOR SATISFACTORY PARENT
        IF(NCMN.GT.2)IKUN0=5
        IF(NCMN.GT.3)IKUN0=8
      ENDIF
      IF(BIC)THEN
        IKUN0=5
        IF(NCMN.GT.2)IKUN0=9
        IF(NCMN.GT.3)IKUN0=12
      ENDIF
C
C INITIALIZE NL FOR OUTER ELECTRON STABILIZATION TO STANDARD ORDER+MXORB
C
      J=MXORB0
      DO N=1,NRSLMX
        DO L=1,N
          J=J+1
          IF(J.GT.NDIM26)THEN
            NRSLMX=N-1
            WRITE(6,*)'***WARNING, NRSLMX REDUCED TO:',NRSLMX
            WRITE(6,*)'   INCREASE NDIM26 IF NEED BE'
            GO TO 424
          ENDIF
          QN(J)=N
          QL(J)=L-1
        ENDDO
      ENDDO
C
 424  IF(MXORB0.GT.NLIT)WRITE(6,*)'***WARNING: MIGHT NOT BE ABLE '
     X                  ,' TO DECODE ORBITAL, INCREASE LIT SPEC.'
C
C DUPLICATE OF POSTP BUT GUESS CATCHES NDIM2/5 SET DIFFERENT
C
      IF(NBINM.GT.NDIM2)THEN
        WRITE(6,847)NBINM
 847    FORMAT(/' INCREASE NDIM2 TO AT LEAST',I5)
        STOP 'ERROR: INCREASE NDIM2'
      ENDIF
      IF(NBINR.GT.NDIM5)THEN
        WRITE(6,848)NBINR
 848    FORMAT(/' INCREASE NDIM5 TO AT LEAST',I5)
        STOP 'ERROR: INCREASE NDIM5'
      ENDIF
C
C PACK
C
      IBNMX=NBINRM*NBINM
      IF(IBNMX.GT.NDIM0)THEN
        NBINRM=NDIM0/NBINM
        IF(NBINRM.LT.NBINM)THEN
          IBNMX=NBINM*NBINM
          WRITE(6,711)IBNMX
 711      FORMAT(' DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM0 TO'
     X           ,' AT LEAST',I4)
          STOP 'ERROR: DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM0'
        ELSE
          WRITE(6,721)NBINRM+1,IBNMX
 721      FORMAT(' REDUCING NTAR2 TO',I5,' BECAUSE OF DIMENSIONS.'
     X          /' INCREASE NDIM0 TO',I5,' IF YOU *REALLY* NEED NTAR2.')
        ENDIF
      ENDIF
C
      NBINC=IABS(NECOR)
      NBINP=MAX(NBINC,NBINM)
      IF(.NOT.BHYBRD)NBINP=MAX(NBINP,NBINRM)
        if(nbinp.gt.9999)then
          write(6,*)
     x         'Target indexing needs adjusting for ntar.gt.9999'
          stop 'Target indexing needs adjusting for ntar.gt.9999'
        endif
C
      BPRNT2=IPRINT.GE.2
      BPRNT1=IPRINT.GE.1
      BPRNT0=IPRINT.GE.0
      BPRTM1=IPRINT.GE.-1
      BPRTM2=IPRINT.GE.-2
C
C********
C ZEROISE
C********
C
      IF(TEMP(1).GT.ZERO)THEN
        DO J=1,JTEMP
          TJ=SQRT(TEMP(J))
          COEF(J)=2.0707D-16/(TJ*TEMP(J))
          TEMP(J)=TEMP(J)/1.5789D5
        ENDDO
        ITST=JTEMP/2+1
      ELSE
        IF(BLSOLD)THEN
          ITST=5
          IF(JTEMP.LE.0)THEN
            JTEMP=JTHTLS
          ELSE
            JTEMP=MIN(JTEMP,JTHTLS)
          ENDIF
        ELSE
          ITST=11
          IF(JTEMP.LE.0)THEN
            JTEMP=JTHTIC
            IF(IREL.GT.0)JTEMP=JTEMP-2
          ELSE
            JTEMP=MIN(JTEMP,JTHTIC)
          ENDIF
        ENDIF
        IF(NCMN.EQ.1)ITST=(3*ITST)/2
      ENDIF
C
C AUGERS
C
      IF(NBIN0.LE.0)THEN
        N=0
        L0=0
        J1=1
        DO L=1,NBINM0                 !INITIAL STATES
          IMX(L,L)=0
          IF(.NOT.BHYBRD)J1=L+1
          DO J=J1,NBINRM              !ALLOW SPACE FOR ALL
            IMX(J,L)=0
            N=N+1
            IF(N.GT.NDIM0)THEN        !CASE NBINM0.GT.NBINM...
              IMX(L,J)=0
              N=N-1
            ELSE
              IMX(L,J)=N
            ENDIF
          ENDDO
          DO J=1,NDIM13
            AANLJ(J,L)=ZERO
          ENDDO
C
          DO M=1,NBINRM               !FINAL PARENTS
            L0=L0+1
            DO I=1,NDIM27
              AAN(L0,I)=ZERO
            ENDDO
            IF(NLMAX.GT.0)THEN
              DO IN=1,NDIM27
                DO I=1,NDIM31
                  AANL(L0,I,IN)=ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IAAMX=N
C
C DR
C
        L0=0
        DO L=1,NBINM                  !INITIAL STATES
          DO M=1,NBINRM               !FINAL PARENTS
            L0=L0+1
            DO K=1,NSYM               !2 SPIN SYSTEMS
              DO I=1,NDIM27
                DO J=1,JTEMP
                  BN(J,I,K,L0)=ZERO
                ENDDO
              ENDDO
            ENDDO
            IF(NLMAX.GT.0)THEN
              DO NS=1,NSYM
                DO IN=1,NDIM27
                  DO K=1,NDIM31
                    DO J=1,JTEMP
                      BNL(J,K,IN,NS,L0)=ZERO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
C
          IF(BRSLF)THEN
            DO I=1,NDIM17
              DO J=1,JTEMP
                AN(J,I,L)=ZERO
              ENDDO
            ENDDO
          ENDIF
          DO I=1,NDIM27
            DO J=1,JTEMP
              ALFN(J,I,L)=ZERO
            ENDDO
          ENDDO
          DO J=1,JTEMP
            ALF(J,L)=ZERO
          ENDDO
        ENDDO
      ENDIF
C
      IF(NRB.LE.0)THEN
        DO L=1,NDIM10
          IMAP(L)=0
        ENDDO
      ENDIF
C
      IF(NBIN0.NE.0)THEN
c
        blog=ebin(1).eq.zero.and.nint(ebin(2)/(ebin(3)-ebin(2))).ne.1
        KMAX=NBINRM
        IF(NBIN0.GT.0)KMAX=1
C
        DO L=1,NBINM
          DO I=1,NDIM27
            TNU(I,L)=ZERO
          ENDDO
        ENDDO
        DO L=1,NBINM
          DO N=1,NBIN1
            DO I=1,NDIM27
              UB(I,N,L)=ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C*********************************************************************
C LOOK FOR ANY TARGET SET-UP: o_str AND/OR LEVELS (REQUIRED FOR BPART)
C*********************************************************************
C
      MR=70
      MRU=MR
C
C CHECK FOR TARGET CF DEF. IN o_str CASE OF HYBRID RESOLUTION
C (ALSO USE CF'S ONLY WHEN INCONSISTENT IN CHANNEL DEF.)
C WE NOW REQUIRE COUPLING SCHEMES TO MATCH (o_str and oca/ls/ic).
C
      NCFR=0
      BFORM=.FALSE.
      FILNAM='o_str'
      INQUIRE(FILE=FILNAM,EXIST=EX)
      IF(EX)THEN
        OPEN(MR,FILE=FILNAM)
        READ(MR,38,END=322)MDUM1,MDUM2
        BFORM=.TRUE.
        READ(MR,'(3X,A1)',END=322)LAB1
        IF(LAB1.EQ.'C')THEN
          F101='(I3,12X,I2,6X,I2,4X,50(I3,I2))'
        ELSE
          F101='(I5,10X,I2,6X,I2,4X,50(I3,I2))'
        ENDIF
        BACKSPACE(MR)
      ENDIF
 322  IF(.NOT.BFORM)THEN
        FILNAM='ou_str'
        INQUIRE(FILE=FILNAM,EXIST=EX)
        IF(EX)THEN
          OPEN(MRU,FILE=FILNAM,FORM='UNFORMATTED')
          READ(MRU,END=1000)MDUM1,MDUM2
        ELSE                                           !EXIT STAGE RIGHT
          IF(BHYBRD)THEN                               !NEEDED
            WRITE(6,*)'NO TARGET DATA ON ',FILNAM,'-REQUIRED FOR HYBRID'
            STOP 'ERROR: NO TARGET DATA ON FILE - REQUIRED FOR HYBRID'
          ELSE                     !ASSUME NOT NEEDED - NO NECOR+BUNDLED
            if(dee*epart.ne.zero)then
              write(6,*)'no target data on ',filnam
     x        ,'- required for damped partitioned'
              stop 
     x 'error: no target data on file - required for damped partitioned'
            else
              GO TO 200            !OR FOR CASE OF BAD LABEL (HOPEFULLY)
            endif
          ENDIF
        ENDIF
      ENDIF
      IF(BHYBRD.AND.MDUM2.GE.0)THEN                      !SUSPICIOUS...
        WRITE(6,*)'TARGET DATA HAS LV.GE.0???'
        STOP 'ERROR: TARGET DATA HAS LV.GE.0 - INVALID FOR HYBRID'
      ENDIF
C
      IF(BFORM)READ(MR,F101)NCFR,NZOLD,NEOLD             !CF HEADER
      IF(.NOT.BFORM)READ(MRU)NCFR,NZOLD,NEOLD
      NEOLD=NEOLD+1
      IF(NCFR.GT.NDIM14)THEN
        WRITE(6,*)'*** INCREASE NDIM14 CONFIGS FOR HYBRID TO',NCFR
        STOP 'ERROR: INCREASE NDIM14'
      ENDIF
      IF(BHYBRD)THEN
        IF(NCFR.GT.NBINRM)THEN
          WRITE(6,*)'*** WARNING: YOUR TARGET ',FILNAM,' FILE HAS',NCFR,
     X    ' CONFIGS BUT YOU HAVE ONLY REQUESTED',NBINRM,' FINAL PARENTS'
        ELSEIF(NCFR.LT.NBINRM)THEN!ASSUME TARGET COORECT,QUIETLY REDUCE
          NBINRM=NCFR
          NBINR=NBINRM+1
        ENDIF
        if(ncfr.gt.9999)then
          write(6,*)
     x         'Hybrid target indexing needs adjusting for ncfr.gt.9999'
          stop 'Hybrid target indexing needs adjusting for ncfr.gt.9999'
        endif
        IF(LCP(1).EQ.0)THEN
          WRITE(6,*)'*** WARNING: MAY NEED TARGET CF NOS IN adasin FOR'
     X             ,' HYBRID OPERATION'
C          STOP'*** ERROR: NEED TARGET CF NOS IN adasin FOR HYBRID MODE'
        ENDIF
      ENDIF
C
      DO I=1,NCFR                                        !CFGS
        IF(BFORM)READ(MR,179,END=1002)IN,NGR,MA0,MB0
     X                              ,(QS0(J),QL0(J),J=1,10)
        IF(.NOT.BFORM)READ(MRU,END=1002)IN,NGR,MA0,MB0
     X                              ,(QS0(J),QL0(J),J=1,10)
        IF(IN.LT.0)THEN
          IF(BHYBRD)THEN
            STOP 'ERROR INVALID TARGET FILE - REQUIRED FOR HYBRID'
          ELSE
            NCFR=0                                !HOPE WE DON'T NEED IT
            GO TO 200
          ENDIF
        ELSE
          JJH(I)=NGR
          WNH(I)=ZERO                             !INITIALIZE
          jkh(i)=0                                !checksum
        ENDIF
        DO 166 J=1,10
          QSH(I,J)=MBLNK
          IF(QL0(J).EQ.MBLNK)GO TO 166
          LMH(I)=J
          M=MOD(QS0(J),50)
          IF(M.GT.0)QSH(I,J)=LIT(M)
          DO K=1,NLIT
            IF(QL0(J).EQ.LIT(K))GO TO 199
          ENDDO
          QLH(I,J)=0
          GO TO 166
 199      QLH(I,J)=K
 166    ENDDO
      ENDDO
C
      IF(BHYBRD)THEN
        IF(BFORM)THEN                          !SKIP AA HEADERS
          DO I=1,3
            READ(MR,103,END=1002)
          ENDDO
        ELSEIF(.NOT.BFORM)THEN                 !THERE SHOULD BE NO RATES
          READ(MRU,END=1002)
          READ(MRU,END=1002)
        ENDIF
        IF(BFORM)READ(MR,121,END=1002)NENG,ECORE
        IF(.NOT.BFORM)READ(MRU,END=1002)NENG,ECORE
C
        IF(BFORM)READ(MR,106,END=1002)MTEST4
        IF(.NOT.BFORM)READ(MRU,END=1002)MTEST4
 106    FORMAT(21X,A4)
        IF(BFORM)BACKSPACE(MR)
        IF(.NOT.BFORM)BACKSPACE(MRU)
C
        BCAH=MTEST4.EQ.MBLNK4
        if(bcah.and.neng.ne.ncfr)stop 'confusion, is o_str CA or not?'
        IF(BCAH.NEQV.BCA)THEN
          WRITE(6,375)BCAH,BCA
 375    FORMAT('*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY BETWEEN'
     X        ,' o_str AND o1 ETC FILES; ONE IS CA THE OTHER NOT:',2L3)
          STOP '*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY'
        ENDIF
C
        IF(BFORM)READ(MR,105,END=1002)MTEST4
        IF(.NOT.BFORM)READ(MRU,END=1002)MTEST4
C
        BTEST=MTEST4.NE.MBLNK4
        IF(BTEST.NEQV.BIC)THEN
          WRITE(6,376)BTEST,BIC
 376    FORMAT('*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY BETWEEN'
     X        ,' o_str AND o1 ETC FILES; ONE IS IC THE OTHER NOT:',2L3)
          STOP '*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY'
        ENDIF
        IF(.NOT.BTEST.NEQV.BLS)THEN
          BTEST=.NOT.BTEST
          WRITE(6,377)BTEST,BLS
 377    FORMAT('*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY BETWEEN'
     X        ,' o_str AND o1 ETC FILES; ONE IS LS THE OTHER NOT:',2L3)
          STOP '*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY'
        ENDIF
C
        DO I=1,NENG                            !GET ENERGIES
          IF(BFORM)READ(MR,123,END=1002)MDUM,MT,IS,IL,IJ,IC,E
          IF(.NOT.BFORM)READ(MRU,END=1002)MDUM,MT,IS,IL,IJ,IC,E
          IS=IABS(IS)                       !CASE CORR.
          IF(BCAH)THEN
            IF(BFORM)THEN
              IS=100000*MT+IS               !AS 2I5
              IS=IS*(2*IL+1)                !UNLIKELY CASE 1 TERM PER CF
c              if(is.ne.jjh(ic))then
c               write(0,*)'i=',i,' is=',is,' ic=',ic,' jjh(ic)=',jjh(ic)
c               stop 'is.ne.jjh test'
c              endif
              JJH(IC)=IS                    !AS ONLY I5 IN CF READ
            ENDIF
            JKH(I)=IC                       !EO->SO
            WNH(IC)=-E-ECORE
          ELSE
            IF(BTEST)THEN
              JW=IJ+1
            ELSE
              JW=IS*(2*IL+1)
            ENDIF
            WNH(IC)=WNH(IC)+JW*E
            jkh(ic)=jkh(ic)+jw              !checksum
          ENDIF
        ENDDO
        IF(.NOT.BCAH)THEN                   !NEED EO INDEX
          DO I=1,NCFR
            WNH(I)=-WNH(I)/jkh(i)-ECORE
c            if(jkh(i).ne.JJH(I))then        !only if no corr.
c              write(6,*)'stat weight checksum failure:',i,jkh(i),jjh(i)
c              stop 'stat weight checksum failure'
c            endif
            JJH(I)=jkh(i)                   !CASE CORR.
          ENDDO
          CALL HPSRTI(NCFR,WNH,JKH)
        ENDIF
        I1=NDIM14
        IF(BLSNEW)I1=I1*NSSYM
        IF(I1.GT.NDIM17)THEN
         WRITE(6,*)'CA/LS HYBRID INDEX ARRAY TOO SMALL, ',
     X             'INCREASE NDIM17 TO',I1
         STOP 'ERROR: CA/LS HYBRID INDEX ARRAY TOO SMALL'
        ENDIF
        DO I=0,I1                           !INITIALIZE CONFIG FLAG
          JVR(I)=0
        ENDDO
        NAUTY=0
      ENDIF
C
 200  IF(EX)CLOSE(MR)                       !MRU=MR
C
C CHECK ANY SET-UP FOR BERIT-WIGNER REDISTRIBUTION: TARGET PARTITIONING
C
      if(bpart)then
c
      tolbi=tolb0
      if(tolbi.lt.zero)tolbi=max(1.d-5,-1.d-7*nzold*(nzold-neold)**2)
c
      ice=0                      !(i)=0
      iepart=0                   !(i)=0
      iwpart=0                   !(0,i)=0
      eparti(0)=-tolbi
      icf0=0
c
c check for any partition of target 
c     (use levels file only as ic always tractable for target)
c
      inquire (file='LEVELS',exist=ex)
c
      if(ex.and.nparti*epart.ne.zero)then
        open(9,file='LEVELS')
        read(9,*)                                    !header
c
        i0=0
        do i=1,nparti
          read(9,993,end=306)iw,idum,idum,idum,icfd,idum,e
          if(icfd.eq.0)go to 306                      !terminator
          if(icfd.gt.ndim4)stop 'increase ndim4 dimension'
          if(ice(icfd).eq.0)then
            icf0=icf0+1
            ice(icfd)=icf0                             !maps so->eo
            if(w0(icf0).eq.zero)w0(icf0)=jjh(icfd)     !icf0->icfd if so
          endif
          do ic=1,icf0
            iwpart(i,ic)=iwpart(i0,ic)        !ic=ice(ic) if so<-eo
          enddo
          icfd=ice(icfd)                               !eo<-so
          iwpart(i,icfd)=iwpart(i0,icfd)+iw+1
          eparti(i)=e+tolbi
          i0=i
          if(bprnt0)write(6,157)i,icfd,eparti(i),iwpart(i,icfd)
        enddo
c
        do i=i0+1,ndim13
          read(9,993,end=305)iw,idum,idum,idum,icfd,idum,e
          if(icfd.eq.0)go to 305                       !terminator
          if(icfd.gt.ndim4)stop 'increase ndim4 dimension'
          if(ice(icfd).eq.0)then
            icf0=icf0+1
            ice(icfd)=icf0                             !maps so->eo
            if(w0(icf0).eq.zero)w0(icf0)=jjh(icfd)     !icf0->icfd if so
          endif
          icfd=ice(icfd)                               !so<-eo
          if(e.gt.eparti(i0))then
            ip=i0+1
            eparti(ip)=eparti(i0)*frake
            do ic=1,icf0
              iwpart(ip,ic)=iwpart(i0,ic)        !ic=ice(ic) if so<-eo
            enddo
            iwpart(ip,icfd)=iwpart(ip,icfd)+iw+1
            if(bprnt1)write(6,157)i-1,i0,eparti(i0)
     x                          ,(iwpart(i0,ic),ic=1,icf0)
            i0=ip
          else
            iwpart(i0,icfd)=iwpart(i0,icfd)+iw+1
          endif
        enddo
c
 305    continue
        if(bprnt1)write(6,157)i-1,i0,eparti(i0)
     x                      ,(iwpart(i0,ic),ic=1,icf0)
c
 306    continue
        nparti=i0
c
        if(epart.gt.zero)then                    !map bin e onto part e
          if(nbin.gt.ndim13)then
            write(6,*) 'too many bin energies for partition,'
     x                ,' increase ndim13 to',nbin
            stop'too many bin energies for partition, increase ndim13'
          endif
          do l=1,lmax
            i0=0
            do i=1,nbin
              e=ebin(i)+e1c(l)
  55          if(e.gt.eparti(i0))then
                i0=i0+1
                if(i0.le.nparti)go to 55
              endif
              i0=i0-1
              iepart(i,l)=i0
            enddo
          enddo
          if(bprnt2)then
            do i=1,nbin
              write(6,158)i,iepart(i,l),eparti(i0)  !-wnp(ic) don't know
     x                  ,(iwpart(i0,ic),ebin(i)-zero,ic=1,icf0)
            enddo
          endif
        endif
      endif
c
      if(ex)close(9)
c
      endif
C
C*****************************************
C INITIALIZE REGULAR UNITS FOR DR (NL, LV)
C*****************************************
C
      BFORM=.FALSE.
      IF(BLS)O='ols'
      IF(BCA)O='oca'
      IF(BIC)O='oic'
      IFILE=1
cpar      ifile=ifile+iam                                           !par
      IC1=IFILE/10
      IC2=IFILE-10*IC1
      IC0=ICHAR('0')
      IC1=IC1+IC0
      IC2=IC2+IC0
      IF(IFILE.LT.10)THEN
        O1='o'//CHAR(IC2)
        O1U='o'//CHAR(IC2)//'u'
      ELSE
        O1='o'//CHAR(IC1)//CHAR(IC2)
        O1U='o'//CHAR(IC1)//CHAR(IC2)//'u'
      ENDIF
C
 330  FILNAM=O1
      INQUIRE(FILE=FILNAM,EXIST=EX)
      IF(EX)THEN
        OPEN(MR,FILE=FILNAM)
        READ(MR,38,END=332)MDUM1,MDUM2
        BFORM=.TRUE.
        READ(MR,'(3X,A1)',END=332)LAB1
        IF(LAB1.EQ.'C')THEN
          F101='(I3,12X,I2,6X,I2,4X,50(I3,I2))'
        ELSE
          F101='(I5,10X,I2,6X,I2,4X,50(I3,I2))'
        ENDIF
        BACKSPACE(MR)
        BACKSPACE(MR)
      ENDIF
 332  IF(.NOT.BFORM)THEN
        FILNAM=O1U
        INQUIRE(FILE=FILNAM,EXIST=EX)
        IF(EX)THEN
          OPEN(MRU,FILE=FILNAM,FORM='UNFORMATTED')
        ELSE
          IF(IFILE.GT.0)THEN                   !CHECK FOR PARALLEL FILES
            IFILE=0
cpar            ifile=ifile+iam                                     !par
            IC1=IFILE/10
            IC2=IFILE-10*IC1
            IC0=ICHAR('0')
            IC1=IC1+IC0
            IC2=IC2+IC0
            O1=O//CHAR(IC1)//CHAR(IC2)
            O1U=O//'u'//CHAR(IC1)//CHAR(IC2)
cpar            ifile=-ifile                                        !par
            GO TO 330
          ELSE                                 !EXIT STAGE RIGHT
            WRITE(6,*)'NO RATE INPUT DATA ON FILE o1 OR o1u!!!'
            STOP 'ERROR: NO RATE INPUT DATA ON FILE o1 OR o1u!!!'
          ENDIF
        ENDIF
      ENDIF
C
cpar      ifile=-1               !flag one file per proc            !par
C
C***********
C INITIALIZE
C***********
C
      BBIN=EI(1).GE.ZERO
      BPASS1=.TRUE.
      IF(JCFJ.GT.0.and.iam0.eq.0)THEN
        LVH=999
      ELSE
        LVH=-999
        JCFJ=999999
      ENDIF
      bfirst=.true.
      NCFT=0
      NENG=0
      TOLN=ZERO
      TOLI0=TOLI
      TOLI00=TOLI
      IF(TOLI.LT.ZERO)TOLI=ZERO
      INIT=0
C
      IB0=MAX(NRSLMX,NMN0)
      IF(IB0.GT.NDIM27)THEN
        WRITE(6,399)
        STOP 'ERROR: SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27'
      ENDIF
      IB00=IB0
      DO I=1,IB0
        IBN(I)=I
        RWT(I)=DONE
      ENDDO
      JCFRB=IABS(JCFR)
C
      NRSOL=0
      LV00=-1
C
C ENTRY POINT FOR NEW UNIT
C
 331  NV00=0
      NV0=100000
      LV0=-1
      DO I=1,NDIM26
        MXOCC(I)=0
      ENDDO
      IF(TOLI00.LT.ZERO)TOLI0=TOLI00
      MXORB=min(30,MXORB0)
      MXORBR=0
      BEQN=BHYBRD.AND.NCMX.GT.NCMN
C
 310  NV=0
      IF(BFORM)READ(MR,38,END=1000)NV,LV
      IF(.NOT.BFORM)READ(MRU,END=1000)NV,LV
  38  FORMAT(5X,I5,5X,I5)
C
      IF(NV.EQ.0.AND.NV0.EQ.100000)THEN
        NV=INR1
        LV=-1
      ENDIF
      IF(LV.GT.LMAXZ-1.AND.LV.NE.999)THEN
        WRITE(6,*)'INCREASE PARAMETER LMAXZ TO: ',LV+1
        STOP 'ERROR: INCREASE PARAMETER LMAXZ'
      ENDIF
      IF(LV.LT.0.AND.NR1.LT.0)INIT=INR1
      IF(NV.GT.0.AND.INR1.NE.999.AND.LV.GE.0)THEN
        IF(NV.LT.INR1)THEN
          WRITE(6,39)NV,NR1
  39      FORMAT(' ERROR IN CROSSJ: NV MUST BE .GE. ABS(NR1):',2I6)
          STOP ' ERROR IN CROSSJ: NV MUST BE .GE. ABS(NR1):'
        ENDIF
        IF(NV0.EQ.100000.AND.NR1.GT.0.AND.NV.GT.NR1.AND.NBIN0.LE.0)THEN
          IF(BLS)WRITE(6,388)NR1,NV-1,NR1
 388      FORMAT(' NOTE: TO OBTAIN NON-HYDROGENIC ENERGY LEVELS FOR N='
     X    ,I2,' TO',I2,' RE-RUN AUTOSTRUCTURE FROM NMIN=',I2)
          IF(BIC)WRITE(6,389)NR1,NV-1,NR1
 389      FORMAT(' NOTE: TO OBTAIN TERM LABELS AND NON-HYDROGENIC'
     X    ,' ENERGY LEVELS FOR N='
     X    ,I2,' TO',I2,' RE-RUN AUTOSTRUCTURE FROM NMIN=',I2)
        ENDIF
      ENDIF
      IF(LV.LT.0.AND.LV00.GE.0)THEN
        WRITE(6,*)'***ERROR: RE-ORDER INPUT FILES on(u) ETC SO THAT'
     X           ,'   EQUIVALENT ELECTRON FILES COME FIRST***'
        STOP '***ERROR: RE-ORDER INPUT FILES on(u)'
      ENDIF
C
      BNOT=.NOT.BRAD.OR.LV.LT.0.OR.LV.EQ.999
      TV=NV*NV
      TV3=TV*NV
      BCFA=IABS(JCFA).NE.999999     !.AND.LV.GE.0
      BCFM=JCFR.LT.0                !.AND.LV.GE.0
      BCFP=JCFR.GT.0.AND.JCFR.LE.100.AND.LV.GE.0
      IF(LV0.LT.0.OR.LV.LT.0)THEN
        BUNIT=.TRUE.                               !not used...
        NRSOLZ=0
        WNP0=-DONE
        KFPM=0
        IF(INIT.EQ.0)EIONMN=ZERO
        IF(BBIN)EI(1)=DONE
      ENDIF
      IF(NV00.EQ.0)NV00=NV
      IF(NV0.EQ.100000)GO TO 70
      EX=.FALSE.
      BPASS1=.FALSE.
      IF(NV.GT.0.AND.LV.EQ.LV0)GO TO 37
      IF(TOLI0.LT.ZERO.AND.LV.GE.0.AND.LV0.GE.0)TOLI0=-TOLI0
C
  91  IF(NV.EQ.0)GO TO 1000
  70  IF(LV.GT.LCUT.AND.LV.NE.999)GO TO 1000
C
C**************
C START A NEW L
C**************
C
      LV0=LV
      NV0=NV-1
      NV00=NV
      LV00=LV0-1
      NMN=MAX(NV,NMN0)
      IF(LV.GE.0)THEN
        MXLL=MAX(MXLL,LV+1)                  !+1 AS LV IS UPPER
      ELSE
        MXLL=MAX(MXLL,NV-1)
      ENDIF
      nvx=-1
C
C FLAG RESOLVED PARENTS NEEDED
      BRSLP1=LV.GE.0.AND.NR1.NE.0.AND..NOT.BRSLE
      BRSLP2=LV.LT.0.AND.NR1.LT.0
      BEQN=BEQN.AND.LV.GE.0           !START RYD RAD AT UPPER-N
C
      IF(NBIN0.LE.0)THEN
        DO M=1,NDIM10         !INITIALIZE FOR T OMITTED BY ANY AS BUNDLE
          JK(M)=0
        ENDDO
      ENDIF
C
      IF(NRB.GT.0)THEN
        DO L=1,NDIM13
          IMAP(L)=0
        ENDDO
        READ(98,*,END=37)NMAP,LMAP,MMAP
        IF((-1)**MMAP.LT.0)THEN
          WRITE(6,*)
     X   'EXPECTING MMAP EVEN TO RELABEL SWAPPED PAIRS OF LEVELS'
          STOP 'INPUT ERROR IN IMAP FILE?'
        ENDIF
        IF(LMAP.EQ.LV)THEN
          DO N=1,MMAP
            READ(98,*)I,IMAP(I)
          ENDDO
        ELSEIF(LMAP.GT.LV)THEN
          BACKSPACE(98)
        ENDIF
      ENDIF
C
C**************
C START A NEW N
C**************
C
  37  DO L=1,NBINM
        TC(L)=ZERO
      ENDDO
      BINT=.FALSE.
      IF(NV.LE.NCUT.AND.NV.GE.NMN.AND.LV.GE.LMN)GO TO 85
      BINT=.TRUE.
      IF(LV.LT.LCUT.OR.NV.LT.NMN)GO TO 75
      LV=LV+1
      GO TO 91
  85  CONTINUE
      IF(NBIN0.LE.0)THEN
        DO I=1,NDIM25
          RPS(I)=ZERO
          RMS(I)=ZERO
          RPSL(I)=ZERO
          RMSL(I)=ZERO
        ENDDO
      ENDIF
      DO I=1,NDIM27
        RSUM(I)=ZERO
        RSUMC(I)=ZERO
      ENDDO
      DO I=1,IB0
        IF(NV.EQ.IBN(I))GO TO 60
      ENDDO
      IF(NV.LT.IBN(IB00))THEN
        WRITE(6,*)'NV MUST HAVE NO GAPS UP TO NRSLMX=',IBN(IB00)
        STOP 'ERROR: NV MUST HAVE NO GAPS UP TO NRSLMX'
      ENDIF
      IF(NV.LT.IBN(IB0))THEN
        WRITE(6,391)NV
 391    FORMAT('NV=',I5,' NOT FOUND - LIKELY CAUSE LV .GE. LAST ',
     X  'SEQUENTIAL N, REDUCE LMAX OR RERUN AUTOSTRUCTURE')
        STOP 'ERROR: NV NOT FOUND'
      ENDIF
      IB0=IB0+1
      I=IB0
      IF(I.GT.NDIM27)THEN
        WRITE(6,399)
 399    FORMAT(' SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27')
        STOP 'ERROR: SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27'
      ENDIF
      IBN(I)=NV
      RWT(I)=EMN3(IBN(I-1))-EMN3(IBN(I))     !RADIATIVE WEIGHTING FACTOR
      RWT(I)=RWT(I)*TV3
      if(irwt.gt.0)rwt(i)=DONE
  60  IB=I
  75  NV0=NV
c
c just scale all rydberg data
c
      if(nv.gt.nxtrp)then
c
        if(nvx.le.0)then                     !we have nothing to extrap.
          write(6,*)'***extrapolation error: set nxtrp .ge.',nv
          stop '***extrapolation error: increase nxtrp'
        endif
c
        do i=1,mxorb
          if(qn(i).eq.nvx)qn(i)=nv
        enddo
c
        t1=nvx
        t2=nv
        tx0=(t1/t2)**3
c
        if(bic)then !use kappa av. as we don't know which level is which
c          kappa=lv-1/(lv+1)
c          kappa=-lv-1
          kappa=-1
        else                                 !non-rel energies
          kappa=0
        endif
c
        tex=qdt(qd0,nz0,ne,nv,lv,kappa)
        tex=tex-qdt(qd0,nz0,ne,nvx,lv,kappa)
c        tex=dz/(nvx*nvx)-dz/(nv*nv)
c
c        write(0,*)nv,tex
c                                                        !autoionization
        tx=tx0
        do i=1,numa
          mn=jta(i)
          if(mn.ne.0)then
            if(mn.lt.0)then           !cf no. for hybrid case- TBD in AS
              jca=-mn
            else
              jca=-lcf(jk(mn))
            endif
            mn=qlb(jca,lmx(jca))
            if(qn(mn).ne.nv)then
c              if(nv.eq.nxtrp+1)write(6,*)i
              aa(i)=aa(i)*tx
              ec(i)=ec(i)+tex
            endif
          else              !can't tell, assume no core rearrangement...
            aa(i)=aa(i)*tx
            ec(i)=ec(i)+tex
          endif
        enddo
c                                                              !energies
        do i=1,neng
          ica=lcf(i)
          if(ica.gt.0)then
            mn=qlb(ica,lmx(ica))
            if(qn(mn).eq.nv)then
c              if(nv.eq.nxtrp+1)write(6,*)i
              energ(i)=energ(i)+tex
            endif
          endif
        enddo
c                                                             !radiation
        tx=tx0
        te3=1
        do i=1,numr
          mn=jtr(i)
          if(mn.ne.0)then
            if(mn.lt.0)then                  !cf no. for hybrid case
              jcr=-mn
            else
              j0=jk(mn)
              jcr=lcf(j0)
              i0=jk(itr(i))
              icr=lcf(i0)
              mn=qlb(icr,lmx(icr))
              if(qn(mn).eq.nv)then           !not strictly necess.
                tx=tx0
                te3=energ(i0)-energ(j0)
                te3=te3/(te3-tex)
                te3=te3*te3*te3
              else                           !since ar(i) not used
                tx=1
                te3=1
              endif
            endif
            mn=qlb(jcr,lmx(jcr))
            if(qn(mn).ne.nv)then
              ar(i)=ar(i)*tx*te3
c              if(nv.eq.nxtrp+1)write(6,*)i,te3
            else                             !core radiation
              eatom(i)=eatom(i)+tex
c              if(nv.eq.nxtrp+1)write(6,*)-i
            endif
          endif
        enddo
c
c now "skip" reads
c                                                              !orbitals
        if(bform)read(mr,*,end=1002)
        if(.not.bform)read(mru,end=1002)
c                                                               !configs
        do i=1,ncf
          if(bform)read(mr,179,end=1002)
          if(.not.bform)read(mru,end=1002)
        enddo
c                                                        !autoionization
        if(bform)read(mr,103,end=1002)
        if(.not.bform)read(mru,end=1002)
        if(bform)read(mr,103,end=1002)
c
 9111   if(bform)read(mr,112,end=1002)i1,i2
        if(.not.bform)read(mru,end=1002)i1,i2
        if(i2.eq.0)go to 9113
        go to 9111
c                                                              !energies
 9113   if(bform)read(mr,121,end=1002)
        if(.not.bform)read(mru,end=1002)
        if(bform)read(mr,105,end=1002)
        if(.not.bform)read(mru,end=1002)
c
        do i=1,neng
          if(bform)read(mr,123,end=1002)    !i0,id,id,id,id,ic,e
          if(.not.bform)read(mru,end=1002)  !i0,id,id,id,id,ic,e
c          if(ic.gt.0)write(6,'(I5,2F15.6)')i,e,energ(jk(i0))
        enddo
c                                                             !radiation
        if(bform)read(mr,104,end=1002)
        if(.not.bform)read(mru,end=1002)
        if(bform)read(mr,103,end=1002)
c
 9131   if(bform)read(mr,132,end=1002) i1,i2
        if(.not.bform)read(mru,end=1002)i1,i2
        if(i2.eq.0)go to 9133
        go to 9131
c
 9133   nrr=nv
        nvx=nv
        go to 124
      else
        nvx=nv
      endif
C
C************************************
C READ HEADER, AND MAYBE ORBITAL CODE
C************************************
C
      DO I=1,MXORB
        QND(I)=0
      ENDDO
C
      IF(BFORM)THEN
 299    READ(MR,F101,END=1002) NCFD,NZ0D,NED,(QND(I),QLD(I),I=1,MXORB)
        if(kfpm.eq.0.and.qnd(mxorb).ne.0.and.mxorb.lt.mxorb0)then
          mxorb=mxorb+20
          mxorb=min(mxorb,mxorb0)
          backspace(mr)
          go to 299
        endif
      ELSE
       READ(MRU,END=1002,ERR=300)NCFD,NZ0D,NED,(QND(I),QLD(I),I=1,MXORB)
       GO TO 302
 300   IF(EX)THEN                                !START OF A FILE
         REWIND(MRU)                             !SO REWIND
         READ(MRU)
         READ(MRU)
       ELSE
         STOP 'UNABLE TO READ ORBITAL HEADER...' !SHOULD NOT GET HERE
       ENDIF                                 
      ENDIF
C
 302  NCF=NCFD                     !NOT EOF SO SAFE TO RELABEL
      IF(NCF.EQ.0)then
        write(6,*)'ncf=0, why??'
        stop 'ncf=0, why??'
cRETURN
      endif
c
      NZ0=NZ0D
      NE=NED
C
      DO I=1,MXORB
        QN(I)=QND(I)
        QL(I)=QLD(I)
      ENDDO
C
      IF(NZOLD.NE.0.AND.NZ0.NE.NZOLD)THEN
        WRITE(6,*)'*** ERROR: DIFFERENT ELEMENTS ON TWO FILES, NZ='
     X            ,NZOLD,NZ0
        STOP '*** ERROR: DIFFERENT ELEMENTS ON TWO FILES'
      ENDIF
      NZOLD=NZ0
C
      IF(NEOLD.NE.0.AND.NE.NE.NEOLD)THEN
        if(lv00.ge.0)then
          WRITE(6,*)'*** ERROR: DIFFERENT IONS ON TWO FILES, NE='
     X              ,NEOLD,NE
          STOP '*** ERROR: DIFFERENT IONS ON TWO FILES'
        else
          if(ne.ne.neold-1.or.ne.ne.neold+1)then
            WRITE(6,*)'*** WARNING: DIFFERENT CHARGES ON TWO FILES, NE='
     X                ,NEOLD,NE
          endif
        endif
      ENDIF
      NEOLD=NE
C
      DO I=1,MXORB                 !SHORT ORBITAL LIST
        IF(QN(I).LE.0)GO TO 301
      ENDDO
      I=MXORB+1
 301  MXORB=I-1
C
      IF(NCF.GT.NDIM14)THEN
        WRITE(6,136)NCF
 136  FORMAT(' DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM14 TO',I5)
        STOP 'ERROR: DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM14'
      ENDIF
C
      NZ=NZ0-NE+1
      DZ=NZ*NZ
      DEN=QDT(QD0,NZ0,NE,NV,LV,0)
      IF(IFLAGB.EQ.0)THEN
        IFLAGB=1
C        TOLB0=TOLB
        IF(TOLB0.LE.ZERO)TOLB=MAX(1.5D-7,1.0D-9*DZ*NZ)
        TOLBE=MAX(TOLBE,TOLB)
        IF(BFORM)TOLBE=MAX(TOLBE,2.D-6)
      ENDIF
C
      IF(TEMP(1).LE.ZERO)THEN
        DO J=1,JTEMP
          IF(BLSOLD)THEN
            TEMP(J)=DZ*THTLS(J)
          ELSE
            TEMP(J)=DZ*THTIC(J)
          ENDIF
          TJ=SQRT(TEMP(J))
          COEF(J)=2.0707D-16/(TJ*TEMP(J))
          TEMP(J)=TEMP(J)/1.5789D5
        ENDDO
      ENDIF
C
      IF(BHYBRD)THEN
        NS=1
        IF(BLSNEW)NS=NSSYM
        DO N=1,NS
          DO I=1,NCF
            KAUTY(I,N)=-2
          ENDDO
        ENDDO
      ENDIF
C
C************************
C READ CONFIGURATION DATA
C************************
C
      NCFT0=NCFT
      DO 102 I=1,NCF
C
      IF(BFORM)READ(MR,179,END=1002)NII(I),NGR,MA0,MB0,(QS0(J),QL0(J)
     X,J=1,10)
      IF(.NOT.BFORM)READ(MRU,END=1002)NII(I),NGR,MA0,MB0,(QS0(J),QL0(J)
     X,J=1,10)
 179  FORMAT(2I5,2X,I3,I2,1X,10(I2,A1))
C
      IN=IABS(NII(I))
      NG(IN)=NGR                         !not used...
C
C DECODE CONFIGURATIONS:
C   LMX(I) IS THE NO OF DISTINCT OPEN-SHELL ORBITALS IN CONFIG I.
C   QSB(I,J) IS THE OCCUPATION NO OF ORBITAL J IN CONFIG I.
C   QLB(I,J) IS THE ORBITAL NO OF ORBITAL J IN CONFIG I, J=1,LMX(I).
C   QS0,QL0 CONTAIN EISSNER SPECIFICATION OF CONFIG TO BE DECODED.
C   ICF(I) IS THE CONFIG NO OF CONFIG I IN THE MASTER LIST.
C
      IF(LV00.NE.LV0)ICF(I)=0
C
      DO 16 J=1,10
        QSB(I,J)=MBLNK
        IF(QL0(J).EQ.MBLNK)GO TO 16
        LMX(I)=J
        M=MOD(QS0(J),50)
        QMB(I,J)=M
        mlast=m
        IF(M.GT.0)QSB(I,J)=LIT(M)
        DO K=1,NLIT
          IF(QL0(J).EQ.LIT(K))GO TO 19
        ENDDO
        QLB(I,J)=0
        GO TO 16
  19    QLB(I,J)=K
        MXORBR=MAX(MXORBR,K)
        IF(NII(I).LT.0)THEN
          IF(IMATCH.eq.0.OR.IMATCH.EQ.-I)                !so IMATCH.lT.0
     X    MXOCC(K)=MAX(M,MXOCC(K))
        ELSE
          IF(IMATCH.EQ.-I)MXOCC(K)=MAX(M,MXOCC(K))
        ENDIF
  16  CONTINUE
C
      J=LMX(I)
      IF(NII(I).LT.0)THEN                                     !CONTINUUM
        J1M=J-1
        LMX(I)=J1M
        QSB(I,J)=MBLNK
        if(bpart.and.kfpm.eq.0)then                 !map to target o_str
          do n=1,ncfr
            if(j1m.ne.lmh(n))go to 225
            do j=1,j1m
              if(qsb(i,j).ne.qsh(n,j))go to 225
              if(qlb(i,j).ne.qlh(n,j))go to 225
            enddo
            k=ice(n)                                   !eo<-so, else k=n
            icfi(i)=k
            if(bprnt0)write(6,*)i,n,k
            go to 102
 225        continue
          enddo
          icfi(i)=0
        endif
      ELSE                                                        !BOUND
     X IF(LV00.NE.LV0)THEN
        MM=QLB(I,J)
        IF(MM.GT.MXORB.OR.MM.EQ.0)THEN
          WRITE(6,*)'***ERROR, CF=',NII(I),' USES ORBITAL NO=',MM
     X     ,' WHICH IS NOT DEFINED IN ORBITAL HEADER!!'
          STOP'***ERROR, NEED ORBITAL NOT DEFINED IN HEADER!!'
        ENDIF
        IF(IMATCH.EQ.-I.AND.(QN(MM).NE.NV.OR.LV.LT.0))THEN
          WRITE(6,*)'***ERROR, IMATCH=',IMATCH
     X                      ,' MUST BE A CONTINUUM/RYDBERG CF'
          STOP '***ERROR, IMATCH MUST BE A CONTINUUM/RYDBERG CF'
        ENDIF
C
        IF(BHYBRD)THEN                                 !JUST MATCH BY CF
          J1=J
          IF(QSB(I,J1).NE.LIT(1))THEN
            J1M=J1
          ELSE
            J1M=J1-1
          ENDIF
          DO N0=1,NCFR
            N=JKH(N0)                                  !N=SYMMETRY ORDER
            IF(J1M.NE.LMH(N))GO TO 227
            DO J=1,J1-1
              IF(QSB(I,J).NE.QSH(N,J))GO TO 227
              IF(QLB(I,J).NE.QLH(N,J))GO TO 227
            ENDDO
            IF(QSB(I,J1).NE.LIT(1))THEN
              IF(QLB(I,J1).NE.QLH(N,J1M))GO TO 227
            ENDIF
            ITARH(I)=N0                                 !N0=ENERGY ORDER
C         write(6,*)i,itarh(i)
            GO TO 226
 227        CONTINUE
          ENDDO
          ITARH(I)=9999
        ENDIF
C
C*******************
C SET-UP MASTER LIST
C*******************
C
 226    IF(QN(MM).NE.NV.OR.LV.LT.0)THEN   !SUITABLE
         IF(NBIN0.LE.0)THEN
          if(qn(mm).eq.ncmx.and.mlast.gt.1)beqn=.false.  !may need tweak
          DO 151 N=1,NCFT
            IF(LMX(I).NE.LMT(N))GO TO 151
            DO J=1,LMX(I)
              IF(QSB(I,J).NE.QST(N,J))GO TO 151
              IF(QLB(I,J).NE.QLT(N,J))GO TO 151
            ENDDO
            ICF(I)=N              !OLD MASTER
            GO TO 102
 151      CONTINUE
C NO MATCH, ADD TO MASTER LIST
          NCFT=NCFT+1
          IF(NCFT.GT.NDIM30)THEN
            WRITE(6,*)'*** INCREASE NUMBER OF MASTER CONFIGS NDIM30'
            STOP 'ERROR: INCREASE DIMENSION OF NDIM30'
          ENDIF
          ICF(I)=-NCFT            !NEW
c          DO J=1,NDIM32
c            ICQT(NCFT,J)=0
c          ENDDO
          LMT(NCFT)=LMX(I)
          DO J=1,LMX(I)
            QST(NCFT,J)=QSB(I,J)
            QLT(NCFT,J)=QLB(I,J)
          ENDDO
         ELSE
          ICF(I)=1
         ENDIF
C ATTEMPT TO DETERMINE MAX CORE N
          NCMX0=MAX(NCMX0,QN(MM))
        ELSE
          M=MAX(1,LMX(I)-1)
          M=QLB(I,M)
          NCMX0=MAX(NCMX0,QN(M))   !ONLY APPROX
        ENDIF
      ENDIF
C
 102  CONTINUE                     !END LOOP OVER CONFIG INPUT
C
      IF(NCMX0.GT.IBN(IB00).AND.INIT.EQ.0)THEN
        WRITE(6,*)' INCREASE NRSLMX TO AT LEAST',NCMX0
        STOP 'ERROR: INCREASE NRSLMX'
      ENDIF
C
C*****************************************
C DETERMINE IF CONFIG CONTRIBUTES TO COREX
C*****************************************
C
      IF(LV00.NE.LV0.or.LVH.LT.0)THEN
        IF(BPASS1.AND.LV.GE.0.AND.LVH.GT.0)THEN
          LVH=-1
          WRITE(6,89)
  89    FORMAT(//' *** FIRST NL-BLOCK CONTAINS ANY MASTER CONTRIBUTION')
          WRITE(0,89)
          IF(JCFJ.NE.999999)THEN
            WRITE(6,*)'(UNLESS EXCLUDED BY JCFJ)'
            WRITE(0,*)'(UNLESS EXCLUDED BY JCFJ)'
          ENDIF
        ELSE
          LVH=LV
        ENDIF
        DO I=1,NCF
          LCA(I)=0
        ENDDO
       IF(LCMN.LT.0)THEN                          !N-N'
        DO 120 I=1,NCF
          J=LMX(I)
          IF(NII(I).LT.0)THEN
            IF(NCMN.LT.10)THEN
              LCA(I)=-1              !TAG WRONG CONT FOR CORE EXCITATION
              DO L=1,J
                K=QLB(I,L)
                IF(QN(K).EQ.NCMN)LCA(I)=0                      !TAG GOOD
              ENDDO
            ENDIF
            GO TO 120
          ENDIF
          MM=QLB(I,J)
          IF(QN(MM).NE.NV.AND.LVH.GE.0.OR.I.GT.JCFJ)GO TO 120
          IF(QSB(I,J).EQ.LIT(1))J=J-1
          IF(NCMX.GT.0.AND.QN(QLB(I,J)).NE.NCMX)GO TO 120
          LM=-1
          IF(LVH.LT.0.AND.QN(MM).NE.NV)LM=-2
          IF(NCMN.LT.10)THEN
            M00=MB0
            DO L=1,J
              M00=M00+1
              IF(QLB(I,L).GT.M00)THEN
          	DO M=M00,QLB(I,L)-1
          	  IF(MXOCC(M).NE.0)THEN
          	    IF(QN(M).EQ.NCMN)LCA(I)=LM
          	    GO TO 120
          	  ENDIF
          	ENDDO
          	M00=QLB(I,L)
              ENDIF
              IF(QSB(I,L).NE.LIT(MXOCC(QLB(I,L))))THEN
          	 IF(QN(QLB(I,L)).EQ.NCMN)LCA(I)=LM
          	 GO TO 120
              ENDIF
            ENDDO
            IF(QN(QLB(I,J)).EQ.NCMN)LCA(I)=LM
          ELSE
            LCA(I)=LM	  !ALL ALLOWED BY JCFJ AND NCMX
          ENDIF
  120   ENDDO
       ELSE                                       !NL-N'L'
        DO N=1,MXORBR
          NOCC1(N)=0
        ENDDO
        IM=IMATCH
        JM=LMX(IM)
        IF(NII(IM).GT.0.AND.QSB(IM,JM).EQ.LIT(1))JM=JM-1
        DO J=1,JM
          K=QLB(IM,J)
          NOCC1(K)=QMB(IM,J)
        ENDDO
        if(lcmn.lt.10)then                       !nl-n'l'
          i1=1
          itwo=2
        else                                     !n-n', all l,l'
          i1=99
          itwo=99
        endif
        DO 119 I=1,NCF
          J=LMX(I)
          IF(NII(I).LT.0)THEN
            if(lcmn.lt.10)then
              DO L=1,J
                IF(QLB(I,L).NE.QLB(IM,L).OR.QSB(I,L).NE.QSB(IM,L))
     X             LCA(I)=-1         !TAG WRONG CONT FOR CORE EXCITATION
              ENDDO
            else
              LCA(I)=-1              !TAG WRONG CONT FOR CORE EXCITATION
              DO L=1,J
                K=QLB(I,L)
                IF(QN(K).EQ.NCMN)LCA(I)=0                      !TAG GOOD
              ENDDO
            endif
            GO TO 119
          ENDIF
          MM=QLB(I,J)
          IF(QN(MM).NE.NV.AND.LVH.GE.0.OR.I.GT.JCFJ)GO TO 119
c       if(i.eq.22)then
c         write(*,*)'Hello World'
c       endif
          DO N=1,MXORBR
            NOCC(N)=0
          ENDDO
          DO L=1,J
            K=QLB(I,L)
            NOCC(K)=QMB(I,L)
          ENDDO
          if(qn(mm).ne.nv)then
            do n=k,1,-1
              if(nocc(n).gt.nocc1(n))go to 117
            enddo
            stop 'core-excitation mis-match...?'
  117       k=n
          endif
          NOCC(K)=NOCC(K)-1
          IDIFF=0
          DO N=1,MXORBR
            IF(NOCC1(N).GT.NOCC(N))THEN
              IF(NOCC1(N)-NOCC(N).GT.I1)GO TO 119
              IF(QN(N).NE.NCMN)GO TO 119
              IF(LCMN.LT.10.AND.QL(N).NE.LCMN)GO TO 119
              IDIFF=IDIFF+1
            ELSEIF(NOCC1(N).LT.NOCC(N))THEN
              IF(NOCC(N)-NOCC1(N).GT.I1)GO TO 119
              IF(QN(N).NE.NCMX)GO TO 119
              IF(LCMX.GE.0.AND.QL(N).NE.LCMX)GO TO 119
              IDIFF=IDIFF+1
            ENDIF
          ENDDO
          LM=-1
          IF(LVH.LT.0.AND.QN(MM).NE.NV)LM=-2
          if(idiff.eq.0.and.ncmn.ne.ncmx)lm=0
          IF(IDIFF.LE.ITWO)LCA(I)=LM                   !WINNER
  119   ENDDO
       ENDIF
       IF(LVH.LT.0)THEN
         WRITE(6,*)' '
         WRITE(6,*)'AUTOIONIZING (MASTER) CONFIGURATIONS CONTRIBUTING:'
         DO I=1,NCF
           IF(LCA(I).EQ.-2)WRITE(6,*)I
         ENDDO
         WRITE(6,*)' '
       ENDIF
      ENDIF
C
C**************************
C READ AUTOIONIZATION RATES
C**************************
C
      BFAST=.NOT.BFORM.AND.NLCOR.EQ.0.AND..NOT.BCFM.AND.ACOR.LT.ZERO
     X      .AND.JCFJ.EQ.999999.AND..NOT.BCFA
C
      IF(BFORM)READ(MR,103,END=1002)
      IF(.NOT.BFORM)READ(MRU,END=1002)NZTEST,NDUME
      IF(BFORM)READ(MR,103,END=1002)
 103  FORMAT(A1)
      I=0
 111  I=I+1
C
      IF(BFORM)READ(MR,112,END=1002)I1,I2,IWA,JCA,I3,T1,T2,EION
      IF(.NOT.BFORM)READ(MRU,END=1002)I1,I2,IWA,JCA,I3,T1,T2,EION
 112  FORMAT(5I5,5X,1PE15.5,2(0PF15.6))
C
      IF(I2.EQ.0) GO TO 113
      IF(I.LT.NDIM12)THEN
        ICA=I1                               !ICA(I)
        ITA(I)=I2
        if(bpart)jcai(i)=jca
        IF(I3.EQ.0)THEN              !A.S. TBD (currently JCA=0 if I3=0)
          JTA(I)=JCA                         !JCA IS <0 
        ELSE
          JTA(I)=I3
        ENDIF
        AA(I)=T1
        EC(I)=T2-E1C(1)
C
        AA(I)=ABS(AA(I))
        IF(BFAST)GO TO 111
C
        BUNA=BUNA.OR.I3.EQ.0                 !BUNDLED AUGERS
C
        I=I-1
C****
        IF(ICA.GT.JCFJ)GO TO 111             !ICA(I+1)
        IF(JCFA.GT.0)THEN
          IF(-JCA.GT.JCFA)GO TO 111          !JCA(I+1)
        ELSEIF(JCFA.LT.0.and.lv.ge.0)THEN
          mn=qlb(-jca,lmx(-jca))
          if(qn(mn).eq.nv)go to 111
        ENDIF
C****
        IF(BCFM.AND.JCFRB.NE.ICA)GO TO 111   !ICA(I+1)
        I=I+1
        IF(NLCOR.LE.0)THEN
          AA(I)=AA(I)*ABS(ACOR)
        ELSE
          IF(LV+1.GT.0.AND.LV.LT.NDIM25)AA(I)=AA(I)*ACORL(LV+1)
        ENDIF
      ENDIF
C
      GO TO 111
C
 113  NUMA=I-1
C
      IF(.NOT.BINT.AND.NUMA.GE.NDIM12) THEN
        WRITE(6,73)NUMA
  73    FORMAT(' SR.CROSSJ: NUMBER OF AUTOIONIZATION RATES EXCEEDS '
     X        ,'STORAGE, INCREASE NDIM12 TO',I10)
        STOP   'ERROR: SR.CROSSJ: NUMBER OF AUTOIONIZATION RATES EXCEEDS
     X STORAGE, INCREASE NDIM12'
      ELSE
        NUMAX=MAX(NUMAX,NUMA)
      ENDIF
C
      IF(BUNA.AND..NOT.BHYBRD.AND.NBIN0.LE.0)THEN
        IFLAGR=-1
        GO TO 1000
      ENDIF
C
      IF(NECOR.LE.0)EION=EION+E1C(1)         !.AND.LV.GE.0
      EIONMN=MIN(EIONMN,EION)
C
C**************
C READ ENERGIES
C**************
C
      IF(NENG.GT.0)NENG0=NENG
      IF(BFORM)READ(MR,121,END=1002) NENG,ECORE
      IF(.NOT.BFORM)READ(MRU,END=1002) NENG,ECORE
 121  FORMAT(10X,I5,45X,F15.6)
C
      IF(NENG.EQ.0)THEN
        DO I=1,NDIM10
          JFIRST(I)=0
          JLAST(I)=-1
        ENDDO
        GO TO 124
      ENDIF
      ECORE0=ECORE
      DEN0=DEN
C
      IF(BFORM)READ(MR,105,END=1002)MTEST4
      IF(.NOT.BFORM)READ(MRU,END=1002)MTEST4
 105  FORMAT(26X,A4)
C
      BTEST=MTEST4.NE.MBLNK4                               !IC=TRUE
      IF(BLS.AND.BTEST)THEN
        IF(BCA)THEN
          WRITE(6,371)
 371      FORMAT(' RUN INITIALIZED FOR CA BUT LS/IC DATA FOUND')
          STOP 'ERROR: RUN INITIALIZED FOR CA BUT LS/IC DATA FOUND'
        ELSE
          WRITE(6,370)
 370      FORMAT(' RUN INITIALIZED FOR LS BUT IC DATA FOUND')
          STOP 'ERROR: RUN INITIALIZED FOR LS BUT IC DATA FOUND'
        ENDIF
      ENDIF
      IF(BIC.AND..NOT.BTEST)THEN
        WRITE(6,374)
 374    FORMAT(' RUN INITIALIZED FOR IC BUT LS DATA FOUND')
        STOP 'ERROR: RUN INITIALIZED FOR IC BUT LS DATA FOUND'
      ENDIF
C
      IF(NENG.GT.NDIM13)THEN
        WRITE(6,369)NENG
 369  FORMAT(' NUMBER OF ENERGIES EXCEEDS STORAGE,INCREASE NDIM13 TO'
     X      ,I6)
        STOP 'ERROR: NUMBER OF ENERGIES EXCEEDS STORAGE,INCREASE NDIM13'
      ENDIF
C
      IF(INR1.EQ.999)THEN
        NR1=NCMX0+1
        INR1=NR1
        WRITE(6,*)' '
        WRITE(6,*)'*** NR1 RESET TO:',NR1
      ENDIF
      IF(BEQN.AND.NR1.NE.NCMX0)THEN
        NR1=NCMX0
        WRITE(6,*)' '
        WRITE(6,*)'*** NR1 RESET TO:',NR1
      ENDIF
      IF(INR1.LE.0)THEN
        NR1=0
        BRAD=.FALSE.
        BEQN=.FALSE.
      ENDIF
C
      NRR=NV
      NAUTO=0
      IRAD=NENG+1
      MFLAG=0
      KFLAG=NENG
      NS=1
      LCF0=0
      ltest=0
      iflagw=iflagw0
C
      DO 122 I=1,NENG
C
      IF(BFORM)READ(MR,123,END=1002)IK(I),IT(I),SS(I),LL(I),JJ(I),LCF(I)
     X,ENERG(I)
      IF(.NOT.BFORM)READ(MRU,END=1002)IK(I),IT(I),SS(I),LL(I),JJ(I),LCF(
     XI),ENERG(I)
 123  FORMAT(5X,6I5,F15.6)
c      write(6,123)ik(i),it(i),ss(i),ll(i),jj(i),lcf(i),energ(i)
C
      M=IK(I)
      M=IABS(M)
      IF(M.LE.NDIM10)THEN
        JK(M)=I
        ITAG(M)=0
        JFIRST(M)=0
        JLAST(M)=-1
        KFIRST(M)=0
        KLAST(M)=-1
        IF(IMAP(M).NE.0)THEN
          IF(NMAP.GT.0)THEN
            IF(NV.GE.NMAP)LCF(I)=IMAP(M)
          ELSE
            IF(NV.EQ.-NMAP)LCF(I)=IMAP(M)
          ENDIF
        ENDIF
      ENDIF
      MFLAG=MAX(MFLAG,M)
      K=IABS(LCF(I))
      IF(IRAD.GT.NENG.AND.LCF(I).LT.0)IRAD=I
C
      IF(BLS)THEN
        ltest=ltest+ll(i)
        IF(BCA)THEN
          ILSJ(I)=0
          IF(BFORM)THEN
            SS(I)=100000*IT(I)+SS(I)           !as write i10 read as 2i5
          ENDIF
          IW8=SS(I)
        ELSE
          ILSJ(I)=10000*IABS(SS(I))+100*LL(I)
          IF(SS(I).LT.0)ILSJ(I)=ILSJ(I)+1
          IW8=IABS(SS(I))*(2*LL(I)+1)
        ENDIF
      ENDIF
      IF(BIC)THEN
        ILSJ(I)=100*IABS(JJ(I))
        IF(SS(I).LT.0)ILSJ(I)=ILSJ(I)+1
        IW8=IABS(JJ(I))+1
      ENDIF
C
      IF(NECOR.LT.0.AND.LCF(I).LT.0.AND.LV.GE.0)ENERG(I)=ENERG(I)+E1C(1)
      TE=ECORE+ENERG(I)
C
      IF(BHYBRD)THEN               !FLAG SUITABLE CFGS WHICH STRADDLE IP
        IF(LCF(I).GT.0.AND.ICF(K).EQ.0)THEN                     !RYDBERG
          IF(BLSNEW)THEN
            IABSSS=IABS(SS(I))
            NS=(IABSSS+1)/2
          ENDIF
          IF(IRAD.GT.NENG)THEN                         !STILL TRUE BOUND
            IF(KAUTY(K,NS).LT.0)KAUTY(K,NS)=-1
          ELSE                                         !WE ARE ABOVE IP
            IF(KAUTY(K,NS).EQ.-1)KAUTY(K,NS)=0         !STRADDLES
          ENDIF
        ENDIF
      ENDIF
C
C**************************
C INDEX AUTOIONIZING STATES (AND AVOID DOUBLE COUNTING)
C**************************
C
      IF(LCF(I).GT.0.AND.LCA(K).NE.0.AND.TE.GE.EIONMN-TOLE-TOLN)THEN
        IF(TE.LT.EIONMN+TOLE)IRAD=MAX(I,IRAD)
        IF(IK(I).GT.0)THEN                              !NOT CORRELATION
          NAUTO=NAUTO+1
          IAUTO(NAUTO)=M        !=IABS(IK(I))
          LCA(K)=IABS(LCA(K))
        ENDIF
        GO TO 122
      ENDIF
C
C********************************
C SET-UP TARGET BINS AND INDEXING (ONLY DONE FOR A NEW UNIT).
C********************************
C
      IF(KFPM.GT.NBINP) GO TO 122
      IF(LCF(I).GT.0.OR.IK(I).LT.0)GO TO 122          !BOUND/CORRELATION
      IF(SS(I).EQ.0)GO TO 122                     !UNWANTED BUNDLED CONT
      J1=LMX(K)
      M=QLB(K,J1)
      IF(QN(M).EQ.NV.AND.LV.GE.0)GO TO 122          !CORE RE-ARRANGEMENT
      IF(EIONMN.GE.ZERO)EIONMN=TE
C
      IF(ENERG(I).GT.(WNP0+TOLB))THEN                 !NEW TARGET/PARENT
C
        IF(LCF0.LT.0)THEN                           !FINALIZE OLD WEIGHT
          J1=LMX(-LCF0)+1
          M=QLB(-LCF0,J1)
          if(.not.bca)IWP0=IWP0/(4*QL(M)+2)         !as CA omits already
          IF(IWT(KFPM).EQ.0)THEN
            if(iflagw.eq.0)then
              IWT(KFPM)=IWP0
            else
              write(6,*)'Cannot determine target stat.weight internally'
              write(6,*)'Try reducing NTAR2 to ',kfpm-1
              stop 'Cannot determine target stat. weight internally'
            endif
          ELSEIF(IWP0.NE.IWT(KFPM).and.iflagw.eq.0)THEN
            WRITE(6,*)
            WRITE(6,*) 'TARGET STAT. WEIGHT MIS-MATCH?'
            WRITE(6,*)LCF0,KFPM,IWP0,IWT(KFPM)
            NFLAG2=kfpm-1
            write(6,*)'SHOULD REDUCE NTAR2 TO ',NFLAG2
            WRITE(6,*)
            STOP 'TARGET STAT. WEIGHT MIS-MATCH?'
c            iflagw=1
          ENDIF
        ENDIF
        if(iflagw.eq.0)LCF0=LCF(I)                    !CAN CHECK WEIGHTS
        IWP0=0
C
        WNP0=ENERG(I)
        KFPM=KFPM+1
C
        IF(BBIN)THEN
          EII(KFPM)=TE
C
          IF(KFPM.EQ.2)TOLN=MAX(TOLN,4.D-4*DZ-(EII(2)-EII(1)))
C
          IF(KFPM.LE.NBINM.AND.TOLR.GT.ZERO)THEN  !SET FINAL METAS RANGE
            T=EII(KFPM)-EII(1)
            IF(TOLR.LT.T)TOLR=T
          ENDIF
CKFPM.LE.-NECOR.or.kfpm.le.nbinrm
          IF(lv.ge.-1)THEN  	                  !CHECK TARGET ENERGIES
            T=TE-EII(1)
            T0=E1C(KFPM)
            IF(KFPM.EQ.1)THEN
              IF(BPRNT0)WRITE(6,372)TOLBE
 372          FORMAT(/3X,'IE',10X,'E(N)',14X,'E(N+1)',2X,'TOLB='
     X              ,1PE10.3)
              T0=T0-E1C(1)
            ENDIF
            MMM=MBLNK
            IF(ABS(T-T0).GT.TOLBE)THEN
              MMM=MSTAR
              if(kfpm.le.-necor)then
                DO KK=KFPM+1,NBINC
                  ECORI(KFPM,KK)=ECORI(KFPM,KK)-(T-T0)
                ENDDO
                DO KK=1,KFPM-1
                  ECORI(KK,KFPM)=ECORI(KK,KFPM)+(T-T0)
                ENDDO
              else
                if(e1c(kfpm).eq.zero)then
                  mmm=mblnk
                  iflage=iflage-1
                endif
              endif
              IFLAGE=IFLAGE+1
              E1C(KFPM)=T
            ENDIF
            IF(BPRNT0)                                 !.OR.MMM.NE.MBLNK
     X         WRITE(6,373)KFPM,T0,T,MMM
 373           FORMAT(I5,2F18.8,3X,A4)
          ENDIF
        ENDIF
C
        IF(KFPM.LE.NBINP)THEN
          WNP(KFPM)=-TE
          DO J=1,10
            QSP(KFPM,J)=QSB(K,J)
            QLP(KFPM,J)=QLB(K,J)
          ENDDO
          LMP(KFPM)=LMX(K)
c          lcp(kfpm)=k                                   !for test print
        ELSEIF(BBIN.AND..NOT.BHYBRD)THEN
          IF(NR1.NE.0.AND.NBIN0.LE.0)THEN
            T=EII(KFPM)+QDT(QD0,NZ0,NE,INR1,0,0)
            IF(T.LT.EIONMN)THEN      !OMIT +TOLR SINCE NOT META BY THEN
              WRITE(6,744) T,EIONMN,(EII(J),J=1,KFPM)
 744  FORMAT(' STRONG WARNING IN SR.CROSSJ, THERE MAYBE A FINAL-STATE'
     X,' WITH PARENT NOT SPECIFIED BY NTAR2'//' INCREASE NTAR2 AND/OR'
     X,' NDIM0'//'EBOUND=',F12.3,5X,'EION=',F12.3,' EBIN='/
     X (10F12.6))
c      STOP 'ERROR: INCREASE NTAR2 AND/OR NDIM0'
            ENDIF
          ENDIF
        ENDIF
C
      ELSEIF(KFPM.LE.NBINP)THEN                    !CHECK LABELLING
C
        IF(LMP(KFPM).NE.LMX(K))GO TO 381
        DO J=1,LMX(K)
          IF(QSP(KFPM,J).NE.QSB(K,J))GO TO 381
          IF(QLP(KFPM,J).NE.QLB(K,J))GO TO 381
        ENDDO
        GO TO 382                                  !WE ARE GOOD
c
C FLAG
 381    IF(NCFR.EQ.0)THEN
          WRITE(6,*)'LEVEL ',I,
     X      ' HAS AN INCONSISTENT TARGET CONFIG LABEL='
     X      ,K,', REDUCE NTAR2 TO',KFPM-1
         IF(iflagw.eq.0)NFLAG2=KFPM-1
         iflagw=iflagw+1
c        write(6,*)i,kfpm,k
c        write(6,383)((qsp(kfpm,j),qlp(kfpm,j)),j=1,lmp(kfpm)-1)
c        write(6,383)((qsb(k,j),qlb(k,j)),j=1,lmx(k)-1)
c 383    format(10(a2,i2,3x))
c        STOP ' INCONSISTENT TARGET CONFIG LABEL, REDUCE NTAR2'
        ELSE
          KFLAG=MIN(KFLAG,KFPM-1)
        ENDIF
C
      ENDIF
C
 382  CONTINUE
C
      IF(LCF0.EQ.LCF(I))IWP0=IWP0+IW8
      WNP0=ENERG(I)           !ALLOW FOR ANY DRIFT OF CONTINUUM ENERGIES
C
 122  CONTINUE                  !END ENERGY LOOP READ
C
      IF(MFLAG.GT.NDIM10)THEN
        WRITE(6,368)MFLAG
 368  FORMAT(' NUMBER OF LEVELS EXCEEDS STORAGE,INCREASE NDIM10 TO',I6)
        STOP 'ERROR: NUMBER OF LEVELS EXCEEDS STORAGE,INCREASE NDIM10'
      ELSE
        MENG=MFLAG              !FOR CORRELATION LABELS
      ENDIF
C
C*****************************************
C END ENERGY READ AND INDEXING BY SYMMETRY
C*****************************************
C
      if(bca)then
        if(ltest.ne.0)then
          write(6,371)
          STOP 'ERROR: RUN INITIALIZED FOR CA BUT LS DATA FOUND'
        endif
      elseif(bls)then
        if(ltest.eq.0)then        !simple LS case may have all s-states,
          write(6,361)            !so just warn
 361  format('warning, run initialized for LS but data found may be CA')
        endif
      endif
C
C DETERMINE SAFE TOLI
C
      BRFRST=TOLI0.LT.ZERO.AND.LV.GE.0
      DO I=1,NENG
      K=IABS(LCF(I))
      IF(BRFRST.AND.LCF(I).GT.0.AND.QN(QLB(K,LMX(K))).EQ.NV)THEN
        TE=ENERG(I)+ECORE
        TOLI=TE-DEN-EIONMN
c        write(37,*)nv,lv,k,toli
        BRFRST=.FALSE.
        IF(TOLI.LT.ZERO.OR.KFPM.LT.2)THEN
          TOLI=ZERO                 !LOWER ONLY?
        ELSE
          IF(BLS)THEN
            IE=2
 366        DE=EII(IE)-EII(IE-1)
            IF(TOLI.GT.DE*0.5D0)THEN     !NOT GROUND
              TOLI=TE-DEN-EIONMN-EII(IE)+EII(1)
c              write(38,*)nv,lv,k,i,ie,toli
              IF(TOLI.LT.ZERO)TOLI=ZERO
              IE=IE+1
              IF(IE.LE.KFPM)GO TO 366
            ENDIF
          ENDIF
          IF(BIC)THEN
            IE=2
            IL=ILVTM(IE)
            IF(IL.EQ.0)THEN
              TOLI=ZERO
              GO TO 367
            ENDIF
            ILM=ILVTM(IE-1)
 365        DE=EII(IL)-EII(ILM)
            IF(TOLI.GT.DE*0.5D0)THEN     !NOT GROUND
              TOLI=TE-DEN-EIONMN-EII(IL)+EII(1)
c              write(38,*)nv,lv,k,i,ie,il,toli
              IF(TOLI.LT.ZERO)TOLI=ZERO
              IE=IE+1
              ILM=IL
              IL=ILVTM(IE)
              IF(IL.GT.0.AND.IL.LE.KFPM)GO TO 365
            ENDIF
          ENDIF
        ENDIF
        GO TO 367
      ENDIF
      ENDDO
C
 367  CONTINUE
C
C AS BUNDLED CASE MAY NOT HAVE ALL TARGETS, SO:
C
      KHOLD=KFPM
      IF(KFLAG.LT.NBINP)THEN
        IF(LCP(1).GT.0)THEN
          WRITE(6,*)
          WRITE(6,*)'COMMENT: USING TARGET CONFIG LABELS FROM FILE '
     X,'o_str'                    !or ou_str...
C          KFPM=KFLAG             !ASSUME TARGET CFS CORRECT TO HERE
          KFPM=0                  !MIGHT AS WELL REPLACE ALL...
        ELSE
          WRITE(6,*)'*** ERROR: NEED TARGET CF NOS IN adasin TO USE '
     X              ,'o_str'
          STOP '*** ERROR: NEED TARGET CF NOS IN adasin - see adasout'
        ENDIF
      ENDIF
C
      IF(KFPM.LT.NBINP)THEN
        DO KK=KFPM+1,NBINP
          K=LCP(KK)                                        !TARGET CF NO
          IF(K.GT.NCFR.or.k.eq.0)GO TO 275          !for ioldw=1 (no CF)
C          IF(BHYBRD)THEN                         !NEED ENERGY FOR ADF09
            TE=E1C(KK)+EIONMN
            IF(BBIN)EII(KK)=TE
            WNP(KK)=-TE
C          ENDIF
          DO J=1,10
            QSP(KK,J)=QSH(K,J)
            QLP(KK,J)=QLH(K,J)
          ENDDO
          LMP(KK)=LMH(K)
        ENDDO
        IF(BHYBRD)THEN
          IF(KHOLD.LT.NBINM)GO TO 274                  !WE NEED INITIALS
          KFPM=NBINP
          IF(BBIN.AND.E1C(NBINP+1).GT.ZERO)THEN
            KFPM=KFPM+1
            EII(KFPM)=E1C(KFPM)+EIONMN
          ENDIF
        ELSE
          KFPM=KHOLD
        ENDIF
        GO TO 279
C
 274    WRITE(6,1378)NBINM,KFPM,TOLB
1378  FORMAT(//' *** ERROR: YOU HAVE REQUESTED NTAR1=',I2,' PARENTS BUT'
     X,' ONLY ',I2,' CAN BE DETERMINED FROM YOUR AUTOSTRUCTURE DATA.'/
     X/' *** POSSIBLE CAUSES:'//
     X'1.   THERE ARE NOT NTAR1 TARGET CONTINUA - '
     X,'CHECK YOUR AUTOSTRUCTURE TARGET.'/
     X'2.   CHECK AUTOSTRUCTURE NMETAR/J, TERMS/LEVELS *AND* EMXLS/IC.'/
     X'3.   LEVEL SPLITTING IS .LT. 1.D-6 BUT O1 IS BEING USED - '
     X,'SWITCHED TO UNFORMATTED DATA ON O1U.'/
     X'4.   TOLB IS TOO LARGE - TRY SETTING IT TO LESS THAN HALF OF '
     X,'THE SMALLEST LEVEL SPLITTING: TOLB=',1PE9.2)
        IF(LV.EQ.-1)WRITE(6,1379)
1379  FORMAT(//' *** MOST LIKELY CAUSE:'/
     X/'     SEPARATION BETWEEN TARGET NTAR1'
     X,' AND NTAR1+1 IS COMPARABLE WITH ENERGY DIFFERENCE'/
     X,'     BETWEEN USING N- AND (N+1)-ELECTRON ORBITAL POTENTIALS'/
     x/' *** WORKAROUND (IN AUTOSTRUCTURE):'//
     X'     EITHER SET EIXMLS/IC MANUALLY/EXPLICITLY IN SMINIM'/
     X'     OR INCREASE NMETAR/J TO REACH STATE WELL ABOVE LOWER ONE(S)'
     X)
      STOP'ERROR: NUMBER OF PARENTS REQUESTED EXCEEDS NO. FOUND ON FILE'
C
 275    WRITE(6,378)NBINP,KK-1,TOLB
 378    FORMAT(' ERROR: YOU HAVE REQUESTED NTAR2=',I2,' PARENTS BUT'
     X,' ONLY ',I2,' CAN BE DETERMINED FROM YOUR AUTOSTRUCTURE DATA.'/
     X/' POSSIBLE CAUSES:'/
     X'1.   THERE ARE NOT NTAR2 TARGET CONTINUA - '
     X,'CHECK YOUR AUTOSTRUCTURE TARGET.'/
     X'2.   LEVEL SPLITTING IS .LT. 1.D-6 BUT O1 IS BEING USED - '
     X,'SWITCHED TO UNFORMATTED DATA ON O1U.'/
     X'3.   TOLB IS TOO LARGE - TRY SETTING IT TO LESS THAN HALF OF '
     X,'THE SMALLEST LEVEL SPLITTING: TOLB=',1PE9.2)
      STOP'ERROR: NUMBER OF PARENTS REQUESTED EXCEEDS NO. FOUND ON FILE'
      ENDIF
C
 279  IF(BBIN.AND.EI(1).GE.ZERO)THEN
        IF(KFPM.LE.NBINP)EII(KFPM+1)=0.8D0*EII(NBINP)
        EI(1)=1.2D0*EII(1)
        DO M=1,NBINP
          EI(M+1)=EII(M+1)-0.5D0*(EII(M+1)-EII(M))
          if(bprnt0)write(6,327)m+1,eii(m+1)
 327  format(i5,f20.8)
        ENDDO
      ENDIF
      KFPM=NENG
C
      IF(BPASS1)THEN
      IF(BPRTM2)WRITE(6,90)
  90    FORMAT(//4X,'N  L',6X,'CROSS(MB)',4X,'CROSS(MB)',4X,'CROSS(MB)',
     X  4X,'CROSS(MB)',4X,'CROSS(MB)',4X,'CROSS(MB)',4X,'CROSS(MB)',4X
     X  ,'CROSS(MB)',4X,'CROSS(MB)')
        IF(BPRNT0)WRITE(6,33)
  33    FORMAT(3X,'CF',4X,'J',3X,'IP',' WI','  JP',' WJ',8X,'EC(J)',10X
     X  ,'SUMAN',10X,'SUMAD',10X,'SUMRN',10X,'SUMRD',8X,'CROSS(MB)')
      ENDIF
C
      IF(BPRNT0.AND.NV.GE.0)WRITE(6,34)NV,LV
  34  FORMAT(I5,I3)
C
C**************
C DETERMINE NR2
C**************
C
      IF(TOLR.EQ.ZERO)THEN
        IF(NR2.LT.0)THEN
          IF(BLSOLD)THEN
            IF(NBINRM.GE.2)THEN
              DE=EII(2)-EII(1)
              NR2=INT(NZ/SQRT(DE))
              NR2=MIN(NR2,100)
              NR2=MAX(NR2,15)
            ELSE
              NR2=100
            ENDIF
            WRITE(6,*)' '
            WRITE(6,*)'*** NR2 RESET TO:',NR2
          ELSE
            NR2=9999
          ENDIF
        ENDIF
C
      ELSE
        IF(INIT.EQ.0)THEN
          NR2=NV
        ELSE
          IF(NR2.LE.0)NR2=9999
        ENDIF
      ENDIF
C
C*****************************************************
C SET-UP RESOLVED TERMS/LEVELS AND DUMMY FOR PARENTAGE
C*****************************************************
C
      IF(BRSLP.AND.NR1.NE.0.AND.NRSOLZ.EQ.0)THEN   
C                          .AND.NBINP.AND.BPASS1 !drop to re-map
C
       DO KK=1,NBINP
C
C TERMS
C
        IF(BLS)THEN
C RESOLVED                                         nbin0.le.0
          if(bfirst.and.(lv.ge.0.or.init.ne.0).and.BRSLF     )then !need
c            IF(KK.GT.NBINRM)GO TO 475
            ISP=1
            IF(IWS(KK).GT.1)ISP=3
            NUP=NRSLMX
            IF(NCUT.GT.0)NUP=MIN(NCUT,NRSLMX)
            DO N=INR1,NUP
              IN=(N*(N-1))/2 + MXORB0
              LUP=N
              IF(LCUT.GE.0)LUP=MIN(LCUT+1,N)
              DO L=1,LUP
           	E=-WNP(KK)+QDT(QD0,NZ0,NE,N-1,L-1,0)  !ERR ON CAUTIOUS
           	IF(E.GT.EIONMN+TOLR)GO TO 471
           	E=-WNP(KK)+QDT(QD0,NZ0,NE,N,L-1,0)  !MAYBE AUTOIONIZING
           	IO=IN+L
           	L1=IABS(L-1-IWL(KK))
           	L2=IABS(L-1+IWL(KK))
           	DO LT=L1,L2
           	  DO IS=1,ISP,2
           	    NRSOL=NRSOL+1
           	    IF(NRSOL.LE.NDIM17)THEN
           	      ITARR(NRSOL)=KK
           	      QNV(NRSOL)=0
           	      SSR(NRSOL)=IWS(KK)+2-IS
           	      LLR(NRSOL)=LT
           	      JJR(NRSOL)=0
           	      WNR(NRSOL)=E
           	      LMR(NRSOL)=LMP(KK)+1
c     write(6,*)nrsol,kk,lcp(kk),l-1,ssr(nrsol),lt,lmr(nrsol),wnr(nrsol)
           	      DO J=1,10
           		QSR(NRSOL,J)=QSP(KK,J)
           		QLR(NRSOL,J)=QLP(KK,J)
           	      ENDDO
           	      J=LMR(NRSOL)
           	      QSR(NRSOL,J)=LIT(1)
           	      QLR(NRSOL,J)=IO
           	    ENDIF
           	  ENDDO
           	ENDDO
              ENDDO
 471          IF(L.EQ.1)GO TO 470
            ENDDO
 470        CONTINUE
          endif
c dummy
          ISP=1
          IF(IWS(KK).GT.1)ISP=3
          L00=0
          LUPZ=LMAXZ
          IF(LCUT.GE.0)LUPZ=MIN(LCUT+1,LMAXZ)
          DO L=1,LUPZ			   !DUMMY
            IF(INIT.NE.0)L00=L
            L1=IABS(L-1-IWL(KK))
            L2=IABS(L-1+IWL(KK))
            DO LT=L1,L2
              DO IS=1,ISP,2
        	NRSOLZ=NRSOLZ+1
        	IF(NRSOLZ.LE.NDIM18)THEN
        	  ITARZ(NRSOLZ)=KK
        	  SSZ(NRSOLZ)=IWS(KK)+2-IS
        	  LLZ(NRSOLZ)=LT
        	  JJZ(NRSOLZ)=MAX(INIT,L00)
        	  LMZ(NRSOLZ)=LMP(KK)+1
c        write(6,*)nrsolZ,kk,999,l-1,ssZ(nrsolZ),lt,lmZ(nrsolZ)
        	  DO J=1,10
        	    QSZ(NRSOLZ,J)=QSP(KK,J)
        	    QLZ(NRSOLZ,J)=QLP(KK,J)
        	  ENDDO
        	  J=LMZ(NRSOLZ)
        	  QSZ(NRSOLZ,J)=LIT(1)
        	  QLZ(NRSOLZ,J)=L-1		  !JUST L
        	ENDIF
              ENDDO
            ENDDO
          ENDDO
C
        ENDIF
C
C LEVELS
C
        IF(BIC)THEN
C RESOLVED                                         nbin0.le.0
          if(bfirst.and.(lv.ge.0.or.init.ne.0).and.BRSLF     )then !need
C     	    IF(KK.GT.NBINRM)GO TO 475
     	    NUP=NRSLMX
     	    IF(NCUT.GT.0)NUP=MIN(NCUT,NRSLMX)
     	    DO N=INR1,NUP
     	      IN=(N*(N-1))/2 + MXORB0
     	      LUP=N
     	      IF(LCUT.GE.0)LUP=MIN(LCUT+1,N)
     	      DO L=1,LUP
     	  	E=-WNP(KK)+QDT(QD0,NZ0,NE,N-1,L-1,0)  !ERR ON CAUTIOUS
     	  	IF(E.GT.EIONMN+TOLR)GO TO 473
     	  	E=-WNP(KK)+QDT(QD0,NZ0,NE,N,L-1,0)  !MAYBE AUTOIONIZING
     	  	IO=IN+L
     	  	JV1=IABS(2*L-3)
     	  	JV2=IABS(2*L-1)
     	  	DO JVT=JV1,JV2,2
     	  	  JMIN=IABS(JVT-IWJ(KK))
     	  	  JMAX=IABS(JVT+IWJ(KK))
     	  	  DO JT=JMIN,JMAX,2
     	  	    NRSOL=NRSOL+1
     	  	    IF(NRSOL.LE.NDIM17)THEN
     	  	      ITARR(NRSOL)=KK
     	  	      QNV(NRSOL)=0
     	  	      SSR(NRSOL)=0
     	  	      LLR(NRSOL)=0
     	  	      JJR(NRSOL)=JT 
     	  	      JVR(NRSOL)=JVT		   
     	  	      WNR(NRSOL)=E
     	  	      LMR(NRSOL)=LMP(KK)+1
     	  	      DO J=1,10
     	  		QSR(NRSOL,J)=QSP(KK,J)
     	  		QLR(NRSOL,J)=QLP(KK,J)
     	  	      ENDDO
     	  	      J=LMR(NRSOL)
     	  	      QSR(NRSOL,J)=LIT(1)
     	  	      QLR(NRSOL,J)=IO
     	  	    ENDIF
     	  	  ENDDO
     	  	ENDDO
     	      ENDDO
 473          IF(L.EQ.1)GO TO 472
     	    ENDDO
 472        CONTINUE
          endif
C
          L00=0
          LUPZ=LMAXZ
          IF(LCUT.GE.0)LUPZ=MIN(LCUT+1,LMAXZ)
          DO L=1,LUPZ                      !DUMMY
            IF(INIT.NE.0)L00=L
            JV1=IABS(2*L-3)
            JV2=IABS(2*L-1)
            DO JVT=JV1,JV2,2
              JMIN=IABS(JVT-IWJ(KK))
              JMAX=IABS(JVT+IWJ(KK))
              DO JT=JMIN,JMAX,2
                NRSOLZ=NRSOLZ+1
                IF(NRSOLZ.LE.NDIM18)THEN
                  ITARZ(NRSOLZ)=KK
                  SSZ(NRSOLZ)=MAX(INIT,L00)
                  LLZ(NRSOLZ)=0
                  JJZ(NRSOLZ)=JT 
                  LMZ(NRSOLZ)=LMP(KK)+1
c        write(6,*)nrsolZ,kk,999,l-1,jvt,jjz(nrsolz)
                  DO J=1,10
                    QSZ(NRSOLZ,J)=QSP(KK,J)
                    QLZ(NRSOLZ,J)=QLP(KK,J)
                  ENDDO
                  J=LMZ(NRSOLZ)
                  QSZ(NRSOLZ,J)=LIT(1)
                  QLZ(NRSOLZ,J)=L-1           !JUST L
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
       ENDDO
C 475   CONTINUE
       IF(NRSOL.GT.NDIM17)THEN
         WRITE(6,379)NRSOL
 379     FORMAT(' SR.CROSSJ: INCREASE NDIM17 TO:',I6)
         STOP 'ERROR:  DIMENSION ERROR: INCREASE NDIM17'
       ENDIF
       IF(NRSOLZ.GT.NDIM18)THEN
         WRITE(6,380)NRSOLZ
 380     FORMAT(' SR.CROSSJ: INCREASE NDIM18 TO:',I6)
         STOP 'ERROR: DIMENSION ERROR: INCREASE NDIM18'
       ENDIF
C
      ENDIF
c
      if(bfirst.and.(lv.ge.0.or.init.ne.0).and.nbin0.le.0)
     x                                      bfirst=.false. !table set-up
C
C************************************************
C SET-UP INDEXING FOR TERMS WITHIN MASTER CONFIGS
C************************************************
C
      IF(LV00.NE.LV0.AND.BRSLF)THEN
C
C QTI(L) IS THE NO OF TERMS IN MASTER CONFIG L.
C QTT(M) IS THE POSITION OF CURRENT TERM M WITHIN ITS MASTER CONFIG.
C THE POSITION OF A TERM WITHIN A CONFIG IN THE MASTERLIST NEVER CHANGES
C
        DO L=1,NCFT
          QTI(L)=0
        ENDDO
C
        DO M=1,MENG
          I=JK(M)
          IF(I.GT.0)THEN
            K=IABS(LCF(I))
            IF(ICF(K).NE.0)THEN
              L=IABS(ICF(K))
              QTI(L)=QTI(L)+1
              QTT(M)=QTI(L)
            ENDIF
          ENDIF
        ENDDO
C DIM TEST
        DO L=NCFT0+1,NCFT
          IF(QTI(L).GT.NDIM32)THEN
            WRITE(6,*)'DIMENSION EXCEEDED IN CROSSJ: INCREASE NDIM32 TO'
     X                ,QTI(L)
            STOP 'ERROR: DIMENSION EXCEEDED IN CROSSJ: INCREASE NDIM32'
          ENDIF
        ENDDO
C
C DETERMINE SYMMETRY GROUPS WITHIN A MASTER CONFIG L
C NGG(L) IS NO OF GROUPS FOR CONFIG L
C QTE(L,N) IS A TEMP HOLD OF THE SYMMETRY OF GROUP N
C QTTG(L,QTT) IS THE GROUP NO OF A TERM WITHIN THE CONFIG
C
        DO L=NCFT0+1,NCFT
          NGG(L)=0
          QTE(L,0)=-1
        ENDDO
C
        DO M=1,MENG
          I=JK(M)
          IF(I.GT.0)THEN
            K=IABS(LCF(I))
            IF(ICF(K).LT.0)THEN           !SET UP GROUPS FOR NEW CONFIGS
              L=-ICF(K)
              N=NGG(L)
              IF(QTE(L,N).NE.ILSJ(I))THEN                     !NEW GROUP
                N=N+1
                NGG(L)=N
                IF(N.LE.NDIM6)QTE(L,N)=ILSJ(I)
              ENDIF
              QTTG(L,QTT(M))=NGG(L)
            ENDIF
          ENDIF
        ENDDO
C DIM TEST
        DO L=NCFT0+1,NCFT
          IF(NGG(L).GT.NDIM6)THEN
            WRITE(6,*)'DIMENSION EXCEEDED IN CROSSJ: INCREASE NDIM6 TO'
     X                ,NGG(L)
            STOP 'ERROR: IMENSION EXCEEDED IN CROSSJ: INCREASE NDIM6'
          ENDIF
        ENDDO
C
C DETERMINE ENERGY ORDER POSITION OF TERM M WITHIN ITS SYMMETRY GROUP:
C QTTE(M)
C
        DO L=1,NCFT
          DO N=1,NGG(L)
            QTE(L,N)=0
          ENDDO
        ENDDO
C
        DO I=1,NENG
          K=IABS(LCF(I))
          IF(ICF(K).NE.0)THEN
            L=IABS(ICF(K))
            M=IABS(IK(I))
            N=QTTG(L,QTT(M))
            QTE(L,N)=QTE(L,N)+1
            QTTE(M)=QTE(L,N)
          ENDIF
        ENDDO
C DIM TEST 
        DO L=NCFT0+1,NCFT
          DO N=1,NGG(L)
            IF(QTE(L,N).GT.NDIM19)THEN
            WRITE(6,*)'DIMENSION EXCEEDED IN CROSSJ: INCREASE NDIM19 TO'
     X                ,QTE(L,N)
             STOP 'ERROR: DIMENSION EXCEEDED IN CROSSJ: INCREASE NDIM19'
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      LV00=LV0
C
C***********************************************
C DETERMINE PARENTS AND EXISTING RESOLVED STATES
C***********************************************
C
c      write(777,*)nv,lv,e00,wnr(1)
C
      NS=1
      N0=0
      IFIRST=1
      DO 36 I=1,NENG         !36
C
        IF(NBIN0.LE.0)IRSOL(I)=0
        KK=IABS(IK(I))
        ITAR(KK)=-1
        K=IABS(LCF(I))
        J1=LMX(K)
        M=QLB(K,J1)
C
        IF(LCF(I).LT.0)THEN               !CONTINUUM
C
          IF(QN(M).EQ.NV.AND.LV.GE.0)GO TO 36       !CORE RE-ARRANGEMENT
          TE=ENERG(I)+ECORE
          DO M=IFIRST,NBINP
            IF(TE.GE.EI(M).AND.TE.LT.EI(M+1))THEN
              ITAR(KK)=M
              IF(M.GT.IFIRST)IFIRST=M   
              GO TO 36
            ENDIF
          ENDDO
C
        ELSE                              !DISCRETE
C
          IF(QN(M).EQ.NV.AND.LV.GE.0)THEN
            TOLII=TOLI
          ELSE
            TOLII=ZERO
          ENDIF
C
          IABSSS=IABS(SS(I))
          IABSJJ=IABS(JJ(I))
C HYBRID
          L=0
          N=0
          IF(BHYBRD)THEN
C
            L=ICF(K)
C
            TE=ENERG(I)+ECORE
            IF(ITARH(K).EQ.9999)THEN
              IF(TE-TOLII.LT.EIONMN+TOLR.and.l.eq.0)THEN !.or.lcmn.lt.0)
                WRITE(6,775)I
                STOP 'ERROR: PARENT NOT FOUND: SEE adasout FOR CAUSES'
              ELSE
                if(l.ne.0)then      !.and.lcmn.ge.0    !relax for master
                  bflagp=.true.
                  if(bprnt0)
     x               WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
                else
                  WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
                endif
              ENDIF
            ENDIF
            LTEST=IABS(L)
            IF(BLSNEW)THEN
              NS=(IABSSS+1)/2
              LTEST=1000*NS+LTEST
            ENDIF
            IF(BLS)JW=IABSSS*(2*LL(I)+1)
            IF(BIC)JW=IABSJJ+1
            N=JVR(LTEST)
            IF(N.NE.0)THEN                    !MASTER (NOT "FIRST" TIME)
              IF(N.GT.0)THEN
                IRSOL(I)=-1
                IF(ITARH(K).GT.NBINM0.AND.TE.GE.EIONMN)GO TO 36
                IF(L.LT.0)THEN            !MASTER, BUT STILL FIRST BLOCK
                  IF(BLS)SSR(N)=SSR(N)+JW
                  IF(BIC)JJR(N)=JJR(N)+JW
                  WNR(N)=WNR(N)+JW*TE                    !-TOLII=0
                ENDIF
                IRSOL(I)=N        !FOR POSS FUTURE CODE OF AUGER BREAKUP
              ELSE
                IRSOL(I)=-1                              !NOT WANTED...
              ENDIF
              GO TO 36
            ELSEIF(L.LT.0)THEN                           !NEW MASTER
              NP=ITARH(K)
            ELSEIF(KAUTY(K,NS).EQ.0.and..not.BINT)THEN   !NEW STRADDLER
              NAUTY=NAUTY+1
              IF(NAUTY.GT.NDIM38)THEN
                WRITE(6,*)'TOO MANY STRADDLING CFGS, INCREASE NDIM38'
                STOP 'TOO MANY STRADDLING CFGS, INCREASE NDIM38'
              ENDIF
              KAUTY(K,NS)=NAUTY
              IAUTY(NAUTY)=10000*ITARH(K)+20*IB+LV !NO CHECK IF SAFE!
              IF(BLSNEW)IAUTY(NAUTY)=IAUTY(NAUTY)+1000000*NS
              JWRN(NAUTY)=JW                        !SINCE MUST BE BOUND
              ERN(NAUTY)=JW*TE
              if(te.gt.eionmn)stop 'hybrid: first level auto...'
              JWRD(NAUTY)=0
              ERD(NAUTY)=ZERO
            ELSEIF(KAUTY(K,NS).GT.0)THEN                  !OLD STRADDLER
              NY=KAUTY(K,NS)
              IF(TE.LT.EIONMN)THEN
                ERN(NY)=ERN(NY)+TE*JW
                JWRN(NY)=JWRN(NY)+JW
              ELSE
                ERD(NY)=ERD(NY)+TE*JW
                JWRD(NY)=JWRD(NY)+JW
              ENDIF
            ENDIF
C
C CA: SET-UP ALL PARENTS AND MASTER AS POTENTIAL FINAL RECOMBINED
C
          ELSEIF(BCA.AND.NBIN0.LE.0)THEN               !JUST MATCH BY CF
C
            L=ICF(K)
            IF(L.GT.0)THEN                      !MASTER (NOT FIRST TIME)
              N=JVR(L)
              IF(N.GT.0)THEN
                ITAR(KK)=ITARR(N)
                IRSOL(I)=N
              ELSE
                IRSOL(I)=-1                              !NOT WANTED...
                ITAR(KK)=9999
              ENDIF
            ELSE                                         !ALL
              TE=ENERG(I)+ECORE
              IF(QSB(K,J1).NE.LIT(1))THEN
                J1M=J1
              ELSE
                J1M=J1-1
              ENDIF
              DO NP=1,NBINRM
                IF(J1M.NE.LMP(NP))GO TO 228
                DO J=1,J1-1
                  IF(QSB(K,J).NE.QSP(NP,J))GO TO 228
                  IF(QLB(K,J).NE.QLP(NP,J))GO TO 228
                ENDDO
                IF(QSB(K,J1).NE.LIT(1))THEN
                  IF(QLB(K,J1).NE.QLP(NP,J1M))GO TO 228
                ENDIF
                GO TO 360                              !WE HAVE A WINNER
  228           CONTINUE
              ENDDO
              IF(TE-TOLII.LT.EIONMN+TOLR.and.l.eq.0)THEN !.or.lcmn.lt.0)
                WRITE(6,775)I
                STOP 'ERROR: PARENT NOT FOUND: SEE adasout FOR CAUSES'
              ELSE
                IF(NCFR.GT.0)                    !SHOULD BE ABLE TO FIND
     X             WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR I=',I
                if(l.ne.0)then       !.and.lcmn.ge.0   !relax for master
                  bflagp=.true.
                  if(bprnt0)
     x               WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
                else
                  WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
                endif
              ENDIF
              NP=9999
  360         ITAR(KK)=NP
c      write(6,*)kk,np
            ENDIF
          ENDIF
C
C HYBRID/CA: SET (NEW) FINAL RECOMBINED
C
          IF(L.LT.0.AND.N.EQ.0)THEN
            L=-L
            IF(BLSNEW)L=1000*((IABSSS+1)/2)+L       !HYBRID LS: CF+TOT S
            JVR(L)=-1
            IRSOL(I)=-1                               !SAFE SINCE MASTER
            IF(NP.GT.NBINM0.AND.TE.GE.EIONMN)GO TO 36          !-TOLII=0
C            IF(TE-TOLII.GE.EIONMN)GO TO 36              !DROP ALL AUTO?
            NRSOL=NRSOL+1
            JVR(L)=NRSOL                            !MAP MASTER TO NRSOL
            QNV(NRSOL)=-QN(M)
            QLV(NRSOL)=QL(M)
            IRSOL(I)=NRSOL
            ITARR(NRSOL)=NP
            LLR(NRSOL)=0                              !CASE CA
            IF(BLS)THEN
              JW=IABSSS*(2*LL(I)+1)
              SSR(NRSOL)=JW                           !WEIGHT
              JJR(NRSOL)=-IABSSS                      !FLAG EXISTS
            ENDIF
            IF(BIC)THEN
              JW=IABSJJ+1
              JJR(NRSOL)=JW                           !WEIGHT
              SSR(NRSOL)=-1                           !FLAG EXISTS
            ENDIF
            IF(BHYBRD)TE=TE*JW
            WNR(NRSOL)=TE                             !-TOLII=0
            LMR(NRSOL)=LMX(K)
            DO J=1,10
              QSR(NRSOL,J)=QSB(K,J)
              QLR(NRSOL,J)=QLB(K,J)
            ENDDO
c      write(6,*)L,NRSOL,kk,n
          ENDIF
          if(bca)GO TO 36   !have all, so (INCASE NECOR, THEN BRSLP=.T.)
C
C GENERAL SET-UP OF PARENTS AND RESOLVED FINAL RECOMBINED STATES
C
          IF(.NOT.BRSLP)GO TO 36
C
          L=ICF(K)
          IF(L.GT.0)THEN                        !MASTER (NOT FIRST TIME)
            IF(BRSLF)then
c              N0=ICQT(L,QTT(KK))           !MATCH BY POSITION IN CONFIG
              N0=ICQTG(L,QTTG(L,QTT(KK)),QTTE(KK))!BY EGY ORDER IN GROUP
            endif
C
            IF(N0.GT.0)THEN
              if(itarr(n0).le.0)then
                write(6,*)'nv=',nv,' lv=',lv,' itar=',itarr(n0)
                stop 'master list itarr le 0 ??'
              endif
              ITAR(KK)=ITARR(N0)
              IRSOL(I)=N0
c test energy match
c              tee=abs(wnr(n0)-wnr(1)-energ(i)-ecore+wnr(1))
c              if(tee.gt.0.1*abs(wnr(n0)-wnr(1)))then
c                write(777,*)i,lcf(i),kk,n,energ(i)+ecore-wnr(1)
c     x          ,wnr(n0)-wnr(1),tee
c              endif
            ELSE
              IF(BRSLF)IRSOL(I)=-1                      !NOT WANTED
              ITAR(KK)=9999
c             write(0,*)i,kk,n0
            ENDIF
            GO TO 36
          ELSE
            IF(L.LT.0.AND.BRSLF)then       !INITIALIZE INCASE NOT WANTED
              L=-L
c              ICQT(L,QTT(KK))=0            !MATCH BY POSITION IN CONFIG
              ICQTG(L,QTTG(L,QTT(KK)),QTTE(KK))=0 !BY EGY ORDER IN GROUP
            endif
C
C USE DUMMY SET-UP TO DETERMINE PARENT
C
            IF(BRSLP1.AND.QN(M).EQ.NV.OR.      !RYDBERG
     X         BRSLP2.AND.QN(M).GE.INR1)THEN
              TOLII=TOLI
              J1M=LMX(K)-1
              DO 128 N=1,NRSOLZ
C
                IF(BLS)THEN
                  ijz=JJZ(N)
                  IF(INIT.EQ.0.AND.ijz.NE.0)GO TO 128       !ALREADY GOT
                  IF(INIT.NE.0.AND.ijz.NE.QN(M))GO TO 128       !WRONG N
                  IF(IABSSS.NE.SSZ(N))GO TO 128
                  IF(LL(I).NE.LLZ(N))GO TO 128
                ENDIF
                IF(BIC)THEN
                  isz=SSZ(N)                           !for compiler bug
                  IF(INIT.EQ.0.AND.isz.NE.0)GO TO 128       !ALREADY GOT
                  IF(INIT.NE.0.AND.isz.NE.QN(M))GO TO 128       !WRONG N
                  IF(JJ(I).NE.JJZ(N))GO TO 128
                ENDIF
C
                IF(LMX(K).NE.LMZ(N))GO TO 128
                IF(QL(QLB(K,J1)).NE.QLZ(N,J1))GO TO 128
                DO J=1,J1M
                  IF(QSB(K,J).NE.QSZ(N,J))GO TO 128
                  I1=QLB(K,J)
                  I2=QLZ(N,J)
                  IF(QL(I1).NE.QL(I2))GO TO 128
                  IF(QN(I1).NE.QN(I2))GO TO 128
                ENDDO
C
                ITAR(KK)=ITARZ(N)
                IF(BLS)JJZ(N)=QN(M)+1        !TAG
                IF(BIC)SSZ(N)=QN(M)+1        !TAG
                GO TO 36
 128          CONTINUE
C
              TE=ENERG(I)+ECORE
              IF(TE-TOLII.LT.EIONMN+TOLR.and.BRSLF)THEN
                IF(INIT.EQ.0.AND.NRB.GE.0)THEN
                  WRITE(6,775)I
 775              FORMAT(' ERROR IN CROSSJ: THERE IS A RESOLVED FINAL',
     X          ' STATE I=',I5,', WITH A PARENT NOT SPECIFIED BY NTAR2.'
     X            //' POSSIBLE CAUSES:'
     X            //'   1/ DIMENSIONS, INCREASE NTAR2 AND/OR NDIM0.'
     X            //'   2/ CHECK ORDER OF TARGET SYMMETRIES IN adasin.'
     X            //'   3/ ELSE SEND DEBUGE FILE TO NRB.')
                  OPEN(99,FILE='DEBUGE')
                  WRITE(99,776)NV,LV,I
 776              FORMAT('  NV=',I5,'  LV=',I5,'  IFAIL=',I5) 
                  IF(BLS)WRITE(99,777)
 777        FORMAT('    I    T      2S+1    L        CF     (EI-E1)/RY'
     X            ,'  TAR')
                  IF(BIC)WRITE(99,778) 
 778        FORMAT('    K   LV    T 2S+1    L   2J   CF     (EK-E1)/RY'
     X            ,'  TAR')
C
                  DO J=1,NENG
                    kk=iabs(ik(j))
                    WRITE(99,779)J,IK(J),IT(J),SS(J),LL(J),JJ(J)
     X                          ,LCF(J),ENERG(J),itar(kk)
 779                FORMAT(7I5,F15.6,i5)
                  ENDDO
                  CLOSE(99)
                  STOP 'ERROR: PARENT NOT FOUND: SEE adasout FOR CAUSES'
                ELSE
                  WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR I=',I
                ENDIF
              ENDIF
c              write(0,*)i,kk
              ITAR(KK)=9999
              IF(NRB.GE.0)GO TO 36
            ENDIF
          ENDIF
c
          if(nbin0.gt.0)stop 'nbin0.gt.0'                 !remove evntly
c            itar(kk)=9999
c            go to 36
c          endif
C
C.....ELSE BY ENERGY ALONE: CASE MASTER FIRST TIME THRU OR LV=-1
          TOLII=ZERO
C
C TAKE-OFF ONE-ELECTRON ENERGY 
C ONE-BODY REL TERMS DON'T SEEM TO HELP....SO COMMENT-OUT FOR SPEED
C
          TE0=ENERG(I)+ECORE
          IF(BLS)TE=TE0-QDT(QD0,NZ0,NE,QN(M),QL(M),0)
          IF(BIC)TE=TE0-QDT(QD0,NZ0,NE,QN(M),QL(M),0)    !,-1) KAPPA AVE
CREL      IF(BIC)TE=TE0-QDT(QD0,NZ0,NE,QN(M),QL(M),-QL(M)-1)!TRY J=L+0.5
          IKUN=1
          DO J=1,NBINRM
            IF(TE.GE.EI(J).AND.TE.LT.EI(J+1))THEN
              ITEST=J
CREL              IF(BIC)THEN                      !COMPARE WITH J=L-0.5
CREL                KAPPA=QL(M)-1/(QL(M)+1)
CREL                TE=TE0-QDT(QD0,NZ0,NE,QN(M),QL(M),KAPPA) 
CREL                DO L=1,NBINRM
CREL                  IF(TE.GE.EI(L).AND.TE.LT.EI(L+1))THEN
CREL                    IF(ABS(TE-EII(J)).LT.ABS(TE-EII(L)))GO TO 40
CREL                    ITEST=L
CREL                    GO TO 40
CREL                  ENDIF
CREL                ENDDO
CREL              ENDIF
CREL  40          CONTINUE
c      write(6,*)-nv,lv,i,itest
              ITAR(KK)=ITEST
              GO TO 715
            ENDIF
          ENDDO
C
          IF(NR1.EQ.0)THEN
            IF(TE0-TOLII.LT.EIONMN+TOLR)THEN
              WRITE(6,775)I
              STOP 'ERROR: INCREASE NTAR2 AND/OR NDIM0'
            ENDIF
CELSE MASTER SO DO NOT NEED PARENT
          ENDIF
C
          ITAR(KK)=9999
c      write(6,*)-nv,lv,i,te,ei(1),ei(nbinrm)
          GO TO 36
C
C CHECK LS/J CONSISTANCY WITH PARENT (WE CAN SKIP THIS...)
C
 715      J=ITEST
          icount=0                                     !for compiler bug
C
c      if(i.eq.59)write(0,*)'Hello world 64' !,itest
C SKIP CONFIGURATION CHECK IF NOT GOOD QUANTUM NUMBER
C
c 717      if(itest.ne.j)then
c             write(0,*)-nv,lv,i,itest,j
c          endif
 717      CONTINUE
c       if(i.eq.59)write(0,*)'****Hello world 62',itest,j,imv
          IF(NRB.GE.0)THEN
            IF(QSB(K,J1).NE.LIT(1))THEN
              J1M=J1
            ELSE
              J1M=J1-1
            ENDIF
            IF(J1M.NE.LMP(J))GO TO 718
            DO L=1,J1-1
              IF(QSB(K,L).NE.QSP(J,L))GO TO 718
              IF(QLB(K,L).NE.QLP(J,L))GO TO 718
            ENDDO
            IF(QSB(K,J1).NE.LIT(1))THEN
              IF(QLB(K,J1).NE.QLP(J,J1M))GO TO 718
            ENDIF
          ENDIF
          IF(QN(M).EQ.NV.AND.LV.GE.0)THEN               !INCASE NR1.EQ.0
            TOLII=TOLI
          ELSE
            TOLII=ZERO
          ENDIF
c
          IF(BLS)THEN
            IF(IABSSS-1.NE.IWS(J).AND.IABSSS+1.NE.IWS(J))GO TO 718
            IF(LL(I).GT.IWL(J)+QL(M).OR.LL(I).LT.
     X         IABS(QL(M)-IWL(J)))GO TO 718
          ENDIF
          IF(BIC)THEN
            IF(IABSJJ.GT.IWJ(J)+2*QL(M)+1)GO TO 718
            JMIN=MIN(IABS(IABS(2*QL(M)-1)-IWJ(J))          !USE QL SINCE
     X              ,IABS(IABS(IWJ(J)-1)-2*QL(M)))        !DON'T HAVE JV
            IF(IABSJJ.LT.JMIN)GO TO 718
          ENDIF
c      if(itest.ne.j)write(6,*)-nv,lv,i,itest,j
          ITEST=J
          ITAR(KK)=ITEST
          GO TO 36
c 718            if(i.eq.59)write(0,*)'Hello world 64',itest,j

 718      IF(ITEST-J)763,764,765
 765      J=ITEST-IMV
          IF(J.GT.0.AND.J.LE.NBINRM)GO TO 717
 763      IKUN=IKUN+1
          IF(IKUN.GT.IKUN0)THEN
            IF(TE0-TOLII.GE.EIONMN+TOLR.OR.QN(M).NE.NV.OR.LV.LT.0)THEN
c      write(0,*)i,kk
              ITAR(KK)=9999
              GO TO 36
            ENDIF
            WRITE(6,745)I,ITEST,J
 745  FORMAT(' ERROR IN SR.CROSSJ, RESOLVED PARENT NOT FOUND, I=',3I4/
     X,' TRY INCREASING NTAR2...?')
            STOP 'ERROR IN SR.CROSSJ, RESOLVED PARENT NOT FOUND'
          ENDIF
 764      IMV=-1*IKUN
          if(icount.gt.3*ikun0)go to 763               !for compiler bug
c       if(i.eq.59)write(0,*)'****Hello world 64',itest,j,imv,icount
C         IF(TE.GT.EII(J))IMV=1
          IF(ITEST.LE.1*IKUN)IMV=1*IKUN
          icount=icount+1                              !for compiler bug
          J=ITEST+IMV
c       if(i.eq.59)write(0,*)'****Hello world 65',itest,j,imv,icount
          IF(J.LT.1.OR.J.GT.NBINRM)GO TO 763
c       if(i.eq.59)write(0,*)'****Hello world 66',itest,j,imv,nbinrm
          GO TO 717
        ENDIF
C
  36  CONTINUE                                    !END ENERGY LOOP
c                write(6,*)i,kk,itar(kk),qn(m),ql(m),energ(i)
c                enddo
C
      DO I=1,NCF
        ICF(I)=IABS(ICF(I))
      ENDDO
C
      IF(.NOT.BRSLP)GO TO 100                     !SKIP RESOLVED PARENTS
C
C RESET (NOTE, DOES NOT WORK FOR CASE INIT.NE.0 - NEEDS RYD L)
C
      IF(BLS)THEN
        DO N=1,NRSOLZ
          JJZ(N)=0
        ENDDO
      ENDIF
      IF(BIC)THEN
        DO N=1,NRSOLZ
          SSZ(N)=0
        ENDDO
      ENDIF
C
      IF(.NOT.BRSLF)GO TO 100                      !SKIP RESOLVED STATES
C
C****************************************
C SET-UP INDEXING FOR NEW RESOLVED STATES
C****************************************
C
      DO 127 I=1,NENG
C
      IF(LCF(I).LT.0)GO TO 127
      IF(IRSOL(I).NE.0)GO TO 127             !MASTER, BUT NOT FIRST TIME
      IF(IK(I).LE.0)GO TO 127
      TE=ENERG(I)+ECORE
      IF(TE-TOLI.GE.EIONMN+TOLR)GO TO 705
C
      KK=IABS(IK(I))
      K=LCF(I)
      J1=LMX(K)
      M=QLB(K,J1)
      IF(QN(M).LE.NRSLMX)THEN
        ITEST=0
        IF(QN(M).EQ.NV.AND.LV.GE.0)THEN
          TOLII=TOLI
          IF(NR1.NE.0)THEN    !ALREADY SET-UP
            ITEST=ITAR(KK)
            GO TO 707
          ENDIF
          IF(ITAR(KK).GT.NBINM.AND.TE-TOLII.GE.EIONMN)GO TO 127
          NRSOL=NRSOL+1
          IF(NRSOL.GT.NDIM17)THEN
            WRITE(6,379)
            STOP'ERROR: SR.CROSSJ TOO MANY RESOLVED RESONANCES, INCREASE
     X NDIM17'
          ENDIF
          QNV(NRSOL)=NV
          QLV(NRSOL)=QL(M)
          GO TO 130
        ELSE
          TOLII=ZERO
        ENDIF
C
        IF(BPASS1)GO TO 709      !FIRST BLOCK, SO CANNOT EXIST
C
  707   IF(NR1.NE.0.AND.QN(M).GE.INR1)THEN     !CORE RAD SO MATCH LEVELS
C
          if(nrsol.eq.0)then
            write(0,*)'resolved state table not set-up!'
            stop 'resolved state table not set-up!'
          endif
c
          DO 126 N=1,NRSOL
C
            IF(BLS)THEN
              IF(JJR(N).NE.0)GO TO 126                  !ALREADY MATCHED
              IF(IABS(SS(I)).NE.SSR(N))GO TO 126
              IF(LL(I).NE.LLR(N))GO TO 126
            ENDIF
C
            IF(BIC)THEN
              IF(SSR(N).NE.0)GO TO 126                  !ALREADY MATCHED
              IF(JJ(I).NE.JJR(N))GO TO 126
            ENDIF
C
            IF(LMX(K).NE.LMR(N))GO TO 126
            ITARRN=IABS(ITARR(N))
            DO J=1,J1
              IF(QSB(K,J).NE.QSR(N,J))GO TO 126
              IF(ITARRN.GT.0)THEN
                I1=QLB(K,J)
                I2=QLR(N,J)
                IF(QN(I1).NE.QN(I2))GO TO 126
                IF(QL(I1).NE.QL(I2))GO TO 126
              ELSE
                IF(QLB(K,J).NE.QLR(N,J))GO TO 126
              ENDIF
            ENDDO
C
            IF(ITARRN.NE.ITEST) THEN
c              write(99,*)nv,lv,n,i,ik(i),itarr(n),itest
              ITAR(KK)=ITARRN                            !CORRECT PARENT
            ENDIF
C
            IF(BLS)JJR(N)=-1          !TAG
            IF(BIC)THEN
              SSR(N)=IABS(SS(I))      !TAG
              LLR(N)=LL(I)
            ENDIF
C
            T0=ENERG(I)+ECORE-TOLII
            IF(TE-TOLII.GE.EIONMN.AND.ITARRN.GT.NBINM)ITARR(N)=-ITARRN
            IRSOL(I)=N
            IF(ITARRN.EQ.0)STOP 'ERROR: ITARR EQ 0 ????'
            WNR(N)=T0+TOLII    !omit from print as maybe large, but safe
c              write(6,*)nv,lv,n,i,ik(i),itarr(n),wnr(n),tolii
            GO TO 127
 126      CONTINUE
C
C
c          write(99,*)nv,lv,i,ss(i),ll(i),jv(i),jj(i),itest
          IF(TE-TOLII.LT.EIONMN)THEN
            if(nrb.lt.0)go to 709
            WRITE(6,387)I
 387    FORMAT(' ERROR: UNABLE TO FIND PRE-EXISTING RESOLVED VALENCE'
     X ,' STATE',I6,' TRY INCREASING NTAR2....?')
        STOP 'ERROR: UNABLE TO FIND PRE-EXISTING RESOLVED VALENCE STATE'
          ELSE
            WRITE(6,386)I
 386  FORMAT(' WARNING: UNABLE TO FIND PRE-EXISTING AUTOIONIZING'
     X,' (METASTABLE?) FINAL STATE',I6)
            GO TO 127
          ENDIF
C
        ELSE                                      !VALENCE RAD INTO CORE
C
C TEST
        IF(ICF(K).EQ.0)THEN
          WRITE(6,*)'MASTER LIST PROBLEM FOR CF=',K
          STOP 'ERROR: MASTER LIST PROBLEM?'
        ENDIF
C
C MAY NOT BE IN FIRST BLOCK, NO REASON TO STOP
C        STOP 'ERROR: UNABLE TO FIND PRE-EXISTING RESOLVED CORE STATE'
C
        ENDIF
C
C DOES NOT EXIST SO EXTEND LIST
C
 709    IF(ITAR(KK).GT.NBINM.AND.TE-TOLII.GE.EIONMN)GO TO 127
        NRSOL=NRSOL+1
        QNV(NRSOL)=-QN(M)
        QLV(NRSOL)=QL(M)
 130    IF(NRSOL.GT.NDIM17)THEN
          WRITE(6,379)
          STOP
     X  'ERROR: SR.CROSSJ TOO MANY RESOLVED RESONANCES, INCREASE NDIM17'
        ENDIF
C
        L=ICF(K)
        IF(L.GT.0)ICQT(L,QTT(KK))=NRSOL                          !MASTER
        IF(L.GT.0)ICQTG(L,QTTG(L,QTT(KK)),QTTE(KK))=NRSOL        !MASTER
C
        IRSOL(I)=NRSOL
c        write(6,*)nv,lv,i,itar(kk),nrsol
        ITARR(NRSOL)=ITAR(KK)
        SSR(NRSOL)=IABS(SS(I))     !TAG
        LLR(NRSOL)=LL(I)
        JVR(NRSOL)=-999            !FORCE SKIP OF W IF NOT ALREADY DONE
        IF(BLS)JJR(NRSOL)=-1       !TAG
        IF(BIC)JJR(NRSOL)=JJ(I)
        WNR(NRSOL)=ENERG(I)+ECORE-TOLII
        LMR(NRSOL)=LMX(K)
        DO J=1,10
          QSR(NRSOL,J)=QSB(K,J)
          QLR(NRSOL,J)=QLB(K,J)
        ENDDO
        IF(ITEST.GT.0)THEN
          QLR(NRSOL,LMR(NRSOL))=(NV*(NV-1))/2+QL(M)+1+MXORB0
          if(lv.ne.ql(m).and.lv.ne.999)stop 'lv .ne. ql'
        ENDIF
      ENDIF
C
 127  CONTINUE
C
C**************************************************************
C DETERMINE VALENCE J FOR RESOLVED OUTER ELECTRON STABILIZATION
C**************************************************************
C
 705  IF(BRAD.AND.BIC.AND.LV.GE.0)THEN
        DO 818 I=1,NENG
          JV(I)=-1
          IF(LCF(I).LT.0)GO TO 818                       !SKIP CONTINUUM
          K=LCF(I)
          J1=LMX(K)
          M=QLB(K,J1)
          KK=IABS(IK(I))
C
C ONLY NEED VALENCE STATES
C
          IF(QN(M).EQ.NV.AND.QN(M).GE.INR1.AND.ITAR(KK).LT.9999)THEN
            LVV=QL(M)                  !CASE LV=999
C
C SEE IF THE LEVEL HAS OCCURRED BEFORE (SAME PARENT, LV AND JJ),
C IF SO THEN IT IS (PROBABLY) 2*LV-1 AND SO WE NOW HAVE 2*LV+1  
C
            DO 817 M=1,I-1
              IF(LCF(M).NE.LCF(I))GO TO 817
              IF(JJ(M).NE.JJ(I))GO TO 817
              IF(ITAR(KK).NE.ITAR(IABS(IK(M))))GO TO 817
C
C ASSUME THIS M IS 2*LV-1 AND SO
C
              JV(I)=2*LVV+1
c        write(0,*)i,m,jv(i),'+'
              GO TO 818
  817       CONTINUE
C
C NOT OCCURRED BEFORE AND SO ASSUME 2*LV-1 UNLESS IT CANNOT FORM THE JJ
C 
            JMIN=IABS(IWJ(ITAR(KK))-2*LVV+1)
            JMAX=IABS(IWJ(ITAR(KK))+2*LVV-1)
            IF(IABS(JJ(I)).LT.JMIN.OR.IABS(JJ(I)).GT.JMAX)THEN
              JV(I)=2*LVV+1
            ELSE
              JV(I)=IABS(2*LVV-1)
            ENDIF
c        write(0,*)i,jv(i)
          ENDIF
  818   CONTINUE
      ENDIF
C
C*********************
C READ RADIATIVE RATES
C*********************
C                                        !HYBRID AND NBIN0.GT.0 RE-ENTRY
 100  BFAST=.NOT.BFORM.AND.RMIN.LT.ZERO.AND..NOT.BRADBF 
     X      .AND.ABS(RABS-DONE).LT.TINY  !ALLOW TEST FOR UNRESOLVED DATA
C
      IF(BFORM)READ(MR,104,END=1002)NZTEST
      IF(.NOT.BFORM)READ(MRU,END=1002)NZTEST,NDUME
 104  FORMAT(66X,I2)
C
      IF(.NOT.BRADBF)IRAD=NENG
      IF(NZTEST.LT.1)THEN
        I=1
      ELSE
C
        IF(BFORM)READ(MR,103,END=1002)
C
        I=0
 131    I=I+1
        IF(BFORM)READ(MR,132,END=1002) ICR,I1,IWR,I2,I3,JWR,T1,DEL,T2
        IF(.NOT.BFORM)READ(MRU,END=1002)ICR,I1,IWR,I2,I3,JWR,T1,DEL,T2
 132    FORMAT(6I5,1PE15.5,2(0PF15.6))
C
        IF(I1.EQ.0)GO TO 133
        IF(I.LT.NDIM7)THEN
          ITR(I)=I1
          JCR=I2
          IF(I3.EQ.0)THEN                                !A.S. IS HYBRID
            JTR(I)=-I2                                 !MUST HAVE CF NO.
          ELSE
            JTR(I)=I3
          ENDIF
          AR(I)=T1
          EATOM(I)=T2-EIONMN           !NOW RELATIVE, AVOID CANCELLATION
C
          AR(I)=ABS(AR(I))
          IF(BFAST)GO TO 131
C
          BUNR1=BUNR1.OR.I3.EQ.0
          BUNR2=BUNR2.OR.I2.EQ.0
C 
CR          AR(I)=AR(I)*RABS                             !DONE LATER NOW
          IF(AR(I).LT.RMIN)I=I-1
          IF(JTR(I).GT.0)THEN
            IF(JK(JTR(I)).GT.IRAD)I=I-1
          ENDIF
        ENDIF
        GO TO 131
      ENDIF
C
 133  NUMR=I-1
C
      IF(.NOT.BINT.AND.NUMR.GE.NDIM7) THEN
        WRITE(6,74)NUMR
  74    FORMAT(' SR.CROSSJ: NUMBER OF RADIATIVE RATES EXCEEDS STORAGE,'
     X        ,' INCREASE NDIM7 TO',I10)
        STOP        'ERROR: SR.CROSSJ: NUMBER OF RADIATIVE RATES EXCEEDS
     X STORAGE, INCREASE NDIM7'
      ELSE
        NUMRX=MAX(NUMRX,NUMR)
      ENDIF
C
      IF(NBIN0.LE.0)THEN
        IF(BUNR1)THEN
          IF(.NOT.BHYBRD)IFLAGR=1
        ELSEIF(BUNR2)THEN
          IFLAGR=2
        ENDIF
        IF(IFLAGR.GT.0)GO TO 1000
      ENDIF
C
C****************************************************************
C EXIT POINT IF THIS NV,LV IS NOT REQUIRED, GO AND READ NEW NV,LV
c (also, re-entry point for use of scaled rydberg data.)
C****************************************************************
 124  IF(BINT)GO TO 310
C
C
C*********************************************************
C EVALUATE HYDROGENIC RADIATIVE RATES FOR VALENCE ELECTRON
C*********************************************************
C (CAN'T USE WITH LV=999 AS ONLY STORE SINGLE UPPER LV)
C
      RSUMD=ZERO
      IF(.NOT.BRAD)GO TO 160
      IF(LV.LT.0.OR.LV.EQ.999)GO TO 160
      NMIN=MAX0(NR1,LV)
      NMAX=MIN0(NR2,NV-1)
      IF(NMIN.GT.NMAX)GO TO 160
      LP=LV+1
      TL=LV
      TLP=LP
C
      DO I=NR1,IB-1
        N=IBN(I)
        IF(N.GT.NMAX)GO TO 159
        IF(N.GE.NMIN)THEN
          T=N*N
          DE=DZ*(TV-T)/(TV*T)
            CALL DIPOL(-1,N,NV,ZERO,LP,CP,CM,JDUM)
          T1=TLP*CM(LP)*1.0D10**JDUM(LP)
          T2=ZERO
          IF(LV.GT.0)T2=TL*CP(LV)*1.0D10**JDUM(LV)
          T=(T1+T2)/(TL+TLP)
          T0=RWT(IB)*DE**3*2.6775D9/DZ
          T=T*T0
          RSUM(I)=T
          IF(BPRNT2)WRITE(6,870)N,T
 870      FORMAT(I10,1PE15.4)
          T=T/RWT(IB)
          RSUMC(I)=T                                         !UNWEIGHTED
          IF(NBIN0.eq.0)T=T*RWT(I)   !.GE. rwt=1 for >0
          RSUMD=RSUMD+T
          IF(NBIN0.LE.0)THEN     !N=I
            RPS(I)=T1*T0   !LV->LV+1
            RMS(I)=T2*T0   !LV->LV-1
            IF(N.LE.NLMAX)THEN
              RPSL(I)=RPS(I)/(2*LV+1)
              RMSL(I)=RMS(I)/(2*LV+1)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
 159  IF(NMAX.LT.NV-1.AND.NBIN0.GE.0)THEN
        TT=NMAX
        TT=TT**3*(EMN3(NMAX)-EMN3(NV-1))
        RSUMD=RSUMD+T*TT
      ENDIF
C
C N,L SPECIFIC DETECTION PROBABILITIES
C
 160  IF(NFNLMX.GT.0)THEN
        IF(NV.GT.0.AND.NV.LE.NFNLMX)THEN
          IFNL=(NV*(NV-1))/2+LV+1
          FNLV=FNL(IFNL)
C          WRITE (*,*) NV,LV,FNLV
        ELSE
          FNLV=ZERO
        ENDIF
      ELSE
        FNLV=DONE
      ENDIF
C
      IF(NUMA.EQ.0)GO TO 310
      IF(NUMR.EQ.0.AND..NOT.BRAD)GO TO 310
C
C************************************************
C APPLY ENERGY CORRECTIONS TO AUTOIONIZING STATES
C************************************************
C
      TC1=E1C(1)
      IF(NECOR.NE.0)THEN
        IF(LV.GE.0.OR.INIT.NE.0)THEN
          i1=0
          if(-den*.9d0.lt.e1c(1))i1=nbinc
          TC1=ZERO
          J=1
          DO I=1,NUMA
            M=JTA(I)
            IF(M.GT.0)M=ITAR(M)
            IF(M.GT.0.AND.M.LE.NBINC)THEN
              IF(NECOR.LT.0)J=ITAR(ITA(I))
              IF(J.GT.0)THEN
                if(j.le.m-i1.and.necor.lt.0)then
                  write(36,*)nv,lv,j,m,jk(ita(i)),jk(jta(i)),-den
c                  stop 'autoionizing parents inverted'
                endif
                J=MIN(J,NBINC)           !J>M   so uses last specified
CABOVE          EC(I)=EC(I)-E1C(1)
c     write(6,*)'***',nv,lv,i,j,m,jk(ita(i)),jk(jta(i)),ec(i),ecori(m,j)
                EC(I)=EC(I)-ECORI(M,J)
                IF(EC(I).LT.ZERO)THEN
                  AA(I)=ZERO
                  IF(M.EQ.1)THEN
                    IF(ITAG(ITA(I)).GT.0)THEN
                      WRITE(6,149)JK(ITA(I))
  149                 FORMAT(' ERROR: CONFUSION OVER WHETHER ENERGY'
     X      ,' LEVEL',I4,' IS BOUND OR AUTOIONIZING')
      STOP 'ERROR:CONFUSION OVER WHETHER LEVEL IS BOUND OR AUTOIONIZING'
                    ELSE
                      ITAG(ITA(I))=-1                             !BOUND
                    ENDIF
                  ENDIF
                ELSE
                  KK=JK(ITA(I))
                  IF(ENERG(KK)+ECORE0.LT.EIONMN)THEN
                    IF(ITAG(ITA(I)).LT.0)THEN
                      WRITE(6,150)KK
  150                 FORMAT(' ERROR: CONFUSION OVER WHETHER ENERGY'
     X     ,' LEVEL',I4,' IS AUTOIONIZING OR BOUND')
      STOP 'ERROR:CONFUSION OVER WHETHER LEVEL IS AUTOIONIZING OR BOUND'
                    ELSE
                      ITAG(ITA(I))=1                       !AUTOIONIZING
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      IF(TC1.NE.ZERO)THEN
        DO I=1,NUMA
CABOVE          EC(I)=EC(I)-TC1
          IF(EC(I).LT.ZERO)AA(I)=ZERO
        ENDDO
      ELSE
        TC1=E1C(1)
      ENDIF
      E1C(1)=ZERO
C
C
C********************************************
C PROCESS UNSORTED RATES: BEGIN BIG LOOP 1410
C********************************************
C
C FIRST RESCALE N FOR CASE NENG=0
C
      TRR=DONE
      TER=ZERO
      IF(NV.GT.NRR)THEN
        T1=NRR
        T2=NV
        TRR=(T1/T2)**3
        TER=DZ/(NRR*NRR)-DZ/(NV*NV)
      ENDIF
C
C LOOP THRU AUTOIONIZING STATES
C
      DO I=1,NUMA
        IF(JFIRST(ITA(I)).EQ.0)JFIRST(ITA(I))=I
        JLAST(ITA(I))=I
      ENDDO
      IF(NENG.GT.0)THEN
        DO I=1,NUMR
          IF(KFIRST(ITR(I)).EQ.0)KFIRST(ITR(I))=I
          KLAST(ITR(I))=I
        ENDDO
      ENDIF
      JS=9999
C
      DO 1410 IA=1,NAUTO
C
      ITAI=IAUTO(IA)
      IF(JPAR.GT.0.AND.JPAR.NE.ITAR(ITAI))GO TO 1410
C
      IF(RCOR.LT.ZERO)JS=ITAR(ITAI) !USE EXP DE ON AR, AS WELL AS |RCOR|
C
      ITT=JK(ITAI)
      IF(BLS)THEN
        IWAJ=IABS(SS(ITT))*(2*LL(ITT)+1)
        IF(LSPI.GT.0.AND.ILSJ(ITT).NE.LSPI)GO TO 1410
      ENDIF
      IF(BIC)THEN
        IWAJ=IABS(JJ(ITT))+1
        IF(J2PI.GT.0.AND.ILSJ(ITT).NE.J2PI)GO TO 1410
      ENDIF
C
      IF(INIT.NE.0.AND.QN(QLB(LCF(ITT),LMX(LCF(ITT)))).GT.NR2)GO TO 1410
C
      TWAJ=IWAJ
      IF(NNCOR*NV.GT.0)TWAJ=TWAJ*ACORN(NV)
      IF(NLCOR.LT.0.AND.LV+1.GT.0.AND.LV.LT.NDIM25)TWAJ=TWAJ*ACORL(LV+1)
C
      DO M=1,NBINM
        SUMAN(M)=ZERO
      ENDDO
      SUMAD=ZERO
c
      if(bpart)then
        do m=0,ncfr
          sumadi(m)=zero
        enddo
      endif
C
C PERFORM RELEVANT SUMS OF AUGER RATES
C
      DO 141 J=JFIRST(ITAI),JLAST(ITAI)
C
        IF(ITAI.NE.ITA(J))GO TO 141
C
        SUMAD=SUMAD+AA(J)
c
        if(bpart)then
          k=-jcai(j)
          ms=icfi(k)                        !maps (n+1)-cf to target cf
          sumadi(ms)=sumadi(ms)+aa(j)       !so or eo, as iwpart defined
        endif
C
        LS=JTA(J)
        IF(LS.le.0)GO TO 141                             !LOSS TERM ONLY
        LS=ITAR(LS)
        IF(LS.LE.0.OR.LS.GT.NBINM)GO TO 141!INITIAL STATE NOT METASTABLE
        IF(EC(J).LT.EMINC.OR.EC(J).GT.EMAXC)GO TO 141
C
        SUMAN(LS)=SUMAN(LS)+AA(J)
        ECA(LS)=EC(J)+d1m10                 !incase formatted gives zero
C
        JTT=JK(JTA(J))
        IF(LCA(-LCF(JTT)).LT.0)ECA(LS)=-ECA(LS) !TAG WRONG CONT CORE EX.
C
  141 CONTINUE
C
      IF(SUMAD.EQ.ZERO)GO TO 1410
C
      LVV=QL(QLB(LCF(ITT),LMX(LCF(ITT))))
c      if(lv.ge.0.and.lv.ne.999.and.lvv.ne.lv)stop 'lvv error'  !test
      NVV=QN(QLB(LCF(ITT),LMX(LCF(ITT))))
c      if(nvv.ne.nv)stop 'nvv error'    !test
      TVV=NVV
      TVV=TVV*TVV
C
      DEN0=QDT(QD0,NZ0,NE,NVV,LVV,0)
C
      IF(SUMAN(1).NE.ZERO.and.eca(1).gt.zero)THEN
        TJ=EION+ECA(1)-DEN0           !OBS. PARENT ENERGY REL. TO GROUND
c      write(6,*)tj,EION,ECA(1)
      ELSE
        TJ=ENERG(ITT)+ECORE0-DEN0                            !UNADJUSTED
c      write(6,*)-tj,ENERG(ITT),ECORE0
      ENDIF
C
C REPRESENTATIVE AUGER RATES
C (ALSO CHECK CONSISTANCY BETWEEN AUTOIONIZING PARENT AND CONTINUUM)
C
      IF(BHYBRD)THEN
        MS0=ITARH(LCF(ITT))
      ELSE
        MS0=ITAR(ITAI)
      ENDIF
      IF(MS0.GT.0.AND.MS0.LE.NBINRM.AND.LV.LE.LVAMX.AND.NBIN0.LE.0)THEN
        DO LS=1,NBINM0
          IF(MS0.LE.LS.AND.SUMAN(LS).NE.ZERO)THEN
            IF(LV.LT.0)THEN        !PARENTAGE NOT GOOD SO
              MS0=LS
              ITAR(ITAI)=MS0
            ELSE
              WRITE(6,147)MS0,LS,ITT
  147         FORMAT(' ERROR: AUTOIONIZING PARENT',I3,' BELOW CONTINUUM'
     X          ,I3,' FOR LEVEL/TERM=',I6/' TRY INCREASING NTAR2.....?')
c              STOP' ERROR: AUTOIONIZING PARENT BELOW CONTINUUM'
            ENDIF
          ENDIF
          T=LVV
C         IF(LV.LT.0.OR.LV.EQ.999)T=LL(ITT)
          T=2*T+1
          T=T/TVV
          K=IMX(LS,MS0)
          IF(K.GT.0)THEN
            IBB=IB
            IF(INIT.NE.0)IBB=NVV
            AAN(K,IBB)=AAN(K,IBB)+SUMAN(LS)*T
            IF(NVV.LE.NLMAX.AND.LVV.LE.LLMAX)
     X      AANL(K,LVV+1,IBB)=AANL(K,LVV+1,IBB)+SUMAN(LS)
            IF(IRSOL(ITT).GT.0)AANLJ(IRSOL(ITT),LS)=SUMAN(LS)
          ENDIF
        ENDDO
      ENDIF
C
      IF(BNOT.OR.MS0.EQ.9999)MS0=0                         !NO OUTER RAD
C
C ZEROISE FOR SUM OVER RADIATIVE RATES
C
      SUMRN=ZERO
      SUMRD=RSUMD
      SUMRT=ZERO
C
      IF(NBIN0.LE.0)THEN
        DO M=1,NBINRM
          DO K=1,IBN(IB00)
            SUMRNN(K,M)=ZERO
            SUMRN0(K,M)=ZERO
          ENDDO
          SUMRNN(NDIM25,M)=ZERO
          SUMRN0(NDIM25,M)=ZERO
        ENDDO
        DO K=1,IB
          SUMRJ(K)=ZERO
          SUMRJ0(K)=ZERO
        ENDDO        
        IF(NLMAX.GT.0)THEN
          DO M=1,NBINRM
            DO K=1,IBN(IB00)
              ML1=MIN(K,LLMAX+1)
              DO ML=1,ML1
                SUMRNL(ML,K,M)=ZERO
              ENDDO
            ENDDO
            DO ML=1,LLMAX+1
              SUMRNL(ML,NDIM25,M)=ZERO
            ENDDO
          ENDDO
          DO K=1,IB
            SUMRLP(K)=ZERO
            SUMRLM(K)=ZERO
          ENDDO
        ENDIF
      ELSE
        IF(.NOT.BNOT)MS0=1                  !ANY PP RAD ALWAYS IN KMAX=1
        DO K=1,IBN(IB00)
          SUMRN0(K,1)=ZERO
        ENDDO
        SUMRN0(NDIM25,1)=ZERO
        DO K=1,IB
          SUMRJ0(K)=ZERO
        ENDDO                
      ENDIF
C
C PP'D OUTER ELECTRON RADIATN (DETERMINE FINAL STABLE STATES INTERNALLY)
C
      IF(MS0.GT.0)THEN
        IF(NMIN.GT.NMAX)GO TO 333
        TT=TOLR
        IF(MS0.GT.NBINRM0)TT=ZERO
        LVV0=MAX(LVV,1)
C
        DO K=NR1,IB-1
          N=IBN(K)
          IF(N.GT.NMAX)GO TO 333
C
          IF(RSUM(K).GT.ZERO)THEN
            TN=QDT(QD0,NZ0,NE,N,LVV0,0)
            T=TN
            IF(NBIN0.LE.0.AND.N.LE.NLMAX)THEN
              DENP=QDT(QD0,NZ0,NE,N,LVV+1,0)
              IF(LVV.GT.0)THEN
                DENM=QDT(QD0,NZ0,NE,N,LVV-1,0)
                T=DENM
              ENDIF
            ENDIF
            IF((TJ+T).GE.EIONMN+TT)GO TO 333
            SUMRT=SUMRT+RSUM(K)
            IF(TJ+TN.LT.EIONMN)THEN                     !BUNDLE-N ENERGY
              SUMRJ0(K)=RSUMC(K)                             !UNWEIGHTED
              SUMRN=SUMRN+RSUMC(K)
            ENDIF
C
            IF(NBIN0.LE.0)THEN
C PARENT MS0 STILL
              IF(MS0.GT.NBINRM)THEN    !SHOULD'VE BEEN CAUGHT ALREADY 
                WRITE(6,719)ITAI
 719  FORMAT(' WARNING IN SR.CROSSJ, UNRESOLVED RECOMBINATION PARENT'
     X,' IT=',I4,' TRY INCREASING NTAR2/NDIM0')
                STOP ' ERROR: INCREASE NTAR2/NDIM0'
              ENDIF
              if(n.gt.ncmx0)NRX=MAX0(NRX,IABS(MS0))
              IF((TJ+TN).LT.EIONMN+TT)SUMRJ(K)=RSUM(K)  !BUNDLE-N ENERGY
c         write(6,*)ms0,k,tj,t,tj+t,eionmn
              IF(N.LE.NLMAX)THEN
                IF((TJ+DENP).LT.EIONMN+TT)SUMRLP(K)=RPSL(K)
                SUMRLM(K)=RMSL(K)             !e test above at go to 333
              ENDIF
            ENDIF
C
          ENDIF
C
        ENDDO
      ENDIF
C
 333  NRSOL0=0
C
      IF(NUMR.EQ.0)GO TO 421
C
C PERFORM RELEVANT SUMS OF RADIATIVE RATES
C
      DO 243 K=KFIRST(ITAI),KLAST(ITAI)
C
      IF(ITR(K).NE.ITAI)GO TO 243                     !CASE >1 MULTIPOLE
C
      TARK=AR(K)
      MN=JTR(K)
      IF(MN.NE.0)THEN
        IF(MN.LT.0)THEN                                      !GET CF NO.
          JCR=-MN
        ELSE
          JCR=LCF(JK(MN))
        ENDIF
        IF(BHYBRD)THEN                            !FLAGS FOR HYBRID CASE
          MS=ITARH(JCR)
          ITAGK=0
        ELSE
          MS=0
          ITAGK=ITAG(MN)
        ENDIF
        MN=QLB(JCR,LMX(JCR))
        IF(QN(MN).NE.NV)THEN
          TARK=AR(K)*TRR
          TE=TER
          TOLII=ZERO
        ELSE
          TE=ZERO
          TOLII=TOLI
          IF(LV.GE.0)TARK=TARK*RABS
        ENDIF
      ELSE
        ITAGK=0
        TE=ZERO
        TOLII=ZERO
        IF(LV.GE.0)TARK=TARK*RABS
      ENDIF
C
      SUMRD=SUMRD+TARK
C
      IF(BCFP.AND.MN.GT.0.AND.JCFR.ge.JCR)THEN         !TEST MAX CASCADE
        ITAGK=-1
        GO TO 41
      ENDIF
      IF(EATOM(K)+TE-TOLII.GE.TOLR) THEN                        !+EIONMN
c        write(6,*)'skip',jtr(k)
        IF(JCFR.GT.100)SUMRD=SUMRD-TARK
        GO TO 243
      ENDIF
C      IF(EATOM(K)+TE-TOLII.GE.TOLR)GO TO 243 !IF JCFR COMMENTED!+EIONMN
  41  CONTINUE
C
      IF(NBIN0.GT.0.OR.MN.EQ.0)THEN
        MS=1
        IF(MN.GT.0)THEN
          MN=QN(MN)
          N=MN
          IF(MN.EQ.NV)N=NDIM25
        ELSE
          N=NDIM25
        ENDIF
        mn=nv                                      !irwt=1 if nbin0.gt.0
        GO TO 242                                           !BINNED ONLY
      ENDIF
C 
C STORE RESOLVED RADIATION
C
      IF(QN(MN).NE.NV)TARK=TARK*RWT(IB)
      IF(MS.EQ.0)THEN
        M=JK(JTR(K))
        IRSOLM=IRSOL(M)
      ELSE                                      !HYBRID
        L=ICF(JCR)
        IF(BLSNEW)L=1000*((IABS(SS(ITT))+1)/2)+L
        IRSOLM=JVR(L)
      ENDIF
      IF(IRSOLM.GT.0)THEN
        IF(MS.EQ.0)MS=ITARR(IRSOLM)
        IF(MS.LT.0)GO TO 243
        IF(MS.GT.NBINM0.AND.(EATOM(K)+TE-TOLII.GE.ZERO.AND.ITAGK !EIONMN
     X                .EQ.0.OR.ITAGK.GT.0)) GO TO 243
        IF(JS.LE.-NECOR.AND.QN(MN).EQ.NV.AND.LV.GE.0)THEN
          IF(JS.GT.MS)THEN
            TRE=(E1X(JS)-E1X(MS))/(E1C(JS)-E1C(MS))
            TRE=ABS(TRE)**3
c            if(abs(tre-done).gt.0.15)write(6,*)js,ms,tre,e1c(1),e1x(1)
            SUMRD=SUMRD-TARK
            TARK=TARK*TRE
            SUMRD=SUMRD+TARK
c         else
c            write(6,*)nv,lv,js,ms,jk(itai),itai,m,jtr(k)
c            stop 'js.eq.ms'
          ENDIF
        ENDIF
        NRSOL0=NRSOL0+1
        AR0(NRSOL0)=TARK
        IRSOL0(NRSOL0)=IRSOLM
      ELSEIF(IRSOLM.LT.0)THEN                           !UNWANTED MASTER
        GO TO 243
C      ELSE                                             !BUNDLED ONLY
      ENDIF
C
C STORE BUNDLED N/NL RADIATION (CORE AND VALENCE)
C
      ML=QL(MN)+1
      MN=QN(MN)
      N=MN
      IF(MN.EQ.NV)N=NDIM25
      IF(MS.EQ.0)THEN
        MS=ITAR(IABS(IK(M)))                      !IK(M)=JTR(K)
        IF(JS.LE.-NECOR.AND.MN.EQ.NV.AND.LV.GE.0)THEN
          IF(JS.GT.MS)THEN
            TRE=(E1X(JS)-E1X(MS))/(E1C(JS)-E1C(MS))
            TRE=ABS(TRE)**3
c            if(abs(tre-done).gt.0.15)write(6,*)js,ms,tre
            SUMRD=SUMRD-TARK
            TARK=TARK*TRE
            SUMRD=SUMRD+TARK
c          else
c            write(6,*)nv,lv,js,ms,jk(itai),itai,m,jtr(k)
c            stop 'js.eq.ms'
          ENDIF
        ENDIF
      ENDIF
      IF(MS.GT.NBINRM)THEN      !SHOULD'VE ALREADY BEEN CAUGHT
        IF(MN.LT.INR1.OR.LV.LT.0.OR.NR1.EQ.0.AND.BRSLP)THEN
          MS=NBINRM                         !RAD INTO CORE, NO PARENTAGE
        ELSE
          WRITE(6,713)JTR(K)
 713      FORMAT(' ERROR IN SR.CROSSJ, FINAL PARENT NOT FOUND, T=',I4
     X           ,'. INCREASE NTAR2/NDIM0')
            STOP 'ERROR IN SR.CROSSJ, FINAL PARENT NOT FOUND'
        ENDIF
      ENDIF
      BTEST=(EATOM(K)+TE-TOLII.GE.ZERO.AND.ITAGK.EQ.0     !EIONMN
     X                                   .OR.ITAGK.GT.0)
      IF(MS.GT.NBINRM0.AND.BTEST) GO TO 243
C                                        HERE N=NDIM25 OR .LE. IBN(IB00)
      if(n.gt.ibn(ib00).and.n.ne.ndim25)then
        write(6,*)'illegal n-value here:',ib00,ibn(ib00),n
        stop 'illegal n-value here'
      endif
      if(ms.lt.0.or.ms.eq.9999)then
        write(6,714)jtr(k)
 714    format(' error in sr.crossj, final parent not found, t=',i4
     x         ,'. contact nrb')
        write(6,*)itai,itt,k,jtr(k),jcr
        stop 'error in crossj, final parent not found - contact nrb'
      endif
c
      SUMRNN(N,MS)=SUMRNN(N,MS)+TARK                           !WEIGHTED
      IF(MN.LE.NLMAX.AND.QL(MN).LE.LLMAX)
     X                SUMRNL(ML,N,MS)=SUMRNL(ML,N,MS)+TARK
      NRX=MAX0(NRX,MS)
C
 242  BTEST=(EATOM(K)+TE-TOLII.LT.ZERO.AND.ITAGK.EQ.0     !EIONMN
     X                                   .OR.ITAGK.LT.0)
      IF(BTEST)THEN
        TARK0=TARK
        IF(MN.NE.NV)TARK0=TARK0/RWT(IB)                      !UNWEIGHTED
        SUMRN0(N,MS)=SUMRN0(N,MS)+TARK0
        SUMRN=SUMRN+TARK0                           !SUM OVER BOUND ONLY
      ENDIF
C
      SUMRT=SUMRT+TARK                               !SUM OVER ALL FINAL
 243  CONTINUE
C 
C EVALUATE RESOLVED RADIATIVE RATES FOR OUTER ELECTRON STABILIZATION
C
 421  IF(MS0.GT.0.AND.NV.GT.NR1.AND.BRSLF)THEN 
C                                   BUT NRSOL=0 IF LAST TWO FALSE...
      DO 480 N=1,NRSOL      !EVALUATE APPROPRIATE RECOUPLING COEFFICIENT
        IF(ITARR(N).NE.MS0)GO TO 480
        I=LMR(N)
        II=QLR(N,I)
        IF(QN(II).GE.NV)GO TO 480
        IF(QN(II).LT.NR1)GO TO 480
        IF(IABS(LVV-QL(II)).NE.1)GO TO 480
        M=ITT
        IF(BLS)THEN
          IF(IABS(SS(M)).NE.SSR(N))GO TO 480
          IF(IABS(LL(M)-LLR(N)).GT.1)GO TO 480
          IF(LL(M)+LLR(N).EQ.0)GO TO 480
          T0=2*LLR(N)+1
          W=WSQ(LVV,QL(II),LL(M),LLR(N),IWL(ITARR(N)),ISIGN)
     X         /(3*(2*IWL(ITARR(N))+1))
        ENDIF
        IF(BIC)THEN
          IF(IABS(IABS(JJ(M))-IABS(JJR(N))).GT.2)GO TO 480
          IF(IABS(JJ(M))+IABS(JJR(N)).EQ.0)GO TO 480
          IF(IABS(JV(M)-JVR(N)).GT.2)GO TO 480
          T0=IABS(JJR(N))+1
          W=WSQ2(JV(M),JVR(N),IABS(JJ(M)),IABS(JJR(N))
     X     ,IWJ(ITARR(N)),ISIGN)/(3*(IWJ(ITARR(N))+1))
c          write(0,*)w,jv(m),jvr(n),jj(m),jjr(n),iwj(itarr(n))
          T0=T0*(JV(M)+1)*(JVR(N)+1)
          W2=WSQ2(2*LVV,2*QL(II),JV(M),JVR(N),1,ISIGN)/6
c          write(0,*)w2,2*lv,2*ql(ii),jv(m),jvr(n)
          T0=T0*W2
        ENDIF
        IF(QL(II).GT.LVV)THEN
          T=RPS(QN(II))
        ELSE
          T=RMS(QN(II))
        ENDIF
        T=T*T0*W
        NRSOL0=NRSOL0+1
        AR0(NRSOL0)=T
        IRSOL0(NRSOL0)=N
 480  CONTINUE
C
      ENDIF
C
C EVALUATE CROSS SECTIONS
C
      IF(SUMRT*FNLV.EQ.ZERO)GO TO 1410
      CROSSD=SUMAD
      IF(JCFR.LE.200)CROSSD=CROSSD+SUMRD       !ADD-IN RADIATION DAMPING
C
      DO 435 MS=1,NBINM
C
      L=MS
      IF(ECA(L).LT.ZERO)GO TO 435                   !BOUND OR WRONG CONT
C
      CROSSN=TWAJ*SUMAN(L)*FNLV                       !FNLV=EXP SURVIVAL
C
C SUMMED
C
      TI=IWT(L)
      CROSS=CROSSN/(CROSSD*TI)
C
C      IF(CROSS.LE.ZERO)GO TO 435              !WHY NOT?
      IF(NBIN0.GT.0)GO TO 430                               !BINNED ONLY
C
C EVALUATE RATE COEFFICIENTS
C
      MSS=IABS(SS(ITT))
      L0=L
C
      DO M=1,JTEMP
        T=-ECA(L)/TEMP(M)
        IF(T.GT.-75.0D0)THEN
          COFT(M)=CROSS*EXP(T)*COEF(M)
        ELSE
          COFT(M)=ZERO
        ENDIF
      ENDDO
C
C BUNDLED (RECALL MS0=0 IF .NOT.BRAD)
      DO K=1,NBINRM
        IF(BBNFP)L0=(L-1)*NBINRM+K
        NSYS=1
        IF(NSYM.EQ.2.AND.MSS.GT.IWS(L).AND.IWS(L).GT.1        !LS NOT CA
     X      .AND.(.NOT.BRSLF.OR.IWS(K).EQ.IWS(L)))NSYS=2
        DO N=1,IBN(IB00)                                  !N=IBN(N) HERE
          T1=SUMRNN(N,K)
          T2=SUMRN0(N,K)
          DO M=1,JTEMP                                 !BUNDLE-N VALENCE
            BN(M,N,NSYS,L0)=BN(M,N,NSYS,L0)+COFT(M)*T1 !INDEX BY LOWER N
            ALFN(M,IB,L)=ALFN(M,IB,L)+COFT(M)*T2      !INDEX BY UPPER NV
          ENDDO
          IF(N.LE.NLMAX)THEN                          !BUNDLE-NL VALENCE
            M2=MIN(N,LLMAX+1)
            DO ML=1,M2
              T=SUMRNL(ML,N,K)
              DO M=1,JTEMP
                BNL(M,ML,N,NSYS,L0)=BNL(M,ML,N,NSYS,L0)+COFT(M)*T
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        T1=SUMRNN(NDIM25,K)
        T2=SUMRN0(NDIM25,K)
        DO M=1,JTEMP
          BN(M,IB,NSYS,L0)=BN(M,IB,NSYS,L0)+COFT(M)*T1    !BUNDLE-N CORE
          ALFN(M,IB,L)=ALFN(M,IB,L)+COFT(M)*T2
        ENDDO
        IF(NVV.LE.NLMAX.AND.LVV.LE.LLMAX)THEN
          T1=SUMRNL(LVV+1,NDIM25,K)
          DO M=1,JTEMP
            BNL(M,LVV+1,IB,NSYS,L0)=BNL(M,LVV+1,IB,NSYS,L0)
     X                                  +COFT(M)*T1      !BUNDLE-NL CORE
          ENDDO
        ENDIF
C
C ADD-IN HYDROGENIC OUTER
        IF(K.EQ.IABS(MS0))THEN                             !PP'D VALENCE
          DO I=NR1,IB-1                                        !BUNDLE-N
            l00=l0
            N=IBN(I)
            IF(N.LE.NLMAX)THEN                                !BUNDLE-NL
              T0=1
              ML=-1
              lm0=l0
              lp0=l0
              lvp=lvv+1
              lvm=lvv-1
              iflagm=0
              iflagp=0
              IF(BEQN.AND.I.EQ.NCMX0)THEN       !ALLOW FOR EQUIVALENT e-
                JH=JKH(K)
                LH=LMH(JH)
                M=QLH(JH,LH)
                IF(NCMX0.EQ.QN(M).and.qsh(jh,lh).eq.lit(1))THEN
                  ML=4*QL(M)
                  T0=(ML+1)*T0/(ML+2)
                  ML=ML/4
                  mssj=0
                  IF(BLSNEW)mssj=1000*((mss+1)/2)
                  do 177 nc=1,ncft   !try & match to any previous master
                    ncj=nc+mssj
                    if(jvr(ncj).lt.0)go to 177          !as autoionizing
                    lt=lmt(nc)
                    m=qlt(nc,lt)
                    if(qn(m).ne.ncmx0)go to 177
                    np=itarr(jvr(ncj))                     !master parent
                    if(qst(nc,lt).eq.lit(2))then
                      if(ql(m).ne.ml)go to 177
                  if(np.ne.k)then
c               write(0,*)'Hello World 0'
                      go to 177
                  endif
                      if(lvv+1.eq.ml)then
                        if(iflagp.gt.0)stop 'iflagp already set'
                        iflagp=ncj
                        if(iflagm.gt.0)go to 178
                      elseif(lvv-1.eq.ml)then
                        if(iflagm.gt.0)stop 'iflagm already set'
                        iflagm=ncj
                        if(iflagp.gt.0)go to 178
                      endif
                    elseif(qst(nc,lt).eq.lit(1).and.lt.gt.1)then
                      mm=qlt(nc,lt-1)
                      if(qn(mm).ne.ncmx0.or.
     x                   qst(nc,lt-1).ne.lit(1))go to 177
                      if(lt-1.ne.lh)go to 177                !match core
                      do j=1,lh-1
                        if(qsh(jh,j).ne.qst(nc,j))go to 177
                        if(qlh(jh,j).ne.qlt(nc,j))go to 177
                      enddo
                      if(ql(m).eq.ml)then
                        if(lvv+1.eq.ql(mm))then                !reversed
                  if(np.eq.k)then
               write(0,*)'Hello World 1'
                  endif
                          if(bbnfp)lp0=l0-k+np
                          lvp=ml
                          if(iflagp.gt.0)stop 'iflagp already set'
                          iflagp=ncj
                          if(iflagm.gt.0)go to 178
                        elseif(lvv-1.eq.ql(mm))then            !reversed
                  if(np.eq.k)then
               write(0,*)'Hello World 2'
                  endif
                          if(bbnfp)lm0=l0-k+np
                          lvm=ml
                          if(iflagm.gt.0)stop 'iflagm already set'
                          iflagm=ncj
                          if(iflagp.gt.0)go to 178
                        endif
                      elseif(ql(mm).eq.ml)then
                        if(lvv+1.eq.ql(m))then
                  if(np.ne.k)then
               write(0,*)'Hello World 3'
                  endif
                          if(iflagp.gt.0)stop 'iflagp already set'
                          iflagp=ncj
                          if(iflagm.gt.0)go to 178
                        elseif(lvv-1.eq.ql(m))then
                  if(np.ne.k)then
               write(0,*)'Hello World 4'
                  endif
                          if(iflagm.gt.0)stop 'iflagm already set'
                          iflagm=ncj
                          if(iflagp.gt.0)go to 178
                        endif
                      endif
                    endif
  177             enddo
  178             continue
                ENDIF
              ENDIF
C BINDLE-NL
              IF(LVV+1.LT.N.AND.LVV+1.LE.LLMAX)THEN
                T=SUMRLP(I)
                IF(LVV+1.EQ.ML)T=T*T0
                DO M=1,JTEMP
                  BNL(M,lvp+1,I,NSYS,lp0)=BNL(M,lvp+1,I,NSYS,lp0)
     x                                                        +COFT(M)*T
                ENDDO
                if(iflagp.gt.0)then
                  i0=jvr(iflagp)
                  if(i0.le.0)then
                    stop 'iflagp autoionizing....???'
                  endif
                  do m=1,jtemp
                    an(m,i0,l)=an(m,i0,l)+coft(m)*t
                  enddo
                endif
              ENDIF
              IF(LVV.GT.0.AND.LVV-1.LT.N.AND.LVV-1.LE.LLMAX)THEN
                T=SUMRLM(I)
                IF(LVV-1.EQ.ML)T=T*T0
                DO M=1,JTEMP
                  BNL(M,lvm+1,I,NSYS,lm0)=BNL(M,lvm+1,I,NSYS,lm0)
     x                                                        +COFT(M)*T
                ENDDO
                if(iflagm.gt.0)then
                  i0=jvr(iflagm)
                  if(i0.le.0)then
                    stop 'iflagm autoionizing....???'
                  endif
                  do m=1,jtemp
                    an(m,i0,l)=an(m,i0,l)+coft(m)*t
                  enddo
                endif
              ENDIF
              l00=min(lm0,lp0)
            ENDIF
C BUNDLE-N
            T1=SUMRJ(I)
            T2=SUMRJ0(I)
            DO M=1,JTEMP
              BN(M,I,NSYS,l00)=BN(M,I,NSYS,l00)+COFT(M)*T1   !BY LOWER N
              ALFN(M,IB,L)=ALFN(M,IB,L)+COFT(M)*T2    !INDEX BY UPPER NV
            ENDDO
C
          ENDDO                                        !END LOWER-N LOOP
        ENDIF                                          !END PP'D VALENCE
C
      ENDDO                       !END LOOP OVER FINAL PARENTS (BUNDLED)
C
C RESOLVED
      IF(NRSOL0.GT.0)THEN
        DO N=1,NRSOL0
          I0=IRSOL0(N)
          T=AR0(N)
          DO M=1,JTEMP
            AN(M,I0,L)=AN(M,I0,L)+COFT(M)*T
          ENDDO
        ENDDO
      ENDIF
C
C BINNED CROSS SECTIONS
C
 430  IF(CROSS.LE.ZERO)GO TO 435                   !NBIN0 .GT.0 RE-ENTRY
      CROSS=CROSS*CON
C
      IF(NBIN0.NE.0.AND.L.LE.LMAX)THEN
        cross0=cross
        IF(ECA(L).LT.EBIN(1).OR.ECA(L).GE.EBIN(NBIN))GO TO 434
        eh=-done
        DO N=1,NBIN1
          IF(ECA(L).GE.EBIN(N).AND.ECA(L).LT.EBIN(N+1))GO TO 431
        ENDDO
        N=NBIN1
 431    CONTINUE
        ERES=EBIN(N+1)-EBIN(N)
c
        nb=0
        if(bpart)then
          el=eca(l)-epart*10
          eh=eca(l)+epart*10
          el=max(el,ebin(1))
          eh=min(eh,ebin(nbin))
          if(el.gt.eh)go to 434
          nb=n
          eres=1
          t0=cross0*crossd
          ibe=iepart(nb+1,l)
          if(bca)then
            icfp=icf0
            if(bprnt2)write(6,158)0,(ic,iwpart(ibe,ic)/w0(ic),ic=1,icf0)
          else
            icfp=0
            do ic=1,icf0
              if(iwpart(1,ic).ne.0)icfp=ic
              w0(ic)=iwpart(ibe,ic)
            enddo
          endif
          if(bprnt2)write(6,158)ibe,(ic,w0(ic)/iwpart(1,ic),ic=1,icfp)
          do n=1,nbin1
            if(el.ge.ebin(n).and.el.lt.ebin(n+1))go to 432
          enddo
          n=nbin1
        endif
c
 432    CONTINUE
c
        if(bpart)then
          if(epart.gt.zero)then
c
            ie=iepart(n+1,l)
            ta=sumadi(0)                  !unresolved target  !dzero
            do ic=1,icf0
              t=sumadi(ic)*iwpart(ie,ic)
              if(t.ne.zero)then
                if(w0(ic).eq.zero)then
                  stop 'w=0'
                endif
                ta=ta+t/w0(ic)
              endif
            enddo
c
            clor=(epart/6.284)/((ebin(n+1)-eca(l))**2+epart*epart/4)
            cross=t0*clor*eres/(dee*ta*clor+sumrd)
c
            clor=(epart/6.284)/((ebin(n+1)+eca(l))**2+epart*epart/4)
            cross=cross+t0*clor*eres/(dee*ta*clor+sumrd)
c
          else
c
c TBD: distribute over autoionizing levels, of same LS/J/p in nl-block
c if want to use redistributed fluorescence yields (need to precalculate
c SUMAN(NYLD),SUMRD(NYLD))...
c
c            cross=cross*ebin(n+1)       !needed as well to cancel below
c
          endif
        endif
C
        CROSS=CROSS/(ERES*EBIN(N+1))
C
        DO K=1,KMAX
          DO I=1,IBN(IB00)
            T=CROSS*SUMRN0(I,K)
            UB(IB,N,L)=UB(IB,N,L)+T                   !INDEX BY UPPER NV
            TNU(IB,L)=TNU(IB,L)+cross0*SUMRN0(I,K)
          ENDDO
          T=CROSS*SUMRN0(NDIM25,K)
          UB(IB,N,L)=UB(IB,N,L)+T
          TNU(IB,L)=TNU(IB,L)+cross0*SUMRN0(NDIM25,K)
          IF(K.EQ.IABS(MS0))THEN                  !PP'D BUNDLE-N VALENCE
            DO I=NR1,IB-1
              T=CROSS*SUMRJ0(I)
              UB(IB,N,L)=UB(IB,N,L)+T                 !INDEX BY UPPER NV
              TNU(IB,L)=TNU(IB,L)+cross0*SUMRJ0(I)
            ENDDO
          ENDIF     
        ENDDO
c
        if(ebin(n+1).lt.eh)then
          n=n+1
          go to 432
        elseif(nb*epart.ne.zero)then
          eres=ebin(nb+1)-ebin(nb)
          cross=cross0/(eres*ebin(nb+1))
c          n=nb
        endif
c
        CROSS=CROSS*ECA(L)
C
      ENDIF
C
 434  CROSS=CROSS*SUMRN/ECA(L)
      TC(L)=TC(L)+CROSS
C
      IF(BPRNT0)THEN            !n
        WRITE(6,32)LCF(ITT),ITAI,L,IWT(L),MS0,IWAJ,ECA(L),SUMAN(L),SUMAD
     X ,SUMRN,SUMRD,CROSS                         
c                  !ub(ib,n,l)*eres
  32    FORMAT(3I5,I3,I4,I3,6(1PE15.4))
      ENDIF
C
 435  CONTINUE                            !END INITIAL TARGET STATE LOOP
C
 1410 CONTINUE                        !END LOOP OVER AUTOIONIZING STATES
C
C************************************************
C END PROCESSING OF UNSORTED RATES: BIG LOOP 1410
C************************************************
C
      E1C(1)=TC1
      DO L=1,NBINM
        TCN(L)=TC(L)*TV3
      ENDDO
      IF(BPRTM1.AND.NV.GE.0)WRITE(6,35)NV,LV,(TC(L),L=1,NBINM)
  35  FORMAT(I5,I3,2X,9(1PE13.4))
      IF(TC(1).EQ.ZERO)NV00=0
C
C  GO AND READ NEW NL BLOCK
C
      BUNIT=.FALSE.
      GO TO 310
C
C  ABORT
 1002 NV=0
      WRITE(6,1107)
 1107 FORMAT(/' ******WARNING, UNEXPECTED END OF DATA IN SR.CROSSJ !!!!'
     X,/)
C      GO TO 1001
C
 1000 CONTINUE
C
      WRITE(6,*)' '
      IF(LV0.GE.0)THEN
        WRITE(6,*)'RYDBERG AUTOIONIZING CONFIGURATIONS CONTRIBUTING:'
      ELSE
        WRITE(6,*)'AUTOIONIZING CONFIGURATIONS CONTRIBUTING:'
      ENDIF
      DO I=1,NCF
        IF(LCA(I).GT.0)WRITE(6,*)I                    !ANY MASTER=0 HERE
      ENDDO
      WRITE(6,*)' '
C
      IF(BUNA)THEN
        WRITE(6,*)' '
        WRITE(6,*)
     X         'NOTE: AS DATA CONTAINS UNRESOLVED AUTOIONIZATION RATES,'
        WRITE(6,*)'      SET NTAR1.LE.NMETAR/J OF THE AS RUN'
        WRITE(6,*)' '
        IF(IFLAGR.EQ.-1)THEN
          WRITE(6,*)
     X        'ADF09 FILE REQUIRES NTAR2.LT.0 AND THE TARGET o_str FILE'
          WRITE(6,*)' '
        ENDIF
      ENDIF
      IF(BUNR2)THEN
        WRITE(6,*)' '
        WRITE(6,*)'NOTE: AS DATA CONTAINS UNRESOLVED RADIATVE RATES'
        WRITE(6,*)' '
        IF(IFLAGR.EQ.2)THEN
          WRITE(6,*)'NO ADF09 POSSIBLE, TOTALS (NBIN.GT.0) ONLY'
          WRITE(6,*)' '
        ENDIF
      ELSEIF(BUNR1)THEN
        WRITE(6,*)' '
        WRITE(6,*)
     X    'NOTE: AS DATA CONTAINS CONFIGURATION RESOLVED RADIATVE RATES'
        WRITE(6,*)' '
        IF(IFLAGR.EQ.1)THEN
          WRITE(6,*)
     X       'ADF09 FILE REQUIRES NTAR2.LT.0 AND THE TARGET o_str FILE'
          WRITE(6,*)' '
        ENDIF
      ENDIF
      IF(IFLAGR.NE.0)THEN
        WRITE(6,*)'***ERROR: RESOLUTION CONFLICT:'
        STOP '***ERROR: RESOLUTION CONFLICT - SEE adasout'
      ENDIF
C
C
C**********************
C  GO AND READ NEW FILE
C**********************
C
      CLOSE(MR)                !MRU=MR
      IF(IFILE.LT.0)GO TO 1001                                  !NO MORE
C
      IFILE=IFILE+1
      IC1=IFILE/10
      IC2=IFILE-10*IC1
      IC0=ICHAR('0')
      IC1=IC1+IC0
      IC2=IC2+IC0
C
      IF(BFORM)THEN
        IF(O1.EQ.'o1')THEN                                       !SERIAL
          FILNAM='o'//CHAR(IC2)
          IF(IFILE.GE.10)FILNAM='o'//CHAR(IC1)//CHAR(IC2)
        ELSE                                                   !PARALLEL
          FILNAM=O//CHAR(IC1)//CHAR(IC2)
        ENDIF
        INQUIRE(FILE=FILNAM,EXIST=EX)
        IF(EX)OPEN(MR,FILE=FILNAM)
      ELSE
        IF(O1U.EQ.'o1u')THEN                                     !SERIAL
          FILNAM='o'//CHAR(IC2)//'u'
          IF(IFILE.GE.10)FILNAM='o'//CHAR(IC1)//CHAR(IC2)//'u'
        ELSE                                                   !PARALLEL
          FILNAM=O//'u'//CHAR(IC1)//CHAR(IC2)
        ENDIF
        INQUIRE(FILE=FILNAM,EXIST=EX)
        IF(EX)OPEN(MRU,FILE=FILNAM,FORM='UNFORMATTED')
      ENDIF
      IF(EX)GO TO 331
C
 1001 CONTINUE
C
      WRITE(6,1007)NUMAX,NUMRX
 1007 FORMAT(/' MAXIMUM USED NUMA=',I10/' MAXIMUM USED NUMR=',I10)
C
      IF(TOLB.GT.TOLB0)WRITE(6,137)TOLB
 137  FORMAT(/' *** ATTN: TOLB HAS BEEN RESET TO =',1PE10.2,' RYD'/)
C
      IF(IFLAGE.NE.0.and.bprnt0)WRITE(6,1006)IFLAGE
 1006 FORMAT(//'NOTE: ',I4,' UNIT5 TARGET ENERGIES DID NOT MATCH WITH'
     X,' THOSE PRESENT IN THE RATE FILE, SEE ABOVE "***" !'/11X
     X,'ENERGY CORRECTIONS ARE BASED ON THOSE IN THE RATE FILE...')
c
      if(bflagp)then
        if(bprnt0)WRITE(0,*)'*** MASTER PARENT(S) NOT DETERMINABLE...'
        WRITE(6,*)'*** MASTER PARENT(S) NOT DETERMINABLE...'
      endif
C
      IF(NFLAG2.LT.NBINM)THEN                           !PARENT PROBLEMS
        WRITE(6,1008)NFLAG2
 1008   FORMAT('*** PRIOR PROBLEMS DECODING TARGET, TRY AND REDUCE'
     X,' NTAR2 TO LAST SAFE PARENT',I5)
       STOP
     X  'ERROR: PARENT PROBLEMS - SEE adasout FOR DETAILS, REDUCE NTAR2'
      ENDIF
C
C********************************************************
C APPLY RELATIVISTIC (JUTTNER) CORRECTION TO DISTRIBUTION
C********************************************************
C
      IF(IREL.NE.0)THEN
C
        NU=2
        MU=4*NU*NU
C
        DO J=1,JTEMP
C
          THETA=TEMP(J)*DALF/2                            !KT(a.u.)/C**2
          T8=THETA/8
C START-OFF WITH ABRAMOWITZ & STEGUN 9.7.2, AVOIDS SMALL THETA OVERFLOW
          KSUM=-5/LOG10(THETA)
          KSUM=KSUM+1
          KSUM=MIN(KSUM,10)
          T=DONE
          SUM=DONE
          DO K=1,KSUM
            T=T*(MU-(2*K-1)**2)*T8/K
            SUM=SUM+T
c            write(0,*)k,t,sum,done/sum
          ENDDO
C IF NOT CONVERGED, SAFE TO EVALUATE EXPLICITLY NOW AS THETA LARGE
          IF(ABS(T/SUM).GT.D1M4)THEN
c            write(0,*)temp(j),theta,sum,t,t/sum
            TT=DONE/THETA
            FX=SQRT(PI*THETA/2)/(BESSK(2,TT)*EXP(TT))
            FREL(J)=FX
          ELSE
            FREL(J)=DONE/SUM
          ENDIF
c          write(0,*)j,ksum,theta,t/sum,done/sum,fx
C
        ENDDO
C
      ELSE
C
        DO J=1,JTEMP
          FREL(J)=DONE
        ENDDO
C
      ENDIF
C
C
C***********************************
C WRITE DR DATA TO ADAS FORMAT ADF09
C *****            ****        *****
C***********************************
C
      IF(NBIN0.NE.0)THEN                          !BINNED CROSS SECTIONS
        EET(1)=-WNP(1)                            !.LT.0
        DO M=2,LMAX
          EET(M)=-EET(1)-WNP(M)                   !.GT.0
        ENDDO
        IF(NBIN0.GT.0)GO TO 629                   !SKIP ADF09 WRITES
      ENDIF
C
      IF(BLSOLD)LAB4='    '
      IF(BLSNEW)LAB4='/LS/'
      IF(BIC)LAB4='/IC/'
      IF(BCA)LAB4='/CA/'
      IF(NE.LE.NDIM20)THEN
        LAB2=LSQ(NE)
      ELSE
        LAB2='**'
      ENDIF
      WRITE(10,21)LAB2,NZ0,LAB4
      WNP0=WNP(1)*DKCM
      IF(BBNFP.AND..NOT.BHYBRD)THEN
        IF(BCA)THEN
          WRITE(10,253)WNP0,NBINM,NRX
        ELSE
          IF(NRX.LT.100)THEN
            IF(BLSOLD)WRITE(10,202)WNP0,NBINM,NRX
            IF(BLSNEW)WRITE(10,212)WNP0,NBINM,NRX
            IF(BIC)WRITE(10,2332)WNP0,NBINM,NRX
          ELSE
            IF(BLSOLD)WRITE(10,203)WNP0,NBINM,NRX
            IF(BLSNEW)WRITE(10,213)WNP0,NBINM,NRX
            IF(BIC)WRITE(10,2333)WNP0,NBINM,NRX
          ENDIF
          IF(BLSOLD)WRITE(10,22)
          IF(BLSNEW)WRITE(10,23)
          IF(BIC)WRITE(10,233)
        ENDIF
      ELSE
        IF(BCA)THEN
          WRITE(10,254)WNP0,NBINM
        ELSE
          IF(BLSOLD)WRITE(10,222)WNP0,NBINM
          IF(BLSNEW)WRITE(10,223)WNP0,NBINM
          IF(BIC)WRITE(10,224)WNP0,NBINM
        ENDIF
        NRX0=NRX
        NRX=NBINM
      ENDIF
C
C PARENT INDEXING
C
      IF(BLSOLD)THEN
        IWF=28
      ELSE
        IWF=26
      ENDIF
      DO M=1,NRX
        WNP(M)=-WNP(M)*DKCM+WNP0
        IF(BLS)THEN
          TW=IWS(M)*(2*IWL(M)+1)
          TW=(TW-DONE)/DTWO
        ENDIF
        IF(BIC)THEN
          TW=IWJ(M)
          TW=TW/DTWO
        ENDIF
        DO J=1,10
          QS0(J)=MBLNK
          QL0(J)=MBLNK
          IF(J.LE.LMP(M))THEN
            K=QLP(M,J)
            QS0(J)=LIT(QN(K))
            L=MIN0((QL(K)+1),20)
            QL0(J)=LABL(L)
          ENDIF
        ENDDO
        J1=MAX0(5,LMP(M))
        J0=J1-4
        IF(M.GT.NBINM)THEN
          MP=LABL(20)
        ELSE
          MP=MBLNK
        ENDIF
        IF(BLSOLD)
     X  WRITE(10,28)M,(QS0(J),QL0(J),QSP(M,J),J=J0,J1),IWS(M),IWL(M),TW
     X             ,WNP(M),MP
        IF(BIC.OR.BLSNEW)
     X  WRITE(10,26)M,(QS0(J),QL0(J),QSP(M,J),J=J0,J1),IWS(M),IWL(M),TW
     X             ,WNP(M),MP
        IF(BCA)
     X  WRITE(10,27)M,(QS0(J),QL0(J),QSP(M,J),J=J0,J1),TW,WNP(M),MP
      ENDDO
C
      IF(BHYBRD)THEN
        M=JKH(1)
        WNH0=WNH(M)*DKCM
        NRX=NRX0
        WRITE(10,254)WNH0,NRX
        DO M0=1,NRX
          M=JKH(M0)
          WNH(M)=-WNH(M)*DKCM+WNH0
          TW=JJH(M)-1
          TW=TW/DTWO
          DO J=1,10
            QS0(J)=MBLNK
            QL0(J)=MBLNK
            IF(J.LE.LMH(M))THEN
              K=QLH(M,J)
              QS0(J)=LIT(QN(K))
              L=MIN0((QL(K)+1),20)
              QL0(J)=LABL(L)
            ENDIF
          ENDDO
          J1=MAX0(5,LMH(M))
          J0=J1-4
C          IF(M0.GT.NBINM)THEN
C            MP=LABL(20)
C          ELSE
            MP=MBLNK
C          ENDIF
          WRITE(10,27)M0,(QS0(J),QL0(J),QSH(M,J),J=J0,J1),TW,WNH(M),MP
        ENDDO
      ENDIF
C
C RESOLVED INDEXING (AND AUGER RATES)
C
      IF(NRSOL.GT.0)THEN
C
      NSKP=0
      DO M=1,NRSOL
        IRSOL0(M)=M
        IF(BLS.AND.JJR(M).EQ.0)IRSOL0(M)=0 !UNMATCHED SO ASSUME UNWANTED
        IF(BIC.AND.SSR(M).EQ.0)IRSOL0(M)=0 !UNMATCHED SO ASSUME UNWANTED
        IF(BHYBRD)THEN
          IF(BLS)THEN
            WNR(M)=WNR(M)/SSR(M)
            JJR(M)=-JJR(M)
          ENDIF
          IF(BIC)THEN
            WNR(M)=WNR(M)/JJR(M)
            JJR(M)=JJR(M)-1                                  !BACK TO 2J
          ENDIF
        ENDIF
        IF(.NOT.BLSOLD)THEN
          IF(ITARR(M).LE.0)THEN
            IRSOL0(M)=0
          ELSE
            DO L=1,NBINM
              IF(AN(ITST,M,L).GT.ZERO)GO TO 44
            ENDDO
            IRSOL0(M)=0
          ENDIF 
        ENDIF
  44    CONTINUE
        IF(WNR(M).GE.EIONMN)THEN                            !FINAL CATCH
          IF(ITARR(M).GT.NBINM0.and.ITARR(M).LT.9999)IRSOL0(M)=0
        ENDIF
        IF(IRSOL0(M).EQ.0)THEN
          NSKP=NSKP+1
          WNR(M)=ZERO
        ENDIF
      ENDDO
c      write(0,*)'nskp=',nskp
C
      CALL HPSRTI(NRSOL,WNR,JVR(1))
c
c      do m0=1,nrsol
c        m=jvr(m0)
c        write(6,*)m,itarr(m),(qsr(m,i),qlr(m,i),i=1,lmr(m)),wnr(m)
c      enddo
C
      NRSOL=NRSOL-NSKP
      IF(NRSOL.GT.0)THEN
        M=JVR(1)
        WNR0=WNR(M)
      ELSE
        WNR0=ZERO
      ENDIF
      E00=WNR0
      WNR0=-WNR0*DKCM
      IF(BLS)LAB2='LS'
      IF(BCA)LAB2='CA'
      IF(BIC)LAB2='IC'
      IF(BHYBRD)THEN
        IF(BLS)WRITE(10,251)LAB2,WNR0,NRSOL
        IF(BIC)WRITE(10,241)LAB2,WNR0,NRSOL
      ELSE
        IF(BLSOLD)WRITE(10,25)LAB2,WNR0,NRSOL
        IF(BLSNEW)WRITE(10,24)LAB2,WNR0,NRSOL
        IF(BCA)WRITE(10,250)LAB2,WNR0,NRSOL
        IF(BIC)THEN
          IF(NRSOL.LT.10000)THEN
            WRITE(10,2401)LAB2,WNR0,NRSOL
          ELSE
            WRITE(10,2402)LAB2,WNR0,NRSOL
          ENDIF
          WRITE(10,240)
        ENDIF
      ENDIF
C
      DO M0=1,NRSOL
        M=JVR(M0)
        IF(WNR(M).GE.EIONMN)THEN
          MP=LABL(20)
          IAANLJ=ITARR(M)-1
          IAANLJ=MIN(IAANLJ,NBINM0)
          IW6=MIN(IAANLJ,6)
        ELSE
          MP=MBLNK
          IAANLJ=0
        ENDIF
        WNR(M)=(WNR(M)-E00)*DKCM
        ISSR=IABS(SSR(M))
        IF(BLS)THEN
          TW=ISSR*(2*LLR(M)+1)
          TW=(TW-DONE)/DTWO
        ENDIF
        IF(BIC)THEN
          TW=IABS(JJR(M))
          TW=TW/DTWO
        ENDIF
        DO 48 J=1,10
          QS0(J)=MBLNK
          QL0(J)=MBLNK
          IF(J-LMR(M))46,47,48
  47      IF(QNV(M).EQ.0)GO TO 46
          K=IABS(QNV(M))
          QS0(J)=LIT(K)
          L=MIN0((QLV(M)+1),20)
          QL0(J)=LABL(L)
          GO TO 48
  46      K=QLR(M,J)
          IF(QN(K).GT.NRSLMX)THEN
            WRITE(6,*)'***ERROR: ORBITAL CONFUSION BETWEEN onu FILES:'
     X                ,' M0,M,N,L='
            WRITE(6,*)M0,M,QN(K),QL(K)
            STOP'***ERROR: ORBITAL CONFUSION BETWEEN onu FILES'
          ENDIF
          QS0(J)=LIT(QN(K))
          L=MIN0((QL(K)+1),20)
          QL0(J)=LABL(L)
  48    CONTINUE
        J1=MAX0(5,LMR(M))
        J0=J1-4
        IF(BLSOLD)THEN
          WRITE(10,28)M0,(QS0(J),QL0(J),QSR(M,J),J=J0,J1),ISSR,LLR(M)
     X               ,TW,WNR(M),MP
        ELSE
          IF(BRSLF)THEN             !IF(.NOT.BHYBRD.AND.(BIC.OR.BLSNEW))
            WRITE(10,29)M0,ITARR(M),(QS0(J),QL0(J),QSR(M,J),J=J0,J1)
     X           ,ISSR,CLIT(LLR(M)),TW,WNR(M),MP,(AANLJ(M,L),L=1,IW6)
          ELSEIF(BLSNEW)THEN        !IF(.NOT.BCA.AND..NOT.BIC)
            WRITE(10,31)M0,ITARR(M),(QS0(J),QL0(J),QSR(M,J),J=J0,J1)
     X                ,CLIT(JJR(M)),TW,WNR(M),MP,(AANLJ(M,L),L=1,IW6)
          ELSE
            WRITE(10,30)M0,ITARR(M),(QS0(J),QL0(J),QSR(M,J),J=J0,J1)
     X                             ,TW,WNR(M),MP,(AANLJ(M,L),L=1,IW6)
          ENDIF
          IF(IAANLJ.GT.IW6)WRITE(10,290)(AANLJ(M,L),L=IW6+1,IAANLJ)
        ENDIF
      ENDDO
C
      ENDIF
C
C BUNDLED DATA
C
      NRX0=NBINRM0    !FINAL RECOMBINED THAT CONTAIN AUTOIONIZING STATES
      J1=NBINM0       !.LE.NBINM
      IF(BHYBRD)THEN
        J2=NRX0
        JT=J1*J2
      ELSE
        IF(J1.EQ.NRX0)J1=J1-1
        J2=NRX0-1
        JT=0
        DO J=1,J1
          DO K=J,J2
            JT=JT+1
          ENDDO
        ENDDO
      ENDIF
C
C BUNDLED-NL INDEXING AND AUGER RATES
C
      IF(NLMAX.GT.0)THEN
        MXLL=MIN(MXLL,LLMAX)
c        mxll=9
        NLREP=0
        DO I=1,IB0
          N=IBN(I)
          IF(N.GT.NLMAX)GO TO 43
          NLREP=NLREP+MIN(N,MXLL+1)
        ENDDO
  43    IBL0=I-1
        IF(JT.EQ.0)THEN
          WRITE(10,14)NLREP
        ELSEIF(JT.GT.8)THEN
          IF(J2.LT.100.AND.J1.LT.10)THEN
            IF(BHYBRD)THEN
              WRITE(10,5342)NLREP,((K,J,K=1,J2),J=1,J1)
            ELSE
              WRITE(10,5342)NLREP,((K+1,J,K=J,J2),J=1,J1)
            ENDIF
          ELSE
            IF(BHYBRD)THEN
              WRITE(10,5343)NLREP,((K,J,K=1,J2),J=1,J1)
            ELSE
              WRITE(10,5343)NLREP,((K+1,J,K=J,J2),J=1,J1)
            ENDIF
          ENDIF
        ELSE
          IF(BHYBRD)THEN
            WRITE(10,5340)NLREP,((K,J,K=1,J2),J=1,J1)
          ELSE
            WRITE(10,5340)NLREP,((K+1,J,K=J,J2),J=1,J1)
          ENDIF
          WRITE(10,5341)
        ENDIF
C
        IL=0
        M1=1
        DO I=1,IBL0
          N=IBN(I)
          L2=MIN(N,MXLL+1)
          DO L1=1,L2
            IL=IL+1
            IF(NRX0.LT.NBINRM.AND.J1.GT.0)THEN   !REMAP
              JX=0
              DO L=1,NBINM0
                IF(.NOT.BHYBRD)M1=L+1
                DO M=M1,NBINRM                   !STORED
                  IF(M.LE.NRX0)THEN              !WANTED
                    JX=JX+1
                    K=IMX(L,M)
                    UB0(JX,1)=AANL(K,L1,I)
                  ENDIF
                ENDDO
              ENDDO
              WRITE(10,17)IL,N,L1-1,(UB0(J,1),J=1,JX)
            ELSE
              IF(J1.GT.0)WRITE(10,17)IL,N,L1-1,(AANL(J,L1,I),J=1,IAAMX)
              IF(J1.EQ.0)WRITE(10,17)IL,N,L1-1
            ENDIF
          ENDDO
        ENDDO
C
        IF(BHYBRD)THEN
          IF(BLSNEW)THEN
            I1=MOD(JJR(1),2)
            WRITE(10,5413)NAUTY
          ELSE
            LAB1=' '
            WRITE(10,5414)NAUTY
          ENDIF
C
          DO I=1,NAUTY
            IF(JWRN(I)*JWRD(I).EQ.0)GO TO 49         !SHIFTED
            T=JWRN(I)
            ERN(I)=ERN(I)/T
            TE1=(ERN(I)-E00)*DKCM
            T=JWRD(I)
            ERD(I)=ERD(I)/T
            TE2=ERD(I)
            TE2=(TE2-E00)*DKCM
            T=(JWRN(I)+JWRD(I))
            T1=JWRN(I)
            T1=T1/T
            T1=MIN(T1,0.999D0)
            T1=MAX(T1,0.001D0)
            T2=JWRD(I)
            T2=T2/T
            T2=MIN(T2,0.999D0)
            T2=MAX(T2,0.001D0)
            IF(BLSNEW)THEN
              NS=IAUTY(I)/1000000
              IAUTY(I)=IAUTY(I)-NS*1000000
              NS=NS+NS-I1
              LAB1=CLIT(NS)
            ENDIF
            IP=IAUTY(I)/10000
            L=IAUTY(I)-IP*10000
            IB=L/20
            N=IBN(IB)
            L=L-IB*20
            WRITE(10,5417)IP,N,L,LAB1,T1,TE1,T2,TE2
 49       ENDDO
        ENDIF
      ENDIF
C
C BUNDLED-N INDEXING AND AUGER RATES
C
      IF(NLMAX.GT.0)THEN
        LAB5='INREP'
      ELSE
        LAB5='IREP '
      ENDIF
      IF(JT.EQ.0)THEN
        WRITE(10,13)IB0,LAB5
      ELSEIF(JT.GT.8)THEN
        IF(J2.LT.100.AND.J1.LT.10)THEN
          IF(BHYBRD)THEN
            WRITE(10,5332)IB0,LAB5,((K,J,K=1,J2),J=1,J1)
          ELSE
            WRITE(10,5332)IB0,LAB5,((K+1,J,K=J,J2),J=1,J1)
          ENDIF
        ELSE
          IF(BHYBRD)THEN
            WRITE(10,5333)IB0,LAB5,((K,J,K=1,J2),J=1,J1)
          ELSE
            WRITE(10,5333)IB0,LAB5,((K+1,J,K=J,J2),J=1,J1)
          ENDIF
        ENDIF
      ELSE
          IF(BHYBRD)THEN
            WRITE(10,5330)IB0,LAB5,((K,J,K=1,J2),J=1,J1)
          ELSE
            WRITE(10,5330)IB0,LAB5,((K+1,J,K=J,J2),J=1,J1)
          ENDIF
        WRITE(10,5331)
      ENDIF
C
      M1=1
      DO I=1,IB0
        IF(NRX0.LT.NBINRM.AND.J1.GT.0)THEN       !REMAP
          JX=0
          DO L=1,NBINM0
            IF(.NOT.BHYBRD)M1=L+1
            DO M=M1,NBINRM                       !STORED
              IF(M.LE.NRX0)THEN                  !WANTED
                JX=JX+1
                K=IMX(L,M)
                UB0(JX,I)=AAN(K,I)
              ENDIF
            ENDDO
          ENDDO
          WRITE(10,15)I,IBN(I),(UB0(J,I),J=1,JX)
        ELSE
          IF(J1.GT.0)WRITE(10,15)I,IBN(I),(AAN(J,I),J=1,IAAMX)
          IF(J1.EQ.0)WRITE(10,15)I,IBN(I)
        ENDIF
      ENDDO
C
      WRITE(10,1)
      DO J=1,JTEMP
        TEMP(J)=TEMP(J)*1.5789D5
      ENDDO
      JTW10=MIN(JTEMP,10)
      JTW11=JTW10+1
      JTW20=MIN(JTEMP,20)
      JTW21=JTW20+1
C
C PARTIAL DR DATA
C
      L0=0
      DO L=1,NBINM                            !LOOP OVER INITIAL PARENTS
C
        WRITE(10,1)
        LP=IWL(L)+1
        IF(BBNFP)THEN
          IF(BCA)THEN
            TW=(IWS(L)-DONE)/DTWO
            WRITE(10,53)L,TW
          ELSEIF(BLS)THEN
            WRITE(10,9)L,IWS(L),LABL(LP),IWS(L)
          ELSEIF(BIC)THEN
            TW=IWJ(L)
            TW=TW/DTWO
            WRITE(10,52)L,IWS(L),LABL(LP),TW
          ENDIF
          M1=1
          M2=NBINRM
        ELSE
          IF(BCA)THEN
            TW=(IWS(L)-DONE)/DTWO
            WRITE(10,4)L,TW
          ELSEIF(BLS)THEN
            NSYSM=2
            IF(IWS(L).EQ.1)NSYSM=1
            WRITE(10,729)L,IWS(L),LABL(LP),IWS(L),NSYSM
          ELSEIF(BIC)THEN
            NSYSM=1
            TW=IWJ(L)
            TW=TW/DTWO
            WRITE(10,8)L,IWS(L),LABL(LP),TW
          ENDIF            
          M1=L
          M2=L
        ENDIF
C
        WRITE(10,12)(TEMP(J),J=1,JTW20)
        IF(JTEMP.GT.JTW20)WRITE(10,6)(TEMP(J),J=JTW21,JTEMP)
C
C RESOLVED DATA
C
        IF(NRSOL.GT.0)THEN
          DO M0=1,NRSOL
            M=JVR(M0)
            IF(AN(ITST,M,L).GT.ZERO)THEN
              WRITE(10,3)M0,(FREL(J)*AN(J,M,L),J=1,JTW10)
              IF(JTEMP.GT.JTW10)
     X        WRITE(10,6)(FREL(J)*AN(J,M,L),J=JTW11,JTEMP)
            ENDIF
          ENDDO
          WRITE(10,1)
        ENDIF
C
C BUNDLED DATA
C
        DO 708 M=M1,M2                          !LOOP OVER FINAL PARENTS
C
          L0=L0+1
          IF(BBNFP)THEN
            IF(M.GT.NRX)GO TO 708
            NSYSM=1
            IF(NSYM.EQ.2)THEN                                 !LS NOT CA
              IF(BRSLF)THEN
                IF(ABS(IWS(L)-IWS(M)).GT.2)GO TO 708
                IF(IWS(L).GT.1.AND.IWS(L).EQ.IWS(M))NSYSM=2
                LM=IWL(M)+1
                IF(M2.LT.100)THEN
                  WRITE(10,716)M,IWS(M),LABL(LM),IWS(M),NSYSM
                ELSE
                  WRITE(10,726)M,IWS(M),LABL(LM),IWS(M),NSYSM
                ENDIF
              ELSE
                M0=JKH(M)
                TW=JJH(M0)-1
                TW=TW/DTWO
                IF(IWS(L).GT.1)NSYSM=2
                WRITE(10,727)M,TW,NSYSM
              ENDIF
            ENDIF
            IF(BCA)THEN
              IF(BHYBRD)THEN
                M0=JKH(M)
                TW=JJH(M0)-1
              ELSE
                TW=IWS(M)-DONE
              ENDIF
              TW=TW/DTWO
              WRITE(10,18)M,TW
            ELSEIF(BIC)THEN
              IF(BRSLF)THEN
                LM=IWL(M)+1
                TW=IWJ(M)
                TW=TW/DTWO
                IF(M2.LT.100)THEN
                  WRITE(10,2)M,IWS(M),LABL(LM),TW
                ELSE
                  WRITE(10,20)M,IWS(M),LABL(LM),TW
                ENDIF
              ELSE
                M0=JKH(M)
                TW=JJH(M0)-1
                TW=TW/DTWO
                WRITE(10,18)M,TW
              ENDIF
            ENDIF
          ENDIF
C
C LOOP-OVER 2 SPIN-SYSTEMS IF LS, ELSE 1.
C
          DO K=1,NSYSM
            IF(NSYM.EQ.2)THEN                                 !LS NOT CA
              ISPTL=IWS(L)+2*K-3
              IF(BRSLF)THEN
                ISPTM=IWS(M)+2*K-3
                ISPTL=MAX0(ISPTL,ISPTM)
              ENDIF
              IF(ISPTL.EQ.0)ISPTL=2
              WRITE(10,10)K,ISPTL
            ENDIF
C
C BUNDLED-NL
C
            IF(NLMAX.GT.0)THEN
              WRITE(10,7)
              IL=0
              DO I=1,IBL0
                N=IBN(I)
                L2=MIN(N,MXLL+1)
                DO L1=1,L2
                  IL=IL+1
                  IF(BNL(ITST,L1,I,K,L0).GT.ZERO)THEN
                    WRITE(10,3)IL,(FREL(J)*BNL(J,L1,I,K,L0),J=1,JTW10)
                    IF(JTEMP.GT.JTW10)
     X              WRITE(10,6)(FREL(J)*BNL(J,L1,I,K,L0),J=JTW11,JTEMP)
                  ENDIF
                ENDDO
              ENDDO
              WRITE(10,1)
            ENDIF
C
C BUNDLED-N 
C
            WRITE(10,11)LAB5
            DO I=1,IB0
              IF(BN(ITST,I,K,L0).GT.ZERO)THEN
                WRITE(10,3)I,(FREL(J)*BN(J,I,K,L0),J=1,JTW10)
                IF(JTEMP.GT.JTW10)
     X          WRITE(10,6)(FREL(J)*BN(J,I,K,L0),J=JTW11,JTEMP)
              ENDIF
            ENDDO
            WRITE(10,1)
          ENDDO
C
 708    CONTINUE         !END LOOP OVER FINAL PARENTS
C
      ENDDO              !END LOOP OVER INITIAL PARENTS
C
C END OF (PARTIAL) WRITES FOR ADF09
C *********************************
C
      IF(BLSOLD)RETURN
C
C NOW SUM TOTAL FROM GROUND+METASTABLES AND WRITE TO UNIT6 
C AND END OF ADF09 (NOT REQUIRED BY ADAS BUT FOR CONVENIENCE OF OTHERS).
C
C FIRST SUM SEQUENTIAL N
C
 629  IFLAG1=0
      IFLAGN=0
      N1=1                                          !NBIN0.GT.0 RE-ENTRY
      DO I=1,IB0
        IF(IBN(I).GT.N1)GO TO 710
        N1=N1+1
C
        IF(NBIN0.LE.0)THEN
          DO L=1,NBINM
            DO J=1,JTEMP
              ALF(J,L)=ALF(J,L)+ALFN(J,I,L)
            ENDDO
          ENDDO
        ENDIF
C
        IF(NBIN0.NE.0)THEN                        !BINNED CROSS SECTIONS
          DO L=1,LMAX
            tc(l)=zero
            IF(UB(I,1,L).NE.ZERO)IFLAG1=IFLAG1+1
            DO N=1,NBIN1
              SBIN(N,L)=SBIN(N,L)+UB(I,N,L)
              tc(l)=tc(l)+ub(i,n,l)
            ENDDO
            ub0(l,i)=tc(l)
          ENDDO
c          write(6,35)ibn(i),lv0,tc(1),tnu(i,1)
        ENDIF
      ENDDO
C
C INTERPOLATE TO NVINT (=100, DEFAULT) OR NCUT, IF SPECIFIED.
C (NB: ALWAYS INTERPOLATES FROM LAST ABOVE, EVEN IF ZERO - THIS MAY GIVE
C      UNEXPECTED RESULTS WHEN STARTING N IS ABOVE NDIM25=150+1.)
C IF A NEW THRESHOLD OPENS-UP HERE THEN RESULTS ARE INACCURATE!
C
  710 CONTINUE
C
      N2=0
      I=I+1
      i11=mod(ib0-i,2)
      IS=I-i11
      DO I=IS,IB0,2
C
        I0=I
        T1=IBN(I-2)
        T2=IBN(I-1)
        T3=IBN(I)
c      write(6,*)i,ibn(i-2),ibn(i-1),ibn(i)
C
        IF(NBIN0.NE.0)THEN 
          V1=T1**3
          V2=T2**3
          V3=T3**3
        ENDIF
C
  385   N1=IBN(I0-2)
        N2=IBN(I0-1)
        TN1=N1*N1
        TN2=N2*N2
        N1=N1+1+i11
c      write(6,*)i0,ibn(i0-2),ibn(i0-1)
C
        DO NN=N1,N2
          IF(NN.GT.NCUT)GO TO 737
C
          IF(NN.GE.NMN)THEN
c      write(6,*)'n=',nn
            TN=NN
            S1=(T2-TN)*(T3-TN)/((T2-T1)*(T3-T1))
            S2=(T1-TN)*(T3-TN)/((T1-T2)*(T3-T2))
            S3=(T1-TN)*(T2-TN)/((T1-T3)*(T2-T3))
C
            IF(NBIN0.LE.0)THEN
              DO L=1,NBINM
                DO J=1,JTEMP
                  TS=S1*ALFN(J,I-2,L)+S2*ALFN(J,I-1,L)+S3*ALFN(J,I,L)
                  IF(TS.GT.DZERO)ALF(J,L)=ALF(J,L)+TS        !CASE LOW-T
                ENDDO
              ENDDO
            ENDIF
C
            IF(NBIN0.NE.0)THEN                    !BINNED CROSS SECTIONS
c              write(6,*)nn
              S1=S1*V1
              S2=S2*V2
              S3=S3*V3
              TNN=NN*NN
              DE=DZ*(TNN-TN2)/(TN2*TNN)
              DO L=1,LMAX
                tc(l)=zero
                IF(TNU(I0-1,L).NE.ZERO)THEN
                  TT=S1*TNU(I-2,L)+S2*TNU(I-1,L)+S3*TNU(I,L)
                  IF(TT.GT.ZERO)THEN
                    TT=TT/(TN*TNN)
                    TT=TT/TNU(I0-1,L)
                    DO N=1,NBIN1
                      IF(UB(I0-1,N,L).GT.ZERO)THEN
                        IF(N.EQ.1)IFLAG1=IFLAG1+1
                        ERES=EBIN(N+1)-EBIN(N)
                        if(blog)then
                          t=ebin(n+1)-eres/2+de
                          if(de.ge.zero)then
                            do k=n,nbin1
                              if(ebin(k+1).ge.t)go to 398
                            enddo
                            k=nbin
                          else
                            do k=n,1,-1
                              if(ebin(k).lt.t)go to 398
                            enddo
                            k=0
                          endif
  398                     continue
c                          write(6,*)n,k,n+int(de/eres)
                        else
                          NDE=INT(DE/ERES)
                          K=N+NDE
                        endif
                        IF(K.GT.0.AND.K.LE.NBIN1)THEN
                          TS=UB(I0-1,N,L)*TT*EBIN(N+1)/EBIN(K+1)
                          ts=ts*eres/(ebin(k+1)-ebin(k))!case non-linear
                          SBIN(K,L)=SBIN(K,L)+TS
                          tc(l)=tc(l)+ts
c                          write(6,35)n,k,ebin(n+1),ebin(k),ebin(k+1)
c     x                              ,eres*UB(I0-1,N,L),ts
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO
              IF(NN.EQ.IBN(I0-1))THEN
c                write(6,35)nn,lv0,tc(1),tnu(i0-1,1)
                DO L=1,NBINM
                  UB0(L,I0-1)=TC(L)
                  IF(UB0(L,I0-1).GT.1.1*UB0(L,I0-2))IFLAGN=N2
c                  write(6,*)nn,ub0(l,i0-1),ibn(i0-2),ub0(l,i0-2)
                ENDDO
c              else
c                write(6,35)nn,lv0,tc(1)            !interpolated values
              ENDIF
            ENDIF
C
          ENDIF
C
        ENDDO
C
        i11=0
        I0=I0+1
        IF((I0-1).EQ.I)GO TO 385
C
        IF(NBIN0.NE.0)THEN
          DE=DZ*(TN2-TN1)/(TN2*TN1)
          NDE=INT(DE/(EBIN(NBIN)-EBIN(NBIN1)))
        ELSE
          NDE=0
        ENDIF
        IF(NDE.EQ.0.AND.N2.GE.NVINT)GO TO 712
C
      ENDDO
      IF(N2*NCUT.GT.N2*N2)THEN
        WRITE(6,215)NCUT,N2
        WRITE(0,215)NCUT,N2
      ENDIF
C
      GO TO 737
C
  712 CONTINUE
C
C SIMPSON'S RULE TO N=IBN(IB0)
C
      IS=I+2
      DO I=IS,IB0,2
C
        T1=IBN(I-2)*IBN(I-2)
        T3=IBN(I)*IBN(I)
        H=(T3-T1)/(T1*T3)
        H=H/DTWELV
        T=IBN(I-2)
        T1=T1*T*H
        T2=IBN(I-1)
        T2=T2*T2*T2
        T2=T2*DFOUR*H
        T=IBN(I)
        T3=T3*T*H
C
        IF(NBIN0.LE.0)THEN
          DO L=1,NBINM
            DO J=1,JTEMP
              ALF(J,L)=ALF(J,L)+T1*ALFN(J,I-2,L)+T2*ALFN(J,I-1,L)
     X                         +T3*ALFN(J,I,L)
            ENDDO
          ENDDO
        ENDIF
C
        IF(NBIN0.NE.0)THEN                        !BINNED CROSS SECTIONS
          DO L=1,LMAX
            IF(UB(I,1,L).NE.ZERO)IFLAG1=IFLAG1+1
            DO N=1,NBIN1
              SBIN(N,L)=SBIN(N,L)+T1*UB(I-2,N,L)+T2*UB(I-1,N,L)
     X  			 +T3*UB(I,N,L)
            ENDDO
          ENDDO
        ENDIF
C
      ENDDO
C
      IF(NBIN0.LE.0)THEN
        DO L=1,NBINM
          DO J=1,JTEMP
            TI=-ALFN(J,IS-2,L)/DTWO
            TF=ZERO
            IF(I-2.EQ.IB0)TF=ALFN(J,IB0,L)/DTWO
            ALF(J,L)=ALF(J,L)+TI+TF
          ENDDO
        ENDDO
      ENDIF
C
      IF(NBIN0.NE.0)THEN                          !BINNED CROSS SECTIONS
        DO L=1,LMAX
          DO N=1,NBIN1
            TI=-UB(IS-2,N,L)/DTWO
            TF=ZERO
            IF(I-2.EQ.IB0)TF=UB(IB0,N,L)/DTWO
            SBIN(N,L)=SBIN(N,L)+TI+TF
          ENDDO
        ENDDO
      ENDIF
C
  737 CONTINUE
C
      IF(IFLAGN.GT.0)THEN
        WRITE(6,734)IFLAGN
        WRITE(0,734)IFLAGN
      ENDIF
      IF(BLOG.AND.IFLAG1.NE.0)THEN
        WRITE(6,733)-EBIN(2)/10
        WRITE(0,733)-EBIN(2)/10
      ENDIF
C
      IF(NBIN0.GT.0)RETURN                   !BINNED CROSS SECTIONS ONLY
C
C WRITE TOTALS TO UNIT6 AND ADF09
C
      WRITE(6,731)
      WRITE(10,731)
      NMXW=MIN(NBINM,10)
C
      DO J=1,JTEMP
        WRITE(6,732)TEMP(J),(FREL(J)*ALF(J,L),L=1,NBINM)
        WRITE(10,732)TEMP(J),(FREL(J)*ALF(J,L),L=1,NMXW)
        TEMP(J)=TEMP(J)/1.5789D5                        !CASE NBIN0.LT.0
      ENDDO 
C
      RETURN
C
C
   1  FORMAT(1X)
   2  FORMAT(85X,'--------------------------'/85X,'PRTF=',I2,2X,
     X'LVLPRT=',' (',I1,A1,1X,F4.1,')')
  18  FORMAT(85X,'--------------------------'/85X,'PRTF=',I2,2X,
     X'CFGPRT=',' (',F7.1,')')
  20  FORMAT(85X,'--------------------------'/85X,'PRTF=',I3,1X,
     X'LVLPRT=',' (',I1,A1,1X,F4.1,')')
   3  FORMAT(I6,5X,1P,10E10.2)
   4  FORMAT(/1X,'-----------------------------------'/
     X1X,'PRT=',I2,2X,'CFGPRT=',' (',F7.1,')'/)
C   5  FORMAT(85X,'--------------------------'/85X,'PRTF=',I2)
   6  FORMAT(11X,1P,10E10.2)
   7  FORMAT(3X,'ILREP'/3X,'-----')
   8  FORMAT(/1X,'-----------------------------------'/
     X1X,'PRT=',I2,2X,'LVLPRT=',' (',I1,A1,1X,F4.1,')'/)
   9  FORMAT(/1X,'--------------------------------'/
     X1X,'PRTI=',I2,2X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2/)
  10  FORMAT(94X,'-----------------'/94X,'SYS=',I2,2X,'SPNSYS=',I2)
  11  FORMAT(3X,A5/3X,'-----')
  12  FORMAT(3X,'INDX TE=',1P,10E10.2/3X,'---- ---',10E10.2)
  13  FORMAT(/3X,'N-SHELL INDEXING & AUGER RATES',23X,'NREP=',I3
     X/3X,'------------------------------'/3X,A5,3X,'N'
     X/3X,'-----',3X,'-')
  14  FORMAT(/3X,'NL-SHELL INDEXING & AUGER RATES',21X,'NLREP=',I3
     X/3X,'-------------------------------'/3X,'ILREP',3X,'N',3X,'L'
     X/3X,'-----',3X,'-',3X,'-')
5413  FORMAT(/3X,'NL-SHELL I.P. STRADDLING CONFGS',21X,'NIPPY=',I3
     X/3X,'-------------------------------'/3X,'INDP',4X,'N',3X,'L',2X
     X,'S',3X,'FWB',9X,'WNB',4X,'FWA',9X,'WNA'/3X,'----',4X,'-',3X,'-'
     X,2X,'-',2X,'----',1X,'------------',2X,'----',1X,'------------')
5414  FORMAT(/3X,'NL-SHELL I.P. STRADDLING CONFGS',21X,'NIPPY=',I3
     X/3X,'-------------------------------'/3X,'INDP',4X,'N',3X,'L'
     X,6X,'FWB',9X,'WNB',4X,'FWA',9X,'WNA'/3X,'----',4X,'-',3X,'-'
     X,5X,'----',1X,'------------',2X,'----',1X,'------------')
  15  FORMAT(I6,I6,9X,1P,8E10.2/(21X,8E10.2))
  17  FORMAT(I6,I6,I4,5X,1P,8E10.2/(21X,8E10.2))
5417  FORMAT(I6,I6,I4,2X,A1,2(2X,F4.3,1X,F12.1))
  21  FORMAT("SEQ='",A2,"'",5X,"NUCCHG=",I2,50X,A4)
 202  FORMAT(/3X,'PARENT TERM INDEXING',12X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2)
 203  FORMAT(/3X,'PARENT TERM INDEXING',12X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I3)
  22  FORMAT(
     X3X,'--------------------'/3X,'INDP',4X,'CODE',17X,'S L   WI',8X,
     X'WNP'/3X,'----',4X,'----',17X,'- -   --',2X,'----------')
 212  FORMAT(/3X,'PARENT TERM INDEXING',17X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2)
 213  FORMAT(/3X,'PARENT TERM INDEXING',17X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I3)
 215  FORMAT(/' *** WARNING:YOUR NMAX=',I4,' DOES NOT MATCH A'
     X,' REPRESENTATIVE N-VALUE.'/5X,'CUT-OFF APPLIED AT N=',I4/)
  23  FORMAT(
     X3X,'--------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',8X,
     X'WNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
  24  FORMAT(/3X,A2,' RESOLVED TERM INDEXING',12X,'BWNR=',F12.1,3X,'NTRM
     X=',I4/3X,'-------------------------'/3X,'INDX',2X,'INDP',3X,'CODE'
     X,17X,'S L   WJ',8X,'WNR',7X,'AUGER RATES:  INDP-INDP'', INDP''=1,'
     X,'...INDP-1'/3X,'----',2X,'----',3X,'----',17X
     X,'- -   --',2X,'----------',6X,11('-'))
  25  FORMAT(/3X,A2,' RESOLVED TERM INDEXING',7X,'BWNR=',F12.1,3X,'NTRM=
     X',I4/3X,'-------------------------'/3X,'INDX',4X,'CODE',17X,'S L
     X WJ',8X,'WNR'/3X,'----',4X,'----',17X,'- -   --',2X,'----------')
 251  FORMAT(/3X,A2,' AVERAGED CONFG INDEXING',11X,'BWNR=',F12.1,3X
     X,'NCFG=',I4/3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'S     WJ',8X,'WNR',7X,"AUGER RATES:  INDP-INDP',"
     X," INDP'=1,...INDP-1"/3X,'----',2X,'----',3X,'----',17X
     X,'- ------',2X,'----------',6X,11('-'))
  26  FORMAT(I6,10X,5(A1,A1,A1,1X),'(',I1,')',I1,'(',F4.1,')',F11.1,A1)
  27  FORMAT(I6,10X,5(A1,A1,A1,1X),'(',F8.1,')',F11.1,A1)
  28  FORMAT(I6,5X,5(A1,A1,A1,1X),'(',I1,')',I1,'(',F4.1,')',F11.1,A1)
  29  FORMAT(I6,I6,4X,5(A1,A1,A1,1X),'(',I1,')',A1,'(',F4.1,')',F11.1,A1
     X,3X,1P,6E10.2)
 290  FORMAT(61X,1P,6E10.2)
  30  FORMAT(I6,I6,4X,5(A1,A1,A1,1X),'(',F8.1,')',F11.1,A1
     X,3X,1P,6E10.2)
  31  FORMAT(I6,I6,4X,5(A1,A1,A1,1X),'(',A1,F7.1,')',F11.1,A1
     X,3X,1P,6E10.2)
  52  FORMAT(/1X,'--------------------------'/
     X1X,'PRTI=',I2,2X,'LVLPRT=',' (',I1,A1,1X,F4.1,')'/)
  53  FORMAT(/1X,'--------------------------'/
     X1X,'PRTI=',I2,2X,'CFGPRT=',' (',F7.1,')'/)
C ,2X,'LVLPRT=',' (',I1,A1,1X,F4.1,')'/)
 157  FORMAT(2I5,1P,E15.4,20(I5))
 158  FORMAT(2I5,1P,E15.4,10(I5,E12.2))
 222  FORMAT(/3X,'PARENT TERM INDEXING',12X,'BWNP=',F12.1,4X,'NPRNT=',I2
     X/3X,'--------------------'/3X,'INDP',4X,'CODE',17X,'S L   WI',8X,
     X'WNP'/3X,'----',4X,'----',17X,'- -   --',2X,'----------')
 223  FORMAT(/3X,'PARENT TERM INDEXING',17X,'BWNP=',F12.1,4X,'NPRNT=',I2
     X/3X,'--------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',8X,
     X'WNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
 224  FORMAT(/3X,'PARENT LEVEL INDEXING',16X,'BWNP=',F12.1,4X,'NPRNT=',
     XI2/3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',
     X8X,'WNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
2332  FORMAT(/3X,'PARENT LEVEL INDEXING',16X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2)
2333  FORMAT(/3X,'PARENT LEVEL INDEXING',16X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I3)
 233  FORMAT(
     X3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',8X,'
     XWNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
 241  FORMAT(/3X,A2,' AVERAGED CONFG INDEXING',11X,'BWNR=',F12.1,3X
     X,'NCFG=',I4/3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'      WJ',8X,'WNR',7X,"AUGER RATES:  INDP-INDP',"
     X," INDP'=1,...INDP-1"/3X,'----',2X,'----',3X,'----',17X
     X,'--------',2X,'----------',6X,11('-'))
2401  FORMAT(/3X,A2,' RESOLVED LEVEL INDEXING',11X,'BWNR=',F12.1,3X
     X,'NLVL=',I4)
2402  FORMAT(/3X,A2,' RESOLVED LEVEL INDEXING',11X,'BWNR=',F12.1,3X
     X,'NLVL=',I5)
 240  FORMAT(3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'S L   WJ',8X,'WNR',7X,"AUGER RATES:  INDP-INDP',"
     X," INDP'=1,...INDP-1"/3X,'----',2X,'----',3X,'----',17X
     X,'- -   --',2X,'----------',6X,11('-'))
 250  FORMAT(/3X,A2,' RESOLVED CONFG INDEXING',11X,'BWNR=',F12.1,3X
     X,'NCFG=',I4/3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'      WJ',8X,'WNR',7X,"AUGER RATES:  INDP-INDP',"
     X," INDP'=1,...INDP-1"/3X,'----',2X,'----',3X,'----',17X
     X,'--------',2X,'----------',6X,11('-'))
 253  FORMAT(/3X,'PARENT CONFG INDEXING',16X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2/
     X3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'      WI',8X,'
     XWNP'/3X,'----',9X,'----',17X,'--------',2X,'----------')
 254  FORMAT(/3X,'PARENT CONFG INDEXING',16X,'BWNP=',F12.1,4X,'NPRNT=',
     XI2/3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'      WI',
     X8X,'WNP'/3X,'----',9X,'----',17X,'--------',2X,'----------')
 5330 FORMAT(/3X,'N-SHELL INDEXING & AUGER RATES',23X,'NREP=',I3
     X/3X,'------------------------------'/3X,A5,3X,'N'
     X,2X,'M''-M = ',8(6X,I2,'-',I1))
 5331 FORMAT(3X,'-----',3X,'-',2X,'----')
 5332 FORMAT(/3X,'N-SHELL INDEXING & AUGER RATES',23X,'NREP=',I3
     X/3X,'------------------------------'/3X,A5,3X,'N'
     X,2X,'M''-M = ',8(6X,I2,'-',I1)/3X,'-----',3X,'-',2X,'----',3X
     X,8(6X,I2,'-',I1)/(21X,8(6X,I2,'-',I1)))
 5333 FORMAT(/3X,'N-SHELL INDEXING & AUGER RATES',23X,'NREP=',I3
     X/3X,'------------------------------'/3X,A5,3X,'N'
     X,2X,'M''-M = ',8(4X,I3,'-',I2)/3X,'-----',3X,'-',2X,'----',3X
     X,8(4X,I3,'-',I2)/(21X,8(4X,I3,'-',I2)))
 5340 FORMAT(/3X,'NL-SHELL INDEXING & AUGER RATES',21X,'NLREP=',I3
     X/3X,'-------------------------------'/3X,'ILREP',3X,'N',3X,'L'
     X,2X,'M''-M = ',2X,8(I2,'-',I1,6X))
 5341 FORMAT(3X,'-----',3X,'-',3X,'-',2X,'----')
 5342 FORMAT(/3X,'NL-SHELL INDEXING & AUGER RATES',21X,'NLREP=',I3
     X/3X,'-------------------------------'/3X,'ILREP',3X,'N',3X,'L'
     X,2X,'M''-M = ',2X,8(I2,'-',I1,6X)/3X,'-----',3X,'-',3X,'-',2X
     X,'----',5X,8(I2,'-',I1,6X)/(27X,8(I2,'-',I1,6X)))
 5343 FORMAT(/3X,'NL-SHELL INDEXING & AUGER RATES',21X,'NLREP=',I3
     X/3X,'-------------------------------'/3X,'ILREP',3X,'N',3X,'L'
     X,2X,'M''-M = ',8(I3,'-',I2,4X)/3X,'-----',3X,'-',3X,'-',2X
     X,'----',3X,8(I3,'-',I2,4X)/(25X,8(I3,'-',I2,4X)))
 716  FORMAT(/1X,'-----------------------------------------'/
     X1X,'PRTF=',I2,2X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2,2X
     X,'NSYS=',I2)
 726  FORMAT(/1X,'-----------------------------------------'/
     X1X,'PRTF=',I3,1X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2,2X
     X,'NSYS=',I2)
 727  FORMAT(/1X,'------------------------------------'/
     X1X,'PRTF=',I2,2X,'CFGPRT=',' (',F7.1,')',2X,'NSYS=',I2)
 729  FORMAT(/1X,'-----------------------------------------'/
     X1X,'PRT=',I2,2X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2,2X
     X,'NSYS=',I2/)
 731  FORMAT(//'    T(K) ',4X,'ALFT( 1)',2X,'ALFT( 2)',2X,'ALFT( 3)'
     X,2X,'ALFT( 4)',2X,'ALFT( 5)',2X,'ALFT( 6)'
     X,2X,'ALFT( 7)',2X,'ALFT( 8)',2X,'ALFT( 9)',2X,'ALFT(10)'
     X/4X,'----',3X,10(2X,'--------'))
 732  FORMAT(1PE10.2,1X,(10E10.2))
 733  FORMAT(/' *** WARNING: UNRESOLVED NEAR-THRESHOLD RESONANCES EXIST'
     X,'! ***'/5X,'TRY REDUCING EWIDTH TO ',1PE8.1,'*UNITS')
 734  FORMAT(/' *** WARNING: NEW SERIES OPENING-UP BELOW N=',I4,
     X' AMIDST NON-SEQUENTIAL N-VALUES, I.E., INTERPOLATION BAD!')
 993  FORMAT(2I2,1X,I3,I2,2I5,3X,F15.8)
C
      END
C
C***********************************************************************
C
       SUBROUTINE DIPOL(JSW,N1,N2,E2,LMAX,CP,CM,JC)
C
C-----------------------------------------------------------------------
C
C  ALAN BURGESS DAMTP CAMBRIDGE, MODS BY NRB.
C
C  SR.DIPOL CALCULATES SQUARES OF HYDROGENIC DIPOLE LENGTH RADIAL MATRIX
C  ELEMENTS FOR BOUND-BOUND OR BOUND-FREE TRANSITIONS.
C
C  BOUND STATES ARE NORMALISED TO UNITY.
C  FREE STATES ARE NORMALISED TO ASYMPTOTIC AMPLITUDE K**(-0.5).
C
C  N.B. DIPOLE ACCELERATION MATRIX ELEMENT = (E12**2/4Z) * DIPOLE LENGTH
C  WHERE E12 = - N1**(-2) + N2**(-2)  FOR BOUND-BOUND
C            = - N1**(-2) + E2        FOR BOUND-FREE
C          Z = REDUCED CHARGE
C  INPUT:
C   FOR BOUND-BOUND,SET JSW=NEGATIVE
C                     N1,N2=PRINCIPAL QUANTUM NUMBERS OF STATES
C                      LMAX=RANGE OF ANGULAR MOMENTUM QUANTUM NUMBERS
C   FOR BOUND-FREE, SET JSW=POSITIVE
C                       N1=BOUND STATE PRINCIPAL QUANTUM NUMBER
C                       E2=FREE STATE ENERGY IN RYDBERGS (=K**2)
C
C  OUTPUT:
C   VECTOR CP(L),L=1,LMAX,CONTAINS SQUARED MATRIX ELEMENTS FOR ANGULAR
C                         MOMENTUM TRANSITIONS FROM L-1 TO L,
C   VECTOR CM(L),L=1,LMAX,CONTAINS SQUARED MATRIX ELEMENTS FOR ANGULAR
C                         MOMENTUM TRANSITIONS FROM L TO L-1,
C               IN BOTH CASES THE TRANSITION IS FROM LOWER TO HIGHER
C               ENERGY, INDEPENDANT OF THE SIGN OF N1-N2 FOR BOUND-BOUND
C               CASES. IF N1=N2 THEN CP(L)=CM(L).
C   VECTOR JC(L),L=1,LMAX WILL USUALLY BE ZERO AND MAY THEN BE IGNORED,
C               BUT FOR EXTREME INPUT VALUES THERE IS POSSIBILITY OF
C               OVER OR UNDERFLOW OF CP(L) OR CM(L),IN WHICH CASE THE
C               OUTPUT VALUES OF CP(L) AND CM(L) SHOULD BE MULTIPLIED
C               BY (1.0D10)**JC(L) TO OBTAIN TRUE VALUES.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
C      PARAMETER (PI=3.14159265359D0)
      PARAMETER (S1=1.0D10)
      PARAMETER (S2=1.0D-10)
      PARAMETER (TEST1=1.0D-20)
      PARAMETER (TEST2=1.0D20)
      PARAMETER (TEST3=0.044D0)
      PARAMETER (TEST4=0.1D0)
      PARAMETER (TEST5=300.0D0)
      PARAMETER (TEST6=1.0D-30)
      PARAMETER (TEST7=1.0D30)
C
      DIMENSION CP(LMAX),CM(LMAX),JC(LMAX)
C
      PI=ACOS(-DONE)
C
      N=N1
      E=E2
      IF(JSW.LE.0)THEN
        EN2=N2
        N3=N2
        IF(N1.EQ.N2)GO TO 59
        IF(N2.LT.N1)THEN
          N=N2
          EN2=N1
          N3=N1
        ENDIF
        E=-DONE/(EN2*EN2)
      ENDIF
C
      EN=N
      ENN=EN*EN
      E1=-DONE/ENN
      JMAX=LMAX
      C1=DONE
      C2=DZERO
      JS=0
      L=N+1
      IF(N.LE.LMAX)THEN
        CP(N)=DONE
        CM(N)=DZERO
        JC(N)=0
        JMAX=N-1
        DO I=L,LMAX
          CP(I)=DZERO
          CM(I)=DZERO
          JC(I)=0
        ENDDO
      ENDIF
C
    9 L=L-1
      IF(L.GT.1)THEN
        EL=L
        ELL=EL*EL
        T1=DONE+ELL*E1
        T2=DONE+ELL*E
        T3=L+L-1
        T4=DONE/(T3+DONE)
        T5=(T3*T1*C2+T2*C1)*T4
        C1=(T1*C2+T3*T2*C1)*T4
        C2=T5
   11   IF(C1*C1.GT.TEST2)THEN
          C1=S2*C1
          C2=S2*C2
          JS=JS+1
          GO TO 11
        ENDIF
        IF(L.LE.LMAX+1)THEN
          CP(L-1)=C1
          CM(L-1)=C2
          JC(L-1)=JS
        ENDIF
        GO TO 9
      ENDIF
C
      JS=0
      T=4
      T=DONE/(T*EN*ENN)
      IF(JSW.LE.0)THEN                          !JSW.LT.0
        ENN2=EN2*EN2
        T1=4
        T1=T1*ENN*ENN2/(ENN2-ENN)
        T1=T1*T1
        T=T*T1*T1/(EN2*ENN2)
        IF(N3.LE.30)THEN
          T=T*((EN2-EN)/(EN2+EN))**(N3+N3)
          GO TO 34
        ENDIF
        E21=E/E1
        IF(E21.LE.TEST4)THEN
          T2=DZERO
          DO J=1,11
            T3=2*(11-J)+1
            T2=DONE/T3+T2*E21
          ENDDO
          T2=T2+T2
        ELSE
          T3=EN/EN2
          T2=LOG((DONE+T3)/(DONE-T3))/T3
        ENDIF
        T2=T2+T2
        T1=T1*EXP(-T2)
C
      ELSE                                      !JSW.GT.0
C
        T1=4
        T1=T1*ENN/(DONE+ENN*E)
        T1=T1*T1
        T=T*T1*T1
        IF(E.LT.TEST3)THEN
          T3=2
          T=T*(PI/T3)
        ELSE
          T4=SQRT(E)
          IF(T4.LE.TEST5)THEN
            T3=(PI+PI)/T4
            T3=DONE-EXP(-T3)
            T3=DONE/T3
          ELSE
            T4=PI/T4
            T3=3
            T3=(DONE+T4+T4*T4/T3)/(T4+T4)
          ENDIF
          T2=2
          T=T*(PI*T3/T2)
        ENDIF
C
        T4=ENN*E
        IF(T4.LE.TEST4)THEN
          T2=DZERO
          DO J=1,11
            T3=2*(11-J)+1
            T2=DONE/T3-T2*T4
          ENDDO
        ELSE
          T3=SQRT(T4)
          T2=ATAN(T3)/T3
        ENDIF
        T2=T2+T2
        T2=T2+T2
        T1=T1*EXP(-T2)
      ENDIF
C                                               !ALL JSW
   34 DO J=1,N
        TJ=J+J
        T2=TJ*(TJ-DONE)
        T2=T2*T2
        T=T*T1/T2
   35   IF(T.LE.TEST1)THEN
          T=T*S1
          JS=JS-1
          GO TO 35
        ENDIF
   37   IF(T.GE.TEST2)THEN
          T=T*S2
          JS=JS+1
          GO TO 37
        ENDIF
      ENDDO
      J=0
C
   40 J=J+1
      IF(J.LE.JMAX)THEN
        TJ=J
        TJ=TJ*TJ
        T1=DONE+TJ*E1
        T2=DONE+TJ*E
        T3=CP(J)
        T3=T2*T*T3*T3
        T4=CM(J)
        T4=T1*T*T4*T4
        L1=JC(J)+JC(J)+JS
C
   42   IF(L1.LT.0)THEN
          IF(T4.GT.TEST6)THEN
            L1=L1+1
            T3=T3*S2
            T4=T4*S2
            GO TO 42
          ENDIF
        ELSEIF(L1.GT.0)THEN
          IF(T3.LT.TEST7)THEN
            L1=L1-1
            T3=T3*S1
            T4=T4*S1
            GO TO 42
          ENDIF
        ENDIF
C
        CP(J)=T3
        CM(J)=T4
        JC(J)=L1
        T=T*T1*T2
   48   IF(T.GT.TEST2)THEN
          T=T*S2
          JS=JS+1
          GO TO 48
        ENDIF
        GO TO 40
      ENDIF
C
      IF(N.LE.LMAX)THEN
        T2=DONE+ENN*E
        T3=CP(N)
        T3=T2*T*T3*T3
        L1=JC(N)+JC(N)+JS
C
   52   IF(L1.LT.0)THEN
          IF(T3.GT.TEST6)THEN
            L1=L1+1
            T3=T3*S2
            GO TO 52
          ENDIF
        ELSEIF(L1.GT.0)THEN
          IF(T3.LT.TEST7)THEN
            L1=L1-1
            T3=T3*S1
            GO TO 52
          ENDIF
        ENDIF
C
        CP(N)=T3
        JC(N)=L1
      ENDIF
C
      RETURN
C
   59 JMAX=LMAX
      IF(N.LE.LMAX)THEN
        DO L=N,LMAX
          CP(L)=DZERO
          CM(L)=DZERO
          JC(L)=0
        ENDDO
        JMAX=N-1
      ENDIF
      T1=9
      T2=4
      T3=(T1/T2)
      T1=EN2*EN2
      T2=T1*T3
      DO J=1,JMAX
        TJ=J
        JC(J)=0
        T=T2*(T1-TJ*TJ)
        CP(J)=T
        CM(J)=T
      ENDDO
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION EMN3(NV)
C
C-----------------------------------------------------------------------
C
C NRB:
C EVALUATE EULER-MACLAURIN SUM N=1 TO INFINITY OF 1/(N+NV)**3
C RELATIVE ACCURACY OF 1.D-3 AT NV=2, 1.D-6 AT NV=4, 1.D-10 AT NV=11.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0)
C
      T2=NV
      T2=T2*T2
      T=T2
      FI=ONE/(2*T)
      T=T*NV
      F0=-ONE/(2*T)
      T=T*NV
      F1=ONE/(4*T)
      T=T*T2
      F3=-ONE/(12*T)
      T=T*T2
      F5=ONE/(12*T)
      T=T*T2
      F7=-ONE*3/(20*T)
C
      EMN3=FI+F0+F1+F3+F5+F7
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE HPSRTI (N,A,IP)
C
C-----------------------------------------------------------------------
C
C IMPLICIT HEAPSORT  BY *MAGNITUDE* OF
C INPUT:  VECTOR A, LENGTH N.
C OUTPUT: DOWN-ORDERED POINTER IN IP, A IS UNCHANGED.
C        (UP-ORDERED CAN BE OBTAINED BY CHANGING .LT. TO .GT. AS BELOW).
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(*),IP(*)
C
      DO I=1,N
        IP(I)=I
      ENDDO
C
      IF(N.LT.2)RETURN
C
      L=N/2+1
      IT=N
C
  1   IF(L.GT.1)THEN
        L=L-1
        IPT=IP(L)
      ELSE
        IPT=IP(IT)
        IP(IT)=IP(1)
        IT=IT-1
        IF(IT.EQ.1)THEN
          IP(1)=IPT
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
C
  2   IF(J.LE.IT)THEN
        IF(J.LT.IT)THEN
          IF(abs(A(IP(J+1))).lT.abs(A(IP(J))))J=J+1  !.lt. down, .gt .up
        ENDIF
        IF(abs(A(IP(J))).lT.abs(A(IPT)))THEN         !.lt. down, .gt .up
          IP(I)=IP(J)
          I=J
          J=J+J
        ELSE
          J=IT+1
        ENDIF
        GO TO 2
      ENDIF
      IP(I)=IPT
      GO TO 1
C
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION QDT(QD,NZ0,NE,N,L,KAPPA)
C
C-----------------------------------------------------------------------
C
C NRB:
C EVALUATES ONE-ELECTRON ENERGY WITH NON-ZERO QUANTUM DEFECT
C
C   : QD0, UNIVERSAL QUANTUM DEFECT GIVEN BY
C         QD0*(NE**1.67-1)/(Z0**.67*Z**.33*(1+L**3))
C         CURRENT VALUE IN FUNCTION QDT IS QD0=0.182
C
C KAPPA= 0 NON-RELATIVISTIC
C      =-1 KAPPA-AVERAGE RELATIVISTIC
C      = L RELATIVISTIC FOR J=L-0.5
C      =-L-1 RELATIVISTIC FOR J=L+0.5
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
      PARAMETER (QD0=0.182D0)
      PARAMETER (ALF2=5.325D-5)
C
      COMMON /QDTS/QDTS(0:30),NQDT
C
      IF(N.LE.0)THEN
        QD=ZERO
        QDT=ZERO
        RETURN
      ENDIF
C
      TZ0=NZ0
      NZ=NZ0-NE+1
      TZ=NZ
      IF(L.LT.0.OR.NE.LE.1)THEN
        QD=ZERO
      ELSE
        IF(NQDT.GT.L)THEN
          QD=QDTS(L)
        ELSE
          TL=L**3+1
          TE=NE
          QD=QD0*(TE**1.667D0-DONE)/(TZ0**0.667D0*TZ**0.333D0*TL)
        ENDIF
      ENDIF
      TN=N
      T3=TN*TN*TN
      TN=TN-QD
C
      QDT=-(TZ/TN)**2
      IF(KAPPA.EQ.0)RETURN
C
      IF(KAPPA.EQ.-1)THEN
        ESO=ZERO
      ELSE
        IF(L.EQ.0)THEN
          WRITE(6,*)'*** FN.QDT ERROR: L=0 FOR KAPPA.NE.-1'
          ESO=ZERO
        ELSE
          ESO=ALF2*NZ**4*KAPPA/(T3*L*(L+1)*(2*L+1))
        ENDIF
      ENDIF
C
      IF(L.EQ.0)THEN
        ED=ALF2*NZ**4/T3
      ELSE
        ED=ZERO
      ENDIF
C
      EM=-ALF2*NZ**4*(4*N/(L+DONE/DTWO)-3)/(4*N*T3)
C
      QDT=QDT+ESO+ED+EM
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION WSQ(A,B,C,D,F,ISIGN)
C
C-----------------------------------------------------------------------
C
C  CALCULATES 3*(2*F+1)*W(A,B,C,D,1,F)**2 WHERE W IS RACAH COEFFICIENT.
C NRB:
C  ISIGN IS THE SIGN OF W BEFORE SQUARING.
C
C-----------------------------------------------------------------------
C
      IMPLICIT INTEGER(A-N)
      IMPLICIT REAL*8 (W)
C
C
      ISIGN=(-1)**(F-A-C)
C
      IF(B-A)10,20,30
C
   10 IF(D-C)11,12,13
   20 IF(D-C)21,22,23
   30 IF(D-C)31,32,33
C
   11 WSQ=.75*DBLE((2*F+1)*(B+D+F+2)*(B+D+F+3)*(B+D-F+1)*(B+D-F+2))/
     +    DBLE((B+1)*(2*B+1)*(2*B+3)*(D+1)*(2*D+1)*(2*D+3))
      RETURN
   12 WSQ=.75*DBLE((2*F+1)*(B+D+F+2)*(-B+D+F)*(B-D+F+1)*(B+D-F+1))/
     +    DBLE((B+1)*(2*B+1)*(2*B+3)*D*(D+1)*(2*D+1))
      ISIGN=-ISIGN
      RETURN
   13 WSQ=.75*DBLE((2*F+1)*(-B+D+F-1)*(-B+D+F)*(B-D+F+1)*(B-D+F+2))/
     +    DBLE((B+1)*(2*B+1)*(2*B+3)*D*(2*D-1)*(2*D+1))
      RETURN
C
   21 WSQ=.75*DBLE((2*F+1)*(B+D+F+2)*(-B+D+F+1)*(B-D+F)*(B+D-F+1))/
     +    DBLE(B*(B+1)*(2*B+1)*(D+1)*(2*D+1)*(2*D+3))
      ISIGN=-ISIGN
      RETURN
   22 G=-(B*(B+1)+D*(D+1)-F*(F+1))
      IF(G.LT.0)ISIGN=-ISIGN
      WSQ=.75*DBLE((2*F+1)*G*G)/             !(2*F+1) IS NOT SQUARED-NRB
     +    DBLE(B*(B+1)*(2*B+1)*D*(D+1)*(2*D+1))
      RETURN
   23 WSQ=.75*DBLE((2*F+1)*(B+D+F+1)*(-B+D+F)*(B-D+F+1)*(B+D-F))/
     +    DBLE(B*(B+1)*(2*B+1)*D*(2*D-1)*(2*D+1))
      RETURN
C
   31 WSQ=.75*DBLE((2*F+1)*(-B+D+F+1)*(-B+D+F+2)*(B-D+F-1)*(B-D+F))/
     +    DBLE(B*(2*B-1)*(2*B+1)*(D+1)*(2*D+1)*(2*D+3))
      RETURN
   32 WSQ=.75*DBLE((2*F+1)*(B+D+F+1)*(-B+D+F+1)*(B-D+F)*(B+D-F))/
     +    DBLE(B*(2*B-1)*(2*B+1)*D*(D+1)*(2*D+1))
      RETURN
   33 WSQ=.75*DBLE((2*F+1)*(B+D+F)*(B+D+F+1)*(B+D-F-1)*(B+D-F))/
     +    DBLE(B*(2*B-1)*(2*B+1)*D*(2*D-1)*(2*D+1))
      RETURN
C
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION WSQ2(A,B,C,D,F,ISIGN)
C
C-----------------------------------------------------------------------
C
C  CALCULATES 3*(2F+1)*W(A,B,C,D,1,F)**2 WHERE W IS RACAH COEFFICIENT.
C NRB:
C  ISIGN IS THE SIGN OF W BEFORE SQUARING,
C*** AND *** A,B,C,D,E,F ARE INPUT TWICE THEIR ACTUAL VALUE.
C
C-----------------------------------------------------------------------
C
      IMPLICIT INTEGER(A-N)
      IMPLICIT REAL*8 (W)
C
C
      ISIGN=(-1)**((F-A-C)/2)
C
      IF(B-A)10,20,30
C
   10 IF(D-C)11,12,13
   20 IF(D-C)21,22,23
   30 IF(D-C)31,32,33
C
   11 WSQ2=.75*DBLE((F+1)*(B+D+F+4)*(B+D+F+6)*(B+D-F+2)*(B+D-F+4))/
     +    DBLE((B+2)*(B+1)*(B+3)*(D+2)*(D+1)*(D+3))/4.
      RETURN
   12 WSQ2=.75*DBLE((F+1)*(B+D+F+4)*(-B+D+F)*(B-D+F+2)*(B+D-F+2))/
     +    DBLE((B+2)*(B+1)*(B+3)*D*(D+2)*(D+1))/2.
      ISIGN=-ISIGN
      RETURN
   13 WSQ2=.75*DBLE((F+1)*(-B+D+F-2)*(-B+D+F)*(B-D+F+2)*(B-D+F+4))/
     +    DBLE((B+2)*(B+1)*(B+3)*D*(D-1)*(D+1))/4.
      RETURN
C
   21 WSQ2=.75*DBLE((F+1)*(B+D+F+4)*(-B+D+F+2)*(B-D+F)*(B+D-F+2))/
     +    DBLE(B*(B+2)*(B+1)*(D+2)*(D+1)*(D+3))/2.
      ISIGN=-ISIGN
      RETURN
   22 G=-(B*(B+2)+D*(D+2)-F*(F+2))
      IF(G.LT.0)ISIGN=-ISIGN
      WSQ2=.75*DBLE((F+1)*G*G)/                !(F+1) IS NOT SQUARED-NRB
     +    DBLE(B*(B+2)*(B+1)*D*(D+2)*(D+1))
      RETURN
   23 WSQ2=.75*DBLE((F+1)*(B+D+F+2)*(-B+D+F)*(B-D+F+2)*(B+D-F))/
     +    DBLE(B*(B+2)*(B+1)*D*(D-1)*(D+1))/2.
      RETURN
C
   31 WSQ2=.75*DBLE((F+1)*(-B+D+F+2)*(-B+D+F+4)*(B-D+F-2)*(B-D+F))/
     +    DBLE(B*(B-1)*(B+1)*(D+2)*(D+1)*(D+3))/4
      RETURN
   32 WSQ2=.75*DBLE((F+1)*(B+D+F+2)*(-B+D+F+2)*(B-D+F)*(B+D-F))/
     +    DBLE(B*(B-1)*(B+1)*D*(D+2)*(D+1))/2.
      RETURN
   33 WSQ2=.75*DBLE((F+1)*(B+D+F)*(B+D+F+2)*(B+D-F-2)*(B+D-F))/
     +    DBLE(B*(B-1)*(B+1)*D*(D-1)*(D+1))/4.
      RETURN
C
      END
C
C***********************************************************************
C
C-----------------------------------------------------------------------
C
C A COLLECTION OF INTEGER ORDER MODIFIED BESSEL FUNCTION ROUTINES.
C
C-----------------------------------------------------------------------
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSK(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION K(N,X) GENERATED BY UP RECURRENCE RELATION:
C
C K(N+1,X)=(2N/X)*K(N,X)+K(N-1,X)                          (A&S: 9.6.26)
C
C-----------------------------------------------------------------------
C
      IF(N.LT.2)STOP 'BAD ARGUMENT N IN BESSK'
C
      TOX=2.0D0/X
      BKM=BESSK0(X)
      BK=BESSK1(X)
      DO J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
      ENDDO
      BESSK=BK
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSK0(X)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION K(0,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.5 AND 9.8.6
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     X                      0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     X                -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
C
      IF (X.LE.2.0D0) THEN
        Y=X*X/4.0D0
        BESSK0=(-LOG(X/2.0D0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*
     X         (P6+Y*P7))))))
      ELSE
        Y=(2.0D0/X)
        BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X          Q7))))))
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSK1(X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION K(1,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.7 AND 9.8.8
C
C-----------------------------------------------------------------------
C
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     X       -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     X                0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
C
      IF (X.LE.2.0D0) THEN
        Y=X*X/4.0D0
        BESSK1=(LOG(X/2.0D0)*BESSI1(X))+(1.0D0/X)*(P1+Y*(P2+Y*(P3+Y*
     X         (P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0D0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X          Q7))))))
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSI(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION I(N,X) GENERATED BY DOWN RECURRENCE RELATION:
C
C I(N+1,X)=-(2N/X)*I(N,X)+I(N-1,X)                         (A&S: 9.6.26)
C
C-----------------------------------------------------------------------
C
      PARAMETER (BIGNO=1.0D10)
      PARAMETER (BIGNI=1.0D-10)
C
      PARAMETER (IACC=40)
C
      IF (N.LT.2)STOP 'BAD ARGUMENT N IN BESSI'
C
      IF (X.EQ.0.0D0) THEN
        BESSI=0.0D0
      ELSE
        TOX=2.0D0/ABS(X)
        BIP=0.0D0
        BI=1.0D0
        BESSI=0.0D0
        M=2*((N+INT(SQRT(DFLOAT(IACC*N)))))
        DO J=M,1,-1
          BIM=BIP+DFLOAT(J)*TOX*BI
          BIP=BI
          BI=BIM
          IF (ABS(BI).GT.BIGNO) THEN
            BESSI=BESSI*BIGNI
            BI=BI*BIGNI
            BIP=BIP*BIGNI
          ENDIF
          IF (J.EQ.N) BESSI=BIP
        ENDDO
        BESSI=BESSI*BESSI0(X)/BI
        IF (X.LT.0.0D0.AND.MOD(N,2).EQ.1) BESSI=-BESSI
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSI0(X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION I(0,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.1 AND 9.8.2
C
C-----------------------------------------------------------------------
C
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     X        1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     X       0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     X      0.2635537D-1,-0.1647633D-1,0.392377D-2/
C
      IF (ABS(X).LT.3.75D0)THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X         (Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSI1(X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION I(1,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.3 AND 9.8.4
C
C-----------------------------------------------------------------------
C
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     X         0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     X       -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     X                   -0.2895312D-1,0.1787654D-1,-0.420059D-2/
C
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X         (Q7+Y*(Q8+Y*Q9))))))))
        IF(X.LT.0.0D0)BESSI1=-BESSI1
      ENDIF
C
      RETURN
      END

