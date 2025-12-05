$PROBLEM    CeftriaxPro Base Model Development
$INPUT      ID TIME DOSE=AMT RATE CONC=DV CONC_ID BW SEX ALB SCR eGFR
            AGE ETHNI PR DVID CMT MDV EVID
$DATA      bs_pr1_1646.dta IGNORE=@
$SUBROUTINE ADVAN13 ;; general nonlinear model
            TOL=9
$MODEL      NCOMP=2 ;; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ;; first (central) compartment
            COMP=(PERIPH) ;; second (peripheral) compartment
$PK
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

EGFRCL = THETA(1)

; exponents for covariates
TVCL = THETA(2)
  CL = TVCL * (eGFR/90)**EGFRCL * EXP(ETA(1))
TVV1 = THETA(3)
  V1 = TVV1 * EXP(ETA(2))
TVV2 = THETA(4)
  V2 = TVV2
 TVQ = THETA(5)
   Q = TVQ

 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
  S1 = V1 ;; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2

; Binding parameters
TVKD = THETA(6)
  KD = TVKD
 TVN = THETA(7)
   N = TVN


$THETA  1.58922 ; ; EGFRCL
 (0,8.31463) ; _; CL
 (0,56.0892) ; __; V1
 (0,43.9815) ; ___; V2
 (0,20.5868) ; ____; Q
 (0,0.0413634) ; _____; KD
 (0,0.840144) ; ______; N
; chech phase 1 trial, SAME order as PK
$OMEGA  0.099116  ;    _______  ; not yet known, outcome picked for next models, put a relatively small number
 0.27993  ;   ________  ; not yet known, outcome picked for next models
 ;; variance covariance matrix for interindividual variability
$SIGMA  0.046373  ;  _________  ; EPS(1), proportional error, not yet known, outcome picked for next models
;; variance covariance matrix for residual error
$DES DADT(1) = -K10*A(1) -K12*A(1) +K21*A(2)  ;; ODE for central    compartment
     DADT(2) =            K12*A(1) -K21*A(2)  ;; ODE for peripheral compartment

$ERROR
IPRED_UNBOUND = F
IF (CONC_ID.EQ.2) IPRED = IPRED_UNBOUND ;; unbound
IF (CONC_ID.EQ.1) IPRED = IPRED_UNBOUND * (1 + (N * ALB)/(KD + IPRED_UNBOUND))

Y = IPRED + IPRED*ERR(1)
W = SQRT(IPRED**2*SIGMA(1,1))
IRES  = CONC - IPRED
IWRES = IRES / W

$ESTIMATION METHOD=1 INTERACTION ;; FOCE-I
            MAXEVAL=9999 SIG=3 PRINT=5

