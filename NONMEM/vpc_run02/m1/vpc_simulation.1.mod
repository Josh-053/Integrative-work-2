;; 1. Based on: run01
;; 2. Description: CeftriaxPro Base Model - Total Plasma Conc.
;; x1. TEAM2
;; 2025-11-30
$PROBLEM    CeftriaxPro Base Model Development
$INPUT      ID TIME DOSE=AMT RATE CONC=DV CONC_ID BW SEX ALB SCR eGFR
            AGE ETHNI PR DVID CMT MDV EVID
$DATA      ../../dataset.csv IGNORE=@
$SUBROUTINE ADVAN13 ;; general nonlinear model
            TOL=9
$MODEL      NCOMP=2 ;; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ;; first (central) compartment
            COMP=(PERIPH) ;; second (peripheral) compartment
$PK
TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1))
TVV1 = THETA(2)
  V1 = TVV1 * EXP(ETA(2))
TVV2 = THETA(3)
  V2 = TVV2
 TVQ = THETA(4)
   Q = TVQ
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
  S1 = V1 ;; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2

; Binding parameters
TVKD = THETA(5)
  KD = TVKD
 TVN = THETA(6)
   N = TVN

$THETA  (0,7.23757) ; ; CL
 (0,57.2322) ; _; V1
 (0,48.7082) ; __; V2
 (0,20.3404) ; ___; Q
 (0,0.0428922) ; ____; KD
 (0,0.85138) ; _____; N
; chech phase 1 trial, SAME order as PK
$OMEGA  0.211723  ;     ______  ; not yet known, outcome picked for next models, put a relatively small number
 0.270926  ;    _______  ; not yet known, outcome picked for next models
 ;; variance covariance matrix for interindividual variability
$SIGMA  0.0566946  ;   ________  ; EPS(1), proportional error, not yet known, outcome picked for next models
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

$SIMULATION (1879376191) ONLYSIMULATION NSUBPROBLEMS=200
$TABLE      ID MDV TIME CONC ONEHEADER NOPRINT NOAPPEND IDFORMAT=I
            FILE=npctab.dta

