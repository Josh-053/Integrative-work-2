;; 1. Based on: run01
;; 2. Description: CeftriaxPro Base Model - Total Plasma Conc.
;; x1. TEAM2
;; 2025-11-30

$PROBLEM CeftriaxPro Base Model Development

$INPUT ID TIME DOSE=AMT RATE CONC=DV CONC_ID BW SEX ALB SCR eGFR AGE ETHNI PR DVID CMT MDV EVID

$DATA
dataset.csv
IGNORE=@

$SUBROUTINES
ADVAN13 ;; general nonlinear model
TOL=9

$MODEL
NCOMP=2                      ;; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ;; first (central) compartment
COMP=(PERIPH)                ;; second (peripheral) compartment

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

$THETA; chech phase 1 trial, SAME order as PK
(0, 7)  ;; CL
(0, 40) ;; V1
(0, 20) ;; V2
(0, 30) ;; Q
(0, 0.1) ;; KD
(0, 1) ;; N

$OMEGA ;; variance covariance matrix for interindividual variability
0.05;; not yet known, outcome picked for next models, put a relatively small number
0.05;; not yet known, outcome picked for next models

$SIGMA;; variance covariance matrix for residual error
0.1 ;; EPS(1), proportional error, not yet known, outcome picked for next models

$DES DADT(1) = -K10*A(1) -K12*A(1) +K21*A(2)  ;; ODE for central    compartment
     DADT(2) =            K12*A(1) -K21*A(2)  ;; ODE for peripheral compartment

$ERROR
IPRED_UNBOUND = F
IF (CONC_ID.EQ.2) THEN
   IPRED = IPRED_UNBOUND
ELSE
   IPRED = IPRED_UNBOUND * (1 + (N * ALB)/(KD + IPRED_UNBOUND))
ENDIF


Y = IPRED + IPRED*ERR(1)
W = SQRT(IPRED**2*SIGMA(1,1))
IRES  = CONC - IPRED
IWRES = IRES / W

$EST
METHOD=1 INTERACTION;; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL ;; new 

$TABLE ;; output table for standard outcomes
ID TIME DV EVID PRED CWRES MDV IPRED WRES RES CWRES NOPRINT ONEHEADER FILE=run02_sdtab

$TABLE ;; output table for PK parameters
ID CL V1 V2 Q K10 K12 K21 KD N NOPRINT NOAPPEND ONEHEADER FILE=run02_patab

$TABLE ;; output table for categorical covariates
ID SEX ETHNI NOPRINT NOAPPEND ONEHEADER FILE=run02_catab

$TABLE ;; output table for continuous covariates
ID BW ALB SCR eGFR AGE NOPRINT NOAPPEND ONEHEADER FILE=run02_cotab

