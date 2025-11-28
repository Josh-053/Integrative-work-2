;; 1. Based on: run00
;; 2. Description: CeftriaxPro Base Model - One Compartment
;; x1. TEAM2
;; 2025-11-24

$PROBLEM CeftriaxPro Base Model Development

$INPUT ID TIME DOSE=AMT RATE CONC=DV CONC_ID BW SEX ALB SCR eGFR AGE ETHNI PR DVID CMT MDV EVID

$DATA
dataset.csv
IGNORE=@
IGNORE=(CONC_ID.EQ.1) ;; ignores rows showing total plasma concentration

$SUBROUTINES
ADVAN1 ;; 1-compartment model
TRANS2 ;; CL, V

$PK
TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1))
 TVV = THETA(2)
   V = TVV * EXP(ETA(2))
   K = CL/V
  S2 = V/1                ;; scale prediction based on DOSE (mmol) and DV (mmol/L)

$THETA; explore dataset & literature, SAME order as PK
(0, 1)  ;; based on phase 1 trial
(0, 70) ;; based on phase 1 trial

$OMEGA ;; variance covariance matrix for interindividual variability
1;; not yet known, outcome picked for next models, put a relatively small number
1;; not yet known, outcome picked for next models

$SIGMA;; variance covariance matrix for residual error
1 ;; EPS(1) for proportional error, not yet known, outcome picked for next models
1 ;; EPS(2) for additive error, not yet known, outcome picked for next models

$ERROR
IPRED = F
Y = IPRED + IPRED*ERR(1) + ERR(2)
W = SQRT((IPRED*SQRT(SIGMA(1,1)))**2 + (SQRT(SIGMA(2,2)))**2)
IRES  = CONC - IPRED
IWRES = IRES / W

$EST
METHOD=1 INTERACTION;; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COV PRINT=E;; second derivatives using -2log(likelihood), blank means sandwich method

$TABLE ;; output table for standard outcomes
ID TIME DV EVID PRED IPRED WRES RES CWRES NOPRINT ONEHEADER FILE=run00_sdtab

$TABLE ;; output table for PK parameters
ID CL V K NOPRINT NOAPPEND ONEHEADER FILE=run00_patab

$TABLE ;; output table for categorical covariates
ID SEX ETHNI NOPRINT NOAPPEND ONEHEADER FILE=run00_catab

$TABLE ;; output table for continuous covariates
ID BW ALB SCR eGFR AGE NOPRINT NOAPPEND ONEHEADER FILE=run00_cotab

