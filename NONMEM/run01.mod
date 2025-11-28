;; 1. Based on: run00
;; 2. Description: CeftriaxPro Base Model - Two Compartment
;; x1. TEAM2
;; 2025-11-24

$PROBLEM CeftriaxPro Base Model Development

$INPUT ID TIME DOSE=AMT RATE CONC=DV CONC_ID BW SEX ALB SCR eGFR AGE ETHNI PR DVID EVID

$DATA
dataset.csv
IGNORE=@
IGNORE=(CONC_ID.EQ.1) ;; ignores rows showing total plasma concentration

$SUBROUTINES
ADVAN3 ;; 2-compartment model
TRANS4 ;; CL, V1, V2, Q

$PK
TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1))
TVV1 = THETA(2)
  V1 = TVV1 * EXP(ETA(2))
TVV2 = THETA(3)
  V2 = TVV2 * EXP(ETA(3))
 TVQ = THETA(4)
   Q = TVQ  * EXP(ETA(4))
   K = CL/V1
 K12 = Q/V1
 K21 = Q/V2
  S2 = V1/1                ;; scale prediction based on DOSE (mmol) and DV (mmol/L)

$ERROR
IPRED=F
Y=IPRED*(1+EPS(1)) ;;proportional = EPS(1)
W = SQRT((IPRED*SQRT(SIGMA(1,1)))**2)
IRES  = CONC - IPRED
IWRES = IRES / W

$THETA; explore dataset & literature, SAME order as PK
(0, 1.16, 3.26)  ;; From Telles et al. (1), max range doubled
(0, 10.04, 22.52);; From Telles et al. (1), max range doubled
(0, 9.54, 38.16) ;; From Simon et al. (2), upper bound = initial value x4 to account for ...
(0, 10.8, 44.32) ;; ... physiological differences

$OMEGA;; variance covariance matrix for interindividual variability
0.00938 ;; taken from run00
0.00306 ;; taken from run00
1
1

$SIGMA;; variance covariance matrix for residual error
0.1

$EST
METHOD=1 INTERACTION;; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5
NOABORT
POSTHOC

$COV ;; second derivatives using -2log(likelihood), blank means sandwich method

$TABLE;; output table
ID TIME DV DVID EVID PRED IPRED IWRES IRES WRES RES CWRES CL V1 V2 Q K K12 K21 SEX ETHNI BW ALB SCR eGFR AGE NOPRINT ONEHEADER FILE=run01

