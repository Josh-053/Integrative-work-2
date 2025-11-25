;; 1. Based on: run01
;; 2. Description: CeftriaxPro Base Model, edited theta ranges
;; x1. TEAM2
;; 2025-11-24

$PROBLEM CeftriaxPro Base Model Development

$INPUT ID TIME DOSE=AMT RATE CONC=DV CONC_ID=DROP BW SEX ALB SCR eGFR AGE ETHNI PR DVID EVID

$DATA
dataset.csv
IGNORE=@

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
Y=F*(1+EPS(1))+EPS(2);; combined: proportional = EPS(1), additive = EPS(2)

$THETA; explore dataset & literature, SAME order as PK
(0, 1.16, 3.26)  ;; From Telles et al. (1)
(0, 10.04, 22.52);; From Telles et al. (1)
(0, 1, 99)
(0, 1, 99)

$OMEGA;; variance covariance matrix for interindividual variability
0.1
0.1
0.1
0.1

$SIGMA;; variance covariance matrix for residual error
0.1
0.1

$EST;; estimation step
METHOD=1 ;; FOCE
INTERACTION
MAXEVALS=9999
PRINT=5

$COV ;; second derivatives using -2log(likelihood), blank means sandwich method

$TABLE;; output table
ID TIME DV DVID EVID PRED WRES RES CWRES CL V1 V2 Q K K12 K21 SEX ETHNI BW ALB SCR eGFR AGE NOPRINT ONEHEADER FILE=run02


