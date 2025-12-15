$PARAM
TVCL = 8.31,      // L/h
TVV1 = 56.1,      // L 
TVV2 = 44,        // 1/h
TVQ  = 20.6
TVKd = 0.0414
TVn  = 0.84

$CMT CENT PERIPH  // 2-compartmental model

$INPUT eGFR = 100 // covariate which has effect on PK parameters

$MAIN
double eGFREffCL = 1.59 ; // eGFR effect on CL

double CL  = TVCL * pow((eGFR/90),eGFREffCL) * exp(ETA(1)); // power function for scaling
double V1  = TVV1 * exp(ETA(2)) ;
double V2  = TVV2  ;
double Q   = TVQ   ;
double Kd  = TVKd  ;
double n   = TVn   ;

double k10 = CL/V1 ;
double k12 = Q/V1  ;
double k21 = Q/V2  ;

$OMEGA
0.0991 // CL
0.28   // V1

$SIGMA
0.0464 // proportional

$ODE
dxdt_CENT   =  - (k10) * CENT - (k12) * CENT + (k21) * PERIPH ;
dxdt_PERIPH =                   (k12) * CENT - (k21) * PERIPH ;

$TABLE 
double IPRED = CENT/V1              ; 
double DV    = IPRED * (1+EPS(1))   ;
double eGFRsim = eGFR               ;

$CAPTURE
IPRED DV eGFRsim
