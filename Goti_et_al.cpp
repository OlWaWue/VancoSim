[PROB]

# Vancomycin popPK-Model

- Author "Oliver"
- Date: "`r date()`"
- Source: "Goti et al. Ther Drug Monit (2018) 40:212–221 "
- URL: "https://doi.org/10.1097/FTD.0000000000000490"

[SET] delta = 0.5, end=24

[THETA] @annotated
4.5 : THETA1 -  clearance (L/h)
58.4 : THETA2 - volume of central compartment (L)
38.4: THETA3 - volume of peripheral compartment (L)
6.5 : THETA4 - intercompartmental clearance (L/h)

0.8 : THETA5 - CrCl on CL
0.7 : THETA6 - DIAL on CL
0.5 : THETA7 - DIAL on Vc


[PARAM] @annotated
DIAL : 0  : 1 - yes, 0 - no
CrCL: 62  : Creatinineclearance (mL/min)
WT : 70 : body weight (kg)

ETA1: 0 : ETA for CL used in MAP
ETA2: 0 : ETA for Vc used in MAP
ETA3: 0 : ETA for Vp used in MAP


[CMT] @annotated
PER: peripheral compartment (ng)

[INIT] @annotated
CENT: 0 : central compartment (ng)

[OMEGA] @annotated @block
ECL: 0.1584                       : ETA on CL
EVC: 0.0000 0.6659                : ETA on Vc
EVP: 0.0000 0.0000 0.3260         : ETA on Vp


[SIGMA] @name SGMA @annotated
PROP: 0.0515 : proportional error
ADD : 11.56 : mg/L

[MAIN]
double TVCL = THETA1*pow( (CrCL/120),THETA5) * pow(THETA6,DIAL);
double TVVc = THETA2*(WT/70) * pow(THETA7,DIAL);

double CL_ind = TVCL*exp(ECL+ETA1);
double V1_ind = TVVc * exp(EVC+ETA2);
double V2_ind = THETA3 * exp(EVP+ETA3);
double Q_ind = THETA4;

[ODE]

dxdt_CENT = - (CL_ind/V1_ind) * CENT + (Q_ind/V2_ind) * PER - (Q_ind/V1_ind) * CENT;
dxdt_PER  = - (Q_ind/V2_ind) * PER + (Q_ind/V1_ind) * CENT;

[TABLE]
capture IPRED = CENT/V1_ind;
double CP = IPRED * (1+PROP)+ADD;

int i = 0;
while(CP < 0) {
  if(++i > 100) {
    mrg::report("Problem simulating positive CP");
        break;
      }
      simeps();
      CP = IPRED * (1+PROP)+ADD;
    }
  
[CAPTURE] @annotated
CP : Plasma concentration (mg/L)
CL_ind : individual clearance (L/h)
V1_ind : individual central volume (L)
V2_ind : individual peripheral volume (L)
Q_ind : individual Q (L/h)
