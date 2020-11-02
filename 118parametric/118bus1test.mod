## The unit commitment problem for the six-bus system##



### MASTER PROBLEM ###
param nCutOpt >= 0 integer;
# the number of cuts
param NUNITS;            # Number of units

param NSEGMENTS;         # Number of segments of energy blocks  second-stage?

param NHOURS;            # Number of hours for each day (at least 24 hours)

param N_DAYS;

param N_GEN;    # number of generators we consider; no use?
param N_SCEN;   # number of scenarios          
    
 
param NLINES;            # number of lines      
param NBUSES;            # number of buses

param NWINDS;            # number of wind farms



set I := 1..NUNITS;

set K := 1..NSEGMENTS;
 # second-stage
set T := 1..NHOURS;
set TT := 1..(NHOURS*N_DAYS);

set S := 1..N_SCEN;
   
# set RK := 1..RNSEGMENTS;
  # reserved segments
# set ZN := 1..NZONES;

  # zone 
set L := 1..NLINES;         # transmission          
set B := 1..NBUSES;

set WG := 1..NWINDS;        # number of wind farms 
set day := 1..N_DAYS;

# param ISTATE {I};       # Number of hours on (+) or off (-) before hour 1 


 
param MDT {I};          # Minimum down time 

param MUT {I};          # Minimum up time
 
param INOFF {I};        # Number of hours that the unit has been off before hour 1

param MBOFF {I};        # Number of hours that the unit must be off from hour 1

param MBON {I};         # Number of hours that the unit must be on from hour 1

param UO {I};           # Initial state
 
 
# appearing in the cuts

param PR {I,K};         # Size of an energy block

 param PMAX {I};         # Maximum power output
 
param PMIN {I} ;        # Minimum power output  
 
param RD {I};           # Ramp-down limit
 
param RSD {I};          # Ramp shut-down limit 

param RSU {I};          # Ramp startup limit

param RU {I};           # Ramp-up limit

param D{T};            # Demand
param Dofdays{TT};      # total demand of days

param LOADDIS {B};                  # load bus distribution factor 

param LOAD {b in B, t in T} = D[t] * LOADDIS[b];   # nodal load

param LINECAP {L};                  # branch flow capacity
param WC{WG,T,S};  # wind actual capacity of a scenario;
param WCofdays{WG,TT,S};  # total wind actual capacity of a scenario;

# main obj cost coefficients
param NLC {I};          # min output cost  
param SUC {I};         # Start up cost 

param SDC {I};          # Shut-down cost 

param SF {L,B};                     # shift factor   

# define the duals
param gen1dual{I,T,S,1..nCutOpt}  default 0;
param delta1dual{I,K,T,S,1..nCutOpt}  default 0;
param wcap1dual{WG,T,S,1..nCutOpt} default 0;
param flowconstdual{T,S,1..nCutOpt} default 0;
param DCfldual{L,T,S,1..nCutOpt} default 0;
 param linecap1dual{L,T,S,1..nCutOpt} default 0;
 param  linecap2dual{L,T,S,1..nCutOpt} default 0;   
 param OutUpdual{I,T,S,1..nCutOpt} default 0; 
  param OutLowdual{I,T,S,1..nCutOpt} default 0;  
  param rampup1dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn1dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;     
   param rampup2dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn2dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;   
  param rampup3dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn3dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;   
  param rampup4dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn4dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;   
  param rampup5dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn5dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;   
  param rampup6dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn6dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;   
  param rampup7dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;  
  param rampdn7dual{I,2..NHOURS,S,S,1..nCutOpt} default 0;   






## define parameters
# the first three sets of parameters are binary parameters determining on-off statues of a generator.
var u {I,T} binary ;    #on or off

var y {I,T} binary ;   #shut down (SD) or not

var z {I,T} binary ;    #start up (SU) or not

var Min_Stage2_Profit ;

var netinj {B,T,S};
# objective function
minimize Expected_Profit:
   sum {i in I, t in T} SUC[i]*z[i,t]+ sum {i in I, t in T} SDC [i]*y[i,t]+
   sum {i in I, t in T} NLC[i]*u[i,t] + Min_Stage2_Profit;






# constraints:
################################
# Key constrains of the three-variable formulation of UC

# corresponding to equation 4: determinting statues at time t of generator i.
 # for t>=2 determine the status of a generator;
subject to 
status1 {i in I, t in T: t>=2}: u[i,t] - u[i,t-1] = z[i,t]- y[i,t];
# determine the on\off status of a geenrator at time UO[i] is the initial status of the ith generator.  
subject to statusI {i in I}: u[i,1]-UO[i] = z[i,1] - y[i,1]; 

#start-up constraints for a generator at t>=2
subject to startup1 {i in I, t in T: t>=2}: 1- u[i,t-1] >= z[i,t] ;
# start-up constraint for a generator at t=1
subject to startupI {i in I}: 1- UO[i] >= z[i,1];
# shut-down constraints for a generator at t>=2
subject to shutdown1 {i in I, t in T: t>=2}: u[i,t-1] >= y[i,t]; 

# shut-down constraints for a generator at t=1
subject to shutdownI {i in I}: UO[i] >= y[i,1] ;

# Min-up/down constraints Equation 5-6
subject to  upI{i in I, t in 1..MBON[i]}: u[i,t] = 1; 
# by MBON definition, this constraint means initially must be on for MBON hours

subject to downI
{i in I, t in 1..MBOFF[i]}: u[i,t] = 0; 
# by MBOFF definition, this constraint means initially must be on for MBOFF hours


subject to startU {i in I, t in MUT[i]..NHOURS}:  u[i,t] >= sum{ tp in (t-MUT[i]+1)..t: tp >= 1} z[i,tp];
#if the unit i is on at time t, it cannot be started up more than once since t-MUT[i]+1.


subject to startD {i in I, t in MUT[i]..NHOURS, tp in (t-MUT[i]+1)..t: tp >= 1}: u[i,t] >= z[i,tp];
    #if the unit i is off at time t, then it should not be turned on since t-MUT[i]+1.

    # Otherwise it violates the minimum up requirement.

subject to DownU {i in I, t in MDT[i]..NHOURS, tp in (t-MDT[i]+1)..t: tp >= 1}: 1- u[i,t] >= y[i,tp]; 
# if the unit i is on at time t, then it should not been shut down since t-MDT[i]+1.

subject to DownD{i in I, t in MDT[i]..NHOURS}:  1- u[i,t] >= sum{ tp in (t-MDT[i]+1)..t: tp >= 1} y[i,tp];  
 #if the unit i is off at time t, it cannot be shut down more than once since t-MDT[i]+1.


# optimality cuts

subj to Cut_Defn {ncut in 1..nCutOpt: ncut>=2}:
#   Min_Stage2_Profit >= 
#      sum {ng in I, t in T, ns in 1..N_SCEN} low_price[ng,t,ns,k] * (lout[ng]*Y[ng,t]) +
#      sum {ng in 1..N_GEN, t in 1..T, ns in 1..N_SCEN} up_price[ng,t,ns,k] * (-uout[ng]*Y[ng,t]) +
#      sum {t in 1..T, ns in 1..N_SCEN}
#         bal_price[t,ns,k] *(Demand[ns,t]) ;

Min_Stage2_Profit >= 1/N_SCEN*(
                     sum{i in I, t in T, s in S} gen1dual[i,t,s,ncut]*(PMIN[i]*u[i,t])+
                     sum{i in I, k in K, t in T, s in S} delta1dual[i,k,t,s,ncut]*(-PR[i,k]) +
                     sum{w in WG, t in T,s in S} wcap1dual[w,t,s,ncut]*(- WC[w,t,s])+
                     sum{t in T, s in S} flowconstdual[t,s,ncut]*(sum{b in B} LOAD[b,t])+
                     sum{l in L, t in T, s in S} DCfldual[l,t,s,ncut]*(sum{b in B} SF[l,b] * netinj[b,t,s])+
                     sum{l in L, t in T, s in S} linecap1dual[l,t,s,ncut]*(-LINECAP [l])+
                     sum{l in L, t in T, s in S} linecap2dual[l,t,s,ncut]*(-LINECAP [l])+   
                     sum{i in I, t in T, s in S} OutUpdual[i,t,s,ncut]*(-PMAX[i]*u[i,t])+  
                     sum{i in I, t in T, s in S} OutLowdual[i,t,s,ncut]*(PMIN[i]*u[i,t]) );  
                  #   sum{i in I, t in T, s in S: t>=2} 
                  #  rampupdual[i,t,s,ncut]*(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+  
       #          + sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==0} rampup1dual[i,t,s,ss,nCutOpt] *(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+ 
 # sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMAX[i]} rampup2dual[i,t,s,ss,nCutOpt] *(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+ 
 #  sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMIN[i]} rampup3dual[i,t,s,ss,nCutOpt]*(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+ 
 # sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMIN[i]} rampup4dual[i,t,s,ss,nCutOpt] *(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+ 
 # sum  {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMAX[i]} rampup5dual[i,t,s,ss,nCutOpt] *(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+ 
 #  sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==0} rampup6dual[i,t,s,ss,nCutOpt] *(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+  
 #  sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==0,p[i,t-1,ss]==0} rampup7dual[i,t,s,ss,nCutOpt] *(-RU[i]*u[i,t-1] - RSU[i]*z[i,t])+ 
 # sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==0} rampdn1dual[i,t,s,ss,nCutOpt] *(- RD[i]*u[i,t] - RSD[i]*y[i,t])+
 # sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMAX[i]} rampdn2dual[i,t,s,ss,nCutOpt] *(- RD[i]*u[i,t] - RSD[i]*y[i,t])+
 #  sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMIN[i]} rampdn3dual[i,t,s,ss,nCutOpt] *(- RD[i]*u[i,t] - RSD[i]*y[i,t])+
 #   sum   {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMIN[i]} rampdn5dual[i,t,s,ss,nCutOpt]*(- RD[i]*u[i,t] - RSD[i]*y[i,t])+
 #   sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMAX[i]} rampdn4dual[i,t,s,ss,nCutOpt] *(- RD[i]*u[i,t] - RSD[i]*y[i,t])+
 #  sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==0} rampdn6dual[i,t,s,ss,nCutOpt] *(- RD[i]*u[i,t] - RSD[i]*y[i,t])+
 #  sum {i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==0,p[i,t-1,ss]==0} rampdn7dual[i,t,s,ss,nCutOpt] *(- RD[i]*u[i,t] - RSD[i]*y[i,t]) );


             #     sum{i in I, t in T, s in S: t>=2} rampdndual[i,t,s,ncut]*(- RD[i]*u[i,t] - RSD[i]*y[i,t])   );                 
              




# ----------------------------------------
### SUBPROBLEMS ###
# parameter definition
 param MC {I,K};         # Marginal cost of a block 
# param MCL {I};         # Marginal cost of a block 
param CENS;              # Cost of unserved energy
param CWC;  # wind curtailment cost

param GENBUS {I};                   # generator bus location id            

param GENMAP {i in I, b in B} = if b==GENBUS[i] then 1 else 0; # generator bus map, 1 indictes generator on that bus

param LINEADM {L};                  # branch admittance
 #inv of xij in the document
param LINEBUSFROM {L};              # branch from bus

param LINEBUSTO {L};                # branch to bus



param line_map_aux1 {l in L, b in B} = if b==LINEBUSFROM[l] then 1 else 0;

param line_map_aux2 {l in L, b in B} = if b==LINEBUSTO[l] then -1 else 0;

param LINEMAP {l in L, b in B} = line_map_aux1[l,b] + line_map_aux2[l,b];  # line map, 1 indicates from bus, -1 indicates to bus

#param BUSREFID;                     # reference bus ID

#param PI = 4 * atan(1);             # constant pi for power angle


param WINDBUS {WG};                 # wind bus location id

param WINDMAP {w in WG, b in B} = if b==WINDBUS[w] then 1 else 0;



#param BMATRIX {1..(NBUSES-1),1..(NBUSES-1)};     # B matrix #reference bus is the last one

 
# variables
var p {I,T,S} >= 0;             
#power output get dispatched for generation unit i at time t under scenario s.

var delta {I,K,T,S} >= 0;         
#power output at block k for generation unit i at time t under scenario s.



var ens {B,T,S} >= 0;         # amount of energy not served.
var dd {B,T,S} >= 0;
var wg {WG,T,S} >= 0;             #wind generation at time t under scenario s. 
# var theta{B,T,S};  #theta factors #theta{B}=0 reference bus
var pf {l in L,T,S};                        # DC power flow 



# objective 
 minimize Exp_Stage2_Profit{s in S}: sum{i in I, t in T, k in K} MC[i,k]*delta[i,k,t,s] + sum{b in B, t in T} ens [b,t,s]*CENS+sum{w in WG, t in T} (WC[w,t,s]-wg[w,t,s])*CWC;
#minimize Exp_Stage2_Profit{s in S}: sum{i in I, t in T, k in K} MC[i,k]*delta[i,k,t,s] + sum{b in B, t in T} ens [b,t,s]*CENS+sum{w in WG, t in T} (-wg[w,t,s])*CWC;
#minimize Exp_Stage2_Profit{s in S}: sum{i in I, t in T} MCL[i]*(p[i,t,s]-PMIN[i]*u[i,t]) + sum{b in B, t in T} ens [b,t,s]*CENS+sum{w in WG, t in T} (WC[w,t,s]-wg[w,t,s])*CWC;
# operational cost unserved cost
                           



##### constraints:
## linearlization constraint 
subject to gen1 {i in I, t in T, s in S}: p[i,t,s] -sum{k in K} delta[i,k,t,s]= PMIN[i]*u[i,t];


subject to delta1 {i in I, k in K, t in T, s in S}: -delta[i,k,t,s] >= -PR[i,k];



## wind power capacity constraints:
subject to wcap{w in WG, t in T, s in S}:  -wg[w,t,s] >= - WC[w,t,s];


# nodal power balance, Equation~(7).
subject to  flowconst{t in T, s in S}: sum{i in I} p[i,t,s] +sum{w in WG} wg[w,t,s] + sum{b in B} ens[b,t,s]= sum{b in B} LOAD[b,t];

# DC flow constant Equation~(8).
#subject to DCfl{l in L, t in T, s in S}: pf[l,t,s] = sum{b in B} SF[l,b] * (sum{i in I} GENMAP[i,b]*p[i,t,s]+sum{w in WG} WINDMAP[w,b]*wg[w,t,s]+ens[b,t,s]-LOAD[b,t]);
#subject to theB{t in T, s in S}: theta[NBUSES,t,s]=0;
# solve theta
# subject to dcflow{l in L, t in T, s in S}: pf[l,t,s] = LINEADM[l] * sum {b in B: LINEMAP[l,b] <> 0} theta[b,t,s]*LINEMAP[l,b];
subject to dem00 {b in B, t in T, s in S}: dd[b,t,s] + ens[b,t,s] = LOAD[b,t];

subject to inj {b in B, t in T, s in S}: netinj[b,t,s] = sum{i in I} GENMAP[i,b] * p[i,t,s] + sum{w in WG} WINDMAP[w,b] * wg[w,t,s] - dd[b,t,s];
subject to DCfl {l in L, t in T, s in S}: pf[l,t,s] = sum{b in B} SF[l,b] * netinj[b,t,s];


# calculate the actual DC flow per line
#subject to linecap1{l in L, t in T, s in S}: pf[l,t,s] <=LINECAP [l];
# subject to linecap2{l in L, t in T, s in S}: pf[l,t,s] >=-LINECAP [l];

subject to linecap1{l in L, t in T, s in S}:  -pf[l,t,s]>=-LINECAP [l]; 
subject to linecap2{l in L, t in T, s in S}:  pf[l,t,s]>=-LINECAP [l]; 



 






# Output conditions Equation(9)

subject to OutUp{i in I, t in T, s in S}:  -p[i,t,s] >= -PMAX[i]*u[i,t] ;



subject to OutLow{i in I, t in T, s in S}:  p[i,t,s] >= PMIN[i]*u[i,t] ;




# Ramp up and down Equation(10)
# Ramp up

#subject to rampup1{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==0}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];

#subject to rampup2{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMAX[i]}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];

#subject to rampup3{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMIN[i]}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];
#subject to rampup4{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMIN[i]}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];

#subject to rampup5{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMAX[i]}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];

#subject to rampup6{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==0}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];

#subject to rampup7{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==0,p[i,t-1,ss]==0}: -p[i,t,s]+p[i,t-1,ss] >= -RU[i]*u[i,t-1] - RSU[i]*z[i,t];




# Ramp down


#subject to rampdn1{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==0}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];

#subject to rampdn2{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMAX[i]}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];

#subject to rampdn3{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMAX[i],p[i,t-1,ss]==PMIN[i]}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];

#subject to rampdn4{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMAX[i]}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];

#subject to rampdn5{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==PMIN[i]}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];

#subject to rampdn6{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==PMIN[i],p[i,t-1,ss]==0}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];

#subject to rampdn7{i in I, t in T, s in S, ss in S: t>=2, p[i,t,s]==0,p[i,t-1,ss]==0}:  -p[i,t-1,s]+ p[i,t,s] >=  - RD[i]*u[i,t] - RSD[i]*y[i,t];




#constrampup1 {i in I, t in T, s in S: t>=2}: pbar[i,t,s] <= p[i,t-1,s] + RU[i]*u[i,t-1] + RSU[i]*z[i,t];

#constrampup2 {i in I, s in S}: pbar[i,1,s] <= PO[i] + RU[i]*UO[i] + RSU[i]*z[i,1];



#constrampdn1 {i in I, t in T, s in S: t>=2}: p[i,t-1,s] <=  p[i,t,s]+ RD[i]*u[i,t] + RSD[i]*y[i,t];

#constrampdn2 {i in I, s in S}:  PO[i] <= p[i,1,s] + RD[i]*u[i,1] + RSD[i]*y[i,1];

#constrampdn3 {i in I, t in T, s in S: t>=2}: pbar[i,t-1,s] <= PMAX[i]*u[i,t] + RSD[i]*y[i,t];







#############################################################################################################


# Evaluation
### SUBPROBLEMS ###
param N_SCENT;   # number of scenarios for evaluation
set SS:= 1..N_SCENT;
     

param WCT{WG,T,SS};  # wind actual capacity of a scenario;


 
# variables
var pT {I,T,SS} >= 0;             
#power output get dispatched for generation unit i at time t under scenario s.

var deltaT {I,K,T,SS} >= 0;         
#power output at block k for generation unit i at time t under scenario s.



var ensT {B,T,SS} >= 0;         # amount of energy not served.
var ddT {B,T,SS} >= 0;
var wgT {WG,T,SS} >= 0;             #wind generation at time t under scenario s. 
# var thetaT{B,T,SS};  #theta factors #theta{B}=0 reference bus
var pfT {l in L,T,SS};                        # DC power flow 
var netinjT {B,T,SS};



# objective 
 minimize Exp_Stage2_ProfitT{s in SS}: sum{i in I, t in T, k in K} MC[i,k]*deltaT[i,k,t,s] + sum{b in B, t in T} ensT [b,t,s]*CENS+sum{w in WG, t in T} (WCT[w,t,s]-wgT[w,t,s])*CWC;
#minimize Exp_Stage2_Profit{s in SS}: sum{i in I, t in T, k in K} MC[i,k]*delta[i,k,t,s] + sum{b in B, t in T} ens [b,t,s]*CENS+sum{w in WG, t in T} (-wg[w,t,s])*CWC;
#minimize Exp_Stage2_Profit{s in SS}: sum{i in I, t in T} MCL[i]*(p[i,t,s]-PMIN[i]*u[i,t]) + sum{b in B, t in T} ens [b,t,s]*CENS+sum{w in WG, t in T} (WCT[w,t,s]-wg[w,t,s])*CWC;
# operational cost unserved cost
                           



##### constraints:
## linearlization constraint 
subject to gen1T {i in I, t in T, s in SS}: pT[i,t,s] -sum{k in K} deltaT[i,k,t,s]= PMIN[i]*u[i,t];


subject to delta1T {i in I, k in K, t in T, s in SS}: -deltaT[i,k,t,s] >= -PR[i,k];



## wind power capacity constraints:
subject to wcapT{w in WG, t in T, s in SS}:  -wgT[w,t,s] >= - WCT[w,t,s];


# nodal power balance, Equation~(7).
subject to  flowconstT{t in T, s in SS}: sum{i in I} pT[i,t,s] +sum{w in WG} wgT[w,t,s] + sum{b in B} ensT[b,t,s]= sum{b in B} LOAD[b,t];

# DC flow constant Equation~(8).
#subject to DCflT{l in L,t in T, s in SS}: pfT[l,t,s] = sum{b in B} SF[l,b] * (sum{i in I} GENMAP[i,b]*pT[i,t,s]+sum{w in WG} WINDMAP[w,b]*wgT[w,t,s]+ensT[b,t,s]-LOAD[b,t]);
# subject to theBT{t in T, s in SS}: thetaT[NBUSES,t,s]=0;
# solve theta
# subject to dcflowT{l in L, t in T, s in SS}: pf[l,t,s] = LINEADM[l] * sum {b in B: LINEMAP[l,b] <> 0} thetaT[b,t,s]*LINEMAP[l,b];
subject to dem00T {b in B, t in T, s in SS}: ddT[b,t,s] + ensT[b,t,s] = LOAD[b,t];

subject to injT {b in B, t in T, s in SS}: netinjT[b,t,s] = sum{i in I} GENMAP[i,b] * pT[i,t,s] + sum{w in WG} WINDMAP[w,b] * wgT[w,t,s] - ddT[b,t,s];
subject to DCflT {l in L, t in T, s in SS}: pfT[l,t,s] = sum{b in B} SF[l,b] * netinjT[b,t,s];


# calculate the actual DC flow per line
subject to linecap1T{l in L, t in T, s in SS}: pfT[l,t,s] <=LINECAP [l];
subject to linecap2T{l in L, t in T, s in SS}: pfT[l,t,s] >=-LINECAP [l];

# subject to linecap1TT{l in L, t in T, s in SS}:  -LINEADM[l] * sum {b in B: LINEMAP[l,b] <> 0} thetaT[b,t,s]*LINEMAP[l,b]>=-LINECAP [l]; 
# subject to linecap2TT{l in L, t in T, s in SS}:  LINEADM[l] * sum {b in B: LINEMAP[l,b] <> 0} thetaT[b,t,s]*LINEMAP[l,b]>=-LINECAP [l]; 



 






# Output conditions Equation(9)

subject to OutUpT{i in I, t in T, s in SS}:  -pT[i,t,s] >= -PMAX[i]*u[i,t] ;



subject to OutLowT{i in I, t in T, s in SS}:  pT[i,t,s] >= PMIN[i]*u[i,t] ;









