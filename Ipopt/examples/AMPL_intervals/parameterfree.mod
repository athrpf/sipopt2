# ===================================================================
# dynamic optimization formulation of the hicks-ray reactor
# model declaration
# victor m zavala  march 2006
# adapted for asNMPC by Hans Pirnay 2009, 2011
# ===================================================================

# define indexes and general variables

param nfe >= 1 integer       ;
param ncp >= 1 integer       ;
param nint >= 1 integer ;

# define mathematical model parameters

param xAL_init    ;
param xAU_init;
param xBL_init;
param xBU_init;
param r1        ;
param r2        ;
param r3        ;

# define dimensions for all indexed variables

set fe := 1..nfe ;  # number of finite elements
set cp := 1..ncp ;  # number of collocation points
set int := 1..nint ; # number of intervallizations

param a{cp,cp} ;    # collocation matrix
param h{fe}    ;    # finite element length

# define the decision variables

var xAL {int,fe,cp}  >= 0 <=1          ;
var xAU {int,fe,cp}  >= 0 <=1          ;
var xBL {int,fe,cp}  >= 0 <=1          ;
var xBU {int,fe,cp}  >= 0 <=1          ;
var xCU {int} ;
var time  >= 1e-3 := 1         ; # control - must be the same for all intervallizations

var p1L{int};
var p1U{int};
var p2L{int};
var p2U{int};

# states first order derivatives
var xALdot{k in int, i in fe, j in cp} = -p1U[k]*xAU[k,i,j];
var xAUdot{k in int, i in fe, j in cp} = -p1L[k]*xAL[k,i,j];
var xBLdot{k in int, i in fe, j in cp} = p1L[k]*xAL[k,i,j] - p2U[k]*xBU[k,i,j];
var xBUdot{k in int, i in fe, j in cp} = p1U[k]*xAU[k,i,j] - p2L[k]*xBL[k,i,j];

var xBmin;


# collocation equations
fecolcxAL{kk in int, i in fe diff{1},j in cp}: xAL[kk,i,j] = xAL[kk,i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*xALdot[kk,i,k];
fecolcxAU{kk in int, i in fe diff{1},j in cp}: xAU[kk,i,j] = xAU[kk,i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*xAUdot[kk,i,k];
fecolcxBL{kk in int, i in fe diff{1},j in cp}: xBL[kk,i,j] = xBL[kk,i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*xBLdot[kk,i,k];
fecolcxBU{kk in int, i in fe diff{1},j in cp}: xBU[kk,i,j] = xBU[kk,i-1,ncp]+time*h[i]*sum{k in cp} a[k,j]*xBUdot[kk,i,k];


fecolxAL0{kk in int, i in 1..1,j in cp}:       xAL[kk,i,j] = xAL_init+time*h[i]*sum{k in cp} a[k,j]*xALdot[kk,i,k];
fecolxAU0{kk in int, i in 1..1,j in cp}:       xAU[kk,i,j] = xAU_init+time*h[i]*sum{k in cp} a[k,j]*xAUdot[kk,i,k];
fecolxBL0{kk in int, i in 1..1,j in cp}:       xBL[kk,i,j] = xBL_init+time*h[i]*sum{k in cp} a[k,j]*xBLdot[kk,i,k];
fecolxBU0{kk in int, i in 1..1,j in cp}:       xBU[kk,i,j] = xBU_init+time*h[i]*sum{k in cp} a[k,j]*xBUdot[kk,i,k];

# state constraint
xCU_def{k in int} : xCU[k] = 1 - xAL[k,nfe,ncp] - xBL[k,nfe,ncp];
xCmax{k in int} : xCU[k] <= 0.2;

# objective function...

xBmin_eqn{kk in int}: xBmin <= xBL[kk, nfe, ncp];

minimize cost: -xBmin;

#-- end of the hicks.mod file --
