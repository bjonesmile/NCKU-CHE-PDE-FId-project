$ontext
An insulated bar with an initial temperature distribution at t=0, having ends
that are subsequently maintain at temperature which may be function of time.
PDE:
T = (theta - theta0)/(theta - theta0), tau = alpha*t/L^2, X = x/L
dT/dtau = d^2T/dX^2, for 0<x<1, 0<t<T
I.C.: T(X,0)=f(x), 0<=x<=1
B.C.: T(0,tau)=g0(tau), 0<t<=1 : T = 0, for 0<=X<=1
      T(1,tau)=g1(tau), 0<t<=1 : T = 1, for X=0 and X=1
$offtext

Set T /t1*t350/ ;
Set X /x1*x20/ ;

*Determination of zone for temperature distribute equation
Set region(X,T) network of grid point;
    region(X,T) = yes;
    region(X,T)$(ord(X)=1) = no;
    region(X,T)$(ord(X)=card(X)) = no;
    region(X,T)$(ord(T)=1) = no;

Set for_sum(X,T) ;
    for_sum(X,T) = no ;
    for_sum(X,T)$(ord(T)=card(T)) = yes;

* Parameters
scalar dt  step space in x direction  ;
scalar dx  step space in y direction  ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load dt dx
$gdxin

parameter r kinematic viscosity ;
          r = dt/(dx*dx) ;

variables
          obj        objective variable
          u(X,T)     temperature distribution ;
          
* Variable bounds and initialization
    u.l(X,T) = 0;
    
    u.up(X,T) = 1;
    u.lo(X,T) = 0;
* Boundary conditions
    u.fx('x1',T) = 1;
    u.fx('x20 ',T) = 1;
    u.fx('x1','t1') = 0;
    u.fx('x20','t1') = 0;
    
Equations
          eobj       objective function
          For_u(X,T) ;
    For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= r*u(X-1,T)+(1-2*r)*u(X,T)+r*u(X+1,T) ;
    
    eobj .. obj =e= sum((X,T)$(for_sum(X,T)), (u(X,T)));
    
Model temp_slab /all/;

*$onecho >bench.opt
*  solvers conopt knitro minos snopt
*$offecho

*temp_distribtion1D.optfile=1;
*option nlp=bench;

Solve temp_slab using nlp maximizing obj;
execute_unload "temp_slab.gdx", u, region, T, X, r;