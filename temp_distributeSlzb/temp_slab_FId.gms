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
parameter init_state(X) slab initial temperature state ;
parameter bc_l(T) left boundary condition of temperature ;
parameter bc_r(T) right boundary condition of temperature ;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load dt dx init_state bc_r bc_l
$gdxin

parameter r kinematic viscosity ;
          r = dt/(dx*dx) ;

positive variable
          FId ;
variables
          obj        objective variable
          u(X,T)     temperature distribution ;
          
* Variable bounds and initialization
    u.l(X,T) = 0;
    u.lo(X,T) = 0;
    u.fx(X,'t1') = init_state(X) ;
    u.fx(X,'t2') = init_state(X) ;
execute_unload "init_state_load.gdx", u;
    
    
* Boundary conditions
    u.fx('x1',T) = bc_l(T) ;
    u.fx('x20 ',T) = bc_r(T) ;
* Uncertain init state

    
Equations
          eobj       objective function
          For_u(X,T) ;
    For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= r*u(X-1,T)+(1-2*r)*u(X,T)+r*u(X+1,T)+u(X,T) ;
    
    eobj .. obj =e= sum((X,T)$(for_sum(X,T)), (u(X,T)));
    
Model temp_slab /all/;

*$onecho >bench.opt
*  solvers conopt knitro minos snopt
*$offecho

*temp_slab_FId.optfile=1;
*option nlp=bench;

Solve temp_slab using nlp maximizing obj;
execute_unload "temp_slab.gdx", u, region, T, X, r;