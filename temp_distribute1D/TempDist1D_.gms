$ontext
An insulated bar with an initial temperature distribution at t=0, having ends
that are subsequently maintain at temperature which may be function of time.
PDE:
du/dt = d^2u/dx^2, for 0<x<1, 0<t<T
I.C.: u(x,0)=f(x), 0<=x<=1
B.C.: u(0,t)=g0(t), 0<t<=T
      u(1,t)=g1(t), 0<t<=T
$offtext

Set T /t1*t7/ ;
Set X /x1*x6/ ;
*Determination of zone for temperature distribute equation
Set region(X,T) network of grid point;
    region(X,T) = yes;
    region(X,T)$(ord(X)=1) = no;
    region(X,T)$(ord(X)=card(X)) = no;
    region(X,T)$(ord(T)=1) = no;
    

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
    u.l(X,T)$(region(X,T)) = 0;
    
    u.up(X,T) = 100;
    u.lo(X,T) = 0;
*Boundary conditions
    u.fx('x1',T) = 100;
    u.fx('x6 ',T) = 100;
    u.fx('x1','t1') = 0;
    u.fx('x6','t1') = 0;
    
Equations
          eobj       objective function
          For_u(X,T) ;
    For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= r*u(X-1,T)+(1-2*r)*u(X,T)+r*u(X+1,T) ;
    
    eobj .. obj =e= u('x3','t7')+u('x4','t7');
    
Model temp_distribtion1D /all/;

*$onecho >bench.opt
*  solvers conopt knitro minos snopt
*$offecho

*temp_distribtion1D.optfile=1;
*option nlp=bench;

Solve temp_distribtion1D using nlp maximizing obj;
execute_unload "temp_distribtion1D_result.gdx", u, region, T, X, r;