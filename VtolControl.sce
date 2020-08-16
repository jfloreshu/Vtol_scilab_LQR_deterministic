//Parameters of the Plane
exec('Vtol_Parameters.sce', -1);

//States matrix
Ap=[0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   0 0 -g -c/m 0 0;
   0 0 0 0 -c/m 0;
   0 0 0 0 0 0];

Bp=[0 0;
   0 0;
   0 0;
   1/m 0;
   0 1/m;
   r/J 0];

Cp=[1 0 0 0 0 0;
   0 1 0 0 0 0];

Dp=[0 0;
   0 0];

sys = syslin('c',Ap,Bp,Cp,Dp);

//checking controllability and observability
[i,j] = size(Ap);
// e=[B, AB, A^2 B,..., A^(n-1) B]
e = cont_mat(sys.A,sys.B);
rankC=rank(e);
if i == rankC then
    disp('Continuous System is Controllable');
end

// o=[C; CA; CA^2;...; CA^(n-1) ]
o = obsv_mat(sys.A, sys.C);
rankO=rank(o);
if j == rankO then
    disp('Continuous System is Observable');
end

tranM=ss2tf(sys); // Matriz de transferencia
disp('Matriz de Transferencia',tranM);

tfc11 = tranM(1,1);
tfc22 = tranM(2,2); 

/* Plot singular values of LTI the model */
tr = trzeros(sys)
w = logspace(-3,3);
sv = svplot(sys,w);
scf(1);
plot2d("ln", w, 20*log(sv')/log(10),leg="Input 1@Input 2")
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");

//Obtenciion de los polors  zeros del modelo de software
scf(2);
plzr(sys);
xtitle("Poles and zeros plot of system","Real", "Imaginarie");

//autovalores
evals=spec(Ap);
///////////////////////////////////---------------------------

//Augment Plant with Integrators at Plant Input
[ns,nc]=size(Bp); //ns= number of inputs; nc=number of controls

Ai=[Ap             Bp;
    0*ones(nc,ns) 0*ones(nc,nc)];



Bi=[0*ones(ns,nc); eye(nc,nc)];
    
Ci=[Cp 0*ones(nc,nc)];

Di=0*ones(nc,nc);

I=eye(nc,nc);

sysi=syslin('c',Ai,Bi,Ci,Di);
//View of the singular values of plant with integrator and the poles and
//zeros
/* Plot singular values of LTI the model */
tri = trzeros(sysi)
w = logspace(-3,3);
svi = svplot(sysi,w);
scf(3);
plot2d("ln", w, 20*log(svi')/log(10),leg="Input 1@Input 2")
xgrid(12)
xtitle("Singular values plot of Plant with Integrator","Frequency (rad/s)", "Amplitude (dB)");

//Obtenciion de los polors  zeros del modelo de software
scf(4);
plzr(sys);
xtitle("Design Plant with integrator:poles and zeros","Real", "Imaginarie");



//lqr controller calculation
//We use the ricatti equation for calculate de gain of the lqr controller
//for this we have  A'*X+X*A-X*B*X+C=0 for function X=riccati(A,B,C,'c','eigen')
C=0.8*eye(8,8);        //State Weighting Matrix
rho=1e-0;       //Cheap control recovery parameter 
                //The smaller the parameter, the better the recovery.
R = rho*eye(nc,nc);//Control Weigthing Matrix


//now we calculate B
B=Bi*inv(R)*Bi';

A=Ai;

//Solv the ricatti equation
X=riccati(A,B,C,'c','eigen');

//the value of the gain G
G=inv(R)*Bi'*X; //<--this value is important mtfk


//---------------------------------------------------------------------

//computing H

//The gain H
H=(ppol(Ai',Ci',[-10,-11,-12,-13,-14,-15,-16,-17]))';  

