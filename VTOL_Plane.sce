//Parameters of the Plane
exec('Vtol_Parameters.sce', -1);

//States matrix
A=[0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   0 0 -g -c/m 0 0;
   0 0 0 0 -c/m 0;
   0 0 0 0 0 0];

B=[0 0;
   0 0;
   0 0;
   1/m 0;
   0 1/m;
   r/J 0];

C=[1 0 0 0 0 0;
   0 1 0 0 0 0];

D=[0 0;
   0 0];

sys = syslin('c',A,B,C,D);

//checking controllability and observability
[i,j] = size(A);
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

//Linear Quadratic Regulator
Q=eye(6,6); R=eye(2,2);

[K, X]=lqr(sys,Q,R);

sysn=syslin('c',sys.A-sys.B*K,sys.B,sys.C,zeros(2,2));

transMn=ss2tf(sysn);

tf11n=transMn(1,1); tf22n=transMn(2,2);
ss11n=tf2ss(tf11n);
ss22n=tf2ss(tf22n);

K1=1/(ss11n.D-ss11n.C*inv(ss11n.A)*ss11n.B);
K2=1/(ss22n.D-ss22n.C*inv(ss22n.A)*ss22n.B);

Pro=[K1 0; 0 K2];

sysnn=syslin('c',sys.A-sys.B*K,sys.B*Pro,sys.C,zeros(2,2));
f=0.01; t=0:0.1:180;

xideal=2*sin(2*%pi*f*t)+3;
yideal=2*sin(2*%pi*f*t+%pi/2)+3;
U=[xideal];
s=poly(0,'s')
deff('[u1,u2]=mymacro(t)',['u1=2*sin(2*%pi*f*t)+3'; 'u2=2*sin(2*%pi*f*t+%pi/2)+3'])
//[yr xr]=csim(mymacro,t,sysnn);



