function xdot=Vtol_NonLinear(u1,u2,u3,u4,u5,u6,u7,u8)
// This is nonlinear Pendulum Model

// Load the parameters
exec('Vtol_Parameters.sce', -1);

// state variables
x=u1;		
y=u2;
theta=u3;
varx=u4;
vary=u5;
vartheta=u6;


// control variables
u_1=u7;	// fuerza F1
u_2=u8; // u2=F2-mg


// Modelo MoDiCA-X
//Estado 1-posicion eje x
e1dot= x;

// estado 2- elevacion eje y
e2dot= y;

//estado 3- angulo de inclinacion
e3dot= theta;

//estado 4- variacion x, velocidad lineal horizontal
e4dot = -g*sin(theta)-(c/m)*varx+(1/m)*u_1*cos(theta)-(1/m)*u_2*sin(theta);

//estadp 5- variacion de y velocidad lineal vertical
e5dot = g*(cos(theta)-1)-(c/m)*vary+(1/m)*u_1*sin(theta)+(1/m)*u_2*cos(theta);

//estadp 6- variacion de theta velocidad anguloar
e6dot = (r/J)*u_1;

//Salida xdot
xdot =[e1dot;e2dot;e3dot;e4dot;e5dot;e6dot]; 

endfunction
