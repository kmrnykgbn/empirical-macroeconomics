/*
 * Core Macroeconomics I HW4 3-2
 * RBC model
 */

var Y C A K L G I;
varexo e z;

parameters alpha beta rho_a rho_g delta psi eta omega Astar Gstar;
alpha = 0.333;
beta  = 0.99;
psi   = 4;
eta = 1;
rho_a   = 0.95;
rho_g   = 0.95;
delta = 0.025;
omega = 0.2;

// calculate steady state
Astar=1;
K_L = ((1/alpha)*((1/beta)+delta-1))^(1/(alpha-1));
Y_L = Astar*(K_L)^alpha;
G_L = omega*Y_L;
C_L=Astar*(K_L)^alpha-delta*K_L-G_L;
I_L=Astar*(K_L)^alpha-C_L-G_L;
Lstar=((1-alpha)/psi)^(1/(1+eta))*C_L^(-1/(1+eta))*Astar*(K_L)^(alpha/(1+eta));
Ystar=Y_L*Lstar;
Gstar=G_L*Lstar;
Kstar=K_L*Lstar;
Cstar=C_L*Lstar;
Istar=I_L*Lstar;

model;
exp(C(+1))/exp(C) = beta*(alpha*exp(A(+1))*(exp(K(+1))/exp(L(+1)))^(alpha-1)+1-delta);
psi*exp(C)*exp(L)^eta = (1-alpha)*exp(A)*(exp(K)/exp(L))^alpha;
exp(Y)=exp(A)*exp(K)^alpha*exp(L)^(1-alpha);
exp(I)=exp(Y)-exp(C)-exp(G);
exp(K)= exp(Y(-1)) + (1-delta)*exp(K(-1)) - exp(C(-1)) - exp(G(-1));
A = (1-rho_a)*(log(Astar)) + rho_a*A(-1) + e;
G = (1-rho_g)*(log(Gstar)) + rho_g*G(-1) + z;
end;

initval;
Y=log(Ystar);
C=log(Cstar);
A=log(Astar);
K=log(Kstar);
L=log(Lstar);
G=log(Gstar);
I=log(Istar);
end;

// Data generation
shocks;
var e; stderr 0.01;
var z; stderr 0.01;
end;

// check steady state value(not log)
exp_val = [Ystar;Cstar;Astar;Kstar;Lstar;Gstar;Istar]

steady;

check;

stoch_simul (order =1, simul_replic =1, periods =1000) C L;
save simuldata;

estimated_params;
psi, gamma_pdf, 4, 0.1;
eta, gamma_pdf, 1, 0.1;
rho_a, beta_pdf, 0.95, 0.02;
rho_g, beta_pdf, 0.95, 0.02;
stderr e, inv_gamma_pdf, 0.01, inf;
stderr z, inv_gamma_pdf, 0.01, inf;
end;

varobs C L;
estimation(order=1, datafile='simuldata.mat', nobs=1000,
mh_replic=2000, mh_nblocks=2, mh_jscale=0.8, mode_check);






