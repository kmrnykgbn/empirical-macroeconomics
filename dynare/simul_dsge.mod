/*
 * Core Macroeconomics I HW4 4-4
 * DSGE model
 */

var Y C A L m i pi w diff_m;
varexo em ea;

parameters beta xi eta zeta psi_l psi_m nu rho_a rho_m pi_star;
beta  = 0.99;
xi   = 1;
eta = 1;
zeta = 1;
psi_l   = 3.5;
psi_m   = 1;
nu  = 11;
rho_a   = 0.9;
rho_m   = 0.5;
pi_star = 0.01;

model;
diff_m = m-m(-1);
exp(w) = psi_l*(1-exp(L))^(-eta)*exp(C)^xi;
psi_m*exp(m)^(-zeta) = (exp(i)/(1+exp(i)))*exp(C)^(-xi);
exp(C)^xi=((1+exp(pi(+1)))/(beta*(1+exp(i))))*exp(C(+1))^xi;
exp(Y)=exp(A)*exp(L);
exp(C)=exp(Y);
exp(w)=((nu-1)/nu)*exp(A);
diff_m = (1-rho_m)*log(pi_star)-(pi-rho_m*pi(-1))+rho_m*diff_m(-1) + em;
A = rho_a*A(-1) + ea;
end;

// calc steady state
A_star=1;
w_star = ((nu-1)/nu)*A_star;
i_star = (1/beta)*(1+pi_star)-1;
L_star = ((A_star*psi_l)/w_star+1)^(-1);
Y_star = A_star*L_star;
C_star = Y_star;
m_star = ((psi_m*(1+i_star)*C_star)/i_star)^(1/zeta);

initval;
A=log(A_star);
w=log(w_star);
i=log(i_star);
pi=log(pi_star);
m=log(m_star);
C=log(C_star);
L=log(L_star);
Y=log(Y_star);
end;

shocks;
var em; stderr 0.01;
end;

// check steady state value(not log)
[Y_star;C_star;A_star;L_star;m_star;i_star;pi_star;w_star]

steady;

check;

stoch_simul (order =1, irf =50);

