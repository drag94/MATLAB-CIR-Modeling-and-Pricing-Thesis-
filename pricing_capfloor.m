%%%
%% Caplet black-sholes

r0 = 0.0551;
k=0.2703;
mu=0.0434;
sigma=0.0272;
PS=0.954;
PT=0.896;
CST=0.9946;
BST=0.8852
quantity=1.06;
K=1
rho=(2*sqrt((k^2)+(2*sigma^2)))/((sigma^2)*(exp(sqrt((k^2)+(2*sigma^2)))-1));
delta=(k+sqrt((k^2)+(2*sigma^2)))/(sigma^2);
r1=(log(CST/K)/BST);
x1=2*r1*(rho+delta+BST);
x2=2*r1*(rho+delta);
v=4*k*mu/(sigma^2);
w1=(2*(rho^2)*r0*exp(sqrt((k^2)+(2*sigma^2))))/(rho+delta+BST);
w2=(2*(rho^2)*r0*exp(sqrt((k^2)+(2*sigma^2))))/(rho+delta);
price=((quantity*PT*ncx2cdf(x1,v,w1)-K*PS*ncx2cdf(x2,v,w2))-quantity*PT+K*PS)

%% CAPLET CIR MONTECARLO
function pathFcir = CirForwardSimulation(r0,k,mu,sigma,strike,T_set,T_ex);
clc
close all
clear all
r0=0.0551;
k=0.2703;
mu=0.0434;
sigma=0.0272;
T_set=2;
T_ex=1;
time_step=250;
delta_t=T_set/time_step;
n_sim=1000000;
c=((sigma^2)*(1-exp(-k*delta_t)))/(4*k) ;
path_R = zeros(n_sim,time_step+1);
path_R(:,1) = r0;
gamma=sqrt((k^2)+2*(sigma^2));
strike=0.06;
Z=randn(n_sim,time_step);

lambda=zeros(n_sim,time_step);
d=(4*mu*k)/(sigma^2);
x=chi2rnd(d-1,[n_sim,time_step]);
path_R = zeros(n_sim,time_step);
path_R(:,1) = r0;
zcb_set= zeros(n_sim,time_step*T_ex);
zcb_ex= zeros(n_sim,time_step*T_ex);
B_Tex=zeros(n_sim,time_step*T_ex);
B_Tset=zeros(n_sim,time_step*T_ex);
C_Tex=zeros(n_sim,time_step*T_ex);
C_Tset=zeros(n_sim,time_step*T_ex);
for j=1:n_sim
    for i=1:time_step-1
        lambda(j,i)=path_R(j,i)*(((4*k*exp(-k*delta_t)))/((1-exp(-k*delta_t))*sigma^2));
        path_R(j,i+1)=(((Z(j,i)+sqrt(lambda(j,i)))^2)+x(j,i))*c;
    end
end
for j=1:n_sim
    for i=1:time_step*T_ex
        n(i)=(T_set-(i-1)/time_step);
        m(i)=(T_ex-(i-1)/time_step);
        B_Tex(j,i)= (2*(exp(gamma*m(i))-1))/((k+gamma)*((exp(gamma*m(i))-1))+2*gamma);
        C_Tex(j,i)=((2*gamma*exp((k+gamma)*m(i)/2))/((k+gamma)*((exp(gamma*m(i))-1))+2*gamma))^((2*k*mu)/(sigma^2));
        B_Tset(j,i)= (2*(exp(gamma*n(i))-1))/((k+gamma)*((exp(gamma*n(i))-1))+2*gamma);
        C_Tset(j,i)=((2*gamma*exp((k+gamma)*n(i)/2))/((k+gamma)*((exp(gamma*n(i))-1))+2*gamma))^((2*k*mu)/(sigma^2));
        zcb_ex(1,1)=C_Tex(1,1)*exp(-path_R(1,1)*B_Tex(1,1));
        zcb_set(1,1)=C_Tset(1,1)*exp(-path_R(1,1)*B_Tset(1,1));
    end
end
for g=1:n_sim
    payoff(g)=max(0,(1-1.05*zcb_set(g,time_step)));
end
mean_payoff=mean(payoff(:))
price_capletCIR=1.06*zcb_ex(1,1)*mean_payoff
end