close all; clear all; clc; pause(0.01), randn('seed',sum(clock)), rand('seed',sum(clock)), warning off
TSTART = tic;
%%%% I estimate the initial parameters
%%% CIR Equation
%% dX_t=(theta1-theta2*Xt)dt+theta3*sqrt(X_t)dW_t
%%% Import Data (dati Bot 12 m...serie mensile 1980-2018 presi da 
% http://www.dt.tesoro.it/it/debito_pubblico/dati_statistici/principali_tassi_di_interesse/) %%%%
[bot_12m,textdata,raw] = xlsread('Data.xls','BOT_12M_1980_2018');
date=datenum(textdata(2:end,1),'mm/dd/yyyy');
rr=bot_12m;
%%% sample length
tt=length(rr);
%%% yield  Plots
%%%%% Figure %%%%%%%
s_start   = '01/01/1980';
s_end    = '01/12/2018';
date_find = datenum([s_start; s_end],'dd/mm/yyyy');
figure();
xx=linspace(date_find(1), date_find(2), tt);
subplot(2,1,1),plot(xx,rr,'linewidth',2);
%%% fix limits for y axis for the yield is between -0.6 and 20%
ylim([-1 25]);
datetick('x','yyyy');
xlabel('Year');ylabel('Interest rates in Percentage');
title(' 12 Months Bot Yield 1980-2018','fontsize',10);

%% from October 2015 the interest rates become negative, therefore i
% remove the observations from November 2015 to the end of the time
% series, and i get 430 observations (from 468 total observations)
% (Cir does not admit negative interest rates)

[idx,logic]=find(rr<0);
X_t=rr(1:idx(1)-1);

s_end_positive=textdata(idx(1)-1,1);
t_pos=length(X_t);
date_find_max_2015 = datenum([s_start; s_end_positive],'dd/mm/yyyy');
xx2=linspace(date_find_max_2015(1), date_find_max_2015(2), t_pos);
plot(xx2,X_t);
datetick('x','yyyy');

%%%%%%%%% MAIN EQUATION %%%%%%%%%%%% 
%% dX_t=theta(beta-X_t)dt+sigma*sqrt(X_t)dW_t
%% theta=alpha
%% beta=mean=mi
%% sigma
%% CIR initial parameters estimation
%% Initial values estim. with OLS having closed form solutions
%% Time series of interest rates observations
%% monthly time step
dt=1/12;
dX_t = diff(X_t);
dX_t = dX_t./X_t(2:end).^0.5;
regressors = [dt./X_t(2:end).^0.5, dt*X_t(2:end).^0.5];
drift = regressors\dX_t; 
%% OLS regressors coefficients estimates
%% residuals to compute variance
res = regressors*drift - dX_t;
alpha0 = -drift(2);
mi0 = -drift(1)/drift(2);
sigma0 = sqrt(var(res, 1)/dt);
%%% Actually parameters should be all positive to have convergence
%%% that is not the case for alpha0
InitialParams = [alpha0 mi0 sigma0]; % Vector of initial parameters
%%% alpha0 is negative so convergence is not ensured.
%%% set it to a small positive value
InitialParams = [0.0001 mi0 sigma0]; % Vector of initial parameters

%%
%%%% Optimization settings
options  =  optimset('fminunc');
options  =  optimset(options , 'Algorithm ','interior-point');
options  =  optimset(options , 'TolFun'      , 1e-006);
options  =  optimset(options , 'TolX'        , 1e-006);
options  =  optimset(options , 'TolCon'      , 1e-006);
options  =  optimset(options , 'Display'     , 'iter');
options  =  optimset(options , 'Diagnostics' , 'on');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxIter'     , 1500);
options  =  optimset(options , 'Jacobian'     ,'off');
options  =  optimset(options , 'MeritFunction'     ,'multiobj');
options  =  optimset(options , 'MaxFunEvals' , 5000);
%%
dt
params_optim = fminsearch(@(params)-CIRlog(X_t, dt, ...
params(1), params(2), params(3)),InitialParams, options );



%% MLE method 
alpha = params_optim(1);
mi = params_optim(2);
sigma = params_optim(3);
alpha0=InitialParams(1)
mi0=InitialParams(2)
sigma0=InitialParams(3)

% estimated initial paramter values
[alpha0,mi0,sigma0; alpha ,mi ,sigma ];

%% and optimal parameters
[CIRlog(X_t,dt,alpha0,mi0,sigma0); ...
CIRlog(X_t,dt,alpha,mi,sigma)];

%% PLOTTING 2-D graphs
figure(1);
ngrid1= 100;
 % plot loglike as a function of alpha given optimal mi and sigma 
 subplot(3,1,1);
 Alpha=(alpha/2):alpha/(ngrid1-1):3*alpha/2;
 LnlAlpha=zeros(1,ngrid1);
 for i=1:ngrid1
 LnlAlpha(i)=CIRlog(X_t, dt, Alpha(i), mi ,sigma, options);
 end;
plot(Alpha,LnlAlpha);
xlabel('\alpha');ylabel('lnL');

 % plot ln L as a function of mi given optimal alpha  and sigma
subplot(3,1,2);
Mi=mi/2:mi/(ngrid1-1):3*mi/2;
LnlMi=zeros(1,ngrid1);
for i=1:ngrid1
    LnlMi(i)=CIRlog(X_t, dt, alpha, Mi(i),sigma, options);
end
plot(Mi,LnlMi);
xlabel('\mu');ylabel('lnL');


% plot ln L as a function of sigma given optimal mi  and alpha
 subplot(3,1,3);
  Sigma=sigma /2:sigma/(ngrid1-1):3*sigma/2;
  LnlSigma=zeros(1,ngrid1);
for i=1:ngrid1
 LnlSigma(i)=CIRlog(X_t,dt,alpha ,mi,Sigma(i), options);
end;
 plot(Sigma,LnlSigma);
xlabel('\sigma');ylabel('lnL');

%% PLOTTING 3-D surfaces
figure(2)
ngrid2=100;
% 
subplot(1,1,1)
Lnl=zeros(ngrid2,ngrid2);
x_plot3D_alpha=[0.03:(0.11-0.03)/(ngrid2-1):0.11];
x_plot3D_mi=[1:(4-1)/(ngrid2-1):4];
x_plot3D_sigma=[0.3:(1.2-0.3)/(ngrid2-1):1.2];
%% plot ln L as a function of alpha and mi given optimal sigma
[Alpha,Mi]=meshgrid(x_plot3D_alpha,x_plot3D_mi);
for i=1:ngrid2
    for j=1:ngrid2
        Lnl(i,j)=CIRlog(X_t, dt, Alpha(i,j), Mi(i,j), sigma);
    end
end
surfc(Alpha,Mi,Lnl)
xlabel('\alpha');ylabel('\mu');zlabel('lnL');
%% plot ln L as a function of sigma and mi given optimal alpha
subplot(1,1,1)
Lnl=zeros(ngrid2,ngrid2);
[Mi, Sigma]=meshgrid(x_plot3D_mi,x_plot3D_sigma  );
for i=1:ngrid2 
    for j=1:ngrid2
        Lnl(i,j)=CIRlog(X_t, dt, alpha  ,Mi(i,j), Sigma(i,j));
    end
end
surfc(Mi,Sigma,Lnl)
xlabel('\mu');ylabel('\sigma');zlabel('lnL');
%% plot ln L as a function of alpha and sigma given optimal mi
subplot(1,1,1)
Lnl=zeros(ngrid2,ngrid2);
[Alpha, Sigma]=meshgrid(x_plot3D_alpha,  x_plot3D_sigma);
for i=1:ngrid2
    for j=1:ngrid2
        Lnl(i,j)=CIRlog(X_t, dt, Alpha(i,j), mi  , Sigma(i,j));
    end
end
surfc(Alpha,Sigma,Lnl)
xlabel('\alpha');ylabel('\sigma');zlabel('lnL');
%%

% Empirical Bot 12m
 figure(1);
 T = length(X_t)/12;
 x0 = X_t(1);
subplot(1,1,1),plot(xx2,X_t,'linewidth',2);
datetick('x','yyyy');
xlabel('Year');ylabel('Interest rates in Percentage');
%% Simulation of Exact and RealX_t
%% 10 simulation paths
 Nsteps = length(X_t);
 nsim_10 =10; %numbers of simulation
 [~, xExact] = Exact(alpha,mi,sigma,T,x0,Nsteps,nsim_10);
 cm=colormap(hsv(nsim_10));
 for i=1:nsim_10
     disp(xExact(i));
     disp(cm(i,:));
    plot(xx2,xExact(:,i),'Color',cm(i,:));
    hold on
 end 
 hold on;
 plot(xx2,X_t, 'b','LineWidth',2);
 xExactLegend_10 = [string([1:nsim_10]) + 'xExact' '12 Months Bot Yield 1980-2015'];
legend(xExactLegend_10);
hold off; 
 xlabel('Year');ylabel('Interest Rates in Percentage');
 datetick('x','yyyy');
 title('12 Months Bot Yield and 10 simulated paths 1980-2015');
%i compute the mean of the paths(10)
mean_exact_10=mean(xExact');
figure();
 plot(xx2,mean_exact_10,'k-',xx2,X_t);
xlabel('Year');ylabel('Interest rates in Percentage');
title('12 Months Bot Yield and mean of 10 simulated paths 1980-2015','fontsize',10);
datetick('x','yyyy');
legend('Mean of 10 paths','12 Months Bot Yield');

%% 100 simulation paths
nsim_100=100;
[~, xExact] = Exact(alpha,mi,sigma,T,x0,Nsteps,nsim_100);
for i=1:nsim_100
    disp(xExact(i));
    plot(xx2, xExact(:,i));
    hold on
end
hold on;
plot(xx2,X_t, 'b', 'LineWidth',2);
hold off
xlabel('Year');ylabel('Interest Rates in Percentage');
datetick('x','yyyy');
title('12 Months Bot Yield and 100 simulated paths 1980-2015');
%i compute the mean of the paths(100)
mean_exact_100=nanmean(xExact');
figure();
plot(xx2,mean_exact_100,'k-',xx2,X_t);
xlabel('Year');ylabel('Interest rates in Percentage');
title('12 Months Bot Yield and mean of 100 simulated paths 1980-2015','fontsize',10);
datetick('x','yyyy');
legend('Mean of 100 paths','12 Months Bot Yield');
%% 1000 simulation paths
nsim_1000=1000;
[~, xExact] = Exact(alpha,mi,sigma,T,x0,Nsteps,nsim_1000);
for i=1:nsim_1000
    disp(xExact(i));
    plot(xx2,xExact(:,i));
    hold on
end
hold on
plot(xx2,X_t,'b','LineWidth',2);
datetick('x','yyyy');
title('12 Months Bot Yield and 1000 simulated paths 1980-2015');
%i compute the mean of the paths(1000)
mean_exact_1000=nanmean(xExact');
figure();
plot(xx2,mean_exact_1000,'k-',xx2,X_t);
xlabel('Year');ylabel('Interest rates in Percentage');
title('12 Months Bot Yield and mean of 1000 simulated paths 1980-2015','fontsize',10);
datetick('x','yyyy');
legend('Mean of 1000 paths','12 Months Bot Yield');

    
    