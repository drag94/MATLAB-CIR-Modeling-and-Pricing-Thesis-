close all; clear all; clc; pause(0.01), randn('seed',sum(clock)), rand('seed',sum(clock)), warning off
TSTART = tic;
%%%% I estimate the initial parameters
%%% CIR Equation
%% dX_t=(theta1-theta2*Xt)dt+theta3*sqrt(X_t)dW_t
%%% Import Data (dati Bot 12 m...serie mensile 1980-2018 presi da 
% http://www.dt.tesoro.it/it/debito_pubblico/dati_statistici/principali_tassi_di_interesse/) %%%%
% Data was transformed later. Use Data.xls
[rates,~,raw] = xlsread('Data.xls','BOT_12M_1980_2018');

date_ranges=datenum(raw(2:end,1));

rates_size = length(rates);
start_end_dates_all_rates = [date_ranges(1), date_ranges(end)];
figure();
x_axis_2018 = linspace(start_end_dates_all_rates(1), start_end_dates_all_rates(2), rates_size);

subplot(2,1,1), plot(x_axis_2018, rates, 'linewidth', 2);

%%% Fix limits for y axis for the yield is between -0.6 and 20%
%ylim([-0.1 0.22]);
datetick('x','yyyy');
xlabel('Year'); ylabel('Interest rates in Percentage');
title('12 Months Bot Yield 1980-2018', 'fontsize', 10);

%% from October 2015 the interest rates become negative, therefore i
% remove the observations from November 2015 to the end of the time
% series, and i get 430 observations (from 468 total observations)
% (Cir does not admit negative interest rates)

negative_rates = find(rates < 0);
positive_rates = rates(1:negative_rates(1) - 1);
length_positive_rate = length(positive_rates);
% +1 because we have data with header. So we have to shift with one
% position
end_positive_rates = datenum(raw(length_positive_rate + 1, 1));

x_axis_positive = linspace(start_end_dates_all_rates(1), end_positive_rates, length_positive_rate);
plot(x_axis_positive, positive_rates);
datetick('x','yyyy');
%ylim([-0.05 0.22]);

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
positive_rates_minus_one = positive_rates(1:end-1);
dX_t = diff(positive_rates);
dX_t = dX_t./positive_rates_minus_one.^0.5;
regressors = [dt./positive_rates_minus_one.^0.5, dt*positive_rates_minus_one.^0.5];
drift = regressors\dX_t;
%% OLS regressors coefficients estimates
%% residuals to compute variance
res = regressors*drift - dX_t;
alpha0 = -drift(2);
mi0 = -drift(1)/drift(2);
sigma0 = sqrt(var(res, 1)/dt);

InitialParams = [alpha0 mi0 sigma0] % Vector of initial parameters
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
function_to_optimize = @(params_to_optimize) CIRlog(positive_rates, dt, ...
params_to_optimize(1), params_to_optimize(2), params_to_optimize(3));

[params_optim, Fval] = fminsearch(function_to_optimize, InitialParams, options);

fprintf('\n alpha = %+3.6f\n mu = %+3.6f\n sigma = %+3.6f\n', params_optim(1), params_optim(2), params_optim(3));
fprintf('log-likelihood = %+3.6f\n', Fval/length_positive_rate);

%% MLE method 
alpha_optim = params_optim(1);
mi_optim = params_optim(2);
sigma_optim = params_optim(3);

%% PLOTTING 2-D graphs
figure(1);
ngrid1 = 100;
% plot loglike as a function of alpha given optimal mi and sigma 
subplot(3,1,1);
Alpha = linspace(alpha_optim/2, 2*alpha_optim, ngrid1);
LnlAlpha=zeros(1,ngrid1);
for i=1:ngrid1
    LnlAlpha(i) = CIRlog(positive_rates, dt, Alpha(i), mi_optim , sigma_optim);
end
plot(Alpha,LnlAlpha);
xlabel('\alpha');ylabel('lnL');

 % plot ln L as a function of mi given optimal alpha  and sigma
subplot(3,1,2);
Mi = linspace(mi_optim/2, 2*mi_optim, ngrid1);
LnlMi=zeros(1,ngrid1);
for i=1:ngrid1
    LnlMi(i) = CIRlog(positive_rates, dt, alpha_optim, Mi(i), sigma_optim);
end
plot(Mi,LnlMi);
xlabel('\mu');ylabel('lnL');

% plot ln L as a function of sigma given optimal mi  and alpha
subplot(3,1,3);
Sigma = linspace(sigma_optim/2, 2*sigma_optim, ngrid1);
LnlSigma=zeros(1,ngrid1);
for i=1:ngrid1
    LnlSigma(i) = CIRlog(positive_rates, dt, alpha_optim, mi_optim, Sigma(i));
end
plot(Sigma,LnlSigma);
xlabel('\sigma');ylabel('lnL');

%% PLOTTING 3-D surfaces
ngrid2=100;
Lnl=zeros(ngrid2,ngrid2);
x_plot3D_alpha = Alpha;
x_plot3D_mi = Mi;
x_plot3D_sigma = Sigma;
%% plot ln L as a function of alpha and mi given optimal sigma
[Alpha_3D,Mi_3D]=meshgrid(x_plot3D_alpha,x_plot3D_mi);
for i=1:ngrid2
    for j=1:ngrid2
        Lnl(i,j)=CIRlog(positive_rates, dt, Alpha_3D(i,j), Mi_3D(i,j), sigma_optim);
    end
end
surfc(Alpha_3D,Mi_3D,Lnl)
xlabel('\alpha');ylabel('\mu');zlabel('lnL');
%% plot ln L as a function of sigma and mi given optimal alpha
subplot(1,1,1)
Lnl=zeros(ngrid2,ngrid2);
[Mi_3D, Sigma_3D]=meshgrid(x_plot3D_mi,x_plot3D_sigma);
for i=1:ngrid2 
    for j=1:ngrid2
        Lnl(i,j)=CIRlog(positive_rates, dt, alpha_optim, Mi_3D(i,j), Sigma_3D(i,j));
    end
end
surfc(Mi_3D,Sigma_3D,Lnl)
xlabel('\mu');ylabel('\sigma');zlabel('lnL');

%% plot ln L as a function of alpha and sigma given optimal mi
subplot(1,1,1)
Lnl=zeros(ngrid2,ngrid2);
[Alpha_3D, Sigma_3D]=meshgrid(x_plot3D_alpha, x_plot3D_sigma);
for i=1:ngrid2
    for j=1:ngrid2
        Lnl(i,j)=CIRlog(positive_rates, dt, Alpha_3D(i,j), mi_optim,Sigma_3D(i,j));
    end
end
surfc(Alpha_3D,Sigma_3D,Lnl)
xlabel('\alpha');ylabel('\sigma');zlabel('lnL');
%%

% Empirical Bot 12m
 figure(1);
 T = length(positive_rates)/12;
 x0 = positive_rates(1);
subplot(1,1,1),plot(x_axis_positive,positive_rates,'linewidth',2);
datetick('x','yyyy');
xlabel('Year');ylabel('Interest rates in Percentage');

%% 10 simulation paths
 Nsteps = length(positive_rates);
 nsim_10 =10; %numbers of simulation
 %% Simulation of Exact and RealX_t
 [~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim_10);
 cm=colormap(hsv(nsim_10));
 for i=1:nsim_10
     disp(xExact(i));
     disp(cm(i,:));
    plot(x_axis_positive,xExact(:,i),'Color',cm(i,:));
    hold on
 end 
 hold on;
 plot(x_axis_positive,positive_rates, 'b','LineWidth',2);
 xExactLegend_10 = [string([1:nsim_10]) + 'xExact' '12 Months Bot Yield 1980-2015'];
legend(xExactLegend_10);
hold off; 
 xlabel('Year');ylabel('Interest Rates in Percentage');
 datetick('x','yyyy');
 title('12 Months Bot Yield and 10 simulated paths 1980-2015');
%i compute the mean of the paths(10)
mean_exact_10=mean(xExact');
figure();
 plot(x_axis_positive,mean_exact_10,'k-',x_axis_positive,positive_rates);
xlabel('Year');ylabel('Interest rates in Percentage');
title('12 Months Bot Yield and mean of 10 simulated paths 1980-2015','fontsize',10);
datetick('x','yyyy');
legend('Mean of 10 paths','12 Months Bot Yield');

%% 100 simulation paths
nsim_100=100;
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim_100);
for i=1:nsim_100
    disp(xExact(i));
    plot(x_axis_positive, xExact(:,i));
    hold on
end
hold on;
plot(x_axis_positive,positive_rates, 'b', 'LineWidth',2);
hold off
xlabel('Year');ylabel('Interest Rates in Percentage');
datetick('x','yyyy');
title('12 Months Bot Yield and 100 simulated paths 1980-2015');
%i compute the mean of the paths(100)
mean_exact_100=nanmean(xExact');
figure();
plot(x_axis_positive,mean_exact_100,'k-',x_axis_positive,positive_rates);
xlabel('Year');ylabel('Interest rates in Percentage');
title('12 Months Bot Yield and mean of 100 simulated paths 1980-2015','fontsize',10);
datetick('x','yyyy');
legend('Mean of 100 paths','12 Months Bot Yield');
%% 1000 simulation paths
nsim_1000=1000;
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim_1000);
for i=1:nsim_1000
    disp(xExact(i));
    plot(x_axis_positive,xExact(:,i));
    hold on
end
hold on
plot(x_axis_positive,positive_rates,'b','LineWidth',2);
datetick('x','yyyy');
title('12 Months Bot Yield and 1000 simulated paths 1980-2015');
%i compute the mean of the paths(1000)
mean_exact_1000=nanmean(xExact');
figure();
plot(x_axis_positive,mean_exact_1000,'k-',x_axis_positive,positive_rates);
xlabel('Year');ylabel('Interest rates in Percentage');
title('12 Months Bot Yield and mean of 1000 simulated paths 1980-2015','fontsize',10);
datetick('x','yyyy');
legend('Mean of 1000 paths','12 Months Bot Yield');

    
    