function [alpha_optim, mi_optim, sigma_optim] = CIREstimation(input_data_path, data_frequency, date_format, data_type)

%%%% I estimate the initial parameters
%%% CIR Equation
%% dX_t=(theta1-theta2*Xt)dt+theta3*sqrt(X_t)dW_t
[rates,~,raw] = xlsread(input_data_path);

date_ranges=raw(1:end,1);
rates=rates(1:end)/100;
rates_size = length(rates);
start_end_dates_all_rates = [datenum(date_ranges(1), date_format), datenum(date_ranges(end), date_format)];
figure();
x_axis_2018 = linspace(start_end_dates_all_rates(1), start_end_dates_all_rates(2), rates_size);
plot(x_axis_2018, rates);

%%% Fix limits for y axis for the yield is between -0.6 and 20%
%ylim([-0.1 0.22]);
datetick('x','yyyy');
xlabel('Year'); ylabel('Interest rates in Percentage');
title(data_type, 'fontsize', 10);

%% Euribor 3-month monthly data Jan-1994 Dec-2011
%  Pribor 3-month daily data 3-Jan-1994 4-October-2007

negative_rates = find(rates < 0);
if ~isempty(negative_rates)
    positive_rates = rates(1:negative_rates(1) - 1);
else
    positive_rates = rates;
end
length_positive_rate = length(positive_rates);

end_positive_rates = datenum(date_ranges(length_positive_rate), date_format);
figure();
x_axis_positive = linspace(start_end_dates_all_rates(1), end_positive_rates, length_positive_rate);
plot(x_axis_positive, positive_rates);
datetick('x','yyyy');
title(strcat(data_type, " - only positive rates"), 'fontsize', 10);
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
dt=1/data_frequency;
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
options  =  optimset('fminsearch');
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
%% Here we minimize a function, but the CIRlog is already negative
% This is the same as maximizing the LL
function_to_optimize = @(params_to_optimize) CIRlog(positive_rates, dt, ...
params_to_optimize(1), params_to_optimize(2), params_to_optimize(3));

[params_optim, Fval] = fminsearch(function_to_optimize, InitialParams, options);

fprintf('\nalpha = %+3.6f\nmu = %+3.6f\nsigma = %+3.6f\n', params_optim(1), params_optim(2), params_optim(3));
fprintf('log-likelihood = %+3.6f\n', -Fval/length_positive_rate);

%% MLE method 
alpha_optim = params_optim(1);
mi_optim = params_optim(2);
sigma_optim = params_optim(3);

fprintf('Inequality 2*alpha*mi >= sigma^2 is satisfied: %s\n', string(2*alpha_optim*mi_optim >= sigma_optim^2))

%% PLOT SETTINGS
ngrid = 50;
Nobs=length_positive_rate;
%% PLOTTING 2-D graphs
figure();
visualization_limit = 1.15;
% plot loglike as a function of alpha given optimal mi and sigma 
subplot(3,1,1);
Alpha = linspace(alpha_optim/visualization_limit, visualization_limit*alpha_optim, ngrid);
LnlAlpha=zeros(1,ngrid);
for i=1:ngrid
    LnlAlpha(i) = -CIRlog(positive_rates, dt, Alpha(i), mi_optim , sigma_optim)/Nobs;
end
plot(Alpha,LnlAlpha);
xlabel('\alpha');ylabel('lnL');

 % plot ln L as a function of mi given optimal alpha  and sigma
subplot(3,1,2);
Mi = linspace(mi_optim/visualization_limit, visualization_limit*mi_optim, ngrid);
LnlMi=zeros(1,ngrid);
for i=1:ngrid
    LnlMi(i) = -CIRlog(positive_rates, dt, alpha_optim, Mi(i), sigma_optim)/Nobs;
end
plot(Mi,LnlMi);
xlabel('\mu');ylabel('lnL');

% plot ln L as a function of sigma given optimal mi  and alpha
subplot(3,1,3);
Sigma = linspace(sigma_optim/visualization_limit, visualization_limit*sigma_optim, ngrid);
LnlSigma=zeros(1,ngrid);
for i=1:ngrid
    LnlSigma(i) = -CIRlog(positive_rates, dt, alpha_optim, mi_optim, Sigma(i))/Nobs;
end
plot(Sigma,LnlSigma);
xlabel('\sigma');ylabel('lnL');

%% PLOTTING 3-D surfaces
Lnl=zeros(ngrid,ngrid);
x_plot3D_alpha = Alpha;
x_plot3D_mi = Mi;
x_plot3D_sigma = Sigma;

%% plot ln L as a function of alpha and mi given optimal sigma
figure();
[Alpha_3D,Mi_3D]=meshgrid(x_plot3D_alpha,x_plot3D_mi);
for i=1:ngrid
    for j=1:ngrid
        Lnl(i,j)=-CIRlog(positive_rates, dt, Alpha_3D(i,j), Mi_3D(i,j), sigma_optim)/Nobs;
    end
end
surfc(Alpha_3D,Mi_3D,Lnl)
xlabel('\alpha');ylabel('\mu');zlabel('lnL');
%% plot ln L as a function of sigma and mi given optimal alpha
figure();
Lnl=zeros(ngrid,ngrid);
[Mi_3D, Sigma_3D]=meshgrid(x_plot3D_mi,x_plot3D_sigma);
for i=1:ngrid 
    for j=1:ngrid
        Lnl(i,j)=-CIRlog(positive_rates, dt, alpha_optim, Mi_3D(i,j), Sigma_3D(i,j))/Nobs;
    end
end
surfc(Mi_3D,Sigma_3D,Lnl)
xlabel('\mu');ylabel('\sigma');zlabel('lnL');

%% plot ln L as a function of alpha and sigma given optimal mi
figure();
Lnl=zeros(ngrid,ngrid);
[Alpha_3D, Sigma_3D]=meshgrid(x_plot3D_alpha, x_plot3D_sigma);
for i=1:ngrid
    for j=1:ngrid
        Lnl(i,j)=-CIRlog(positive_rates, dt, Alpha_3D(i,j), mi_optim,Sigma_3D(i,j))/Nobs;
    end
end
surfc(Alpha_3D,Sigma_3D,Lnl)
xlabel('\alpha');ylabel('\sigma');zlabel('lnL');
%% Empirical Bot 12m
T = length(positive_rates)/data_frequency;
x0 = positive_rates(1);

%% 10 simulation paths
Nsteps = length(positive_rates);
nsim = 10; %numbers of simulation
%% Simulation of Exact and RealX_t
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim);
figure();
cm=colormap(hsv(nsim));
for i=1:nsim
 disp(xExact(i));
 disp(cm(i,:));
plot(x_axis_positive,xExact(:,i),'Color',cm(i,:));
hold on
end 
hold on;
plot(x_axis_positive,positive_rates, 'b','LineWidth',2);
xExactLegend_10 = [string([1:nsim]) + 'xExact' data_type];
legend(xExactLegend_10);
hold off; 
xlabel('Year');ylabel('Interest Rates in Percentage');
datetick('x','yyyy');
title(join([data_type; "and 10 simulated paths"], ' '));

%i compute the mean of the paths(10)
figure();
mean_exact_10=mean(xExact');
plot(x_axis_positive,mean_exact_10,'k-',x_axis_positive,positive_rates);
xlabel('Year');ylabel('Interest rates in Percentage');
title(join([data_type; "and mean of 10 simulated paths"], ' '));
datetick('x','yyyy');
legend('Mean of 10 paths', data_type);

%% 100 simulation paths
nsim=100;
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim);
figure();
for i=1:nsim
    disp(xExact(i));
    plot(x_axis_positive, xExact(:,i));
    hold on
end
hold on;
plot(x_axis_positive,positive_rates, 'b', 'LineWidth',2);
hold off
xlabel('Year');ylabel('Interest Rates in Percentage');
datetick('x','yyyy');
title(join([data_type; "and 100 simulated paths"], ' '));

%i compute the mean of the paths(100)
mean_exact_100=mean(xExact');
figure();
plot(x_axis_positive,mean_exact_100,'k-',x_axis_positive,positive_rates);
xlabel('Year');ylabel('Interest rates in Percentage');
title(join([data_type; "and mean of 100 simulated paths"], ' '));
datetick('x','yyyy');
legend('Mean of 100 paths', data_type);

%% 1000 simulation paths
nsim=1000;
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim);
figure();
for i=1:nsim
    disp(xExact(i));
    plot(x_axis_positive,xExact(:,i));
    hold on
end
hold on
plot(x_axis_positive,positive_rates,'b','LineWidth',2);
datetick('x','yyyy');
title(join([data_type; "and 1000 simulated paths"], ' '));

%i compute the mean of the paths(1000)
mean_exact_1000=mean(xExact');
figure();
plot(x_axis_positive,mean_exact_1000,'k-',x_axis_positive,positive_rates);
xlabel('Year');ylabel('Interest rates in Percentage');
title(join([data_type; "and mean of 1000 simulated paths"], ' '));
datetick('x','yyyy');
legend('Mean of 1000 paths', data_type);