% LIBOR https://www.bankofengland.co.uk/statistics/yield-curves
clc; clear all; clear; rng(1);
data_frequency = 2;
data_type = 'Libor Forward Curve 8 JAN 2020';
format long;
input_data_path = 'LIBOR_SPOT_08012020_DAILY_25ANNI.xlsx';
ValuationDate = '08-Jan-2020';
DateCurveSamples = 0.5:0.5:25;
%% 8 Jan 2020 YIELD CURVE
%%
% Adding a date amount of one day to the starting Settle point
Basis = 1;
CurveDates = daysadd(ValuationDate,360*DateCurveSamples,Basis);
[SpotRates,~,~] = xlsread(input_data_path);
negative_rates = ~isempty(find(SpotRates < 0, 1));
SpotRates = SpotRates/100;
figure();
plot(DateCurveSamples,SpotRates)
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xlabel('Maturity'); ylabel('Interest rates in percentage')
title("Spot Rate Curve LIBOR 8 JAN 2020")

% Forward curve plotted
[ForwardRates, CurveDates] = zero2fwd(SpotRates, CurveDates, ... 
ValuationDate);
figure();
plot(DateCurveSamples,ForwardRates)
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xlabel('Maturity'); ylabel('Interest rates in percentage')
title("Forward Rate Curve LIBOR 8 JAN 2020")

%%
positive_rates = ForwardRates;
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
length_positive_rate = length(positive_rates);
fprintf('log-likelihood = %+3.6f\n', -Fval/length_positive_rate);

%% MLE method 
alpha_optim = params_optim(1);
mi_optim = params_optim(2);
sigma_optim = params_optim(3);

fprintf('Inequality 2*alpha*mi >= sigma^2 is satisfied: %s\n', string(2*alpha_optim*mi_optim >= sigma_optim^2))

%% SIMULATION (DO NOT CHANGE)
T = length(positive_rates)/data_frequency;
x0 = positive_rates(1);
Nsteps = length(positive_rates);
x_axis_positive = linspace(0, 25, length_positive_rate);

%% 1 SIMULATION
nsim = 1; %numbers of simulation
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim);
figure();
cm=colormap(hsv(nsim));
for i=1:nsim
 disp(xExact(i));
 disp(cm(i,:));
plot(x_axis_positive,xExact(:,i),'Color',cm(i,:));
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xlabel('Maturity'); ylabel('Interest rates in percentage')
hold on
end 
hold on;
plot(x_axis_positive,positive_rates, 'b','LineWidth',2);
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xExactLegend_1 = ["Simulated Forward Curve" data_type];
leg = legend(xExactLegend_1);
leg.Location = 'southeast';
hold off; 
xlabel('Year');ylabel('Interest Rates in Percentage');
title(join([data_type; "and 1 simulated paths"], ' '));

%% 100 simulation paths
nsim=100;
[~, xExact] = Exact(alpha_optim,mi_optim,sigma_optim,T,x0,Nsteps,nsim);
figure();
for i=1:nsim
    disp(xExact(i));
    plot(x_axis_positive, xExact(:,i));
    curtick = get(gca, 'YTick');
    set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
    hold on
end
hold on;
plot(x_axis_positive,positive_rates, 'b', 'LineWidth',2);
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
hold off
xlabel('Year');ylabel('Interest Rates in Percentage');
title(join([data_type; "and 100 simulations"], ' '));

%i compute the mean of the paths(100)
mean_exact_100=mean(xExact');
figure();
plot(x_axis_positive,mean_exact_100,'k-',x_axis_positive,positive_rates);
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
xlabel('Year');ylabel('Interest rates in Percentage');
title(join([data_type; "and mean of 100 simulations"], ' '));
leg = legend('Mean of 100 paths', data_type);
leg.Location = 'southeast';

%% CAPBYCIR
Strike_CAP = 0.0086;
Strike_FLOOR = 0.008;
Contract_steps = 10; % Since it is every 6 months = 2 every year, this means 10/2 = 5 years
ValuationDateAsNum = datenum(ValuationDate);

Rates = mean_exact_100(1:Contract_steps)';
Dates = cellstr(datestr([ValuationDateAsNum CurveDates(1:Contract_steps)']));
EndDates = Dates(2:end)';
Compounding = 1;
RateSpec = intenvset('ValuationDate', ValuationDate, 'StartDates', ValuationDate, 'EndDates',EndDates,'Rates', Rates, 'Compounding', Compounding);

%%
NumPeriods = length(EndDates);

Settle = ValuationDate;
Maturity = '08-Jan-2025'; 
CapReset = '2';
Basis = 1;
Principal = 100000;
CIRTimeSpec = cirtimespec(ValuationDate, Maturity, NumPeriods, 'Basis', Basis,...
    'Compounding', Compounding);
CIRVolSpec = cirvolspec(sigma_optim, alpha_optim, mi_optim);
CIRT = cirtree(CIRVolSpec, RateSpec, CIRTimeSpec);

[PriceCAP,~] = capbycir(CIRT,Strike_CAP,Settle,Maturity, 'CapReset', CapReset, ...
    'Basis', Basis, 'Principal', Principal) 

[PriceFLOOR,~] = floorbycir(CIRT,Strike_FLOOR,Settle,Maturity, 'FloorReset', CapReset, ...
    'Basis', Basis, 'Principal', Principal) 
%% SWAPTION PARAMS

ExerciseDatesEU = Dates(end-2);
SwapSettlementEU = ExerciseDatesEU;
SwapMaturityEU = Maturity;
Spread = 0.11; % https://www.learningmarkets.com/understanding-the-libor-spread/
SwapReset = 2;
Strike = 0.007;

%% EU
OptionTypeEU = 0;
%PUT
OptSpec = 'put';  
[PriceSwapEuPut,PriceTreeSwapEuPut] = swaptionbycir(CIRT,OptSpec,Strike,ExerciseDatesEU,Spread,SwapSettlementEU,SwapMaturityEU,'SwapReset',SwapReset, ...
'Basis',Basis,'Principal',Principal, 'AmericanOpt', OptionTypeEU)

%CALL
OptSpec = 'call';  
[PriceSwapEuCall,PriceTreeSwapEuCall] = swaptionbycir(CIRT,OptSpec,Strike,ExerciseDatesEU,Spread,SwapSettlementEU,SwapMaturityEU,'SwapReset',SwapReset, ...
'Basis',Basis,'Principal',Principal, 'AmericanOpt', OptionTypeEU)


%% USA
ExerciseDatesUSA = [CIRT.dObs(end-5) CIRT.dObs(end-1)];
SwapSettlementUSA = ExerciseDatesUSA(1);
SwapMaturityUSA = Maturity;

OptionTypeUSA = 1;
%PUT
OptSpec = 'put';  
[PriceSwapUsaPut,PriceTreeSwapUsaPut] = swaptionbycir(CIRT,OptSpec,Strike,ExerciseDatesUSA,Spread,SwapSettlementUSA,...
    SwapMaturityUSA,'SwapReset',SwapReset, ...
'Basis',Basis,'Principal',Principal, 'AmericanOpt', OptionTypeUSA)

%CALL
OptSpec = 'call';  
[PriceSwapUsaCall,PriceTreeSwapUsaCall] = swaptionbycir(CIRT,OptSpec,Strike,ExerciseDatesUSA,Spread,SwapSettlementUSA,...
    SwapMaturityUSA,'SwapReset',SwapReset, ...
'Basis',Basis,'Principal',Principal, 'AmericanOpt', OptionTypeUSA)