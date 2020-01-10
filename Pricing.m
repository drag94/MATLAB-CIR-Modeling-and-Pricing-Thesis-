%% DOC: https://it.mathworks.com/help/fininst/swaptionbyblk.html - https://it.mathworks.com/help/fininst/floorbyblk.html - 
%% https://it.mathworks.com/help/fininst/capbyblk.html
% LIBOR + OIS: https://www.bankofengland.co.uk/statistics/yield-curves
clc; clear all; clear; rng(1);
%% 8 Jan 2020 YIELD CURVE
format long;
ValuationDate = datenum('08-Jan-2020');
input_data_path = 'LIBOR_SPOT_08012020_DAILY_25ANNI.xlsx';
Compounding = 2;
Settle = datestr(addtodate(ValuationDate, 1, 'year'));
DateCurveSamples = 0.5:0.5:25;
%DateCurveSamples = [3/12 6/12 9/12 1:30];

% CAP-FLOOR PARAMS
Maturity = {datestr(addtodate(datenum(Settle), 1, 'year')), datestr(addtodate(datenum(Settle), 2, 'year'))};
Strike = [0.0085;0.0085];
Volatility = 0.12; % From here: http://www.cboe.com/products/vix-index-volatility/volatility-on-interest-rates for LIBOR - VIX
Shift = 0.7;
Principal = 100000;
Reset = 1;

% SWAPTION PARAMS
ExerciseDate = datestr(addtodate(datenum(Settle), 2, 'year'));
SwVolatility = 0.69; % From here: http://www.cboe.com/products/vix-index-volatility/volatility-on-interest-rates for LIBOR - SRVIX
SwMaturity = datestr(addtodate(datenum(Settle), 3, 'year'));
%%
% Adding a date amount of one day to the starting Settle point
Basis = 1;
CurveDates = daysadd(ValuationDate,360*DateCurveSamples,Basis);
[SpotRates,~,~] = xlsread(input_data_path);
negative_rates = ~isempty(find(SpotRates < 0, 1));
SpotRates = SpotRates/100;
figure();
plot(DateCurveSamples,SpotRates)
title("Spot Rate Curve LIBOR")

% Forward curve plotted
[ForwardRates, CurveDates] = zero2fwd(SpotRates, CurveDates, ... 
ValuationDate);
figure();
plot(DateCurveSamples,ForwardRates)
title("Forward Rate Curve LIBOR")

%%

% Compounding set to -1 because it is continuos (stated in the ECB site)
Curve = intenvset('Rates',ForwardRates,'StartDate',Settle,'EndDates',CurveDates,'Compounding',Compounding,'Basis',Basis, ...
    'ValuationDate',ValuationDate);
% If there is at least one negative rate value, apply Shift Black Model
if negative_rates
    disp("Shift is active on model's application");
    [PriceCAP, Caplets] = capbyblk(Curve, Strike, Settle, Maturity, Volatility, 'Shift', Shift, ...
        'Principal', Principal, 'Reset', Reset, 'ValuationDate', ValuationDate, 'Basis', Basis)
else
    [PriceCAP, Caplets] = capbyblk(Curve, Strike, Settle, Maturity, Volatility, ...
        'Principal', Principal, 'Reset', Reset, 'ValuationDate', ValuationDate, 'Basis', Basis)
end

if negative_rates
    disp("Shift is active on model's application");
    [PriceFLOOR, Floorlets] = floorbyblk(Curve, Strike, Settle, Maturity, Volatility, 'Shift', Shift, ...
        'Principal', Principal, 'Reset', Reset, 'ValuationDate', ValuationDate, 'Basis', Basis)
else
    [PriceFLOOR, Floorlets] = floorbyblk(Curve, Strike, Settle, Maturity, Volatility, ...
        'Principal', Principal, 'Reset', Reset, 'ValuationDate', ValuationDate, 'Basis', Basis)
end
%% SWAPTION BLK

OptSpecCall = 'call';
OptSpecPut = 'put';

if negative_rates
    disp("Shift is active on model's application");
    PriceSwCall = swaptionbyblk(Curve,OptSpecCall,Strike,Settle,ExerciseDate, SwMaturity,SwVolatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal,'Shift',Shift)
else 
    PriceSwCall = swaptionbyblk(Curve,OptSpecCall,Strike,Settle,ExerciseDate, SwMaturity,SwVolatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal)
end

if negative_rates
    disp("Shift is active on model's application");
    PriceSwPut = swaptionbyblk(Curve,OptSpecPut,Strike,Settle,ExerciseDate, SwMaturity,SwVolatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal,'Shift',Shift)
else
    PriceSwPut = swaptionbyblk(Curve,OptSpecPut,Strike,Settle,ExerciseDate, SwMaturity,SwVolatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal)
end