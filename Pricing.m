%% DOC: https://it.mathworks.com/help/fininst/swaptionbyblk.html - https://it.mathworks.com/help/fininst/floorbyblk.html - 
% https://it.mathworks.com/help/fininst/capbyblk.html
clc; clear all; clear; rng(1);
%% 8 Jan 2020 YIELD CURVE
format long;
ValuationDate = datenum('08-Jan-2020');
input_data_path = 'YieldsSpotAAABond_08012020.xlsx';
Settle = datestr(addtodate(ValuationDate, 1, 'year'));
DateCurveSamples = [3/12 6/12 9/12 1:30];

% CAP-FLOOR PARAMS
Maturity = {datestr(addtodate(datenum(Settle), 1, 'year')), datestr(addtodate(datenum(Settle), 2, 'year'))};
Strike = 0.04;
Volatility = 0.31; % From here: http://www.cboe.com/delayedquote/advanced-charts?ticker=SRVIX for LIBOR
Shift = 0.7;
Principal = 100000;
Reset = 4;
ExerciseDate = datestr(addtodate(datenum(Settle), 2, 'year'));

% SWAPTION PARAMS
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

% Forward curve plotted
[ForwardRates, CurveDates] = zero2fwd(SpotRates, CurveDates, ... 
ValuationDate);
figure();
plot(DateCurveSamples,ForwardRates)

% Compounding set to -1 because it is continuos (stated in the ECB site)
Curve = intenvset('Rates',SpotRates,'StartDate',ValuationDate,'EndDates',CurveDates,'Compounding',-1,'Basis',Basis, ...
    'ValuationDate',ValuationDate);
if negative_rates
    [PriceCAP, Caplets] = capbyblk(Curve, Strike, Settle, Maturity, Volatility, 'Shift', Shift, ...
        'Principal', Principal, 'Reset', Reset, 'ValuationDate', ValuationDate, 'Basis', Basis)
else
    [PriceCAP, Caplets] = capbyblk(Curve, Strike, Settle, Maturity, Volatility, ...
        'Principal', Principal, 'Reset', Reset, 'ValuationDate', ValuationDate, 'Basis', Basis)
end

if negative_rates
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
    PriceSwCall = swaptionbyblk(Curve,OptSpecCall,Strike,Settle,ExerciseDate, SwMaturity,Volatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal,'Shift',Shift)
else 
    PriceSwCall = swaptionbyblk(Curve,OptSpecCall,Strike,Settle,ExerciseDate, SwMaturity,Volatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal)
end

if negative_rates
    PriceSwPut = swaptionbyblk(Curve,OptSpecPut,Strike,Settle,ExerciseDate, SwMaturity,Volatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal,'Shift',Shift)
else
    PriceSwPut = swaptionbyblk(Curve,OptSpecPut,Strike,Settle,ExerciseDate, SwMaturity,Volatility,'Basis',Basis,'Reset',Reset,...
        'Principal',Principal)
end