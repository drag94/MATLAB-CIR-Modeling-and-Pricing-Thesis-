close all; clear all; clc, rng(1), warning off

frequency_daily = 250;
frequency_monthly = 12;

date_format_pribor = 'dd mmm yyyy'; %Pribor 3M
date_format_euribor = 'yyyymmm'; %Euribor 3M
%%
CIREstimation('Pribor3mdaily.xlsx', frequency_daily, date_format_pribor)

%%
CIREstimation('Eur3Mmonthly.xlsx', frequency_monthly, date_format_euribor)


