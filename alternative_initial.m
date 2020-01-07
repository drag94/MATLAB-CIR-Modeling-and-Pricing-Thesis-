%%% r(0) regredendo i prezzi di una serie di bond a scadeza breve sulla
%%% loro maturity
%%%carico dati di Bot (ZCB italia a breve)%%% 
%%%regresione del loro yield sulla Maturity
%%%% Anno di 50 gg


[beta]=mvregress(YY,XX)
r_0=beta(1)

TT=360

%%%% data odierna
start_date   = '19/12/2019';

%%% dati Bot
s_end1 =  '31/12/2019'
P_1 =99.905
date_find_1 = datenum([start_date; s_end1],'dd/mm/yyyy');
xx_1=[date_find_1(2)-date_find_1(1)]/TT

s_end2 =  '14/01/2020'
P_2 = 100.043
date_find_2 = datenum([start_date; s_end2],'dd/mm/yyyy');
xx_2=[date_find_2(2)-date_find_2(1)]/TT

s_end3 =  '14/02/2020'
P_3 = 100.102
date_find_3 = datenum([start_date; s_end3],'dd/mm/yyyy');
xx_2=linspace(date_find_2(1), date_find_2(2),date_find_2(2)-date_find_2(1));
xx_3=[date_find_3(2)-date_find_3(1)]/TT

s_end4 =  '31/03/2020'
P_4 = 100.072
date_find_4 = datenum([start_date; s_end4],'dd/mm/yyyy');
xx_4=[date_find_4(2)-date_find_4(1)]/TT


s_end5 =  '30/04/2020'
P_5 = 100.125
date_find_5 = datenum([start_date; s_end5],'dd/mm/yyyy');
xx_5=[date_find_5(2)-date_find_5(1)]/TT


s_end6 =  '14/05/2020'
P_6 = 100.066
date_find_6 = datenum([start_date; s_end6],'dd/mm/yyyy');
xx_6=[date_find_6(2)-date_find_6(1)]/TT


s_end7 =  '12/06/2020'
P_7 = 100.073
date_find_7 = datenum([start_date; s_end7],'dd/mm/yyyy');
xx_7=[date_find_7(2)-date_find_7(1)]/TT


%%%% all together
YY=[P_1;P_2; P_3;P_4;P_5;P_6;P_7]
II=[ones(7,1)]
PP=[ xx_1;xx_2;xx_3;xx_4;xx_5;xx_6; xx_7]
XX=[II,PP]

