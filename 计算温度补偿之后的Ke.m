h = 0.6; %mm
L1 = 20;% FBG2到FBG1的距离
L2 = 20; % 自由端到FBG2的距离
L = 60; %固定端到自由端距离
KT1 = 15.446; %15.446 pm/℃
KT2 = 14.219; %14.219 pm/℃;
delta_T = 10; %10
s = [0     1     2     3     4     5     6     7     8     9    10    11    12];

FBG1_lambda0 = 1553.266e3 - 5.69e3;
FBG1_lambda = [1553.266 1553.4135 1553.5497 1553.7267 1553.8665 1553.9978 1554.145 1554.28575 1554.4196 1554.5994 1554.7686 1554.8854 1555.000]*1e3 - 5.69e3;
FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0;

FBG2_lambda0 = 1538.3003e3 - 5.69e3;
FBG2_lambda = [1538.3003 1538.2257 1538.13475 1538.0365 1537.950 1537.86025 1537.7983 1537.7104 1537.6185 1537.549 1537.473125 1537.37367 1537.3164]*1e3 - 5.69e3;
FBG2_delta_lambda = FBG2_lambda0 - FBG2_lambda;

% FBG1_delta_lambda = FBG1_delta_lambda - KT1*delta_T;
% FBG2_delta_lambda = FBG2_delta_lambda - KT2*delta_T;

FBG1_lambda = FBG1_lambda - KT1*delta_T;
FBG2_lambda = FBG2_lambda - KT2*delta_T;

kb1 = polyfit(s,FBG1_lambda,1);
kb2 = polyfit(s,FBG2_lambda,1);

ke1  = 2*L^3*kb1(1)/(3*(L1+L2)*h)/1e6
ke2  = 2*L^3*kb2(1)/(3*L1*h)/1e6

% ke1 = 2*L^3*(FBG1_delta_lambda(end) - KT1*delta_T)/(3*(L1+L2)*s(end)*h)/1e6
% ke2 = 2*L^3*(FBG2_delta_lambda(end) - KT2*delta_T)/(3*L1*s(end)*h)/1e6