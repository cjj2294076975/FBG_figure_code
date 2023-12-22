% 已知温度标定曲线以及曲率表达式求解KT
h = 0.6; %单位: mm
close all;clc;clear
FBG1_K_epsilon = 0.87862;
FBG2_K_epsilon = 0.99464;
% K = (delta_lambda/lambda - FBG1_Kt*delta_T)/(K_epsilon*h)
delta_T = [25 30 35 40 45 50 55.1 60];
FBG1_lambda0 = 1547.5760;
FBG2_lambda0 = 1532.6103;
FBG1_lambda = [1553.229 1553.314 1553.413 1553.462 1553.539 1553.618 1553.699 1553.783] - 5.69; %升温FBG1
FBG2_lambda = [1538.287 1538.373 1538.452 1538.520 1538.587 1538.652 1538.727 1538.793] - 5.69; %升温FBG2
FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0;
FBG2_delta_lambda = FBG2_lambda - FBG2_lambda0;
kb1 = polyfit(FBG1_lambda0*delta_T,FBG1_delta_lambda,1);
kb2 = polyfit(FBG2_lambda0*delta_T,FBG2_delta_lambda,1);
FBG1_Kt = kb1(1)/FBG1_lambda0;
FBG2_Kt = kb2(1)/FBG2_lambda0;
fprintf("FBG1_Kt: %s\n",num2str(FBG1_Kt))
fprintf("FBG2_Kt: %s\n",num2str(FBG2_Kt))

% FBG1_Kt: 6.4491e-09
% FBG2_Kt: 6.0534e-09