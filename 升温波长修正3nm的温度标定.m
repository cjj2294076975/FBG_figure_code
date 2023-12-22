clc;clear;close all;
% 升温波长修正3nm的温度标定
% 原来5.69nm

delta_lambda1 = 3;
delta_lambda2 = 3;

FBG1_lambda0 = 1553.229 - delta_lambda1; %nm T0 = 25;
FBG2_lambda0 = 1538.287 - delta_lambda2; %nm T0 = 25;
T0 = 25;

%% 升温
s_T = [25 30 35 40 45 50 55.1 60]; %升温
TT = linspace(min(s_T),max(s_T),20);
s_FBG1_lambda = [1553.229 1553.314 1553.413 1553.462 1553.539 1553.618 1553.699 1553.783] - delta_lambda1; %升温FBG1
s_FBG2_lambda = [1538.287 1538.373 1538.452 1538.520 1538.587 1538.652 1538.727 1538.793] - delta_lambda2; %升温FBG2

s_kb1 = polyfit(s_T,s_FBG1_lambda,1); %升温FBG1拟合
s_kb2 = polyfit(s_T,s_FBG2_lambda,1); %升温FBG2拟合

figure(1);
FontSize = 10; %设置字体大小
% scatter(s_T,s_FBG1_lambda, 'HandleVisibility', 'off');hold on;
% plot(TT,s_kb1(1)*TT+s_kb1(2),'r');hold on; %拟合FBG1升温数据
plot(TT,s_kb1(1)*TT+s_kb1(2),'r-o');hold on; %拟合FBG1升温数据
xlabel("Temperature/℃");ylabel("Wavelength/nm");grid on;
legend("FBG1 Temperature-rising \it{λ}\rm{_B=0.015446}\it{T}\rm{+1549.8505} \it{R}^2\rm{=0.9973}", 'FontName', 'Times New Roman');

figure(2);
% scatter(s_T,s_FBG2_lambda, 'HandleVisibility', 'off');hold on;
% plot(TT,s_kb2(1)*TT+s_kb2(2),'r');hold on; %拟合FBG2升温数据
plot(TT,s_kb2(1)*TT+s_kb2(2),'r-o');hold on; %拟合FBG2升温数据
xlabel("Temperature/℃");ylabel("Wavelength/nm");grid on;
legend("FBG2 Temperature-rising \it{λ}\rm{_B=0.014219}\it{T}\rm{+1534.9444} \it{R}\rm{^2=0.9978}", 'FontName', 'Times New Roman');

% 输出
fprintf("升温过程 FBG1中心波长:%s;温度灵敏度k:%s nm/℃;初始值b:%s;\n",num2str(FBG1_lambda0),num2str(s_kb1(1)),num2str(s_kb1(2)));
fprintf("升温过程 FBG2中心波长:%s;温度灵敏度k:%s nm/℃;初始值b:%s;\n",num2str(FBG2_lambda0),num2str(s_kb2(1)),num2str(s_kb2(2)));

T = 25;
FBG1_T_lambda0 = vpa(s_kb1(1)*T+s_kb1(2),8)
FBG2_T_lambda0 = vpa(s_kb2(1)*T+s_kb2(2),8)

% corr_matrix = corrcoef(s_T, s_FBG1_lambda);
% s_kb1_R2 = corr_matrix(1, 2)^2
% corr_matrix = corrcoef(s_T, s_FBG2_lambda);
% s_kb2_R2 = corr_matrix(1, 2)^2
