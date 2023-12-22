clc;clear;close all;

% 原来5.69nm
FBG1_lambda0 = 1553.229 - 3; %nm T0 = 25;
FBG2_lambda0 = 1538.287 - 3; %nm T0 = 25;
T0 = 25;

%% 升温
s_T = [25 30 35 40 45 50 55.1 60]; %升温
s_FBG1_lambda = [1553.229 1553.314 1553.413 1553.462 1553.539 1553.618 1553.699 1553.783] - 3; %升温FBG1
s_FBG2_lambda = [1538.287 1538.373 1538.452 1538.520 1538.587 1538.652 1538.727 1538.793] - 3; %升温FBG2

% s_FBG1_delta_lambda = s_FBG1_lambda - FBG1_lambda0;
% s_FBG2_delta_lambda = s_FBG2_lambda - FBG2_lambda0;

s_kb1 = polyfit(s_T,s_FBG1_lambda,1); %升温FBG1拟合
s_kb2 = polyfit(s_T,s_FBG2_lambda,1); %升温FBG2拟合

corr_matrix = corrcoef(s_T, s_FBG1_delta_lambda);
s_kb1_R2 = corr_matrix(1, 2)^2
corr_matrix = corrcoef(s_T, s_FBG2_delta_lambda);
s_kb2_R2 = corr_matrix(1, 2)^2

%% 降温
% j_T = [60 55 50 45 40 35 30 25]; %降温
% j_FBG1_lambda = [1553.783 1553.694 1553.611 1553.528 1553.446 1553.362 1553.286 1553.229] - 5.69; %降温FBG1
% j_FBG2_lambda = [1538.793 1538.715 1538.621 1538.552 1538.479 1538.391 1538.321 1538.287] - 5.69; %降温FBG2
% 
% j_FBG1_delta_lambda = j_FBG1_lambda - FBG1_lambda0;
% j_FBG2_delta_lambda = j_FBG2_lambda - FBG2_lambda0;
% 
% j_kb1 = polyfit(j_T,j_FBG1_delta_lambda,1); %降温FBG1拟合
% j_kb2 = polyfit(j_T,j_FBG2_delta_lambda,1); %降温FBG2拟合
% corr_matrix = corrcoef(j_T, j_FBG1_delta_lambda);
% j_kb1_R2 = corr_matrix(1, 2)^2
% corr_matrix = corrcoef(j_T, j_FBG2_delta_lambda);
% j_kb2_R2 = corr_matrix(1, 2)^2

figure(1);
FontSize = 10; %设置字体大小
TT = linspace(min(s_T),max(s_T),20);
scatter(s_T,s_FBG1_delta_lambda, 'HandleVisibility', 'off');hold on;
plot(TT,s_kb1(1)*TT+s_kb1(2),'r-o');hold on; %拟合FBG1升温数据
% scatter(j_T,j_FBG1_delta_lambda, 'HandleVisibility', 'off');hold on;
% plot(TT,j_kb1(1)*TT+j_kb1(2),'k-+');hold on; %拟合FBG1降温数据
xlabel("Temperature/℃");ylabel("Wavelength/nm");grid on;
% legend("FBG1 Temperature-rising Δλ=0.015446T-0.37851 R^2=0.9973","FBG1 Temperature-reducing Δλ=0.016064T-0.41936 R^2=0.9979", 'FontName', 'Times New Roman');
% legend("FBG1 Temperature-rising Δλ=0.015446T-0.37851 R^2=0.9973", 'FontName', 'Times New Roman');
legend("FBG1 Temperature-rising λ_B=0.015446T+1547.1605 R^2=0.9973", 'FontName', 'Times New Roman');

% 输出
fprintf("升温过程 FBG1中心波长:%s;温度灵敏度k:%s nm/℃;初始值b:%s;\n",num2str(FBG1_lambda0),num2str(s_kb1(1)),num2str(s_kb1(2)));
fprintf("升温过程 FBG2中心波长:%s;温度灵敏度k:%s nm/℃;初始值b:%s;\n",num2str(FBG2_lambda0),num2str(s_kb2(1)),num2str(s_kb2(2)));
% fprintf("降温过程 FBG1中心波长:%s;温度灵敏度k:%s nm/℃;初始值b:%s;\n",num2str(FBG1_lambda0),num2str(j_kb1(1)),num2str(j_kb1(2)));
% fprintf("降温过程 FBG2中心波长:%s;温度灵敏度k:%s nm/℃;初始值b:%s;\n",num2str(FBG2_lambda0),num2str(j_kb2(1)),num2str(j_kb2(2)));

T = 20;
FBG1_T_lambda0 = vpa(s_kb1(1)*T+s_kb1(2),8)
FBG2_T_lambda0 = vpa(s_kb2(1)*T+s_kb2(2),8)

figure(2);
TT = linspace(min(j_T),max(j_T),20);
scatter(s_T,s_FBG1_delta_lambda, 'HandleVisibility', 'off');hold on;
% scatter(j_T,j_FBG1_delta_lambda, 'HandleVisibility', 'off');hold on;
plot(TT,s_kb2(1)*TT+s_kb2(2),'r-o');hold on; %拟合FBG2升温数据
%plot(TT,j_kb2(1)*TT+j_kb2(2),'k-+');hold on; %拟合FBG2降温数据
xlabel("Temperature/℃");ylabel("Wavelength/nm");grid on;
%legend("FBG2 Temperature-rising Δλ=0.014219T-0.3426 R^2=0.9978","FBG2 Temperature-reducing Δλ=0.01494T-0.4021 R^2=0.9932", 'FontName', 'Times New Roman');
%legend("FBG2 Temperature-rising Δλ=0.014219T-0.3426 R^2=0.9978", 'FontName', 'Times New Roman');
legend("FBG2 Temperature-rising λ_B=0.014219T+1532.2544 R^2=0.9978", 'FontName', 'Times New Roman');

% text(25, 0.55, 'FBG2 Temperature-rising Δλ=0.015446T-0.37851 R^2=0.9978', 'FontName', 'Times New Roman', 'FontSize', FontSize);
% text(25, 0.5, 'FBG2 Temperature-reducing Δλ=0.016064T-0.41936 R^2=0.9932', 'FontName', 'Times New Roman', 'FontSize', FontSize);