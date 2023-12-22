clc;clear;close all;
% 8月28日实验
% 双点实验规划
% 曲率线性插值(固定端-FBG传感器)
% 等曲率(FBG传感器-自由端)
% 固定端-20-FBG2-20-FBG1-20-自由端
h = 0.6; %mm
L1 = 20;% FBG2到FBG1的距离
L2 = 20; % 自由端到FBG2的距离
L = 60; %固定端到自由端距离
s0 = 36.60; %初始位移
KT1 = 15.446; %15.446 pm/℃
KT2 = 14.219; %14.219 pm/℃;

disp(['自由端距离FBG1',string(20),'FBG1距离FBG2',string(20),'FBG2距离固定端',string(L1+L2)]);

%% FBG1位移灵敏度系数以及应变灵敏度系数计算
s = s0 - [36.60 35.6 34.6 33.6 32.6 31.6 30.6 29.6 28.6 27.6 26.6 25.6 24.6];
FBG1_lambda0 = 1553.266e3 - 5.69e3;
FBG1_lambda = [1553.266 1553.4135 1553.5497 1553.7267 1553.8665 1553.9978 1554.145 1554.28575 1554.4196 1554.5994 1554.7686 1554.8854 1555.000]*1e3 - 5.69e3;
% FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0 - KT1*delta_T;
FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0;
% FBG1_lambda(8) = 1554.294e3 - 5.69e3; %1554.28575
% FBG1_lambda(2) = 1553.44e3 - 5.69e3; %1554.28575

FBG1_kb_s = polyfit(s,FBG1_delta_lambda,1); %k1: FBG1的位移灵敏度系数
FBG1_fit_s = linspace(min(s),max(s),13); %拟合的位移s
FBG1_s_fit_delta_lambda = FBG1_kb_s(1)*FBG1_fit_s+FBG1_kb_s(2);

disp(['FBG1中心波长',string(FBG1_lambda0/1000),'位移灵敏度系数',string(FBG1_kb_s(1))]);

FBG1_epsilon = 3*(L1+L2)*h*s/(2*L^3)*10^6; %FBG1理论应变
FBG1_kb_epsilon = polyfit(FBG1_epsilon,FBG1_delta_lambda,1);
FBG1_fit_epsilon = linspace(min(FBG1_epsilon),max(FBG1_epsilon),13);
FBG1_e_fit_delta_lambda = FBG1_kb_epsilon(1)*FBG1_fit_epsilon+FBG1_kb_epsilon(2);
FBG1_epsilon_err = max(abs(FBG1_delta_lambda - FBG1_e_fit_delta_lambda))/max(FBG1_delta_lambda); %FBG1应变线性度

disp(['FBG1中心波长',string(FBG1_lambda0/1000),'应变灵敏度系数',string(FBG1_kb_epsilon(1))]);
disp(['FBG1应变关于波长变化量的线性度',string(FBG1_epsilon_err)]);

%% FBG2位移灵敏度系数以及应变灵敏度系数计算
FBG2_lambda0 = 1538.3003e3 - 5.69e3;
%FBG2_lambda = [1538.3003 1538.2257 1538.13475 1538.0365 1537.950 1537.86025 1537.7983 1537.7104 1537.6185 1537.549 1537.473125 1537.37367 1537.3164]*1e3 - 5.69e3+0.0146*2e3;
FBG2_lambda = [1538.3003 1538.2257 1538.13475 1538.0365 1537.950 1537.86025 1537.7983 1537.7104 1537.6185 1537.549 1537.473125 1537.37367 1537.3164]*1e3 - 5.69e3;
FBG2_delta_lambda = FBG2_lambda - FBG2_lambda0;
% FBG2_delta_lambda = FBG2_lambda - FBG2_lambda0 - KT2*delta_T;

FBG2_kb_s = polyfit(s,FBG2_delta_lambda,1); %k1: FBG2的位移灵敏度系数
FBG2_fit_s = linspace(min(s),max(s),13); %拟合的位移s
FBG2_s_fit_delta_lambda = FBG2_kb_s(1)*FBG2_fit_s + FBG2_kb_s(2);
FBG2_s_err = 1 - max(abs(FBG2_delta_lambda - FBG2_s_fit_delta_lambda))/max(FBG2_delta_lambda); %FBG2位移线性度

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'位移灵敏度系数',string(FBG2_kb_s(1))]);

FBG2_epsilon = 3*L2*h*s/(2*L^3)*10^6; %FBG2理论应变
FBG2_kb_epsilon = polyfit(FBG2_epsilon,FBG2_delta_lambda,1);
FBG2_fit_epsilon = linspace(min(FBG2_epsilon),max(FBG2_epsilon),13);
FBG2_e_fit_delta_lambda = FBG2_kb_epsilon(1)*FBG2_fit_epsilon+FBG2_kb_epsilon(2);
FBG2_epsilon_err = max(abs(FBG2_delta_lambda - FBG2_e_fit_delta_lambda))/max(FBG2_delta_lambda); %FBG2应变线性度

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'应变灵敏度系数',string(FBG2_kb_epsilon(1))]);
disp(['FBG2应变关于波长变化量的线性度',string(FBG2_epsilon_err)]);

%%绘图
figure(1); %FBG1 位移 - 波长变化量 - 应变
FontSize = 14;
scatter3(s,FBG1_epsilon,FBG1_delta_lambda);hold on;
plot3(FBG1_fit_s,FBG1_fit_epsilon,FBG1_s_fit_delta_lambda);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_1/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('Δλ_{B1}/pm', 'FontName', 'Times New Roman', 'FontSize', FontSize);grid on;
%text(0, 100, {'FBG1 K_s = 146.436pm/mm','FBG1 K_ε = 0.879pm/με'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);

figure(2); %FBG2 位移 - 波长变化量 - 应变
scatter3(s,FBG2_epsilon,FBG2_delta_lambda);hold on;
plot3(FBG2_fit_s,FBG2_fit_epsilon,FBG2_s_fit_delta_lambda);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_2/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('Δλ_{B2}/pm', 'FontName', 'Times New Roman', 'FontSize', FontSize);grid on;
%text(0, 100, {'FBG2 K_s = 82.887pm/mm','FBG2 K_ε = 0.995pm/με'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
%title('FBG应变-波长变化量');
% legend("FBG1曲线","FBG2曲线");

%无温度补偿

% K1 = 2*FBG1_delta_lambda(point_num)/(h*FBG1_kb_epsilon(1))*1e-6;
% K2 = 2*FBG2_delta_lambda(point_num)/(h*FBG2_kb_epsilon(1))*1e-6;

%有温度补偿
% FBG1_kb_epsilon(1) = 0.7125; %30
% FBG2_kb_epsilon(1) = -0.6995;

% FBG1_kb_epsilon(1) = 0.8670 %0
% FBG2_kb_epsilon(1) = -0.9839

% ks1 = 3e6*h*(L1+L2)*FBG1_kb_epsilon(1)/(2*L^3)
% ks2 = 3e6*h*L1*FBG1_kb_epsilon(1)/(2*L^3)

%% 末端位移
point_num = 5; % 1-13: 2 5 8 12
% delta_T = -2; %ΔT∈[-10-10] 7 11mm -2℃
% delta_T = -1.5; %ΔT∈[-10-10] 1mm    -1.5℃
delta_T = -0.5; %ΔT∈[-10-10] 4mm    ℃
b = 2e-3; %纤芯半径2um

% K1 = 2*FBG1_delta_lambda(point_num)/(h*FBG1_kb_epsilon(1))*1e-6;
% K2 = 2*FBG2_delta_lambda(point_num)/(h*FBG2_kb_epsilon(1))*1e-6;

% K1 = 2*(FBG1_delta_lambda(point_num))/(FBG1_kb_epsilon(1)*(h+2*b))*1e-6;
% K2 = 2*(FBG2_delta_lambda(point_num))/(FBG2_kb_epsilon(1)*(h+2*b))*1e-6;

K1 = 2*(FBG1_delta_lambda(point_num) - 2*KT1*delta_T)/(FBG1_kb_epsilon(1)*h)*1e-6;
K2 = 2*(FBG2_delta_lambda(point_num) - 2*KT2*delta_T)/(FBG2_kb_epsilon(1)*h)*1e-6;

% 曲率线性插值
K = [K1 K1 K2 0]
sss = [0 20 40 60];
fix_point_number = 100;
%[ss,kk] = linear_interf(sss,K,fix_point_number);

% 曲率三次样条插值
[ss,kk] = spline_interf(sss,K,fix_point_number);

% 绘制L关于曲率的曲线
figure(3);
plot(ss,kk);xlabel('s/mm');ylabel('K/mm^{-1}');title('位移关于曲率的关系');
R = 1./kk;
ds = L/fix_point_number; %连续点的间距
theta = kk*ds; %theta = K*ds
% 局部坐标
x = R - R.*cos(theta);
y = zeros(1,fix_point_number);
z = R.*sin(theta);
% 全局坐标
TT = eye(4);
AA = [0 0 0 1]';
for i =1:100
    if kk(i)==0
        x(i) = 0;z(i)=ds;
    end
end
for i =1:100
    if i == 1
        A = [0 0 ds 1]';
    elseif i>1
        dabc = [x(i) y(i) z(i)]';
        rot_T = roty(theta(i)/pi*180);
        T = [rot_T dabc;0 0 0 1];
        TT = TT*T;
        A = TT*[dabc;1];
    end
    AA = [AA A];
end
AA = [roty(90) [0 0 0]';0 0 0 1]*AA;
X = AA(1,:);Y = AA(2,:);Z = AA(3,:);
figure(3);
plot3(X,Y,Z);
xlabel("x/mm");ylabel('y/mm');zlabel('z/mm');title("三维曲线重构");

% 末端误差计算
w1 = s(point_num) % 实际末端偏移量
w2 = abs(Z(end))% 理论末端偏移量
error = abs(w1-w2)/w1*100

% 曲面拟合
mesh_x = [X X X];mesh_x = sort(mesh_x);
mesh_y = [linspace(-12.5,-12.5,101) Y linspace(12.5,12.5,101)];mesh_y = sort(mesh_y);
[mesh_X,mesh_Y] = meshgrid(mesh_x,mesh_y);
mesh_z = [Z Z Z];mesh_z = sort(mesh_z);
mesh_Z = repmat(mesh_z,303,1);
figure(4);
mesh(mesh_X,mesh_Y,mesh_Z);xlabel("x/mm");ylabel('y/mm');zlabel('z/mm');
% axis equal;
%title("悬臂梁三维形态重构");
text(0,0,-w2,num2str(w2));