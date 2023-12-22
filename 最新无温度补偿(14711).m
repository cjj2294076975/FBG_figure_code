clc;clear;close all;
%% 1mm 4mm 7mm 11mm 无温度补偿
h = 0.6; %mm 厚度
L1 = 21; %自由端到FBG2的距离
L2 = 21; %FBG1到FBG2的距离
L3 = 19; %FBG1距离固定端
L = L1+L2+L3; %固定端到自由端距离
disp(['自由端距离FBG1',string(L1),'FBG1距离FBG2',string(L2),'FBG2距离固定端',string(L3)]);

%% FBG1位移灵敏度系数以及应变灵敏度系数计算
s = [0 1 4 7 11];
FBG1_lambda0 = 1553.266e3 - 5.69e3;
FBG1_lambda = [1553.266 1553.4135 1553.8665 1554.28575 1554.8854]*1e3 - 5.69e3;
FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0;

FBG1_kb_s = polyfit(s,FBG1_delta_lambda,1); %FBG1_kb_s(1): FBG1的位移灵敏度系数
FBG1_fit_s = linspace(min(s),max(s),100); %拟合的位移s
FBG1_s_fit_delta_lambda = FBG1_kb_s(1)*FBG1_fit_s+FBG1_kb_s(2);
disp(['FBG1中心波长',string(FBG1_lambda0/1000),'位移灵敏度系数',string(FBG1_kb_s(1))]);

FBG1_epsilon = 3*(L1+L2)*h*s/(2*L^3)*10^6; %FBG1理论应变
FBG1_kb_epsilon = polyfit(FBG1_epsilon,FBG1_delta_lambda,1)
k1 = 3*(L1+L2)*h/(2*L^3)*1e6;
FBG1_k_epsilon = FBG1_kb_s(1)/k1
FBG1_fit_epsilon = linspace(min(FBG1_epsilon),max(FBG1_epsilon),100);
% FBG1_e_fit_delta_lambda = FBG1_kb_epsilon(1)*FBG1_fit_epsilon+FBG1_kb_epsilon(2);
% FBG1_epsilon_err = max(abs(FBG1_delta_lambda - FBG1_e_fit_delta_lambda))/max(FBG1_delta_lambda); %FBG1应变线性度

disp(['FBG1中心波长',string(FBG1_lambda0/1000),'应变灵敏度系数',string(FBG1_kb_epsilon(1))]);
% disp(['FBG1应变关于波长变化量的线性度',string(FBG1_epsilon_err)]);

%% FBG2位移灵敏度系数以及应变灵敏度系数计算
FBG2_lambda0 = 1538.3003e3 - 5.69e3;
FBG2_lambda = [1538.3003 1538.2257 1537.950 1537.7104 1537.37367]*1e3 - 5.69e3;
FBG2_delta_lambda = FBG2_lambda0 - FBG2_lambda;

FBG2_kb_s = polyfit(s,FBG2_delta_lambda,1); %FBG2_kb_s(1): FBG2的位移灵敏度系数
FBG2_fit_s = linspace(min(s),max(s),100); %拟合的位移s
FBG2_s_fit_delta_lambda = FBG2_kb_s(1)*FBG2_fit_s + FBG2_kb_s(2);
% FBG2_s_err = 1 - max(abs(FBG2_delta_lambda - FBG2_s_fit_delta_lambda))/max(FBG2_delta_lambda); %FBG2位移线性度

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'位移灵敏度系数',string(FBG2_kb_s(1))]);

FBG2_epsilon = 3*L2*h*s/(2*L^3)*10^6; %FBG2理论应变 FBG1_epsilon = 3*L1*h*s/(2*L^3)*10^6; %FBG1理论应变
FBG2_kb_epsilon = polyfit(FBG2_epsilon,FBG2_delta_lambda,1);
FBG2_fit_epsilon = linspace(min(FBG2_epsilon),max(FBG2_epsilon),100);
% FBG2_e_fit_delta_lambda = FBG2_kb_epsilon(1)*FBG2_fit_epsilon+FBG2_kb_epsilon(2);
% FBG2_epsilon_err = max(abs(FBG2_delta_lambda - FBG2_e_fit_delta_lambda))/max(FBG2_delta_lambda); %FBG2应变线性度

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'应变灵敏度系数',string(FBG2_kb_epsilon(1))]);
% disp(['FBG2应变关于波长变化量的线性度',string(FBG2_epsilon_err)]);

%%绘图
FontSize = 14;
figure(1); %FBG1 位移 - 波长变化量 - 应变
scatter3(s,FBG1_epsilon,FBG1_delta_lambda/1000);hold on;
plot3(FBG1_fit_s,FBG1_fit_epsilon,FBG1_s_fit_delta_lambda/1000);
text(0, 100, {'k_{s1} = 146.801pm/mm','k_{ε1} = 0.881pm/με'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_1/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('Δλ_{B1}/nm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
grid on;

figure(2); %FBG1 位移 - 波长变化量 - 应变
scatter3(s,FBG2_epsilon,FBG2_delta_lambda/1000);hold on;
plot3(FBG2_fit_s,FBG2_fit_epsilon,FBG2_s_fit_delta_lambda/1000);
text(0, 100, {'k_{s2} = 84.5745pm/mm','k_{ε2} = 1.0149pm/με'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_2/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('Δλ_{B2}/nm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
grid on;

%% 末端位移
b = 0*1e-3; %光纤半径125um
point_num = 5; % 1-5: 0 1 4 7 11（取2-5） 第i个点
K1 = 2*FBG1_delta_lambda(point_num)/((h+2*b)*FBG1_kb_epsilon(1))*1e-6;
K2 = 2*FBG2_delta_lambda(point_num)/((h+2*b)*FBG2_kb_epsilon(1))*1e-6;

% 曲率三次样条插值
K = [K1 K1 K2 0];
sss = [0 L1 L1+L2 L];
fit_point_number = 100;
[ss,kk] = spline_interf(sss,K,fit_point_number);

% 绘制L关于曲率的曲线
figure(3);
plot(ss,kk);xlabel('s/mm');ylabel('K/mm^-1');title('位移关于曲率的关系');
R = 1./kk;
ds = L/fit_point_number; %连续点的间距
theta = kk*ds; %theta = K*ds
% 局部坐标
x = R - R.*cos(theta);
y = zeros(1,fit_point_number);
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
% figure(3);
% plot3(X,Y,Z);
% xlabel("x/mm");ylabel('y/mm');zlabel('z/mm');title("三维曲线重构");

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
%title("悬臂梁三维形态重构");
text(0,0,-w2,num2str(round(w2,2)));