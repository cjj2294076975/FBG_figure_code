clc;clear;close all;
%% 1mm 4mm 8mm 12mm 有温度补偿 初始温度28
correct_lambda = 3e3; %对中心波长修正3000pm
h = 0.6; %mm 厚度
L1 = 20.5; %自由端到FBG2的距离
L2 = 20.5; %FBG1到FBG2的距离
L3 = 19; %FBG1距离固定端
L = L1+L2+L3; %固定端到自由端距离
KT1 = 15.446; %15.446 pm/℃
KT2 = 14.219; %14.219 pm/℃;

disp(['自由端距离FBG1',string(L1),'FBG1距离FBG2',string(L2),'FBG2距离固定端',string(L3)]);

%% FBG1位移灵敏度系数以及应变灵敏度系数计算
s = [0 1 4 8 12];
% FBG1_lambda0 = 1553.266e3 - correct_lambda - delta_T*KT1; %20℃初始中心波长 单位pm
FBG1_lambda0 = 1550.1594e3; %20℃ 3nm修正

FBG1_lambda = [1553.266 1553.4135 1553.8665 1554.4196 1555.0]*1e3 - correct_lambda;
FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0;

FBG1_kb_s = polyfit(s,FBG1_lambda,1); %FBG1_kb_s(1): FBG1的位移灵敏度系数
FBG1_fit_s = linspace(min(s),max(s),100); %拟合的位移s
FBG1_s_fit_lambda = FBG1_kb_s(1)*FBG1_fit_s+FBG1_kb_s(2);
disp(['FBG1中心波长',string(FBG1_lambda0/1000),'位移灵敏度系数',string(FBG1_kb_s(1))]);

FBG1_epsilon = 3*(L1+L2)*h*s/(2*L^3)*10^6; %FBG1理论应变
FBG1_kb_epsilon = polyfit(FBG1_epsilon,FBG1_delta_lambda,1);
FBG1_fit_epsilon = linspace(min(FBG1_epsilon),max(FBG1_epsilon),100);

disp(['FBG1中心波长',string(FBG1_lambda0/1000),'应变灵敏度系数',string(FBG1_kb_epsilon(1))]);
% disp(['FBG1应变关于波长变化量的线性度',string(FBG1_epsilon_err)]);

%% FBG2位移灵敏度系数以及应变灵敏度系数计算
% FBG2_lambda0 = 1538.3003e3 - correct_lambda - delta_T*KT2; %20℃初始中心波长
FBG2_lambda0 = 1535.2288e3;%20℃ 3nm修正
FBG2_lambda = [1538.3003 1538.2257 1537.950 1537.6185 1537.3164]*1e3 - correct_lambda;
FBG2_delta_lambda = FBG2_lambda - FBG2_lambda0;

FBG2_kb_s = polyfit(s,FBG2_lambda,1); %FBG2_kb_s(1): FBG2的位移灵敏度系数
FBG2_fit_s = linspace(min(s),max(s),100); %拟合的位移s
% FBG2_s_fit_delta_lambda = FBG2_kb_s(1)*FBG2_fit_s + FBG2_kb_s(2);
FBG2_s_fit_lambda = FBG2_kb_s(1)*FBG2_fit_s + FBG2_kb_s(2);
% FBG2_s_err = 1 - max(abs(FBG2_delta_lambda - FBG2_s_fit_delta_lambda))/max(FBG2_delta_lambda); %FBG2位移线性度

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'位移灵敏度系数',string(FBG2_kb_s(1))]);

FBG2_epsilon = -3*L2*h*s/(2*L^3)*10^6; %FBG2理论应变 FBG1_epsilon = 3*L1*h*s/(2*L^3)*10^6; %FBG1理论应变
FBG2_kb_epsilon = polyfit(FBG2_epsilon,FBG2_delta_lambda,1);
FBG2_fit_epsilon = linspace(max(FBG2_epsilon),min(FBG2_epsilon),100);
% FBG2_e_fit_delta_lambda = FBG2_kb_epsilon(1)*FBG2_fit_epsilon+FBG2_kb_epsilon(2);
% FBG2_epsilon_err = max(abs(FBG2_delta_lambda - FBG2_e_fit_delta_lambda))/max(FBG2_delta_lambda); %FBG2应变线性度

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'应变灵敏度系数',string(FBG2_kb_epsilon(1))]);
% disp(['FBG2应变关于波长变化量的线性度',string(FBG2_epsilon_err)]);

%%绘图
FontSize = 14;
figure(1); %FBG1 位移 - 波长变化量 - 应变
scatter3(s,FBG1_epsilon,FBG1_lambda/1000);hold on;
plot3(FBG1_fit_s,FBG1_fit_epsilon,FBG1_s_fit_lambda/1000);
text(0,0,1550, {'k_{s1} = 144.083pm/mm, R^2=0.9998','k_{ε1} = 0.843pm/με, R^2=0.9998'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_1/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('λ_{B1}/nm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
zticks([1550,1550.3,1550.6,1550.9,1551.2,1551.5,1551.8,1552.1]);
zlim([1550,1552.1]);
ax = gca;
ax.TickDir = 'in';
grid on;

figure(2); %FBG1 位移 - 波长变化量 - 应变
scatter3(FBG2_epsilon,s,FBG2_lambda/1000);hold on;
plot3(FBG2_fit_epsilon,FBG2_fit_s,FBG2_s_fit_lambda/1000);
text(0,0,1535, {'k_{s2} = -82.84pm/mm, R^2=0.9988','k_{ε2} = 0.970pm/με, R^2=0.9988'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_2/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('λ_{B2}/nm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
% zticks([1531.5,1531.7,1531.9,1532.1,1532.3,1532.5,1532.7]);
% zlim([1531.5,1532.7]);
ax = gca;
ax.TickDir = 'in';
grid on;

%% 末端位移
b = 0*1e-3; %光纤半径125um 暂时不考虑
point_num = 5; % 1-5: 0 1 4 8 12（取2-5） 第i个点
delta_T = 7; %ΔT∈[-10 - 10] ℃

K1 = 2*(FBG1_delta_lambda(point_num) - KT1*delta_T)/(FBG1_kb_epsilon(1)*h)*1e-6;
K2 = -2*(FBG2_delta_lambda(point_num) - KT2*delta_T)/(FBG2_kb_epsilon(1)*h)*1e-6;

% 曲率三次样条插值
K = [K1 K1 K2 0];
sss = [0 L1 L1+L2 L];
fit_point_number = 1000;
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
for i =1:fit_point_number
    if kk(i)==0
        x(i) = 0;z(i)=ds;
    end
end
for i =1:fit_point_number
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
% AA = [roty(90) [0 0 0]';0 0 0 1]*AA;
AA = AA(1:3,:)'*roty(-90);
X = AA(:,1);Y = AA(:,2);Z = AA(:,3);
X = X';Y = Y';Z = Z';
figure(3);
plot3(X,Y,Z);
xlabel("x/mm");ylabel('y/mm');zlabel('z/mm');title("三维曲线重构");

% 末端误差计算
w1 = s(point_num); % 实际末端偏移量
w2 = abs(Z(end));% 理论末端偏移量
error = (w2-w1)/w1*100;
disp([w1,w2,error])

% 曲面拟合
mesh_x = [X';X';X'];
mesh_x = sort(mesh_x);
mesh_y = [linspace(-12.5,-12.5,fit_point_number+1)';Y';linspace(12.5,12.5,fit_point_number+1)'];
mesh_y = sort(mesh_y);
[mesh_X,mesh_Y] = meshgrid(mesh_x,mesh_y);
mesh_z = [Z Z Z];mesh_z = sort(mesh_z,'descend');
mesh_Z = repmat(mesh_z,3*(fit_point_number+1),1);

figure(4);
mesh(mesh_X,mesh_Y,mesh_Z);xlabel("x/mm");ylabel('y/mm');zlabel('z/mm');
zlim([-12,6]);
ax = gca;
ax.TickDir = 'in';
text(60,0,-w2,num2str(round(w2,3)));
