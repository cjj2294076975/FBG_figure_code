clc;clear;close all;

%% 1mm 4mm 8mm 12mm 无温度补偿
correct_lambda = 3e3; %对中心波长修正3000pm
h = 0.6; %厚度
L1 = 20.5; %FBG1到自由端的距离
L2 = 41; %FBG2到自由端的距离
L = 60; %固定端到自由端距离
KT1 = 14.219; %14.219 pm/℃; 已修改 1535
KT2 = 15.446; %15.446 pm/℃; 已修改 1550

Temperature = 20;
FBG1_lambda0 = 1000*(0.014218706020243*Temperature + 1.534944402260315e+03); %Temperature对应的FBG1的中心波长
FBG2_lambda0 = 1000*(0.015445655439910*Temperature + 1.549850491573111e+03); %Temperature对应的FBG2的中心波长

disp(['FBG1距离自由端',string(L1),'FBG2距离自由端',string(L2),'悬臂梁总长',string(L)]);

%% FBG1位移灵敏度系数以及应变灵敏度系数计算 已修改
FBG1_lambda = [1538.3003 1538.2257 1537.950 1537.6185 1537.3164]*1e3 - correct_lambda;%已修改
FBG1_delta_lambda = FBG1_lambda - FBG1_lambda0; %已修改

s = [0 1 4 8 12];
FBG1_kb_s = polyfit(s,FBG1_delta_lambda,1); %FBG2_kb_s(1): FBG1的位移灵敏度系数 已修改
FBG1_fit_s = linspace(min(s),max(s),100); %拟合的位移s 已修改
FBG1_s_fit_delta_lambda = FBG1_kb_s(1)*FBG1_fit_s + FBG1_kb_s(2); %已修改
FBG1_kb_s_corr_matrix = corrcoef(s, FBG1_delta_lambda);%已修改
FBG1_kb_s_R2 = FBG1_kb_s_corr_matrix(1, 2)^2;%已修改
disp(['FBG1中心波长',string(FBG2_lambda0/1000),'位移灵敏度系数',string(FBG1_kb_s(1))]);%已修改

FBG1_epsilon = -3*L1*h*s/(2*L^3)*10^6; %FBG1理论应变
FBG1_kb_epsilon = polyfit(FBG1_epsilon,FBG1_delta_lambda,1);
FBG1_fit_epsilon = linspace(max(FBG1_epsilon),min(FBG1_epsilon),100);
FBG1_kb_epsilon_corr_matrix = corrcoef(FBG1_epsilon, FBG1_delta_lambda);
FBG1_kb_epsilon_R2 = FBG1_kb_s_corr_matrix(1, 2)^2;

disp(['FBG1中心波长',string(FBG1_lambda0/1000),'应变灵敏度系数',string(FBG1_kb_epsilon(1))]);

%% FBG2位移灵敏度系数以及应变灵敏度系数计算

FBG2_lambda = [1553.266 1553.4135 1553.8665 1554.4196 1555.0]*1e3 - correct_lambda;% 已修改
FBG2_delta_lambda = FBG2_lambda - FBG2_lambda0; % 已修改

FBG2_kb_s = polyfit(s,FBG2_delta_lambda,1); %FBG2_kb_s(1): FBG2的位移灵敏度系数 已修改
FBG2_fit_s = linspace(min(s),max(s),100); %拟合的位移s
FBG2_s_fit_delta_lambda = FBG2_kb_s(1)*FBG2_fit_s+FBG2_kb_s(2); % 已修改
disp(['FBG2中心波长',string(FBG2_lambda0/1000),'位移灵敏度系数',string(FBG2_kb_s(1))]);% 已修改
FBG2_kb_s_corr_matrix = corrcoef(s, FBG2_delta_lambda);% 已修改

FBG2_epsilon = 3*L2*h*s/(2*L^3)*10^6; %FBG2理论应变 已修改
FBG2_kb_epsilon = polyfit(FBG2_epsilon,FBG2_delta_lambda,1); %已修改
FBG2_fit_epsilon = linspace(min(FBG2_epsilon),max(FBG2_epsilon),100); %已修改
FBG2_kb_epsilon_corr_matrix = corrcoef(FBG2_epsilon, FBG2_delta_lambda);%已修改
FBG2_kb_epsilon_R2 = FBG2_kb_s_corr_matrix(1, 2)^2;%已修改

disp(['FBG2中心波长',string(FBG2_lambda0/1000),'应变灵敏度系数',string(FBG2_kb_epsilon(1))]);

%% 末端位移 不修改
point_num = 5; % 1-5: 0 1 4 8 12（取2-5） 第i个点
delta_T = 6.5; %ΔT∈[-10 - 10] ℃

K1 = -2*(FBG1_delta_lambda(point_num) - KT1*delta_T)/(FBG1_kb_epsilon(1)*h)*1e-6;
K2 = 2*(FBG2_delta_lambda(point_num) - KT2*delta_T)/(FBG2_kb_epsilon(1)*h)*1e-6;

% 曲率三次样条插值
K = [K2 K2 K1 0]; %固定端 FBG2 FBG1 自由端
sss = [0 L-L2 L-L1 L]; 
fit_point_number = 1000; %拟合点数
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
w2 = round(abs(Z(end)),3);% 理论末端偏移量
error = (w2-w1)/w1*100; %百分比
disp([w1,w2,error])

% 曲面拟合
mesh_x = [X X X];
mesh_x = sort(mesh_x);
mesh_y = [linspace(-12.5,-12.5,fit_point_number+1) Y linspace(12.5,12.5,fit_point_number+1)];
mesh_y = sort(mesh_y);
[mesh_X,mesh_Y] = meshgrid(mesh_x,mesh_y);
mesh_z = [Z Z Z];mesh_z = sort(mesh_z,'descend');
mesh_Z = repmat(mesh_z,3*(fit_point_number+1),1);
figure(4);
mesh(mesh_X,mesh_Y,mesh_Z);xlabel("x/mm");ylabel('y/mm');zlabel('z/mm');
zlim([-12,6]);
% title("悬臂梁三维形态重构");
text(60,0,-w2,num2str(round(w2,3)));
ax = gca;
ax.TickDir = 'in';
