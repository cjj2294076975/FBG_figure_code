close all;clc;
%% 绘制位移s - 应变e - 波长变化量delta_lambda
s = [0     1     2     3     4     5     6     7     8     9    10    11    12];
FBG1_fit_s = [0     1     2     3     4     5     6     7     8     9    10    11    12];
FBG2_fit_s = [0     1     2     3     4     5     6     7     8     9    10    11    12];
FBG1_epsilon = 1.0e+03 *[0    0.1667    0.3333    0.5000    0.6667    0.8333    1.0000    1.1667    1.3333    1.5000 1.6667    1.8333    2.0000];
FBG2_epsilon = 1.0e+03 *[0    0.0833    0.1667    0.2500    0.3333    0.4167    0.5000    0.5833    0.6667    0.7500 0.8333    0.9167    1.0000];
FBG1_delta_lambda = 1.0e+03 *[0    0.1475    0.2837    0.4607    0.6005    0.7318    0.8790    1.0197    1.1536    1.3334 1.5026    1.6194    1.7340];
FBG2_delta_lambda = [0   74.6000  165.5500  263.8000  350.3000  440.0500  502.0000  589.9000  681.8000  751.3000 827.1750  926.6300  983.9000];

FBG1_s_fit_delta_lambda = 1.0e+03 *[0.0034    0.1498    0.2963    0.4427    0.5891    0.7356    0.8820    1.0284    1.1749    1.3213 1.4677    1.6142    1.7606];
FBG2_s_fit_delta_lambda = 1.0e+03 *[0.0071    0.0900    0.1728    0.2557    0.3386    0.4215    0.5044    0.5873    0.6702    0.7530 0.8359    0.9188    1.0017];

FBG1_fit_epsilon = 1.0e+03 *[0    0.1667    0.3333    0.5000    0.6667    0.8333    1.0000    1.1667    1.3333    1.5000 1.6667    1.8333    2.0000];
FBG2_fit_epsilon = 1.0e+03 *[0    0.0833    0.1667    0.2500    0.3333    0.4167    0.5000    0.5833    0.6667    0.7500 0.8333    0.9167    1.0000];

%%绘图
%无温度补偿FBG1
FontSize = 14;
figure(1); %FBG1 位移 - 波长变化量 - 应变
scatter3(s,FBG1_epsilon,FBG1_delta_lambda/1000);hold on;
plot3(FBG1_fit_s,FBG1_fit_epsilon,FBG1_s_fit_delta_lambda/1000);
text(0, 2000, 2,{'k_{s1} = 146.436pm/mm','k_{ε1} = 0.879pm/με'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_1/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('Δλ_{B1}/nm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
grid on;

%有温度补偿FBG1


figure(2); %FBG1 位移 - 波长变化量 - 应变
%无温度补偿FBG2
% s = [12 11 10 9 8 7 6 5 4 3 2 1 0];
scatter3(s,-FBG2_epsilon,-FBG2_delta_lambda/1000);hold on;
plot3(FBG2_fit_s,-FBG2_fit_epsilon,-FBG2_s_fit_delta_lambda/1000);
text(0, 0, 0,{'k_{s2} = -82.887pm/mm','k_{ε2}= -0.995pm/με'}, 'FontName', 'Times New Roman', 'FontSize', FontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);
xlabel('s/mm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
ylabel("ε_2/με", 'FontName', 'Times New Roman', 'FontSize', FontSize);
zlabel('Δλ_{B2}/nm', 'FontName', 'Times New Roman', 'FontSize', FontSize);
grid on;

%有温度补偿FBG2

