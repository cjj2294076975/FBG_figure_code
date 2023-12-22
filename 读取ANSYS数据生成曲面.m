%% 读取ANSYS数据生成曲面
clc;clear;close all;
data = xlsread("C:\Users\Tsai\Desktop\悬臂梁论文\ansys_excel\11.28最新14812\12mm.xlsx");
x = data(:,1)+10;
y = data(:,2);
z = data(:,3);
[m,n] = size(data);
X = [];Y = [];Z = [];
for i=1:m
    if x(i)<0
        X = [X x(i)];Y = [Y y(i)];Z = [Z z(i)];
    end
end
AA = [X;Y;Z];
AA = AA'*roty(-180);
X = AA(:,1)';Y = AA(:,2)';Z = AA(:,3)';
[m,n] = size(X);
% X = sort(X);
X = [X;X;X];
Z = [Z;Z;Z];
Y = [linspace(-12.5,-12.5,n);Y;linspace(12.5,12.5,n)];
mesh(X,Y,Z)
xlabel("x/mm");ylabel("y/mm");zlabel("z/mm");
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
text(60,0,-1,num2str(abs(min(Z(:)))), 'FontName', 'Times New Roman', 'FontSize', 10);