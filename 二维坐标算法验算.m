%% 曲率算法验算
% f(x)=sin(x) x∈[0,pi]
% L = f_0^pi sqrt(1+f(x)'^2)dx
clc;clear;close all;
syms x;
f = x.^2;
theta0 = 180 - atand(subs(diff(f),x,0));
L = int(sqrt(1+diff(f)^2),0,pi);
L = double(L);
X = [];Y = [];
K = abs(diff(f,2))/(1+diff(f,1)^2)^1.5
% K = subs(K,x,pi)
n = 11;
for i =0:n
    syms x;
    i_L = int(sqrt(1+(diff(f,1))^2),0,x); %第i段为长度为iL/n
    xx = solve(i_L-i*L/n==0)
    yy = subs(f, x, xx);
    X = [X xx];Y = [Y yy];
end
KK = [];
syms x;
for i=1:12
    k = subs(K,x,X(i));
    k = double(k);
    KK = [KK k];
end
m = 100;
sss = linspace(0,L,m);
KKK = spline(linspace(0,L,12),KK,sss); %KKK三次样条曲率插值
lK = [];
for i=1:11
    lk = linspace(KK(i),KK(i+1),10);%线性曲率插值
    lK = [lK lk];
end
%% 三次样条曲率插值
% figure(1); %绘制曲率
% scatter(linspace(0,L,12),KK,'k');hold on;plot(sss,KKK,'r');hold on;plot(linspace(0,L,110),lK,'.-b');
% title('曲率');xlabel('弧长');ylabel('曲率');legend('实际曲率','三次样条插值曲率','线性插值曲率');
[~,m] = size(KKK);
ds = L/m;
R = 1./KKK;
theta = KKK*ds;

%% 局部坐标 
x = R.*(1 - cos(theta));
y = R.*sin(theta);
%% 全局坐标
for i =1:m
    if KKK(i)==0
        x(i) = 0;y(i)=ds;
    end
end
AA = [0 0 1]';
T = [cosd(theta0) -sind(theta0) 0;sind(theta0) cosd(theta0) 0; 0 0 1];
for i =1:m
    if i==1
        A = [0 ds 1]';
    else   
        P = [1 0 x(i-1);0 1 y(i-1); 0 0 1];
        RR = [cos(theta(i-1)) -sin(theta(i-1)) 0;sin(theta(i-1)) cos(theta(i-1)) 0; 0 0 1];
        t = P*RR;
        T = T*inv(t);
        A = T*[x(i) y(i) 1]';
    end
    AA = [AA A];
end
figure(2);
XX = AA(1,:);ZZ = AA(2,:);plot(XX,ZZ,'r');hold on;
scatter(X,Y,'k');title('实际值-拟合值误差计算');hold on;

%% 线性曲率插值
[~,m] = size(lK);
ds = L/m;
R = 1./lK;
theta = lK*ds;

%% 局部坐标 
x = R.*(1 - cos(theta));
y = R.*sin(theta);
%% 全局坐标
for i =1:m
    if lK(i)==0
        x(i) = 0;y(i)=ds;
    end
end
AA = [0 0 1]';
T = [cosd(theta0) -sind(theta0) 0;sind(theta0) cosd(theta0) 0; 0 0 1];
for i =1:m
    if i==1
        A = T*[0 ds 1]';
    else   
        P = [1 0 x(i-1);0 1 y(i-1); 0 0 1];
        RR = [cos(theta(i-1)) -sin(theta(i-1)) 0;sin(theta(i-1)) cos(theta(i-1)) 0; 0 0 1];
        T = T*inv(P*RR);
        A = T*[x(i) y(i) 1]';
    end
    AA = [AA A];
end
XX = AA(1,:);ZZ = AA(2,:);plot(XX,ZZ,'.-b');
xlabel("s/mm");ylabel('s/mm');title("曲线重构");legend('控制点','三次样条拟合曲线','线性插值拟合曲线');
