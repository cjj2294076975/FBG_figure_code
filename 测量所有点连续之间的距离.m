X = AA(1,:);Y = AA(2,:);Z = AA(3,:);
s = 0;
[~,n] = size(X);
for i =1:n-1
    ss = sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2+(Z(i+1)-Z(i))^2);
    s = s+ss;
end
disp(s)