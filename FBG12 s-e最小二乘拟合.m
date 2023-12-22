s = [0 1 4 7 11];
FBG1_e = [0 92 310 400 650];
FBG2_e = [0 450 900 1700 2500];
kb1 = polyfit(s,FBG1_e,1)
kb2 = polyfit(s,FBG2_e,1)
ss = linspace(0,11,100);
FBG1_ee = kb1(1)*ss+kb1(2);
FBG2_ee = kb2(1)*ss+kb2(2);
figure(1);
scatter(s,FBG1_e);hold on;
plot(ss,FBG1_ee,'-')
figure(2);
scatter(s,FBG2_e);hold on;
plot(ss,FBG2_ee,'-')