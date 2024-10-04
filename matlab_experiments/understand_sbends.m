Fam='test'; 
L=0.204240; 
Angle=0.010851; 
K1=-0.801418; 

Klist = -10:0.1:10;
R15 = [];
R25 = [];
for index = 1:length(Klist)
  R = findelemm66(atsbend(Fam,L,Angle,Klist(index)));
  R15(index) = R(1,5);
  R25(index) = R(2,5);
end
disp(R);

figure(1);
plot(Klist, R15, '-x');

figure(2);
plot(Klist, R25, '-x');

