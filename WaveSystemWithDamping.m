%Damping

clearvars;
close all;

%Time Parameters
h = 0.01; %Time Step
tfinal = 250; %Final Time

n = round(tfinal/h); % Number of Time Steps
tv = 0:h:tfinal; %Time vector

%Equation Parameters
length = 5; %Length of Box
f_0 = 0; %Forcing coefficient
gamma = 1;
gravity = 1;
j = 1;
%damping = 0.5;

kj = (j * pi)/length;
k2j = kj*2;

gj = gravity + gamma * (kj^2);
g2j = gravity + gamma * (kj^2);

%Tj = tanh(kj*hei);
Tj = sqrt((3*gamma*(kj^2)/gravity)/(1+(gamma*(kj^2)/gravity)));
%T2j = tanh(k2j * hei);
T2j = (2*Tj)/(1+Tj^2);

height = atanh(Tj)/kj;

omega = 2*sqrt(gj*kj*Tj); %Frequency of Forcing

omegaj = sqrt(gj*kj*Tj);
omega2j = sqrt(g2j*k2j*T2j);

damping = 2*gj*kj*Tj; %critical damping value

alpha = 0; %Non-Linearity (0: Off, 1: On)

%Initial Conditions System of 4
aj_0 = 0.01; 
bj_0 = 0;
a2j_0 = 0;
b2j_0 = 0;

%RK-4 Calculations System of 4
aj_values = zeros(1, n+1);
aj_values(1) = aj_0;
bj_values = zeros(1, n+1);
bj_values(1) = bj_0;
a2j_values = zeros(1, n+1);
a2j_values(1) = a2j_0;
b2j_values = zeros(1, n+1);
b2j_values(1) = b2j_0;


for i=1:n %RK-4
%System of Four
    %k1 terms
    k11 = f1(tv(i),aj_values(i), bj_values(i), a2j_values(i), b2j_values(i), kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k21 = f2(tv(i),aj_values(i), bj_values(i), a2j_values(i), b2j_values(i), kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k31 = f3(tv(i),aj_values(i), bj_values(i), a2j_values(i), b2j_values(i), kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k41 = f4(tv(i),aj_values(i), bj_values(i), a2j_values(i), b2j_values(i), kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    
    %k2 terms
    k12 = f1(tv(i),aj_values(i)+0.5*h*k11, bj_values(i)+0.5*h*k21, a2j_values(i)+0.5*h*k31, b2j_values(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k22 = f2(tv(i),aj_values(i)+0.5*h*k11, bj_values(i)+0.5*h*k21, a2j_values(i)+0.5*h*k31, b2j_values(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k32 = f3(tv(i),aj_values(i)+0.5*h*k11, bj_values(i)+0.5*h*k21, a2j_values(i)+0.5*h*k31, b2j_values(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k42 = f4(tv(i),aj_values(i)+0.5*h*k11, bj_values(i)+0.5*h*k21, a2j_values(i)+0.5*h*k31, b2j_values(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    
    %k3 terms
    k13 = f1(tv(i),aj_values(i)+0.5*h*k12, bj_values(i)+0.5*h*k22, a2j_values(i)+0.5*h*k32, b2j_values(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k23 = f2(tv(i),aj_values(i)+0.5*h*k12, bj_values(i)+0.5*h*k22, a2j_values(i)+0.5*h*k32, b2j_values(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k33 = f3(tv(i),aj_values(i)+0.5*h*k12, bj_values(i)+0.5*h*k22, a2j_values(i)+0.5*h*k32, b2j_values(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k43 = f4(tv(i),aj_values(i)+0.5*h*k12, bj_values(i)+0.5*h*k22, a2j_values(i)+0.5*h*k32, b2j_values(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
   
    %k4 terms
    k14 = f1(tv(i),aj_values(i)+h*k13, bj_values(i)+h*k23, a2j_values(i)+h*k33, b2j_values(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k24 = f2(tv(i),aj_values(i)+h*k13, bj_values(i)+h*k23, a2j_values(i)+h*k33, b2j_values(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k34 = f3(tv(i),aj_values(i)+h*k13, bj_values(i)+h*k23, a2j_values(i)+h*k33, b2j_values(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    k44 = f4(tv(i),aj_values(i)+h*k13, bj_values(i)+h*k23, a2j_values(i)+h*k33, b2j_values(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping);
    
    %rk4 formula
    aj_values(i+1) = aj_values(i)+(h/6)*(k11+2*k12+2*k13+k14);
    bj_values(i+1) = bj_values(i)+(h/6)*(k21+2*k22+2*k23+k24);
    a2j_values(i+1) = a2j_values(i)+(h/6)*(k31+2*k32+2*k33+k34);
    b2j_values(i+1) = b2j_values(i)+(h/6)*(k41+2*k42+2*k43+k44);
end

figure;
plot(tv,aj_values,'-k')
xlabel('t'); ylabel('aj'); 
title('RK Solution for aj');

figure;
plot(tv,bj_values,'-k')
xlabel('t'); ylabel('bj'); 
title('RK Solution for bj');

figure;
plot(aj_values,bj_values,'-k')
xlabel('aj'); ylabel('bj');
title('RK solution for bj/aj');
if(alpha == 1)
    figure;
    plot(tv,a2j_values,'-k')
    xlabel('t'); ylabel('a2j'); 
    title('RK Solution for a2j');

    figure;
    plot(tv,b2j_values,'-k')
    xlabel('t'); ylabel('b2j'); 
    title('RK Solution for b2j');
end

function dajdt = f1(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping)
    dajdt = kj*Tj*bj + ((2*(kj^2)*(1-Tj*T2j)*aj*b2j - (kj^2)*(1+(Tj^2))*a2j*bj))*alpha;
end

function dbjdt = f2(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping)
    dbjdt = -(f_0*cos(omega*t)+gj)*aj - (damping * bj) - ((2*(kj^2)*(1-Tj*T2j)*bj*b2j))*alpha;
end

function da2jdt = f3(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping)
    da2jdt = k2j*T2j*b2j + (2*(kj^2)*(1-Tj*T2j)*aj*bj)*alpha;
end

function db2jdt = f4(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f_0, omega, alpha, damping)
    db2jdt = -(f_0*cos(omega*t)+g2j)*a2j - (damping * b2j) + (.5*(kj^2)*(1+(Tj^2))*(bj^2))*alpha;
end
