clearvars;
close all;

% input parameter
h = 0.01; % time step
tfinal = 25; % final time for computation
k = 0.01; % increment for tank length
u = 0.01; % increment for tank width
lfinal = 5;
wfinal = 5;
% the number of time steps
n = round(tfinal/h); % the number of time steps
m = round(lfinal/k); % number of length steps
p = round(wfinal/u); % number of width steps
tv = 0:h:tfinal; % vector for time
% tv = linspace(0,tfinal,n+1);

%Parameters
len = 5;
f0 = 0;
gam = 1;
j = 4;

kj = (j*pi)/len;
k2j = kj*2;

gravity = 1;
gj = gravity + gam * (kj^2);
g2j = gravity + gam * (k2j^2);

%Tj = tanh(kj*hei);
Tj = sqrt((3*gam*(kj^2)/gravity)/(1+(gam*(kj^2)/gravity)));
%T2j = tanh(k2j * hei);
T2j = (2*Tj)/(1+Tj^2);

hei = atanh(Tj)/kj;

omega = 2*sqrt(gj*kj*Tj); %frequency
%omega = 1;

omegaj = sqrt(gj*kj*Tj)
omega2j = sqrt(g2j*k2j*T2j)

%Nonlinearity on/off
alpha = 1;
third = 1;

C1 = (kj^2/(2*omegaj))*((2*g2j+gj+gj*Tj^2)-(2*g2j*Tj*T2j))
C2 = ((kj^2*gj)/(2*omegaj*g2j))*((2*g2j+gj+gj*Tj^2)-(2*g2j*Tj*T2j))



%Initial Conditions for System of 4
x10 = 0.001;% + 0.002i; % initial condition for x1, aj
x20 = 0; % initial condition for x2, bj
x30 = x10 * (1/sqrt(2)); %- 0.002i; % initial condition for x3, a2j
x40 = 0; % initial condition for x4, b2j

%For Fixed Point Sol A2j = (1/sqrt(2)) * (Aj)

%Initial Conditions for Envelope
x50 = real(x10)/2; % Aj
x60 = real(x30)/2; % A2j

x20 = -2*imag(x50)*(gj/omegaj); % initial condition for x2, bj
x40 = -2*imag(x60)*(g2j/omega2j); % initial condition for x4, b2j


% Improved Euler method (RK-4)
x1_num2 = zeros(1,n+1); % initialization of a row vector x1_num
x1_num2(1) = x10;

x2_num2 = zeros(1,n+1); % initialization of a row vector x2_num
x2_num2(1) = x20;

x3_num2 = zeros(1,n+1); % initialization of a row vector x3_num
x3_num2(1) = x30;

x4_num2 = zeros(1,n+1); % initialization of a row vector x4_num
x4_num2(1) = x40;

x5_num2 = zeros(1, n+1);
x5_num2(1) = x50;

x6_num2 = zeros(1, n+1);
x6_num2(1) = x60;

tankl = 0:k:lfinal;
tankw = 0:u:lfinal;
solution = zeros(m+1,n+1);

for i=1:n %RK-4
%System of Four
    %k1 terms
    k11 = f1(tv(i),x1_num2(i), x2_num2(i), x3_num2(i), x4_num2(i), kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k21 = f2(tv(i),x1_num2(i), x2_num2(i), x3_num2(i), x4_num2(i), kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k31 = f3(tv(i),x1_num2(i), x2_num2(i), x3_num2(i), x4_num2(i), kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    k41 = f4(tv(i),x1_num2(i), x2_num2(i), x3_num2(i), x4_num2(i), kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    
    %k2 terms
    k12 = f1(tv(i),x1_num2(i)+0.5*h*k11, x2_num2(i)+0.5*h*k21, x3_num2(i)+0.5*h*k31, x4_num2(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k22 = f2(tv(i),x1_num2(i)+0.5*h*k11, x2_num2(i)+0.5*h*k21, x3_num2(i)+0.5*h*k31, x4_num2(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k32 = f3(tv(i),x1_num2(i)+0.5*h*k11, x2_num2(i)+0.5*h*k21, x3_num2(i)+0.5*h*k31, x4_num2(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    k42 = f4(tv(i),x1_num2(i)+0.5*h*k11, x2_num2(i)+0.5*h*k21, x3_num2(i)+0.5*h*k31, x4_num2(i)+0.5*h*k41, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    
    %k3 terms
    k13 = f1(tv(i),x1_num2(i)+0.5*h*k12, x2_num2(i)+0.5*h*k22, x3_num2(i)+0.5*h*k32, x4_num2(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k23 = f2(tv(i),x1_num2(i)+0.5*h*k12, x2_num2(i)+0.5*h*k22, x3_num2(i)+0.5*h*k32, x4_num2(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k33 = f3(tv(i),x1_num2(i)+0.5*h*k12, x2_num2(i)+0.5*h*k22, x3_num2(i)+0.5*h*k32, x4_num2(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    k43 = f4(tv(i),x1_num2(i)+0.5*h*k12, x2_num2(i)+0.5*h*k22, x3_num2(i)+0.5*h*k32, x4_num2(i)+0.5*h*k42, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
   
    %k4 terms
    k14 = f1(tv(i),x1_num2(i)+h*k13, x2_num2(i)+h*k23, x3_num2(i)+h*k33, x4_num2(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k24 = f2(tv(i),x1_num2(i)+h*k13, x2_num2(i)+h*k23, x3_num2(i)+h*k33, x4_num2(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third);
    k34 = f3(tv(i),x1_num2(i)+h*k13, x2_num2(i)+h*k23, x3_num2(i)+h*k33, x4_num2(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    k44 = f4(tv(i),x1_num2(i)+h*k13, x2_num2(i)+h*k23, x3_num2(i)+h*k33, x4_num2(i)+h*k43, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha);
    
    %rk4 formula
    x1_num2(i+1) = x1_num2(i)+(h/6)*(k11+2*k12+2*k13+k14);
    x2_num2(i+1) = x2_num2(i)+(h/6)*(k21+2*k22+2*k23+k24);
    x3_num2(i+1) = x3_num2(i)+(h/6)*(k31+2*k32+2*k33+k34);
    x4_num2(i+1) = x4_num2(i)+(h/6)*(k41+2*k42+2*k43+k44);
    
    %Surface Solution
    %solution(:,i) = x1_num2(i)*exp(-i*kj*tank);
 
    solution(:,i) = x1_num2(i)*(cos(kj*tankl)-i*sin(kj*tankl));
    
end
% 
myWriter = VideoWriter('Wave');
myWriter.FrameRate = 15;
open(myWriter);


% Plotting the first iteration
% figure;
% g = plot(tankl,solution(:,1),'-b');
% xlabel('Tank Length'); ylabel('Tank Height');
% xlim([0 5])
% ylim([-5 5])
figure;
x0=0;
y0=0;
width=1920;
height=1080;
set(gcf,'position',[x0,y0,width,height])
[X, Y] = meshgrid(tankl, tankw);
Z = x1_num2(1) * (cos(kj*X)-i*sin(kj*X)) + x3_num2(1) * (cos(k2j*X)-i*sin(k2j*X));
s = surf(X,Y,Z);
colormap(winter);
shading interp;
brighten(0);
%xlabel('Tank Length'); ylabel('Tank Width'); zlabel('Tank Height');
xlim([0 5])
ylim([0 5])
zlim([-5 5])
% Iterating through the length of the time array
for k = 2:5:n+1
     % Updating the line
     s.ZData = x1_num2(k) * (cos(kj*X)-i*sin(kj*X)) + x3_num2(k) * (cos(k2j*X)-i*sin(k2j*X));;
     %g.YData = solution(:,k);
     % Updating the title
     title(sprintf('Time: %0.2f sec', tv(k)),...
     'Interpreter','Latex');
     % Delay
     pause(0.001)
     movieVector = getframe(gcf);
     writeVideo(myWriter, movieVector);
end
 
close(myWriter);

function dajdt = f1(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha, third)
    dajdt = kj*Tj*bj + ((2*(kj^2)*(1-Tj*T2j)*aj*b2j - (kj^2)*(1+(Tj^2))*a2j*bj))*alpha - ((Tj*kj^3)*(3-2*Tj*T2j)*(bj*aj^2))*third;
end

function dbjdt = f2(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha,third)
 dbjdt = -(f0*cos(omega*t)+gj)*aj - ((2*(kj^2)*(1-Tj*T2j)*bj*b2j))*alpha+((Tj*kj^3)*(3-2*Tj*T2j)*(aj*bj^2)+(2/3)*(kj^4)*(aj^3))*third;
end

function da2jdt = f3(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha)
 da2jdt = (k2j*T2j*b2j) + (2*(kj^2)*(1-Tj*T2j)*aj*bj)*alpha;
end

function db2jdt = f4(t, aj, bj, a2j, b2j, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha)
 db2jdt = -g2j*a2j + (.5*(kj^2)*(1+(Tj^2))*(bj^2))*alpha;
end

function dAjdt = f5(t, Aj, A2j, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha, C1, C2)
 dAjdt = 1i*A2j*conj(Aj)*C1;
end

function dA2jdt = f6(t, Aj, A2j, kj, k2j, Tj, T2j, gj, g2j, f0, omega, alpha, C1, C2)
 dA2jdt = 1i*(Aj^2)*C2;
end

