clearvars;
close all;
%Real length of box: 12.2 cm
%Water depth needs to be: 
starting_l = pi/0.1;
final_l = pi/10;
%increment = pi/0.1;
%n = (final_l-starting_l)/increment+1;
n=101;
increment = (final_l-starting_l)/(n-1);
gamma = 1;
j = 1;
gravity = 1;

l = starting_l:increment:final_l;
h = zeros(1,n);
hBar = zeros(1,n);
Tj = zeros(1,n);
kj = zeros(1,n);
kjBar = zeros(1,n);
for i=1:n
   
    kj(i) = (j*pi)/l(i);
    kjBar(i) = kj(i)*sqrt(gamma/gravity);
    %Tj(i) = sqrt((3*(gamma/gravity)*kj(i)^2))/((1+((gamma/gravity)*kj(i)^2)));
    Tj(i) = sqrt((3*kjBar(i)^2)/(1+kjBar(i)^2));
    h(i) = atanh(Tj(i))/kj(i);
    hBar(i) = h(i)*sqrt(gravity/gamma);
    
end
   
%figure;
%plot(h,l,'-k')
%xlabel('h'); ylabel('l');
%title('l vs h');

figure;
plot(kjBar,Tj,'-k')
xlabel('kjBar'); ylabel('Tj');
title('kjBar vs Tj');

figure;
plot(hBar,kjBar,'-k')
xlabel('hBar'); ylabel('kjBar');
title('hBar vs kjBar');

figure;
plot(hBar,l)