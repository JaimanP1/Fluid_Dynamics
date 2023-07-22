%envelope equation

clear all;

%% DO NOT CHANGE

C1 = ParametersClass.getConstant1();
C2 = ParametersClass.getConstant2();
freq = ParametersClass.getFrequency();

%%

tic
%initial condition
IC_A = complex(.001,.0);
IC_A2 = complex(.0,-.0);
x0 = [IC_A; IC_A2];
%complex IC?

%compute trajectory
timeVars = ParametersClass.getTimeVars();
dt = timeVars(1);
t_final = timeVars(2);
tspan = 0:dt:t_final;
n = timeVars(3);

%solution vector
Y = zeros(2,n);

%time stamp 17:30 for more efficient way to do this
Y(:,1) = x0;
xin = x0;
for i=1:tspan(end)/dt
    %better way to do the for loop?
    time = i*dt;
    xout = rk4SingleStep(@(t,x)EnvelopeSystem_v3(t,x),dt,time,xin);
    %navigate to the next empty column and input xout
    Y(:,i) = xout;
    %X = [X xout];
    xin = xout;
end

E = C2*(abs(Y(1,:))).^2 + C1*(abs(Y(2,:))).^2;
dA1dt = i*C1*Y(2,:).*conj(Y(1,:));
L = (abs(Y(1,:)).^2).*abs(Y(2,:)).*cos((2*angle(Y(1,:))-angle(Y(2,:))));
%option 2. WORKS

toc

subplot(3,2,1), plot(tspan,abs(Y(1,:)), 'k'), title('modulus A1 vs t');
subplot(3,2,2), plot(tspan,abs(Y(2,:)), 'k'), title('modulus A2 vs t');
subplot(3,2,3), plot(abs(Y(1,:)), abs(dA1dt)), title('mod dA1dt vs mod A1');
%check index thing
subplot(3,2,4), plot(tspan,E, 'k'), title('E vs t');
subplot(3,2,5), plot(tspan,L), title('L vs t')
