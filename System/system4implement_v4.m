clear all; %doing this so can call envelope solutions



%make sure these match with system4 file!!!!!
%find way to consolidate into one function
%using these to compute hamiltonian

%% Calling the Envelope equation
EnvelopeImplement_v4;
%%
tic

k = ParametersClass.getWavenumber();
k2 = ParametersClass.getWavenumber2();
T = ParametersClass.getHyper_tang();
T2 = ParametersClass.getHyper_tang2();
g = ParametersClass.getModified_g();
g2 = ParametersClass.getModified_g2();    
freq = ParametersClass.getFrequency();
freq2 = ParametersClass.getFrequency2();
alpha = ParametersClass.getAlpha();
forcing = ParametersClass.getForcing();


%%

%initial condition
IC_a = 2*real(IC_A);
IC_a2 = 2*real(IC_A2);
IC_b = -2*imag(IC_A)*(g/freq);
IC_b2 = -2*imag(IC_A2)*(g2/freq2);
x0 = [IC_a; IC_a2; IC_b; IC_b2; IC_a; IC_b];
%last two are for linear conditions

%compute trajectory
timeVars = ParametersClass.getTimeVars();
dt = timeVars(1);
t_final = timeVars(2);
tspan = 0:dt:t_final;
n = timeVars(3);

%solution vector
X = zeros(6,n);

%% 


X(:,1) = x0;
xin = x0;
for i=1:tspan(end)/dt
    time = i*dt;
    xout = rk4SingleStep(@(t,x)system4_v6(t,x),dt,time,xin);
    
    X(:,i) = xout;
    
    %X = [X xout];
    %old way added a column after every iteration
    %for t_final = 1000, went from 206 to .75 seconds

    xin = xout;

end

H = .5*(g*(X(1,:).^2) + k*T*(X(3,:).^2)) + .5*(g2*(X(2,:).^2) + k2*T2*(X(4,:).^2)) + alpha*2*(k^2)*(1 - T*T2)*X(1,:).*X(3,:).*X(4,:) - alpha*.5*(k^2)*(1 + (T^2))*(X(3,:).^2).*X(2,:);

toc
%%

figure;
subplot(3,3,1), plot(tspan,X(1,:), 'k'), title('a vs t');
subplot(3,3,2), plot(tspan,X(2,:), 'k'), title('a2 vs t');
subplot(3,3,3), plot(tspan,X(3,:), 'k'), title('b vs t');
subplot(3,3,4), plot(tspan,X(4,:), 'k'), title('b2 vs t');
subplot(3,3,5), plot(tspan,H, 'k'), title('Hamiltonian vs t');
subplot(3,3,6), plot(X(1,:),X(2,:), 'k'), title('a vs a2');
subplot(3,3,7), plot(X(3,:),X(4,:), 'k'), title('b vs b2');
subplot(3,3,8), plot(X(1,:),X(3,:), 'k'), title('a vs b');
subplot(3,3,9), plot(X(2,:),X(4,:), 'k'), title('a2 vs b2');

%%
%finding the first 500 time steps
linear_time = length(tspan)/time * t_final;
figure;
subplot(2,1,1), plot(tspan(1:linear_time),X(1,1:linear_time), 'k'), title('a and linear a vs t')
hold on;
subplot(2,1,1), plot(tspan(1:linear_time),X(5,1:linear_time), 'r');
subplot(2,1,2), plot(tspan(1:linear_time),X(3,1:linear_time), 'k'), title('b and linear b vs t')
hold on; 
subplot(2,1,2), plot(tspan(1:linear_time),X(6,1:linear_time), 'r');
%%
%graphs of envelope and a/b
%vectors must be of same length, make sure tspan is the same
figure;
subplot(2,1,1), plot(tspan,abs(Y(1,:))*2, 'r'), title('mod A1 and a vs t');
hold on;
subplot(2,1,1), plot(tspan,X(1,:), 'k');
subplot(2,1,2), plot(tspan,abs(Y(2,:))*2, 'r'), title('mod A2 adn a2 vs t');
hold on;
subplot(2,1,2), plot(tspan,X(2,:), 'k');
%%



%find ways to print alpha and forcing values

