function dx = system4_v6(t,x)
%system of four nonlinear differential equations


k = ParametersClass.getWavenumber();
k2 = ParametersClass.getWavenumber2();
T = ParametersClass.getHyper_tang();
T2 = ParametersClass.getHyper_tang2();
g = ParametersClass.getModified_g();
g2 = ParametersClass.getModified_g2();    
freq = ParametersClass.getFrequency();
alpha = ParametersClass.getAlpha();
forcing = ParametersClass.getForcing();
gamma = ParametersClass.getGamma();
rho = ParametersClass.getRho();
%using these function calls makes the program ~4x slower for t>1000


dx = [
k*T*x(3) + ((2*(k^2)*(1-T*T2)*x(1)*x(4) - (k^2)*(1+(T^2))*x(2)*x(3)) - ((k^3)*T*(3-(2*T*T2))*x(3)*(x(1)^2)))*alpha;
(k2*T2*x(4)) + (2*(k^2)*(1-T*T2)*x(1)*x(3) + ((k^3)*T*(3-(2*T*T2))*x(1)*(x(3)^2)) + (2/3)*(gamma/rho)*(k^4)*(x(1)^3))*alpha;
-(forcing*cos(2*freq*t)+g)*x(1) - (2*(k^2)*(1-T*T2)*x(3)*x(4))*alpha; 
-g2*x(2) + (.5*(k^2)*(1+(T^2))*(x(3)^2))*alpha;
%linear a
k*T*x(6);
%linear b
-(forcing*cos(2*freq*t)+g)*x(5);
];

end