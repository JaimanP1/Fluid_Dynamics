function dx = EnvelopeSystem_v3(t,x)
%EnvelopeSystem

C1 = ParametersClass.getConstant1();
C2 = ParametersClass.getConstant2();

%%forcing stuff
forcing = ParametersClass.getForcing();
freq = ParametersClass.getFrequency();
freq2 = ParametersClass.getFrequency2();
g = ParametersClass.getModified_g();
g2 = ParametersClass.getModified_g2();
nu = ParametersClass.getNu();
F1 = freq*forcing/(4*g);
F2 = freq2*forcing/(4*g2);
%F2 should be multiplied by focing2 instead of just forcing, but here I am
%assuming these coefficients are the same

%using equations 2.11-2.12, find way to include detuning parameters
dx = [
i*C1*x(2)*conj(x(1)) + (1i*F1*conj(x(1))) - (1/2)*nu*x(1);
i*C2*(x(1))^2 + (1i*F2*conj(x(2))) - (1/2)*nu*x(2);
%second equation is technically using nu2, but since they are coefficients,
%I am working off the assumption that they are the same
];
end