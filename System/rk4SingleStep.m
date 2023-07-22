function xout = rk4SingleStep(f, dt, tk, xk)

f1 = f(tk,xk);
f2 = f(tk+dt/2, xk+(dt/2)*f1);
f3 = f(tk+dt/2, xk+(dt/2)*f2);
f4 = f(tk+dt, xk+dt*f3);

xout = xk + (dt/6)*(f1+2*f2+2*f3+f4);
end


