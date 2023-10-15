classdef ParametersClass
   

    methods(Static)
        %%
        %parameters to change
        
        %controls nonlinearity, set to 0 for linear
        function alpha = getAlpha()
        alpha = 1;
        end
        
        function gravity = getGravity()
        gravity = 1;
        end
        
        %non dimensionalized surface tension parameter
        function gamma = getGamma()
        gamma = 1;
        end
        
        %amplitude of forcing function
        function forcing = getForcing()
        forcing = .01;
        end
        
        %length of tank used, in cm
        function length = getLength()
        length = 12.3;
        end
        
        %frequency
        function mode = getMode()
        mode = 1;
        end

        %density
        function rho = getRho()
        rho = 1;
        end

        %viscosity
        function nu = getNu()
        nu = 0;
        end

        %%
        %dependent parameters, do not change
        
        function wavenumber = getWavenumber()
        j = ParametersClass.getMode();
        l = ParametersClass.getLength();
        wavenumber = j * (pi) / l;
        end
        
        function wavenumber2 = getWavenumber2()
        k = ParametersClass.getWavenumber();
        wavenumber2 = 2 * k;
        end
        
        function hyper_tang = getHyper_tang()
        gamma = ParametersClass.getGamma();
        k = ParametersClass.getWavenumber();
        gravity = ParametersClass.getGravity();
        hyper_tang = sqrt((3*gamma*(k^2)/gravity)/(1+(gamma*(k^2)/gravity)));
        end
        
        function hyper_tang2 = getHyper_tang2()
        T = ParametersClass.getHyper_tang();
        hyper_tang2 = 2*T/(1+T^2);
        end
        
        function height = getHeight()
        T = ParametersClass.getHyper_tang();
        k = ParametersClass.getWavenumber();
        height = atanh(T)/k;
        end
        
        function modified_g = getModified_g()
        gravity = ParametersClass.getGravity();
        gamma = ParametersClass.getGamma();
        k = ParametersClass.getWavenumber();
        modified_g = gravity + gamma * (k^2);
        end
        
        function modified_g2 = getModified_g2()
        gravity = ParametersClass.getGravity();
        gamma = ParametersClass.getGamma();
        k2 = ParametersClass.getWavenumber2();
        modified_g2 = gravity + gamma * (k2^2);
        end
        
        function frequency = getFrequency()
        g = ParametersClass.getModified_g();
        T = ParametersClass.getHyper_tang();
        k = ParametersClass.getWavenumber();
        frequency = sqrt(g*k*T);
        end

        function frequency2 = getFrequency2()
        g2 = ParametersClass.getModified_g2();
        T2 = ParametersClass.getHyper_tang2();
        k2 = ParametersClass.getWavenumber2();
        frequency2 = sqrt(g2*k2*T2);
        end

        function constant1 = getConstant1()
        l = ParametersClass.getLength();
        j = ParametersClass.getMode();
        k = ParametersClass.getWavenumber();
        k2 = ParametersClass.getWavenumber2();
        gamma = ParametersClass.getGamma();
        gravity = ParametersClass.getGravity();
        T = ParametersClass.getHyper_tang();
        T2 = ParametersClass.getHyper_tang2();
        h = ParametersClass.getHeight();
        g = ParametersClass.getModified_g();
        g2 = ParametersClass.getModified_g2();    
        freq = ParametersClass.getFrequency();
        alpha = ParametersClass.getAlpha();
        forcing = ParametersClass.getForcing();
        constant1 = (k^2/(2*freq))*(2*g2+g+g*(T^2)-2*g2*T*T2);
        end

        function constant2 = getConstant2()
        l = ParametersClass.getLength();
        j = ParametersClass.getMode();
        k = ParametersClass.getWavenumber();
        k2 = ParametersClass.getWavenumber2();
        gamma = ParametersClass.getGamma();
        gravity = ParametersClass.getGravity();
        T = ParametersClass.getHyper_tang();
        T2 = ParametersClass.getHyper_tang2();
        h = ParametersClass.getHeight();
        g = ParametersClass.getModified_g();
        g2 = ParametersClass.getModified_g2();    
        freq = ParametersClass.getFrequency();
        alpha = ParametersClass.getAlpha();
        forcing = ParametersClass.getForcing();
        constant2 = (k^2/(2*freq))*(g/g2)*(2*g2+g+g*(T^2)-2*g2*T*T2);
        end

        %Ignore this function
%         function Y = getSystemVar()
%         k = ParametersClass.getWavenumber();
%         k2 = ParametersClass.getWavenumber2();
%         T = ParametersClass.getHyper_tang();
%         T2 = ParametersClass.getHyper_tang2();
%         g = ParametersClass.getModified_g();
%         g2 = ParametersClass.getModified_g2();    
%         freq = ParametersClass.getFrequency();
%         alpha = ParametersClass.getAlpha();
%         forcing = ParametersClass.getForcing();
%         Y = [k, k2, T, T2, g, g2, freq, alpha, forcing];
%         end

        function timeVars = getTimeVars()
        dt = .005;
        t_final = 200;
        n = t_final/dt + 1;
        timeVars = [dt; t_final; n];
        end

    end
end