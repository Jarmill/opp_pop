    function [E_out] = energy_blackbox(x_in)

        if all(x_in==0)
            E_out = Inf*[1; 1];
            return
        end

        af = x_in(1:end-1);
        I0 = x_in(end);

        %black-box the energy
        ah = [0; af; 2*pi];
        da = diff(ah);
        
        %compute the current
        % I0 = 0;
        % I_step = uf.*da;
        % I0_val = cumsum([0; I_step]);
        
        I_s = I0*ones(size(ah));
        
        alpha_val = linspace(0, 2*pi, N);
        I_val = zeros(1, N);
        I_val(1) = I0;
        
        Na = (length(I_s)-1);
        
        %compute the current
        for i = 1:Na
            ucurr = uf(i);
            Iprev = I_s(i);
            dt = da(i);
            I_s(i+1) = ucurr*(1-exp(-kappa*dt))/kappa + Iprev*exp(-kappa*dt);
        end
        
        %compute the energy
        % E_s = zeros(Na, 1)*a(1);
        energy = 0;
        for i = 1:Na
            ucurr = uf(i);
            Iprev = I_s(i);
            dt = da(i);
        
            E0 = (ucurr-kappa*Iprev)*(3*ucurr + kappa*Iprev);
            E_denom = 2*kappa^3;
            E_num_1 = 2*ucurr^2*kappa*dt;
            E_num_2 = exp(-2*kappa*dt)*(ucurr - kappa*Iprev) * ...
                (kappa*Iprev + ucurr*(4*exp(kappa*dt)-1));
            
            energy = energy +  (E_num_1 + E_num_2 - E0)/E_denom;
        end
        
        E_out = [energy; I_end];
    end