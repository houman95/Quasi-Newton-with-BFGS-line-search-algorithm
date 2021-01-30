function [alpha flag] = sectioning(f, gf, c1, c2, aj, bj, epsilon)
    phi_0 = f(0);
	g0 = gf(0);
    flag = 0;
    while(1)
	    Phi_aj = f(aj);
        PhiPrime = gf(aj);
        alpha = choose(aj,bj);
        phi_i = f(alpha);        
        if ((aj-alpha) * PhiPrime) <= epsilon  %exit condition
            %alpha = 0.01;
            %disp('Sectioning Help')
            flag = 1;
            break;
        end
        if (( phi_i > phi_0 + c1 * alpha * g0 ) || ( phi_i >= Phi_aj ))
            bj = alpha;
        else
            g_alpha = gf(alpha);
			if abs(g_alpha) <= -c2*(g0)
               break;
            end            
            if (bj-aj)* g_alpha >=0
                bj = aj;
            end
            aj = alpha;
        end
    end
end