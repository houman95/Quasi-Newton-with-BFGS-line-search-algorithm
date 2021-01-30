function [alphaStar,success] = bracketingphase(f, gf, c1, c2, alpha_1, alpha_max, epsilon)
    success=0;
    alpha_prev = 0;
    phi_prev = f(0);
    alpha = alpha_1;
	phiPrime_0 = gf(0);	
	phiAlpaha_i = f(alpha);
	phi_0 = f(0);
	while(1)
		phiAlpaha_i = f(alpha);
		Phiprimeai = gf(alpha);
		if(alpha==alpha_max)
            disp('Hellp Bracketing')
			break;
		end
		if((phiAlpaha_i>phi_0+c1*alpha*phiPrime_0)||(phiAlpaha_i>=phi_prev))
			a = alpha_prev;
			b = alpha;
			[alphaStar,success] = sectioning(f, gf, c1, c2, a, b, epsilon);
			break;
		end	
		if(abs(Phiprimeai)<=-c2*phiPrime_0)
            alphaStar = alpha;
			break;
		elseif(Phiprimeai>=0)
			a = alpha;
			b = alpha_prev;
			[alphaStar,success] = sectioning(f, gf, c1, c2, a, b, epsilon);
			break;		
		end
		phi_prev = phiAlpaha_i;
		temp = alpha;
		if(2*alpha-alpha_prev<alpha_max)			
			alpha = mean([2*alpha-alpha_prev, min(alpha_max,alpha + 10*(alpha-alpha_prev))]);			
		else
			alpha = alpha_max;			
		end
		alpha_prev = temp;
	end
end