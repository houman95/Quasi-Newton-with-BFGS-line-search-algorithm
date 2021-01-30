function [x_min, f_min, iter] = SD_ls(f, gf, x0, Stop_tol, GSS_tol,varargin)
    iter = 0;
    x = x0
    while(1)
        iter = iter+1;
        P = -gf(x);
        phi = @(a) f(x + a*P);
		derv_phi = @(a)(P'*gf(x+a*P));
		phi_prime_0 = derv_phi(0);
        if(iter==1)
            alpha1 = 1;
            f_now = f(x0);
        else
            alpha1_cand = -2*abs(df)/phi_prime_0;
            alpha1 = min([1,1.01*alpha1_cand]);
        end
        [alpha,success] =  bracketingphase(phi,derv_phi, 0.01, 0.9,alpha1,8,GSS_tol);
        x_new = x + alpha*P;
        if(success || (norm(x_new - x,2)<= Stop_tol))
            break
        end
		f_new = f(x_new); 
		df = f_new - f_now;
        x = x_new;
		f_now = f_new;
		if(iter >500)
		break;
		end
    end
    x_min = x_new;
    f_min = f_new;    
end