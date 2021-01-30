function [xmin, fmin, iter] = Newton_ls(f,gf,Hf,x0,Stop_tol,GSS_tol,varargin)
    iter = 0;
    x = x0;
    while(1)
        iter = iter + 1;
        H = Hf(x);
        try(chol(H));
            
           % fprintf("tried %d",iter);
        catch
            mu=abs(min(eig(H)));
            H = H + (1.0001)*mu*eye(size(H));
            %fprintf("catched %d",iter);
        end        
        p = -inv(H)*gf(x);
        phi = @(a) f(x + a*p);
		derv_phi = @(a)(p'*gf(x+a*p));
		phi_prime_0 = derv_phi(0);
        if(iter==1)
            alpha1 = 1;
            f_now = f(x0);
        else
            alpha1_cand = -2*abs(df)/phi_prime_0;
            alpha1 = min([1,1.01*alpha1_cand]);
        end
		[alpha,success] =  bracketingphase(phi,derv_phi, 0.01, 0.9,alpha1,8,GSS_tol);
%         if success
%             disp('Hii')
%             p = -gf(x);
%         end
        x_new = x + alpha*p;
        if((norm(x_new - x,2)<=Stop_tol))
            break
        end
        if(iter > 100)
            break
        end
        f_new = f(x_new); 
		df = f_new - f_now;
        x = x_new;
		f_now = f_new;

    end
    xmin = x_new;
    fmin = f(xmin);
end