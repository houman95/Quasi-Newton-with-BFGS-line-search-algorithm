function [xmin, fmin, iter] = BFGS(f,gf,x0,Stop_tol,GSS_tol,varargin)
    iter = 0;
    x = x0;
    g = gf(x);
	C0 = eye(length(g));
	p = -C0*g;
    phi = @(a) f(x + a*p);
    derv_phi = @(a)(p'*gf(x+a*p));
    alpha =  bracketingphase(phi,derv_phi, 1e-4, 0.9,1,10,GSS_tol);
    x_new = x + alpha*p;
    g_new = gf(x_new);
    delta = x_new - x;
    gamma = g_new - g;
    C = (gamma'*delta/(gamma'*gamma))*C0;
	x= x_new;
	g = g_new;
    while(1)
        p = -C*g_new;
        phi = @(a) f(x + a*p);
		derv_phi = @(a)(p'*gf(x+a*p));
		alpha =  bracketingphase(phi,derv_phi, 1e-4, 0.9,1,10,GSS_tol);
        x_new = x + alpha*p;
		g_new = gf(x_new);
		delta = x_new - x;
		gamma = g_new - g;
		a = 1/(delta'*gamma);
		C = (C0 - a*delta*gamma')*C*(C0 - a*gamma*delta') + a*delta*delta';
        if(norm(x_new - x,2)<=Stop_tol)
            break
        else
            %disp('error is:');disp(norm(x_new - x,2));
        end
        if(iter < 1000)
            %disp(iter);
        end
        x= x_new;
		g = g_new;
        iter = iter+1;
    end
    xmin = x_new;
    fmin = f(xmin);
end