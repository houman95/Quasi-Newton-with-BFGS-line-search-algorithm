function alpha_j = choose(a,b)
    tau2 = 0.1;
    tau3 = 0.5;
	alpha_j = (a+tau2*(b-a)+b-tau3*(b-a))/2;
end