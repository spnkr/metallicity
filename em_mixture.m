function [w] = em_mixture(x,mu1,var1,mu2,var2)
	
	max_iters = 100;
	j=1;
	
	w_prime = 0.1;
	w = w_prime + 10;
	W = zeros(1,1);
	
	while abs(w-w_prime) > 0.000001 && j < max_iters
		w = w_prime;
		w_prime = mean(expected_w(x,w,mu1,var1,mu2,var2));
		
		W(j) = w;
		j=j+1;
	end
	
	if fig
		plot(W,'.-');
		title(strcat(['w=' num2str(w_prime) ' after ' num2str(j) ' iterations']));
		xlabel('Iteration');
		ylabel('w');
	end
	
	%expected value of w given the data and the initial parameter
	%if pass multiple it will return a matrix of individually processed
	%probabilities
	function p = expected_w(x,w,mu1,var1,mu2,var2)
		p_1 = (w.*normpdf(x,mu1,var1));
		p_mixture = w.*normpdf(x,mu1,var1) + (1-w).*normpdf(x,mu2,var2);
		p = p_1 ./ p_mixture;
	
		

		
		
		