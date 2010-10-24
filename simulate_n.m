function [x,f] = simulate_n(pi,mu,sigma2,n)
	m = length(pi);
	if sum(pi) ~= 1
		error('sum pi must = 1')
	end
	pic = cumsum(pi);
	sigma = sqrt(sigma2);
	
	x = zeros(n,1);
	f = zeros(n,m);
	
	for i=1:n
		r0=rand;
		ndx = sum(pic<=r0)+1;
		x(i) = mu(ndx) + sigma(ndx).*randn(1,1);
	end
	
	for i=1:m
		f(:,i) = normpdf(x,mu(i),sigma(i));
	end
	
	
	
	
end