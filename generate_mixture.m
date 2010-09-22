function [x] = generate_mixture(us1,us2,m,p)
	x = zeros(m,1);

	for i=1:m
		r=rand;
		if r<p
			x(i) = normrnd(us1(1),us1(2));
		else
			x(i) = normrnd(us2(1),us2(2));
		end
	end
