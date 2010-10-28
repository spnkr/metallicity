function [I,S,V,correl,stdev] = fisher_info(pi,f)
	m = length(pi);
	g = m-1;
	I = zeros(g,g);

	denom = (f*pi).^2;
	
	for k=1:g
		for r=1:g
			kk = ((f(:,k)-f(:,m)).*(f(:,r)-f(:,m)))./denom;
			I(k,r) = sum(kk(isfinite(kk)));
		end
	end
	
	
	
	
	
	
	
	
	S = inv(I);
	
	for i=1:m
		S(m,i) = 0;
		S(i,m) = 0;
		
		if m==i
			for j=1:g
				for k=1:g
					S(m,m) = S(m,m) + S(j,k);
				end
			end
		else
			for j=1:g
				S(m,i) = S(m,i) + S(j,i);
				S(i,m) = S(i,m) + S(i,j);
			end
			S(m,i) = -S(m,i);
			S(i,m) = -S(i,m);
		end
	end
	
	V = zeros(m,1);
	for i=1:m
		V(i) = abs(S(i,i));
	end
	
	stdev = sqrt(V);
	
	correl = S;
	for i=1:m
		for j=1:m
			correl(i,j) = correl(i,j)/(stdev(i)*stdev(j));
		end
	end
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
%	n = size(f,1);	
% 	denom = zeros(n,1);
% 	for i=1:n
% 		for j=1:m
% 			denom(i) = denom(i) + (pi(j).*f(i,j));
% 		end
% 	end
% 	denom = denom.^2;

% 	for i=1:n
% 		for k=1:g
% 			for r=1:g
% 				%at i, d2l/dkdr
% 				numer = (f(i,k)-f(i,m))*(f(i,r)-f(i,m));
% 				I(k,r) = I(k,r) + numer/denom(i);
% 			end
% 		end
% 	end
% 	I
	
	