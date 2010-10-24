function [I,S,V,stdev] = fisher_info(pi,f)
	m = length(pi);
	g = m-1;
	
	I = zeros(g,g);
	
	denom = (f*pi).^2;
	
	for i=1:g
		for k=1:g
			ii = f(:,i)./(denom.*f(:,k));
			I(i,k) = sum(ii(isfinite(ii)));
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
		end
	end
	
	V = zeros(m,1);
	for i=1:m
		V(i) = abs(S(i,i));
	end
	
	stdev = sqrt(V);
	
	

	
	
	