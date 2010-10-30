function [S,V,correl,stdev] = bootstrap_covar(mi, pistar, varargin)
	load_args
	
	pistar=pistar';
	B = size(pistar,1);
	m = size(pistar,2);
	
	pimean = sum(pistar)./B;
	
	if 1==11
	S = zeros(m,m);
	for b=1:B
		for i=1:m
			for j=1:m
				S(i,j) = S(i,j) + (pistar(b,i)-pimean(i))*(pistar(b,j)-pimean(j));
			end
		end
	end
	S=S./(B-1);
	end
	
	pimean = repmat(pimean,B,1);
	S=(pistar'-pimean')*((pistar'-pimean')');
	S = S./(B-1);
	
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