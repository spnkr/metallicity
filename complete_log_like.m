function l = complete_log_like(f,p,n,m)
	l = 0;
	for i=1:n
		l0 = 0;
		for j=1:m
			l0 = l0 + p(j)*f(i,j);
		end
		lgl0 = log(l0);
		if ~isfinite(lgl0)
			lgl0=0;
		end
		l = l + lgl0;
	end
end