function [x] = bootstrap_generate_m_n(mi, B)
	doplot=false;

	data = mi.x;
	ds = size(data,1);
	data2=[];
	for i=1:B
		ndx = ceil(rand.*ds);
		x(i,1:2) = data(ndx,1:2);
	end
	
	x = x(x(:,1)>-3,[1 2]);
	
	global im;
	if doplot && im>0
		fig
		scatter(x(:,1),x(:,2),1,'k','filled');
		flabel('Fe/H','\alpha/Fe',[num2str(n) ' m-n sampled data points']);
	end

end