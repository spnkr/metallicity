function dotplot_light(x,y)
	
	mc = [0.7 0.8 .9];
	oc = [0.8 0.9 1];
	
	if 1==1
		%green
		mc = [0.6 0.7 .2];
		oc = [0.7 0.8 .3];
	end
	if 1==11
		%purple
		mc = [0.6 0.2 .8];
		oc = [0.7 0.3 .9];
	end
	if 1==11
		%darkred
		mc = [0.8 0.1 .0];
		oc = [0.9 0.1 .0];
	end
	if 1==11
		%red
		mc = [0.8 0.4 .4];
		oc = [0.9 0.4 .4];
	end
	
	plot(x, y, 'Color', mc, 'Marker', '.', 'LineStyle', 'none')
	plot(x, y, 'Color', oc, 'Marker', 'o', 'LineStyle', 'none')
	