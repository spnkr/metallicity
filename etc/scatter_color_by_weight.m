function scatter_color_by_weight(x,y,w,varargin)
	global im;
	load_args
	ttl = arg('title','');
	xax = arg('x','');
	yax = arg('y','');
	vary_size = arg('vary_size',false);
	dot_size = arg('dot_size',1);
	cmapname = arg('colormap','winter');
	change_back = arg('change_back',NaN);
	
	
	if im > -1
		fg = figure(im);
		im=im+1;
		
		colormap(cmapname);
		curcolmap = colormap;
		curmapsize = size(curcolmap,1);
		minz = (min(w)); maxz = (max(w));
		mapz = round(1 + ((w) - maxz) ./ (minz-maxz) .* (curmapsize-1));
		mapz(~isfinite(mapz))=1;
		
		if vary_size
			mapz_inv = round(1 + ((w) - minz) ./ (maxz-minz) .* (curmapsize-1));
		else
			mapz_inv = dot_size;
		end
		
		curcolmap = colormap;
		curmapsize = size(curcolmap,1);
		
		scatter(x,y,mapz_inv,curcolmap(mapz,:), 'filled');
		
		flabel(xax, yax, ttl);
		if isfinite(change_back)
			whitebg(fg, change_back);
		end
	end
	