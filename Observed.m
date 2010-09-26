classdef Observed < handle
	
	properties
		dpath='data/obsdata2.dat';
		name;
		data;
		x,y,w;
	end
	
	
	methods(Static)
		
	end
	
	
	methods
		function ob = Observed(varargin)
			load_args
			dpath = arg('path',ob.dpath);
			ob.name = arg('name','halo000');
			
			ob.data = load(dpath);
			x = ob.data(:,1);
			y = ob.data(:,2);
			w = ob.data(:,3);
			
			ndx = x>-3;
			ob.x=x(ndx);
			ob.y=y(ndx);
			ob.w=w(ndx);
		end
		
		
		
		function plot(ob,varargin)
			load_args
			
			overlay = arg('overlay',false);
			dot_size = arg('dot_size',5);
			axl = arg('axis',[-3 0 -.5 1]);
			cmp = arg('colormap','winter');
			
			scatter_color_by_weight(ob.x,ob.y,log(ob.w),struct(...
				'overlay',overlay,...
			'title', ob.name, 'x','Fe/H', 'y','\alpha/Fe','colormap',cmp,'dot_size',dot_size));
			axis(axl)
		end

		
		
		
		
	end
	
end