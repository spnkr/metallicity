classdef Model < handle
	
	properties
		data;
		props;
		cache;
		precision;
		stepsize;
		xranges;
		yranges;
	end
	
	
	methods(Static)
		
	end
	
	
	methods
		function mo = Model(varargin)
			load_args
			
			stepsize = arg('step',0.1);
			precision = arg('precision',10);
			mo.precision = precision;
			mo.stepsize = stepsize;
			
			disp(strcat(['Loading with step of ' num2str(stepsize) ' and prec of ' num2str(precision) '...']))
			
			models = load('data/modeldata2.dat');
			xbin = models(:,4);
			ybin = models(:,5);
			
			
			
			mo.data = {};
			mo.props = {};
			mo.cache = {};
			k=1;
			for n=min(ybin):max(ybin)
				for m=min(xbin):max(xbin)
					mo.data{k} = models(models(:,4)==m & models(:,5)==n,[1 2 3]);
					v = models(models(:,4)==m & models(:,5)==n,[6 7 9 8]);
					mo.props{k} = [min(v(:,1)) min(v(:,2)) min(v(:,3)) min(v(:,4))];
					k=k+1;
				end
			end
			
			mo.xranges = {};
			mo.yranges = {};
			
			for k=1:length(mo.data)
				mo.xranges{k} = min(mo.data{k}(:,1)):stepsize:max(mo.data{k}(:,1));
				mo.yranges{k} = min(mo.data{k}(:,1)):stepsize:max(mo.data{k}(:,1));
				
				mo.xranges{k} = round(mo.xranges{k}.*precision)./precision;
				mo.yranges{k} = round(mo.yranges{k}.*precision)./precision;
				
				[X,Y] = meshgrid(mo.xranges{k}, mo.yranges{k});
 				X = round(X.*precision)./precision;
 				Y = round(Y.*precision)./precision;
				
				p = size(X,1);
				q = size(X,2);
				Z = zeros(p*q,3);
				
				d=1;
				for i=1:p
					for j=1:q
						xs = X(i,j);
						ys = Y(i,j);

						a = mo.f_ab_live(k,xs,ys);

						Z(d,:) = [a xs ys];
						d=d+1;
					end
				end
				
				mo.cache{k} = Z;

				k=k+1;
			end
			
			disp('Loaded.')
		end
		
		
		function save(mo,path)
			save(path,'mo');
			disp('Saved');
		end
	
		function out = f_ab(mo,mndx,x,y)
			data = mo.cache{mndx};
			xr = mo.xranges{mndx};
			yr = mo.yranges{mndx};
			
			x = round(x.*mo.precision)./mo.precision;
			y = round(y.*mo.precision)./mo.precision;
			
			xndx = find(xr==x);
			yndx = find(yr==y);
			
			out = data((yndx-1)*length(xr)+xndx,1);
		end
		
		
		
		
		
		function out = f_ab_live(mo,mndx,x,y)
			data = mo.data{mndx};
			
			xa = data(:,1);
			ya = data(:,2);
			
			x_bound = max(xa(xa<=x));
			
			if numel(x_bound)==0
				out=0;
				return;
			end
			
			y_bound = max(ya(ya<=y));
			
			if numel(y_bound)==0
				out=0;
				return;
			end
			
			out = data(xa==x_bound & ya==y_bound,3);
		end

		
		
		
		
	end
	
end