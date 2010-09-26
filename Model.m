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
		
		function mo = generate(stepsize,precision,varargin)
			load_args
			
			do_save = arg('do_save',true);
			save_to = arg('save_to','cache/models_01.mat');
			
			mo = Model(struct('step',stepsize,'precision',precision));
			
			if do_save
				mo.save(save_to);
			end
		end
		
		function mod = load(name)
			load(name);
			mod = mo;
		end
		
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
					mo.data{k} = models(models(:,4)==m & models(:,5)==n,[1 2 3 4 5]);
					mo.data{k}(:,[1 2]) = round(mo.data{k}(:,[1 2]).*mo.precision)./mo.precision;
					v = models(models(:,4)==m & models(:,5)==n,[6 7 9 8]);
					mo.props{k} = [min(v(:,1)) min(v(:,2)) min(v(:,3)) min(v(:,4))];
					k=k+1;
				end
			end
			
			mo.xranges = {};
			mo.yranges = {};
			
			for k=1:length(mo.data)
				mo.xranges{k} = min(mo.data{k}(:,1)):stepsize:max(mo.data{k}(:,1));
				mo.yranges{k} = min(mo.data{k}(:,2)):stepsize:max(mo.data{k}(:,2));
				
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
			show(['Saved to ' path])
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

		
		
		
		function plot(mo,varargin)
			load_args
			
			overlay = arg('overlay',false);
			ob = arg('observed',NaN);
			vx = cell2mat(varargin);
			vx.subplot = true;
			
			global im;
			fg = figure(im);
			clf(fg,'reset');
			im=im+1;
			
			if overlay
				hold on
			end
			
			k=1;
			for i=1:length(mo.cache)
				if ~overlay
					subplot(5,5,k);
				end
				vx.model = k;
				mo.plot_single(vx);
				k=k+1;
			end
			
			
			if overlay
				ob.plot(struct('overlay',true));
				hold off
				title('All models and observed data')
			end
			
		end
		
		function plot_single(mo,varargin)
			load_args
			
			k = arg('model',1);
			step_size = arg('step_size',0.01);
			precision = arg('precision',100);
			subplot = arg('subplot',false);
			overlay = arg('overlay',false);
			
			[X,Y] = meshgrid(min(mo.xranges{k}):step_size:max(mo.xranges{k}), ...
				min(mo.yranges{k}):step_size:max(mo.yranges{k}));
 			X = round(X.*precision)./precision;
 			Y = round(Y.*precision)./precision;
			
			Z = zeros(size(X,1),size(X,2));
			p = size(Z,1);
			q = size(Z,2);
						
			for i=1:p
				for j=1:q
					xs = X(i,j);
					ys = Y(i,j);
					
					a = mo.f_ab(k,xs,ys);
					
					if ~isfinite(a)
						a=0;
					end

					if numel(a)==0
						a=0;
					end

					Z(i,j) = a;
				end
			end


			if ~subplot && ~overlay
				fig
			end
			
			AlphaData = Z;
			AlphaData(AlphaData>0) = 1;
			colormap winter
			surf(X,Y,Z,'AlphaData',AlphaData,'AlphaDataMapping','scaled','FaceAlpha','flat','EdgeAlpha',0,...
				'CData',log(Z));
			flabel('Fe/H','\alpha/Fe','F_{a,b}')
			view(0,90);
			axis([-3 0 -.5 1]);

			
			title(strcat(['' num2str(mo.props{k}(1)) '-' num2str(mo.props{k}(2)) ...
				' Gyr; M=' num2str(mo.props{k}(3))  '-' num2str(mo.props{k}(4)) ]))
			
		end
		
		
		function out = nonzeros(mo,k)
			out = mo.data{k}(mo.data{k}(:,3)>0,1:3);
		end
		
		
		function test_cache(mo,varargin)
			load_args
			
			sigma = 0;
			stepsize = arg('step', mo.stepsize);
			single_model_no = arg('model_no',NaN);
			
			if ~isfinite(single_model_no)
				mns = 1:size(mo.data,2);
			else
				mns = [single_model_no];
			end
			
			for kp=1:length(mns)
				model_no = mns(kp);
				
				[X,Y] = meshgrid(min(mo.data{model_no}(:,1)):stepsize:max(mo.data{model_no}(:,1)), ...
					min(mo.data{model_no}(:,2)):stepsize:max(mo.data{model_no}(:,2)));
				X = round(X.*mo.precision)./mo.precision;
				Y = round(Y.*mo.precision)./mo.precision;

				Z = zeros(size(X,1),size(X,2));
				p = size(Z,1);
				q = size(Z,2);

				for i=1:p
					for j=1:q
						xs = X(i,j);
						ys = Y(i,j);

						a = mo.f_ab(model_no,xs,ys);
						b = mo.f_ab_live(model_no,xs,ys);
						v = norm(a-b);
						sigma = sigma + v;

						if v ~= 0
							err=[a b xs ys]
						end
					end
				end
				
				
				disp(strcat(['diff of ' num2str(sigma) ' for model no ' num2str(model_no)]))
			end
		end
		
		
	end
	
end