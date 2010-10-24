classdef Mixture < handle
	
	properties
		filename;
		
		x;
		f;
		
		true_loglike;
		est_loglike;
		true_pi;
		est_pi;
		est_pi_all;
		est_loglike_all;
		
		em_optimal_iters;
		em_run_iters;
		em_run_seconds;
		em_init_p;
		
		models;
		num_models=0;
		bin_step=0.1;
		xrange;
		yrange;
	end
	
	
	methods(Static)
		function mix = load(filename)
			load(strcat(['cache/' filename '.mat']));
			mix=mi;
		end
		
		function tp = get_true_pi(halonum)
			tp=NaN;
			if halonum==3
				tp = [	14.467074  ;
					7.4354783  ;
					20.991296  ;
					9.2340355  ;
					33.655754  ;
					8.0191336  ;
					2.4001462  ;
					2.5734037  ;
				   0.11883542  ;
				  0.076990272  ;
				   0.35844400  ;
				   0.25313549  ;
				   0.22529346  ;
				  0.048890239  ;
				  0.024226458  ;
				  0.046833900  ];
				tp = tp./100;
			end
		end
	end
	
	
	methods
		%loading
		function mi = Mixture(varargin)
			load_args
			
			mi.filename = arg('save_as',NaN);
			mi.em_optimal_iters = arg('em_optimal_iters',NaN);
			obs_path = arg('obs_path',NaN);
			model_path = arg('model_path',NaN);
			mi.f = arg('f',NaN);
			mi.true_pi = arg('true_pi',NaN);
			mi.bin_step = arg('bin_step',mi.bin_step);
			graph = arg('graph',false);
			if graph
				fig
			end
			
			x = load(obs_path);
			x = x(x(:,1)>-3,[1 2]);
			mi.x=x;
			
 			if ~isfinite(mi.f)
 				models = load(model_path);
				
				xbin = models(:,4);
				ybin = models(:,5);
				x_bin_num = max(xbin)+1;
				y_bin_num = max(ybin)+1;
				
				for i=0:x_bin_num-1
					for j=0:y_bin_num-1
						ndx = models(:,4)==i & models(:,5)==j;
						if sum(models(ndx,3))>0
							mi.num_models = mi.num_models + 1;
						end
					end
				end
				
				xrange = min(models(:,1)):mi.bin_step:max(models(:,1));
				yrange = min(models(:,2)):mi.bin_step:max(models(:,2));
				num_models = mi.num_models;
				
				
				
				xn=length(xrange);
				yn=length(yrange);
				k=0;
				sq = zeros(xn,yn);
				f = zeros(num_models,xrange,yrange);
				indic=zeros(xn*yn,1);
				m = 0;
				for i=0:x_bin_num-1
					for j=0:y_bin_num-1
						m = m + 1;
						idx = ((models(:,4)==j) & (models(:,5)==i));
						distx = models(idx==1,:);
						n = sum(distx(:,3)~=0);
						if (n==0) continue; end
						l = 0; k = k + 1;
						indic(m) = 1;
						for ix = 1:xn,
							for iy = 1:yn,
								l = l + 1;
								sq(ix,iy) = distx(l,3);
								f(k,ix,iy) = distx(l,3);
							end
						end
						if graph
							subplot(5,5,m); mesh(xrange,yrange,sq'); hold on;
						end
					end
				end
				
				if graph
					hold off
				end
				
				f = 100*f;
				fval = EvalDens(x(:,1),x(:,2),f,xrange,yrange,num_models,xn,yn,mi.bin_step/2);
				id = (sum(fval')~=0);
				fval = fval(id,:);
				mi.f = fval;
				
				
				mi.models = models;
				mi.xrange=xrange;
				mi.yrange=yrange;
			end
			
			if isfinite(mi.true_pi)
				mi.true_loglike = mi.complete_log_like(mi.f,mi.true_pi,length(mi.x),mi.num_models);
			end
			
			
			
		end
		
		function save(mi)
			if isfinite(mi.filename)
				save(strcat(['cache/' mi.filename '.mat']),'mi');
				disp(strcat(['Saved results to cache/' mi.filename '.mat']));
			end
		end
		
		function save_as(mi,filename)
			mi.filename = filename;
			mi.save();
		end
		
		%\llp &= \sumn \log \Big( \summ \pi_j \fab(x_i,y_i)  \Big)
		function l = complete_log_like(mi,f,p,n,m)
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
		
		function [p,ll,P,LL] = em(mi,varargin)
			load_args
			car = cell2mat(varargin);
			
			if ~isfinite(arg('max_iters',NaN)) && ~isfinite(arg('max_seconds',NaN)) && isfinite(mi.em_optimal_iters)
				car.max_iters = mi.em_optimal_iters;
			end
			
			car.num_models = mi.num_models;
			car.p_actual = mi.true_pi;
			car.true_loglike = mi.true_loglike;
			
			
			[p,ll,P,LL,init_p,counter,tmr] = em(mi.x,mi.f,car);
			
			mi.est_pi = p;
			mi.est_loglike = ll;
			mi.est_pi_all = P;
			mi.est_loglike_all = LL;
			mi.em_run_seconds = tmr;
			mi.em_run_iters = counter;
			mi.em_init_p = init_p;
			
			mi.save
		end
		
	end
	
end