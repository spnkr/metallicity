classdef Mixture < handle
	
	properties
		filename;
		
		x;
		f;
		
		loglike_true;
		loglike_est;
		pi_true;
		pi_est;
		pi_est_all;
		loglike_est_all;
		
		em_run_iters;
		em_run_seconds;
		em_init_p;
		
		models;
		num_models=0;
		bin_step=0.1;
		xrange;
		yrange;
		
		pval_lrt;
		
		info;
		covar;
		variance;
		stdev;
	end
	
	
	methods(Static)
		function mix = load(filename)
			load(strcat(['cache/' filename '.mat']));
			mix=mi;
		end
		
		function tp = get_pi_true(halonum)
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
			obs_path = arg('obs_path',NaN);
			model_path = arg('model_path',NaN);
			mi.f = arg('f',NaN);
			mi.x = arg('x',NaN);
			mi.pi_true = arg('pi_true',NaN);
			mi.bin_step = arg('bin_step',mi.bin_step);
			graph = arg('graph',false);
			if graph
				fig
			end
			
			if ~isfinite(mi.x)
				x = load(obs_path);
				x = x(x(:,1)>-3,[1 2]);
				mi.x=x;
			end
			
 			if isfinite(mi.f)
				mi.num_models = size(mi.f,2);
			else
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
			
			if isfinite(mi.pi_true)
				mi.loglike_true = mi.complete_log_like(mi.f,mi.pi_true,length(mi.x),mi.num_models);
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
		
		function [p,ll,P,LL] = em(mi,varargin)
			load_args
			car = cell2mat(varargin);
			
			n=arg('n',NaN);
			if isfinite(n)
				mi.x = mi.x(1:n,:);
				mi.f = mi.f(1:n,:);
				if isfinite(mi.pi_true)
					mi.loglike_true = mi.complete_log_like(mi.f,mi.pi_true,length(mi.x),mi.num_models);
				end
			end
			
			car.num_models = mi.num_models;
			car.p_actual = mi.pi_true;
			car.loglike_true = mi.loglike_true;
			
			global im;
			[p,ll,P,LL,init_p,counter,tmr] = em(mi.x,mi.f,car);
			im=im+1;
			
			mi.pi_est = p;
			mi.loglike_est = ll;
			mi.pi_est_all = P;
			mi.loglike_est_all = LL;
			mi.em_run_seconds = tmr;
			mi.em_run_iters = counter;
			mi.em_init_p = init_p;
		end
		
		function update_stats(mi)
			if isfinite(mi.pi_true)
				mi.loglike_true = mi.complete_log_like(mi.f,mi.pi_true,length(mi.x),mi.num_models);
			end
			
			mi.get_pval_lrt();
			
			[I,S,V,stdev] = fisher_info(mi.pi_est,mi.f);
			mi.info = I;
			mi.covar = S;
			mi.variance = V;
			mi.stdev = stdev;
		end
		
		function text_summary(mi)
			if isfinite(mi.pi_true)
				pi_true_diff = 100.*[mi.pi_est mi.pi_true mi.pi_est-mi.pi_true]
			else
				pi = mi.pi_est
			end
		end
		
		
		
		%--plots
		function plot_covar(mi)
			plot_covar_dots(mi);
		end
		
		function plot_stdev(mi)
			fig
			w=(100.*(mi.pi_est./range(mi.pi_est))).^1;
			w=w./range(w);
			scatter(1:mi.num_models,mi.stdev,400.*w,'r','filled')
			flabel('j','\sigma_j','Information based standard deviation of \pi')
		end
		
		
		%--internal
		%\llp &= \sumn \log \Big( \summ \pi_j \fab(x_i,y_i)  \Big)
		function l = complete_log_like(mi,f,p,n,m)
			l = complete_log_like(f,p,n,m);
		end
		
		
		
		function [p,t] = get_pval_lrt(mi)
			t=-2*(mi.loglike_true - mi.loglike_est);
			df=mi.num_models-1;
			
			%chi2inv(.95,15)

			chance_less_than_or_equal_to = 100*chi2cdf(t,df);
			chance_more_extreme = 100-chance_less_than_or_equal_to;
			p = chance_more_extreme;
			mi.pval_lrt = p;
		end
	end
	
end