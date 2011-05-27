classdef Mixture < handle
	
	properties
		filename;
		model_path;
		model_skip_ndx=[];
		
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
		models_sparse;
		num_models=0;
		bin_step=0.1;
		xrange;
		yrange;
		
		pval_lrt;
		zscores;
		
		info;
		covar;
		variance;
		stdev;
		correl;
	end
	
	
	methods(Static)
		function mix = load(varargin)
			filename = strcat(cell2mat(varargin));
			
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
			else
				tp = [1.679 21.38 20.53 3.806 1.002 8.522 21.56 12.29 1.841 0.151 1.662 1.536 1.582 0.945...
					0.053 0.282 0.347 0.542 0.111 0.003 0.045 0.117]';
				tp = tp./100;
			end
		end
	end
	
	
	methods
		function a = grid_size(mi)
			a=mi.bin_step;
		end
		%--loading
		function mi = Mixture(varargin)
			load_args
			
			with_headers = arg('with_headers',false);
			mi.filename = arg('save_as',NaN);
			obs_path = arg('obs_path',NaN);
			model_path = arg('model_path',NaN);
			mi.f = arg('f',NaN);
			x = arg('x',NaN);
			mi.pi_true = arg('pi_true',NaN);
			mi.bin_step = arg('bin_step',mi.bin_step);
			graph = arg('graph',false);
			model_obj = arg('model_obj',NaN);
			
			if graph
				fig
			end
			
			if ~isfinite(x)
				x = load(obs_path);
				x = x(x(:,1)>-3,[1 2]);
			end
			mi.x=x;
			
			
 			if isfinite(mi.f)
				mi.num_models = size(mi.f,2);
			else
				if ~isfinite(model_obj)
					if with_headers
						[h,models] = hdrload(model_path);
					else
						models = load(model_path);
					end
					mi.model_path = model_path;
				
					xbin = models(:,4);
					ybin = models(:,5);
					x_bin_num = max(xbin)+1;
					y_bin_num = max(ybin)+1;
					num_models=0;
					models_sparse=[];
				else
					models = model_obj;
					mi.model_path = model_path;

					xbin = models(:,4);
					ybin = models(:,5);
					x_bin_num = max(xbin)+1;
					y_bin_num = max(ybin)+1;
					num_models=0;
					models_sparse=[];
				end
				
 				
				
				
				kk=1;
				for i=0:x_bin_num-1
					for j=0:y_bin_num-1
						ndx = models(:,4)==i & models(:,5)==j;
						if sum(models(ndx,3)) == 0
							mi.model_skip_ndx(kk) = 0;
						else
							mi.model_skip_ndx(kk) = 1;
						end
						kk=kk+1;
					end
				end
				
				for i=0:x_bin_num-1
					for j=0:y_bin_num-1
						ndx = models(:,4)==i & models(:,5)==j;
						if sum(models(ndx,3))>0
							num_models = num_models + 1;
							nd = models(ndx,1:3);
% 							nd = nd(nd(:,3)>0,:);
							a = size(models_sparse,1);
							b = a + size(nd,1);
							models_sparse(a+1:b,1:4) = [nd num_models.*ones(size(nd,1),1)];
						end
					end
				end
				
				mi.num_models = num_models;
				mi.models_sparse = models_sparse;
				
				xrange = min(models(:,1)):mi.bin_step:max(models(:,1));
				yrange = min(models(:,2)):mi.bin_step:max(models(:,2));
				num_models = mi.num_models;
				
				
				
				xn=length(xrange);
				yn=length(yrange);
				k=0;
				sq = zeros(xn,yn);
				warning off
				f = zeros(num_models,xrange,yrange);
				warning on
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
								try
									sq(ix,iy) = distx(l,3);
									f(k,ix,iy) = distx(l,3);
								catch
									sq(ix,iy) = 0;
									f(k,ix,iy) = 0;
								end
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
				
				if sum(id) ~= size(fval,1)
					fval = fval(id,:);
					mi.x = mi.x(id,:);
					%warning('found a 0 fval!');
				end
				
				
				mi.f = fval;
				
				
				mi.models = models;
				mi.xrange=xrange;
				mi.yrange=yrange;
			end
			
			if isfinite(mi.pi_true)
				try
					mi.loglike_true = mi.complete_log_like(mi.f,mi.pi_true,length(mi.x),mi.num_models);
				catch
					display('failed loglike');
				end
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
		
		
		
		
		
		
		%--processing
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
			if im>0
				im=im+1;
			end
			
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
			
			mi.pval_lrt = mi.get_pval_lrt();
			
			[I,S,V,correl,stdev] = fisher_info(mi.pi_est,mi.f);
			mi.info = I;
			mi.covar = S;
			mi.variance = V;
			mi.stdev = stdev;
			mi.correl = correl;
			
			mi.zscores = mi.get_zscores();
		end
		
		
		
		
		
		
		%--summaries
		function text_summary(mi)
			if isfinite(mi.pi_true)
				'pi_est; pi_true; diff'
				pi_summary = [rounder(100.*mi.pi_est,10000) rounder(100.*mi.pi_true,10000) rounder(100.*abs(mi.pi_true-mi.pi_est),10000)]
			else
				pi_summary = mi.pi_est
			end
		end
		
		
		
		
		
		
		%--plots
		function plot_history(mi)
			plot_formation_history(mi);
		end
		
		function plot_correl(mi,size_is_correl,color_neg)
			plot_correl_dots(mi,mi.correl,size_is_correl,color_neg);
		end
		
		function plot_stdev(mi)
			
			global im;
			aa=figure(im);
			im=im+1;
			clf(aa)
			clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
						'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];

			hold on
			w = 10000.*(mi.pi_est);
			for i=1:mi.num_models
				wi=max(w(i),10);
				scatter(i,mi.stdev(i),wi,clrs(i),'filled')
			end
			flabel('\pi_j','\sigma_j','Information based standard deviation of \pi')
			hold off
		end
		
		function plot_zscores(mi)
			global im;
			aa=figure(im);
			im=im+1;
			clf(aa)
			clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
						'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];

			hold on
			plot([0 mi.num_models], [1.96 1.96], 'k:');
			plot([0 mi.num_models], [-1.96 -1.96], 'k:');
			plot([0 mi.num_models], 2.*[1.96 1.96], 'k--');
			plot([0 mi.num_models], 2.*[-1.96 -1.96], 'k--');
% 			plot([0 mi.num_models], 3.*[1.96 1.96], 'k-');
% 			plot([0 mi.num_models], 3.*[-1.96 -1.96], 'k-');
			
			w = 10000.*(mi.pi_est);
			for i=1:mi.num_models
				wi=max(w(i),10);
				scatter(i,mi.zscores(i),wi,clrs(i),'filled')
			end
			
			flabel('\pi_h','Z-score','Z-scores for \pi, and \pm 1 and 2 \sigma');
			hold off
		end
		
		function plot_info_error_bars(mi,nstdevs)
			global im;
			fg=figure(im);im=im+1;
			clf(fg)

			clrss = {'c', 'm', 'k', 'g'};
			clrg = [0.7 0.4];
			k=1;
			hold on
			for i=-nstdevs:-1
				errorbar(1:mi.num_models,mi.pi_est,-i.*mi.stdev,strcat([clrss{-i} '.']),'LineWidth',2,'Color',clrg(k).*ones(1,3));
				k=k+1;
			end
			
			h1=plot(mi.pi_true,'ko','MarkerSize',15,'LineWidth',2);
			h2=plot(mi.pi_est,'gx','MarkerSize',20,'LineWidth',2,'Color',[10 198 28]./255);
			
			legend([h1 h2],'True \pi', 'Est \pi')
			hold off
			flabel('j','\pi',['Information based error bars \pm 1 & ' num2str(nstdevs) '\sigma - ' mi.filename])
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
		end
		
		function z = get_zscores(mi)
			if ~isfinite(mi.pi_true)
				z=[];
			else
				z = (mi.pi_est-mi.pi_true)./mi.stdev;
			end
		end
	end
	
end