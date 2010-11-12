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
			else
				tp = [1.679 21.38 20.53 3.806 1.002 8.522 21.56 12.29 1.841 0.151 1.662 1.536 1.582 0.945...
					0.053 0.282 0.347 0.542 0.111 0.003 0.045 0.117]';
				tp = tp./100;
			end
		end
	end
	
	
	methods
		%--loading
		function mi = Mixture(varargin)
			load_args
			
			mi.filename = arg('save_as',NaN);
			obs_path = arg('obs_path',NaN);
			model_path = arg('model_path',NaN);
			mi.f = arg('f',NaN);
			x = arg('x',NaN);
			mi.pi_true = arg('pi_true',NaN);
			mi.bin_step = arg('bin_step',mi.bin_step);
			graph = arg('graph',false);
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
 				models = load(model_path);
				mi.model_path = model_path;
				
				xbin = models(:,4);
				ybin = models(:,5);
				x_bin_num = max(xbin)+1;
				y_bin_num = max(ybin)+1;
				num_models=0;
				models_sparse=[];
				
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
				
				if sum(id) ~= size(fval,1)
					fval = fval(id,:);
					mi.x = mi.x(id,:);
					sepr
					sepr
					sepr
					sepr
					sepr
					sepr
					sepr(mi.filename);
					warning('found a 0 fval!');
					sepr
					sepr
					sepr
					sepr
					sepr
					sepr
				end
				
				
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
		
		function [S,V,correl,stdev] = bootstrap_covariance(mi,pistar)
			[S,V,correl,stdev] = bootstrap_covar(mi,pistar);
		end
		
		
		
		
		%--summaries
		function text_summary(mi)
			if isfinite(mi.pi_true)
				'pi_est; pi_true; diff; zscore'
				pi_summary = 100.*[round(10000.*mi.pi_est)./10000 round(10000.*mi.pi_true)./10000 ...
					round(10000.*mi.pi_est-mi.pi_true)./10000 round(1000.*mi.zscores)./(1000*100)]
			else
				pi_summary = mi.pi_est
			end
		end
		
		
		
		
		
		
		%--plots
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
		
		function [conf,conf2] = bootstrap_plot_bars_both(mi,pistar)
			n=size(pistar,2);
			m=size(pistar,1);

			T = sqrt(n).*(pistar - repmat(mi.pi_est,1,n));
			alpha=0.95;
			conf = zeros(m,3);


			for j=1:m
				t = T(j,:);
				t = sort(t);
				lndx = ceil(n*((1-alpha)/2));
				undx = ceil(n - lndx);

				mu=mi.pi_est(j);
				conf(j,1) = mu;
				conf(j,2) = mu - t(undx)/sqrt(n);
				conf(j,3) = mu - t(lndx)/sqrt(n);
			end

			global im;
			fg=figure(im);clf(fg);

			hold on

			for i=1:m
				plot(i.*ones(1,2), conf(i,[2 3]), 'g.-');
			end

			plot(mi.pi_est,'k.')
			plot(mi.pi_true,'rx')
			hold off
			flabel('j','\pi',[num2str(alpha*100) '% Boostrap error bars (n=' num2str(size(pistar,2)) '). center=green, non=red'])




			conf2 = zeros(mi.num_models,3);
			B = size(pistar,2);
			clev=alpha;
			lndx = ceil(B*((1-clev)/2));
			undx = ceil(B - lndx);
			for i=1:mi.num_models
				sv=sort(pistar(i,:));
				mu=mi.pi_est(i);
				conf2(i,1) = mu;
				conf2(i,2) = sv(lndx);
				conf2(i,3) = sv(undx);
			end

			% errorbar(1:mi.num_models,conf(:,1),conf(:,2),conf(:,3),'k.');
			hold on

			for i=1:mi.num_models
				plot(i.*ones(1,2), conf2(i,[2 3]), 'r.-');
			end

			plot(mi.pi_est,'k.')
			plot(mi.pi_true,'rx')
			hold off
		end
		
		function conf = bootstrap_plot_bars(mi,pistar)
			n=size(pistar,2);
			m=size(pistar,1);

			T = sqrt(n).*(pistar - repmat(mi.pi_est,1,n));
			alpha=0.95;
			conf = zeros(m,3);


			for j=1:m
				t = T(j,:);
				t = sort(t);
				lndx = ceil(n*((1-alpha)/2));
				undx = ceil(n - lndx);

				mu=mi.pi_est(j);
				conf(j,1) = mu;
				conf(j,2) = mu - t(undx)/sqrt(n);
				conf(j,3) = mu - t(lndx)/sqrt(n);
			end

			global im;
			fg=figure(im);clf(fg);

			hold on

			for i=1:m
				plot(i.*ones(1,2), conf(i,[2 3]), 'c.-');
			end
			
			%alt method test
			if 1==111
			for j=1:m
				t = T(j,:);
				t = sort(t);
				lndx = ceil(n*((1-alpha)/2));
				undx = ceil(n - lndx);

				mu=mi.pi_est(j);
				conf(j,1) = mu;
				conf(j,2) = mu + t(undx)/sqrt(n);
				conf(j,3) = mu + t(lndx)/sqrt(n);
			end
			for i=1:m
				plot(i.*ones(1,2), conf(i,[2 3]), 'r.-');
			end
			
			end
			%end

			plot(mi.pi_est,'k.')
			plot(mi.pi_true,'rx')
			hold off
			flabel('j','\pi',[num2str(alpha*100) '% Boostrap error bars (n=' num2str(size(pistar,2)) '). x=true, .=est'])
		end

		function plot_bootstrap_covar(mi, S, V, correl, stdev, pistar)
			w = max(mi.pi_est.*2000,10);
			global im;
			
			fg=figure(im);im=im+1;
			clf(fg)
			hold on
			scatter(1:mi.num_models,mi.variance,w,[0 0 0],'filled')
			scatter(1:mi.num_models,V,w,[.5 .1 .3])
			hold off
			flabel('j','Variance',['Filled: Info variance. Open: Boostrap variance (n=' num2str(size(pistar,2)) ')'])
			axis([0 size(pistar,1) min(min(min([mi.variance;V])),-.01) max(max(max([mi.variance;V])),.01)])


			fg=figure(im);im=im+1;
			clf(fg)
			hold on
			scatter(1:mi.num_models,mi.stdev,w,[0 0 0],'filled')
			scatter(1:mi.num_models,stdev,w,[.5 .1 .3])
			hold off
			flabel('j','Standard deviation',['Filled: Info \sigma. Open: Boostrap \sigma (n=' num2str(size(pistar,2)) ')'])
			axis([0 size(pistar,1) min(min(min([mi.stdev;stdev])),-.01) max(max(max([mi.stdev;stdev])),.01)])


			fg=figure(im);im=im+1;
			clf(fg)

			stddelta=mi.stdev-stdev;
			scatter(1:mi.num_models,stddelta,w,[0 0 0],'filled')

			flabel('\pi','Standard deviation delta',['Info v. Boostrap \sigma (n=' num2str(size(pistar,2)) ')'])
			axis([0 size(pistar,1) min(min(stddelta),-.1) max(max(stddelta),.1)])



			fg=figure(im);im=im+1;
			clf(fg)

			w = zeros(mi.num_models,mi.num_models);
			for i=1:mi.num_models
				for j=1:mi.num_models
					w(i,j) = mi.pi_est(i) + mi.pi_est(j);
				end
			end
			w=max(200.*reshape(w,1,size(w,1)^2),10);

			hold on
			stddelta=reshape(mi.correl,1,size(mi.correl,1)^2)-reshape(correl,1,size(correl,1)^2);
			scatter(1:mi.num_models^2,stddelta,w,[0 0 0],'filled')
			hold off
			flabel('j','Correlation delta',['Info v. Boostrap Correlation (n=' num2str(size(pistar,2)) ')'])
			axis([0 length(stddelta) min(min(stddelta),-.1) max(max(stddelta),.1)])





			fg=figure(im);im=im+1;
			clf(fg)
			rc = ceil(sqrt(mi.num_models));
			for i=1:mi.num_models
				subplot(rc,rc,i);
				hist(pistar(i,:))
			end


		end
		
		function plot_info_error_bars(mi,nstdevs)
			global im;
			fg=figure(im);im=im+1;
			clf(fg)

			clrss = {'c', 'm', 'k', 'g'};
			hold on
			for i=-nstdevs:-1
				errorbar(1:mi.num_models,mi.pi_est,-i.*mi.stdev,strcat([clrss{-i} '.']));
			end
			
			plot(mi.pi_true,'rx')
			plot(mi.pi_est,'k.')
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