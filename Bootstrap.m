classdef Bootstrap < handle
	
	properties
		pi;
		pi_true;
		loglike;
		mi_name;
		filename;
		
		correlation,variance,covariance,stdev;
		
		type_string;
	end
	
	
	methods(Static)
		function bo = generate_parametric(filename,mi,base_number,B,N,varargin)
			bo = Bootstrap.generate(filename,false,mi,base_number,B,N,cell2mat(varargin));
			bo.type_string = 'Parametric';
		end
		
		function bo = generate_non_parametric(filename,mi,base_number,B,N,varargin)
			bo = Bootstrap.generate(filename,true,mi,base_number,B,N,cell2mat(varargin));
			bo.type_string = 'Nonparametric';
		end
		
		% private
		function bo = generate(filename,do_m_n,mi,base_number,B,N,varargin)
			file_base = strcat(['temp/' filename '_']);
			load_args
			interactive = arg('interactive',false);
			tic
			global im;
			xim=im;
			if interactive
				im=1;
			else
				im=-1;
			end
			vargs = struct(	'max_iters', arg('max_iters',200),...
							'quick_print', arg('quick_print',200),...
							'interactive',interactive,...
							'init_p',arg('init_p',NaN));
			pmt = bootstrap_generate(file_base,base_number,mi,B,N,do_m_n,vargs);
			im=xim;
			toc
			
			bo = Bootstrap(filename,mi,file_base,(base_number+1):B);
			bo.save
		end
	end
	
	
	methods
		%--loading
		function bo = Bootstrap(filename,mi,mi_name_base,mi_ndx_range)
			bo.filename = filename;
			bo.mi_name = mi.filename;
			
			J = length(mi_ndx_range);
			pi = zeros(mi.num_models,J);
			loglike = zeros(J,1);
			
			for j=1:J
				i = mi_ndx_range(j);
				mi = Mixture.load(strcat([mi_name_base num2str(i)]));
				pi(:,j) = mi.pi_est;
				loglike(j) = mi.loglike_est;
			end
			
			bo.pi = pi;
			bo.loglike = loglike;
			
			bo.update_covariance;
		end
		
		function save(bo)
			if isfinite(bo.filename)
				save(strcat(['cache/' bo.filename '.mat']),'bo');
				disp(strcat(['Saved results to cache/' bo.filename '.mat']));
			end
		end
		
		function save_as(bo,filename)
			bo.filename = filename;
			bo.save();
		end
		
		
		
		
		
		
		
		%--calcs
		function [S,V,correl,stdev] = update_covariance(bo)
			mi = Mixture.load(bo.mi_name);
			[S,V,correl,stdev] = bootstrap_covar(mi,bo.pi);
			bo.variance = V;
			bo.covariance = S;
			bo.correlation = correl;
			bo.stdev = stdev;
		end
		
		
		
		
		
		
		
		
		
		
		
		
		%--plots
		function [conf,conf2] = plot_error_bars_both(bo)
			mi = Mixture.load(bo.mi_name);
			n=size(bo.pi,2);
			m=size(bo.pi,1);

			T = sqrt(n).*(bo.pi - repmat(mi.pi_est,1,n));
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
			flabel('j','\pi',[bo.mi_name ': ' num2str(alpha*100) ...
				'% Boostrap error bars (n=' num2str(size(bo.pi,2)) '). center=green, non=red'])




			conf2 = zeros(mi.num_models,3);
			B = size(bo.pi,2);
			clev=alpha;
			lndx = ceil(B*((1-clev)/2));
			undx = ceil(B - lndx);
			for i=1:mi.num_models
				sv=sort(bo.pi(i,:));
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
		
		function conf = plot_error_bars(bo)
			mi = Mixture.load(bo.mi_name);
			n=size(bo.pi,2);
			m=size(bo.pi,1);

			T = sqrt(n).*(bo.pi - repmat(mi.pi_est,1,n));
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
			im=im+1;
			hold on

			for i=1:m
				h1=plot(i.*ones(1,2), conf(i,[2 3]), 'c.-', 'Color', [.5 .5 .5]);
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

			h2=plot(mi.pi_est,'b*');
			h3=plot(mi.pi_true,'go');
			
			legend([h2 h3],'Est \pi','True \pi')
			
			hold off
			flabel('j','\pi',[bo.mi_name ': ' num2str(alpha*100) '% ' bo.type_string ...
				' boostrap errors. B=' ...
				num2str(size(bo.pi,2)) '.'])
		end

		function plot_covar(bo)
			[S,V,correl,stdev] = covariance(bo);
			mi = Mixture.load(bo.mi_name);
			w = max(mi.pi_est.*2000,10);
			global im;
			
			fg=figure(im);im=im+1;
			clf(fg)
			hold on
			scatter(1:mi.num_models,mi.variance,w,[0 0 0],'filled')
			scatter(1:mi.num_models,V,w,[.5 .1 .3])
			hold off
			flabel('j','Variance',[bo.mi_name ': Filled: Info variance. Open: Boostrap variance (n=' ...
				num2str(size(bo.pi,2)) ')'])
			axis([0 size(bo.pi,1) min(min(min([mi.variance;V])),-.01) max(max(max([mi.variance;V])),.01)])


			fg=figure(im);im=im+1;
			clf(fg)
			hold on
			scatter(1:mi.num_models,mi.stdev,w,[0 0 0],'filled')
			scatter(1:mi.num_models,stdev,w,[.5 .1 .3])
			hold off
			flabel('j','Standard deviation',[bo.mi_name ': Filled: Info \sigma. Open: '...
				'Boostrap \sigma (n=' num2str(size(bo.pi,2)) ')'])
			axis([0 size(bo.pi,1) min(min(min([mi.stdev;stdev])),-.01) max(max(max([mi.stdev;stdev])),.01)])


			fg=figure(im);im=im+1;
			clf(fg)

			stddelta=mi.stdev-stdev;
			scatter(1:mi.num_models,stddelta,w,[0 0 0],'filled')

			flabel('\pi','Standard deviation delta',[bo.mi_name ': Info v. Boostrap \sigma (n=' ...
				num2str(size(bo.pi,2)) ')'])
			axis([0 size(bo.pi,1) min(min(stddelta),-.1) max(max(stddelta),.1)])



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
			flabel('j','Correlation delta',[bo.mi_name ': Info v. Boostrap Correlation (n=' ...
				num2str(size(bo.pi,2)) ')'])
			axis([0 length(stddelta) min(min(stddelta),-.1) max(max(stddelta),.1)])





			fg=figure(im);im=im+1;
			clf(fg)
			rc = ceil(sqrt(mi.num_models));
			for i=1:mi.num_models
				subplot(rc,rc,i);
				hist(bo.pi(i,:))
			end


		end
		
		function plot_data_coverage(bo,bootname)
			mi = Mixture.load(bo.mi_name);
			h=figure(2);
			clf(h);
			subplot(1,2,1);
			fi=1;
			ndxa=1:size(mi.x,1);
			ndx = mi.f(ndxa,fi)>0;
			F=mi.f(ndx,fi);
			scatter(mi.x(ndx,1), mi.x(ndx,2), max(1,F), [.1 .1 .2], 'filled');
			flabel('x','y',['halo3'])
			halo3_F = [min(F) max(F)]
			axis([-3 0 -.3 .6])

			
			mii = Mixture.load(bootname);
			subplot(1,2,2);
			fi=1;
			ndx = mii.f(ndxa,fi)>0;
			F=mii.f(ndx,fi);
			scatter(mii.x(ndx,1), mii.x(ndx,2), max(1,F), [.5 .4 .1], 'filled');
			flabel('x','y',['GEN'])
			gen_F = [min(F) max(F)]
			axis([-3 0 -.3 .6])

		end
		
		function plot_spread(bo)
			global im;
			h=figure(im);
			im=im+1;
			clf(h);
			subplot(1,2,1);
			hold on
			m=size(bo.pi,1);
			maxy=max(.5,max(max(abs(bo.pi))));
			
			for i=1:m
				plot([i i], [0 maxy], 'Color', [.5 .5 .5], 'LineStyle',':')
			end

			h1=plot(1:m,bo.pi(:,1),'k.');
			if size(bo.pi,2) > 1
				plot(1:m,bo.pi(:,2:size(bo.pi,2)),'k.');
			end
			h2=plot(1:m,bo.pi_true,'b*');
			
			legend([h1 h2],'Est \pi','True \pi');
			

			flabel('j','\pi_j',[bo.type_string ' B=' num2str(length(bo.loglike)) '. ' bo.mi_name])
			axis([0 m 0 maxy])
			hold off

			subplot(1,2,2);
			hold on
			m=size(bo.pi,1);
			delt = mean(bo.pi')'-bo.pi_true;
			
			maxy=max(.5,max(max(abs(delt))));
			miny=-maxy;
			
			for i=1:m
				plot([i i], [miny maxy], 'Color', [.5 .5 .5], 'LineStyle',':')
			end

			plot(1:m,delt,'r>')

			flabel('j','\pi_j',['mean(est \pi) - true \pi'])
			axis([0 m miny maxy])

			hold off
		end
	end
	
end