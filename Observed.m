classdef Observed < handle
	
	properties
		dpath='data/obsdata2.dat';
		model_number;
		point_count;
		
		
		name;
		data;
		x,y,lum;
		f_ab;
		num_models;
		
		actual_log_like=NaN;
		p_actual=NaN;
	end
	
	
	methods(Static)
		function [obx, mox] = load(model_no,sample_size)
			mox = Model.load(strcat(['cache/models_' num2str(model_no) '.mat']));
			load(strcat(['cache/observed_' num2str(model_no) '_' num2str(sample_size) 'k.mat']));
			obx=ob;
			disp(strcat(['Loaded model data ' num2str(model_no) ' for n=' num2str(sample_size) 'k']))
		end
		
		function save(ob)
			path = strcat(['cache/observed_' num2str(ob.model_number) '_' num2str(ob.point_count) 'k.mat']);
			save(path,'ob');
			show(['Saved to ' path])
		end
	end
	
	
	methods
		%loading
		function ob = Observed(varargin)
			load_args
			dpath = arg('path',ob.dpath);
			ob.dpath = dpath;
			ob.name = arg('name','halo000');
			data = arg('data',NaN);
			ob.p_actual = arg('p_actual',NaN);
			
			if ~isfinite(data)
				ob.data = load(dpath);
			else
				ob.data = data;
			end
			
			x = ob.data(:,1);
			y = ob.data(:,2);
			lum = ob.data(:,3);
			
			ndx = x>-3;
			ob.x=x(ndx);
			ob.y=y(ndx);
			ob.lum=lum(ndx);
		end
		
		function load_models(ob,mo)
			n = length(ob.x);
			m = length(mo.data);
			
			f_ab = NaN.*zeros(n,m);
			
			for i=1:n
				for j=1:m
					try
					tval = mo.f_ab(j,ob.x(i),ob.y(i));
					catch
						ob.x(i)
						ob.y(i)
						sepr
						tval = mo.f_ab(j,max(ob.x(i)),max(ob.y(i)));
					end
					try
						f_ab(i,j) = tval;
					catch
						[ob.x(i) ob.y(i)]
						tval
						sepr
						f_ab(i,j) = max(tval);
					end
				end
			end
			
			ob.f_ab = f_ab;
			ob.num_models = size(f_ab,2);
			
			
			if isfinite(ob.p_actual)
				n = length(ob.x);
				m = size(ob.f_ab,2);
				ob.actual_log_like = ob.complete_likelihood(ob.p_actual,n,m);
			end
		end
		
		
		
		
		
		
		%---
		function [p,P] = em_bodhi(ob, varargin)
			load_args
			
			max_iters = arg('max_iters',100000);
			min_iters = arg('min_iters',5);
			n = arg('n',length(ob.x));
			min_norm = arg('min_norm',0.0000000001);
			init_str = arg('init','rand(m,1)');
			p = arg('p',NaN);
			interactive = arg('interactive',true);
			max_seconds = arg('max_seconds',NaN);
			baseline_p = arg('baseline_p',NaN);
			
			obs_path = arg('obs_path','');
			model_path = arg('model_path','');
			
			
			P = zeros(16,1);
			norms = [];
			ll = [];
			
			fid = fopen(obs_path,'r'); x = fscanf(fid,'%f',[3,inf]); x= x';
			

			xmin = -3; xmax = 0; ymin = -0.5; ymax = 1;
			idx = ((x(:,1) >= xmin) & (x(:,1) <= xmax) & (x(:,2) >= ymin) & (x(:,2) <= ymax));

			x  = x(idx,:);
			
			fid = fopen(model_path,'r'); data = fscanf(fid,'%f',[9,inf]); data = data';



			% return;
 			r = sqrt(2*0.05*0.05);
 			t = (1/8:1/4:1)'*2*pi; x1 = r*sin(t); x2 = r*cos(t);

 			k = 0; s = 0; 
			for i=0:4,
 				for j =0:4,
 					k = k + 1;
 					idx = ((data(:,4)==j) & (data(:,5)==i) & (data(:,3)~=0));
 					distx = data(idx==1,:); n = length(distx);
 					%subplot(5,5,k); axis([-3 0 -0.5 1]);
 					if (n ~= 0) s = s + 1; end
%  					for l = 1:n,
%  						subplot(5,5,k); fill(distx(l,1)+x1,distx(l,2)+x2,[distx(l,3) 0 0]); hold on;
%  					end
 					%xlabel('[\alpha/Fe]'); ylabel('[Fe/H]'); axis([-3 0 -0.5 1]);
 					% return;
 				end
 			end
 			T = s;
 			xn = 30; yn = 17;
 			xgrid = -2.95:0.1:-0.05;
 			ygrid = -0.65:0.1:0.95;
 

 			k=0; sq = zeros(xn,yn);
 			f = zeros(T,xgrid,ygrid);
 			indic= zeros(25,1); m = 0;
 			for i=0:4,
 				for j =0:4,
 					m = m + 1;
 					idx = ((data(:,4)==j) & (data(:,5)==i));
 					distx = data(idx==1,:); n = sum(distx(:,3)~=0);
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
 					%figure(5); subplot(5,5,m); mesh(xgrid,ygrid,sq'); hold on;
 				end
 			end

			f = 100*f;
			fval = EvalDens(x(:,1),x(:,2),f,xgrid,ygrid,T,xn,yn);
			id = (sum(fval')~=0); fval = fval(id,:); 

			% return
			%The EM algorithm
			u = gamrnd(1/T,1,T,1);
			pi_0 = u/sum(u); 
			pi_1 = ones(T,1)/T;
			pi_0 = pi_1 + 1;
			k = 0;
			N = length(fval);
			w = zeros(N,T);
			
			init_p = pi_1;

			global im;
			if interactive
				tfgr = figure(im);
				clf(tfgr);
			end
			Tmr = tic;
			
			for counter=1:max_iters
				k = k + 1;
				s = 0;
				for i = 1:N,
					w(i,:) = pi_1.*fval(i,:)'; 
					tot = sum(w(i,:));
					w(i,:) = w(i,:)/tot;
					s = s + log(tot);
				end
				pi_0 = pi_1; 
				pi_1 = sum(w)'/N;
				
				tl_a = counter/50;tl_whole = floor(tl_a);tl_part = tl_a-tl_whole;
				if tl_part == 0
					pi_1
				end
				
				P(:,counter) = pi_1;
				norms(counter) = norm(pi_0 - pi_1);
				ll(counter) = s;
				
				
				if interactive
					ob.plot_progress(norms,ll,P,pi_1,N,16,counter,-counter,im,init_p,baseline_p);
				end
				
				if isfinite(max_seconds)
					tmr = toc(Tmr);
					el_sec = tmr;
					if max_seconds < el_sec
						disp('time reached');
						break;
					end
				else
					if norms(counter) < min_norm && counter > min_iters
						%disp('min norm reached; stopping')
						break;
					end
					
					if counter==max_iters
						warning('min norm not reached!')
					end
				end
			end
			tmr = toc(Tmr);
			
			ob.plot_progress(norms,ll,P,pi_1,N,16,counter,tmr,im,init_p,baseline_p);
			
			
			
			p = pi_1;
			
		end
		
		
		%----
		
		
		
		
		%em
		function [p,P,ll,LL] = em(ob, varargin)
			load_args
			
			max_iters = arg('max_iters',100000);
			min_iters = arg('min_iters',5);
			n = arg('n',length(ob.x));
			min_norm = arg('min_norm',0.00000001);
			init_str = arg('init','rand(m,1)');
			p = arg('p',NaN);
			interactive = arg('interactive',true);
			baseline_p = arg('baseline_p',NaN);
			
			max_seconds = arg('max_seconds',NaN);
			
			m = size(ob.f_ab,2);

			if ~isfinite(p)
				p = eval(init_str);
				p(p==0) = max(p);
			end
			
			p = p./sum(p);
			
			init_p = p;
			
			P = NaN.*ones(m,max_iters);
			norms = zeros(1,1);
% 			plike = zeros(1,1);
			ll = zeros(1,1);
			
			w = zeros(n,m);
			
			
			global im;
			if interactive
				tfgr = figure(im);
				clf(tfgr);
			end
			
			Tmr = tic;
			P(:,1) = p;
			for counter=1:max_iters
				p0 = p;
				
				if interactive
					ob.plot_progress(norms,ll,P,p,n,m,counter,-counter,im,init_p,baseline_p);
				end
				
				for j=1:m
					for i=1:n
						p_f_ak_bk = zeros(m,1);

% 						xi = ob.x(i);
% 						yi = ob.y(i);

						for k=1:m
							p_f_ak_bk(k) = p0(k).*ob.f_ab(i,k);
						end
						
						wij = (p0(j).*ob.f_ab(i,j)) ./ sum(p_f_ak_bk);
						if ~isfinite(wij)
							wij = 0;
						end
						w(i,j) = wij;
					end
					
					p(j) = sum(w(:,j))/n;
				end
				
				P(:,counter+1) = p;
				
				%print every 10 p's
				tl_a = counter/10;tl_whole = floor(tl_a);tl_part = tl_a-tl_whole;
				if tl_part == 0
					p
				end
				
				
				norms(counter) = abs(norm(p-p0));
				
				%plike(counter) = ob.partial_likelihood(p,n,m);
				ll(counter) = ob.complete_likelihood(p,n,m);
				
				if isfinite(max_seconds)
					tmr = toc(Tmr);
					el_sec = tmr;
					if max_seconds < el_sec
						disp('time reached');
						break;
					end
				else
					if norms(counter) < min_norm && counter > min_iters
						%disp('min norm reached; stopping')
						break;
					end
					
					if counter==max_iters
						warning('min norm not reached!')
					end
				end
			end
			tmr = toc(Tmr);
			
			
			ob.plot_progress(norms,ll,P,p,n,m,counter,tmr,im,init_p,baseline_p);
			
			
			
			im=im+1;
			
			LL = isfinite(ll);
			ll = ll(length(isfinite(ll)));
			
		end
		
		%fxy=log(\prod\pi*f_a,b)
		function l = partial_likelihood(ob,p,n,m)
			l = 0;
			for j=1:m
				for i=1:n
					l0 = log(p(j)) + log(ob.f_ab(i,j));
					if isfinite(l0) && ~isinf(l0)
						l = l+l0;
					end
				end
			end
		end
		
		%includes z
		function l = complete_likelihood(ob,p,n,m)
			l = 0;
			for i=1:n
				l0 = 0;
				for j=1:m
					l0 = l0 + p(j)*ob.f_ab(i,j);
				end
				lgl0 = log(l0);
				if ~isfinite(lgl0)
					lgl0=0;
				end
				l = l + lgl0;
			end
		end
		
		function plot_progress(ob,norms,ll,P,p,n,m,counter,tmr,im,init_p,baseline_p)
			figure(im);
			spr=3;spc=3;
			pi_colors = 10.*(1+init_p);
			pi_colors = [11.139; 19.786; 18.486; 10.506; 14.662; 13.257; 16.302; 12.303; 15.799; 16.032; 15.999; 14.484; 10.354; 15.138; 14.077;  11.08];
			y_ax_min_span = 1;
			clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
				'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
			
			subplot(spr,spc,1);
			roll_window_size = 50;
			mini_roll_window_size = 10;
			
			norms_pct_change = norms';
			norms_pct_change = (norms_pct_change-circshift(norms_pct_change,1))./circshift(norms_pct_change,1);
			
			plot(norms_pct_change,'k.-');
			kk = length(norms_pct_change);
			if kk>5
				kkdelt = std(norms_pct_change(5:kk));
				axis([5 kk min(norms_pct_change(5:kk))-kkdelt max(norms_pct_change(5:kk))+kkdelt])
			end
			if tmr<0
				flabel('Trial','% change',strcat(['(In progress) ' num2str(counter) ' runs, n='...
					num2str(n)]));
			else
				flabel('Trial','% change',strcat([num2str(counter) ' runs in ' ...
					num2str(tmr) 's, n=' num2str(n)]));
			end
			
			
			if isfinite(baseline_p)
				spbi = subplot(spr,spc,7);
				scatter(1:ob.num_models,(baseline_p-p)./p,20.*(1+init_p),pi_colors)
				flabel('j','% \Delta', '% \Delta v other estimate');
				
				if counter > 20 && max(abs((baseline_p-p)./p)) > 1
					fds=11;
				end
				
				axis(spbi, 'tight')
				currax = axis(spbi);
				if abs(abs(currax(3)) - abs(currax(4))) < y_ax_min_span
					currax(3) = currax(3)-y_ax_min_span;
					currax(4) = currax(4)+y_ax_min_span;
					axis(currax)
				end
				
				y_ax_min_span = 0.5;
 				spbi = subplot(spr,spc,8);
 				scatter(1:ob.num_models,baseline_p-p,20.*(1+init_p),pi_colors,'filled')
 				flabel('j','\Delta', '\Delta v other estimate');
				axis(spbi, 'tight')
				currax = axis(spbi);
				if abs(abs(currax(3)) - abs(currax(4))) < y_ax_min_span
					currax(3) = currax(3)-y_ax_min_span;
					currax(4) = currax(4)+y_ax_min_span;
					axis(currax)
				end
				y_ax_min_span = 1;
				
% 				subplot(spr,spc,7);
% 				scatter(1:ob.num_models,baseline_p,20.*(1+init_p),pi_colors)
% 				hold on
% 				scatter(1:ob.num_models,p,20.*(1+init_p),pi_colors,'filled')
% 				hold off
% 				flabel('j','\pi_j', '\pi vs other estimate');
			end
			
			
		
			
			subplot(spr,spc,3);
			
			P = P(:,1:counter);
			
			plot(1:size(P,2),P(1,:),strcat([clrs(1) '-']));
			for i=2:m
				hold on
				plot(1:size(P,2),P(i,:),strcat([clrs(i) '.-']));
				hold off
			end
			flabel('Trial','\pi_j',strcat(['Weights over time']));
			
			
			
			
			
			if counter>2
				subplot(spr,spc,6);
				
				PPC = P(:,1:counter)';
				ppc0=circshift(PPC,1);
				ppc0 = ppc0';
				PPC = PPC';
				ppc0=(PPC-ppc0)./ppc0;
				ppc0=ppc0(:,1:size(ppc0,2));
				PPC = 100.*ppc0;
				PPC = PPC(:,2:size(PPC,2));

% 				plot(2:size(PPC,2)+1,PPC(1,:),strcat([clrs(1) '-']));
% 				for i=2:m
% 					hold on
% 					plot(2:size(PPC,2)+1,PPC(i,:),strcat([clrs(i) '.-']));
% 					hold off
% 				end
% 				flabel('Trial','% \Delta',strcat(['% change \pi over last \pi']));
				
				subplot(spr,spc,6);
				
				roll_window_size = mini_roll_window_size;
				if counter>roll_window_size
					ppc_cols = counter-roll_window_size:counter-1;
					PPC = PPC(:,ppc_cols);
				else
					ppc_cols = 1:size(PPC,2);
				end
				
				plot(ppc_cols,PPC(1,:),strcat([clrs(1) '.-']));
				for i=2:m
					hold on
					plot(ppc_cols,PPC(i,:),strcat([clrs(i) '.-']));
					hold off
				end
				flabel('Trial','% \Delta',strcat(['Roll ' num2str(roll_window_size) ' % \Delta'...
					'\pi over last \pi']));
			end
			
			
			
			
			
			
			subplot(spr,spc,9);
			scatter(1:ob.num_models,init_p,20.*(1+init_p),pi_colors)
			hold on
			scatter(1:ob.num_models,p,20.*(1+p),pi_colors,'filled')
			hold off
			flabel('j','\pi_j',strcat(['Starting and ending \pi']));
			
			if isfinite(ob.p_actual)
				y_ax_min_span = 0.5;
				ppdiff=ob.p_actual-p;
				spbi=subplot(spr,spc,5);
				scatter(1:ob.num_models,ppdiff,20.*(1+ppdiff),pi_colors,'filled')
				flabel('j','\Delta',strcat(['\Delta v true']));
				axis(spbi, 'tight')
				currax = axis(spbi);
				if abs(abs(currax(3)) - abs(currax(4))) < y_ax_min_span
					currax(3) = currax(3)-y_ax_min_span;
					currax(4) = currax(4)+y_ax_min_span;
					axis(currax)
				end
				y_ax_min_span = 1;
				
				spbi=subplot(spr,spc,4);
				ppdiff=100.*(ob.p_actual-p)./p;
				dotweights = 20.*(1+init_p);
				sx = 1:ob.num_models;
				
				pct_diff_limit=35;
				
				ndx = ppdiff<pct_diff_limit & ppdiff>-pct_diff_limit;
				if sum(ndx)>0
					scatter(sx(ndx),ppdiff(ndx),dotweights(ndx),pi_colors(ndx))
				end
				hold on
				ndx = ppdiff>=pct_diff_limit | ppdiff<=-pct_diff_limit;
				if sum(ndx)>0
					scatter(sx(ndx),pct_diff_limit.*ones(1,sum(ndx)),50.*ones(1,sum(ndx)),pi_colors(ndx),'filled')
				end
				hold off
				
				axis(spbi, 'tight')
				currax = axis(spbi);
				if abs(abs(currax(3)) - abs(currax(4))) < y_ax_min_span
					currax(3) = currax(3)-y_ax_min_span;
					currax(4) = currax(4)+y_ax_min_span;
					axis(currax)
				end
				
				flabel('j','% \Delta', '% \Delta v true');
			end
			
			
			
			subplot(spr,spc,2);
			plot(ll,'b.-')
			hold on
			plot([0 length(ll)], [ob.actual_log_like ob.actual_log_like], 'r--');
			hold off
			flabel('Trial','l(\pi|x,y)',strcat(['l(\theta)=' num2str(ll(length(ll)))]));
			
			
			
		end
		
		function [best_ll, best_p, all_p, all_ll] = em_multi(ob, mo, varargin)
			load_args
			global im;
			xim=im;
			
			tic
			best_p = [];
			best_ll = NaN;
			
			count = arg('count',1);
			max_seconds = arg('max_seconds', NaN);
			max_iters = arg('max_iters', NaN);
			
			
			if isfinite(max_seconds)
				cont_cond_str = strcat([num2str(max_seconds) 's_']);
			else
				cont_cond_str = strcat([num2str(max_iters) 'i_']);
			end
			
			save_path = arg('save',strcat(['cache/em_multi_m' num2str(ob.model_number) '_rand_' cont_cond_str ...
				num2str(count) 'x.mat']));
			
			all_p = zeros(ob.num_models,count);
			all_ll = zeros(count,1);
			for i=1:count
				im=xim;
				[p,P,ll] = ob.em(cell2mat(varargin));
				all_p(:,i) = p;
				all_ll(i) = ll;
				
				if ~isfinite(best_ll) || best_ll < ll
					best_ll = ll;
					best_p = p;
				end
				
				disp(strcat(['Finished run ' num2str(i)]))
			end
			
			p = best_p;
			
			save(save_path,'p','all_p','all_ll','best_ll');
			disp(strcat(['Saved p,best_ll,all_p,all_ll to ' save_path]));
			
			ob.plot_weight_changes(struct('all_p',all_p));
			
			p_pct_change=circshift(all_p,1);
			p_pct_change=(all_p-p_pct_change)./p_pct_change;
			p_pct_change=p_pct_change(2:size(p_pct_change,1),:);
			ob.plot_weight_changes(struct('all_p',p_pct_change,'pct_change',true));
			
			
			print_pi(best_p,best_ll,mo);
			ob.plot_differences(best_p);
			toc
		end
		
		
		
		
		
		
		
		
		
		
		%convergence
		function ap = load_em_pis(ob,varargin)
			load_args
			error('no disk load atm')
			pth = arg('path','cache/all_p_50_full_runs.mat');
			
			load(pth);
			ap = all_p;
		end
		
		function plot_weight_changes(ob,varargin)
			load_args
			
			all_p = arg('all_p',NaN);
			pct_change = arg('pct_change',false);
			
			if ~isfinite(all_p)
				disp('Loading p_all from disk');
				all_p = ob.load_em_pis(cell2mat(varargin));
			end
			
			global im;
			fg=figure(im);
			clf(fg)
			im=im+1;
			
			m = size(all_p,1);
			p = size(all_p,2);

			clrs = rand(m,3);

			num_std_devs = 2;

			lngrd = 0:p+1;
			lnones = ones(1,p+2);

			for i=1:m
				subplot(4,4,i);
				
				if ~pct_change
					if isfinite(ob.p_actual)
						plot(lngrd,ob.p_actual(i).*lnones,'-','LineWidth',10, 'Color', [0.9 .9 0.9])
					end

					hold on
					plot(lngrd,mean(all_p(i,:)).*lnones,'-','Color',clrs(i,:))

					plot(lngrd,(num_std_devs*std(all_p(1,:))+mean(all_p(i,:))).*lnones,':','Color',clrs(i,:))
					plot(lngrd,(-num_std_devs*std(all_p(1,:))+mean(all_p(i,:))).*lnones,':','Color',clrs(i,:))
				end

				plot(all_p(i,:),'.','Color',clrs(i,:))
				if ~pct_change
					hold off
				end
				
				title(strcat(['\pi_{' num2str(i) '}']))
			end
		end
		
		
		
		
		
		
		
		
		%sim
		function [p,P,ll,LL] = simulate(ob, mo, varargin)
			load_args
			
			m = length(mo.data);
			n = arg('sample',100);
			save_data = arg('save_data', strcat(['cache/generated_' num2str(n) '_auto.mat']));
			obs_save = arg('obs_save', 'cache/observed_generated.mat');
			
			tic
			
			
			P = arg('p',ob.p_actual);%[.015 .2 .005 .07 .12 .005 .005 .01 .27 .1 .05 .004 .001 .02 .095 .03];
			if sum(P) ~= 1
				warning('P must sum to 1; renormalizing')
				P = P./sum(P)
			end

			if m ~= length(P)
				error('P and models are different lengths')
			end
			PC = cumsum(P);

			data = zeros(n,5);
			R1 = rand(n,1);
			RX = rand(n,1);
			RY = rand(n,1);

			grid_size = arg('grid_size',.1);

			for i=1:n
				r0 = rand();
				pi_ndx = sum(PC<=r0)+1;

				r1 = R1(i);
				nzd = mo.nonzeros(pi_ndx);
				cnzd = cumsum(nzd(:,3));
				nzda = [nzd cnzd];

				grid_ndx = sum(cnzd<=r1)+1;

				x = nzda(grid_ndx,1);
				y = nzda(grid_ndx,2);

				rx = RX(i)*grid_size;
				ry = RY(i)*grid_size;

				x = x + rx;
				y = y + ry;

				fab = mo.f_ab(pi_ndx,x,y);
				data(i,:) = [x y fab pi_ndx r0];
			end

			data = data(data(:,3)>0,:);


			fig
			subplot(1,3,1);
			scatter(data(:,1),data(:,2),data(:,3)./sum(data(:,3)),'k','filled');
			flabel('Fe/H','\alpha/Fe',[num2str(n) ' generated data points']);

			subplot(1,3,2);
			scatter(data(:,1),data(:,2),10,'r','filled');
			flabel('Fe/H','\alpha/Fe',[num2str(n) ' generated data points']);

			subplot(1,3,3);
			plot(P,'k.')
			flabel('j','\pi_j','True weights');


			save(save_data,'data');


			toc


			ob = Observed(struct('name','generated halo','data',data));
			ob.load_models(mo);
			ob.save(obs_save);



			'doing em'
			[p,all_p,ll,LL] = ob.em(cell2mat(varargin));
			'finished em'



			fig
			subplot(1,2,1)
			plot(P,'k.')
			hold on
			plot(p,'g.')
			hold off
			flabel('j','\pi_j','Actual (black) v predicted (green)');

			subplot(1,2,2)
			plot(P-p,'r.')
			flabel('j','Real - actual','Differences');

			
			
		end
		
		
		function plot_differences(ob,p)
			fig
			subplot(1,2,1);
			plot(ob.p_actual,'k.');
			hold on
			plot(p,'g.');
			hold off
			flabel('j','\pi_j','Actual (black) v predicted (green)');

			subplot(1,2,2);
			plot(ob.p_actual-p,'r.');
			flabel('j','Real - actual','Differences');
		end
		
		
		
		
		
		
		
		%visualization
		function plot(ob,varargin)
			load_args
			
			overlay = arg('overlay',false);
			dot_size = arg('dot_size',5);
			axl = arg('axis',[-3 0 -.5 1]);
			cmp = arg('colormap','winter');
			
			scatter_color_by_weight(ob.x,ob.y,log(ob.lum),struct(...
				'overlay',overlay,...
			'title', ob.name, 'x','Fe/H', 'y','\alpha/Fe','colormap',cmp,'dot_size',dot_size));
			axis(axl)
		end
		
	end
	
end