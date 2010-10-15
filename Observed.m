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
			
			max_seconds = arg('max_seconds',NaN);
			
			m = size(ob.f_ab,2);

			if ~isfinite(p)
				p = eval(init_str);
				p(p==0) = max(p);
			end
			
			p = p./sum(p);
			
			init_p = p;
			
			changing_norms = [];
			
			P = NaN.*ones(m,max_iters);
			norms = zeros(1,1);
% 			plike = zeros(1,1);
			ll = zeros(1,1);
			
			w = zeros(n,m);
			
			global im;
			figure(im);
			
			Tmr = tic;
			P(:,1) = p;
			for counter=1:max_iters
				p0 = p;
				
				if interactive
					norms_pct_change = norms';
					norms_pct_change = (norms_pct_change-circshift(norms_pct_change,1))./circshift(norms_pct_change,1);
			
					ob.plot_progress(norms,norms_pct_change,ll,P,p,n,m,counter,-counter,im,init_p);
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
			
			
			
			norms_pct_change = norms';
			norms_pct_change = (norms_pct_change-circshift(norms_pct_change,1))./circshift(norms_pct_change,1);
			ob.plot_progress(norms,norms_pct_change,ll,P,p,n,m,counter,tmr,im,init_p);
			
			for i=2:length(norms_pct_change)
				if sign(norms_pct_change(i)) ~= sign(norms_pct_change(i-1))
					changing_norms(length(changing_norms)+1) = i;
				end
			end
			if length(changing_norms > 0)
				changing_norms
			end
			
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
		
		function plot_progress(ob,norms,norms_pct_change,ll,P,p,n,m,counter,tmr,im,init_p)
			figure(im);
			spr=3;spc=3;
			
			subplot(spr,spc,1);
			roll_window_size = 50;
			
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
			
			
			
			subplot(spr,spc,6);
			
			if kk>roll_window_size
				norms_pct_change = norms_pct_change(kk-roll_window_size:kk);
			end
			plot(norms_pct_change,'k.-','Color',[.8 .8 .8]);
			
			flabel('Trial','% change',['Rolling ' num2str(roll_window_size) ' norm % change']);
			
			subplot(spr,spc,8);
			kk = length(norms);
			if kk>roll_window_size
				norms = norms(kk-roll_window_size:kk);
			end
			plot(norms,'k.-','Color',[.8 .8 .8]);
			
			flabel('Trial','norm',['Rolling ' num2str(roll_window_size) ' norm']);
			
			
			pi_colors = 10.*(1+init_p);
			
			

			clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
				'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
			
			subplot(spr,spc,3);
			
			P = P(:,1:counter);
			
			plot(1:size(P,2),P(1,:),strcat([clrs(1) '-']));
			for i=2:m
				hold on
				plot(1:size(P,2),P(i,:),strcat([clrs(i) '.-']));
				%scatter(1:size(P,2),P(i,:),5.*(1+P(i,:)),pi_colors(i).*ones(size(P,2),1),'filled')
				hold off
			end
			flabel('Trial','\pi_j',strcat(['Weights over time']));
			
			
			subplot(spr,spc,9);
			if counter>roll_window_size
				PR = P(:,counter-roll_window_size:counter);
			else
				PR = P;
			end
			plot(1:size(PR,2),PR(1,:),strcat([clrs(1) '-']));
			for i=2:m
				hold on
				plot(1:size(PR,2),PR(i,:),strcat([clrs(i) '.-']));
				%scatter(1:size(P,2),P(i,:),5.*(1+P(i,:)),pi_colors(i).*ones(size(P,2),1),'filled')
				hold off
			end
			flabel('Trial','\pi_j',strcat(['Rolling ' num2str(roll_window_size) ' weights over time']));
			
			
% 			subplot(spr,spc,3);
% 			scatter(1:ob.num_models,P(:,1),20.*(1+P(:,1)),10.*(1+P(:,1)))
% 			flabel('j','\pi_j','Starting weights');
% 			
% 			subplot(spr,spc,4);
% 			scatter(1:ob.num_models,p,20.*(1+p),10.*(1+p),'filled')
% 			flabel('j','\pi_j',strcat(['Weights. \pi^{(0)}=' init_str]));
			
			
			
			subplot(spr,spc,4);
			scatter(1:ob.num_models,init_p,20.*(1+init_p),pi_colors)
			hold on
			scatter(1:ob.num_models,p,20.*(1+p),pi_colors,'filled')
			hold off
			flabel('j','\pi_j',strcat(['Starting and ending \pi']));

% 			ppdiff=init_p-p;
% 			subplot(spr,spc,5);
% 			scatter(1:ob.num_models,ppdiff,20.*(1+ppdiff),pi_colors,'filled')
% 			flabel('j','Abs difference',strcat(['\pi^{(0)} vs current \pi']));
			
			if isfinite(ob.p_actual)
				ppdiff=ob.p_actual-p;
				subplot(spr,spc,5);
				scatter(1:ob.num_models,ppdiff,20.*(1+ppdiff),pi_colors,'filled')
				flabel('j','Abs difference',strcat(['Current \pi versus true']));
			end
			
			
% 			if isfinite(ob.p_actual) & size(P,2) > 1
% 				subplot(spr,spc,6);
% 				pdc=p;
% 				pdc = (pdc-P(:,size(P,2)-1))./P(:,size(P,2)-1);
% 				plot(pdc,'b.-')
% % 				hold on
% % % 					pdc = p';
% % % 					pdc = (pdc-ob.p_actual')./ob.p_actual';
% % % 					plot(pdc,'r.-')
% % % 					pdc = p';
% % % 					pdc = (pdc-init_p')./init_p';
% % % 					plot(pdc,'g.-')
% % 				hold off
% 				flabel('j','% change',['% \Delta \pi over last (' num2str(counter) ')']);
% 			end
			
			
			subplot(spr,spc,2);
			plot(ll,'b.-')
			hold on
			plot([0 length(ll)], [ob.actual_log_like ob.actual_log_like], 'r--');
			hold off
			flabel('Trial','l(\pi|x,y)',strcat(['l(\theta)=' num2str(ll(length(ll)))]));
			
			
			llr = ll;
			kk=length(llr);
			if kk>roll_window_size
				llr = llr(kk-roll_window_size:kk);
			end
% 			subplot(spr,spc,8);
% 			plot(llr,'b.-')
% 			flabel('Trial','l(\pi|x,y)',['Rolling ' num2str(roll_window_size) ' l(\theta)']);
			
			llr = ll';
			llr = (llr-circshift(llr,1))./circshift(llr,1);
			kk=length(llr);
			if kk>roll_window_size
				llr = llr(kk-roll_window_size:kk);
			end
			subplot(spr,spc,7);
			plot(llr,'b.-','Color',[.7 .8 .8])
			flabel('Trial','l(\pi|x,y)',['% \Delta l(\theta) roll ' num2str(roll_window_size)]);
			
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