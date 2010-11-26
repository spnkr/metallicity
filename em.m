function [p,ll,P,LL,init_p,counter,tmr] = em(x,f,varargin)
	load_args
	
	max_iters = arg('max_iters',100000);
	min_iters = arg('min_iters',5);
	n = arg('n',length(x));
	min_norm = arg('min_norm',0.0000000001);
	p0_eval = arg('p0_eval',NaN);%rand(num_models,1)
	init_p = arg('init_p',NaN);
	interactive = arg('interactive',true);
	max_seconds = arg('max_seconds',NaN);
	baseline_p = arg('baseline_p',NaN);
	quick_print = arg('quick_print', 50);
	if ~isfinite(quick_print)
		quick_print = max_iters;
	end
	true_loglike = arg('true_loglike',NaN);
	num_models=arg('num_models',NaN);
	p_actual=arg('p_actual',NaN);
	
	%1.235 to 1.236 over past 25 runs is the slowest rate of convergence
	ll_stop_prec=arg('ll_stop_prec',1000);
	ll_stop_lookback=arg('ll_stop_lookback',25);
	
	norms = [];
	LL = [];
	
	% return
	%The EM algorithm
% 	u = gamrnd(1/num_models,1,num_models,1);
% 	pi_0 = u/sum(u); 
% 	pi_0 = pi_1 + 1;
	
	if isfinite(init_p)
		pi_1 = init_p;
	elseif isfinite(p0_eval)
		pi_1 = eval(p0_eval);
	else
		pi_1 = ones(num_models,1)/num_models;
	end
	pi_1 = pi_1./sum(pi_1);
	
	
	k = 0;
	N = length(f);
	w = zeros(N,num_models);

	init_p = pi_1;
	P = zeros(num_models,1);

	global im;
	if interactive && im>0
		tfgr = figure(im);
		clf(tfgr);
	end
	Tmr = tic;

	for counter=1:max_iters
		k = k + 1;
		s = 0;
		for i = 1:N,
			w(i,:) = pi_1.*f(i,:)'; 
			tot = sum(w(i,:));
			w(i,:) = w(i,:)/tot;
			s = s + log(tot);
		end
		pi_0 = pi_1; 
		pi_1 = sum(w)'/N;

		tl_a = counter/quick_print;tl_whole = floor(tl_a);tl_part = tl_a-tl_whole;
		if tl_part == 0
			pi_1.*100
		end
		
		P(:,counter) = pi_1;
		norms(counter) = norm(pi_0 - pi_1);
		ll = s;
		LL(counter) = ll;

		tl_a = counter/quick_print;tl_whole = floor(tl_a);tl_part = tl_a-tl_whole;
		if interactive && (counter < 20 || tl_part==0)
			plot_progress(norms,LL,P,pi_1,N,num_models,counter,-counter,im,init_p,baseline_p,...
				true_loglike,num_models,p_actual,quick_print);
		end
		
		if 1==11
		if counter>min_iters && (...
					sign((norms(counter)-norms(counter-1))/norms(counter-1)) ~= ...
						sign((norms(counter-1)-norms(counter-2))/norms(counter-2))...
				)
			disp(strcat(['norm percent change sign flipped over last; stopping at ' num2str(counter)]));
			break;
		end
		end
		
		
		if counter>min_iters && counter>ll_stop_lookback && round(ll_stop_prec*LL(counter))/ll_stop_prec...
				== round(ll_stop_prec*LL(counter-1))/ll_stop_prec && ...
				round(ll_stop_prec*LL(counter))/ll_stop_prec == ...
				round(ll_stop_prec*LL(counter-ll_stop_lookback))/ll_stop_prec
			disp(strcat(['ll did not change more than ' num2str(ll_stop_prec) ' dec places '...
				'in last ' num2str(ll_stop_lookback) ' iters; stopping at ' num2str(counter)]));
			break;
		end
		
		if isfinite(max_seconds)
			tmr = toc(Tmr);
			el_sec = tmr;
			if max_seconds < el_sec
				disp('time reached');
				break;
			end
		elseif counter==max_iters
			break;
		end
	end
	tmr = toc(Tmr);

	plot_progress(norms,LL,P,pi_1,N,num_models,counter,tmr,im,init_p,baseline_p,true_loglike,...
		num_models,p_actual,quick_print);

	

	p = pi_1;
end
	
	function plot_progress(norms,ll,P,p,n,m,counter,tmr,im,init_p,baseline_p,true_loglike,...
		num_models,p_actual,quick_print)
		if im<=0
			return
		end
		figure(im);
		spr=2;spc=3;
		pi_colors = [11.139; 19.786; 18.486; 10.506; 14.662; 13.257; 16.302; 12.303; 15.799; 16.032; 15.999; 14.484; 10.354; 15.138; 14.077;  11.08; 16.032; 15.999; 14.484; 10.354; 15.138; 14.077;  11.08];

		pi_colors = pi_colors(1:m);

		y_ax_min_span = 1;
		clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
			'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
			'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
			'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
		
		clrint = [  0.23309      0.16965      0.49939      ;
				  0.71771      0.93246      0.97015      ;
				  0.75749      0.40726      0.21711      ;
				  0.45631      0.16391      0.10073      ;
				  0.36744      0.13855      0.86039      ;
				  0.78951     0.094644      0.86793      ;
				  0.45087      0.25999      0.65318      ;
				 0.058676      0.34455      0.64045      ;
				  0.56037      0.93054       0.6123      ;
				 0.067491     0.052629      0.23289      ;
				  0.70662      0.62035      0.75643      ;
				 0.039426      0.51577      0.45674      ;
				  0.54358        0.868      0.88962      ;
				  0.41002      0.27975      0.87272      ;
				  0.50077      0.68458      0.59843      ;
				   0.4534      0.59915     0.096282      ;
				  0.46896      0.31409      0.42626      ;
				  0.90027      0.11488      0.11124      ;
				  0.24429      0.39885      0.79126      ;
				  0.62341      0.11236      0.99919      ;
				  0.34019       0.9606      0.47105      ;
				  0.46092      0.57113      0.54085      ;
				  0.20918      0.24431      0.64035      ;
				  0.98938      0.19832      0.29173      ;
				  0.32373      0.80925      0.43316      ;
				  0.72088      0.65298      0.55861      ;
				  0.32011      0.65065      0.63479      ;
				  0.57024      0.80118      0.10692      ;
				0.0075955      0.36107      0.93856      ;
				  0.23355      0.54492      0.48231      ];
		
		clrint = min(1,clrint.*2.5);
		
		
		clrint = [  1 0 0      ;
				  0 1 0      ;
				  0 0 1      ;
				    .8 .2 0      ;
				  0 .8 .2      ;
				  .2 0 .8      ;
				    .2 .8 0      ;
				  .5 .5 0      ;
				  0 .5 .5      ;
				    .5 0 .5      ;
				  .2 .1 .6      ;
				  .6 .1 .2      ;
				    .2 .1 .6      ;
				  .2 .6 .1     ;
				  .1 .6 .2     ;
				    1 0 0      ;
				  0 1 0      ;
				  0 0 1      ;
				    1 0 0      ;
				  0 1 0      ;
				  0 0 1      ;
				    1 0 0      ;
				  0 1 0      ;
				  0 0 1      ;
				    1 0 0     ];
		
		
		subplot(spr,spc,4);
		roll_window_size = quick_print;

		norms_pct_change = norms';
		norms_pct_change = (norms_pct_change-circshift(norms_pct_change,1))./circshift(norms_pct_change,1);

		plot(norms_pct_change,'k.-');
		kk = length(norms_pct_change);
		if kk>5
			kkdelt = std(norms_pct_change(5:kk));
			if isfinite(kkdelt)
				axis([5 kk min(norms_pct_change(5:kk))-kkdelt max(norms_pct_change(5:kk))+kkdelt])
			end
		end
		if tmr<0
			flabel(['Trial ' num2str(1000*norms_pct_change(length(norms_pct_change)))],'% change',strcat(['(In progress) ' num2str(counter) ' runs, n='...
				num2str(n)]));
		else
			flabel('Trial','% change',strcat([num2str(counter) ' runs in ' ...
				num2str(tmr) 's, n=' num2str(n)]));
		end


		if isfinite(baseline_p)
			spbi = subplot(spr,spc,7);
			scatter(1:num_models,(baseline_p-p)./p,20.*(1+init_p),pi_colors)
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
			scatter(1:num_models,baseline_p-p,20.*(1+init_p),pi_colors,'filled')
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

		%plot(1:size(P,2),P(1,:),strcat([clrs(1) '.-']));
		for i=1:m
			hold on
			plot(0:size(P,2)-1,P(i,:),'Color',clrint(i,:));
			hold off
		end
		flabel('Trial','\pi_j',strcat(['Weights over time']));
		axis([0 max(1,size(P,2)-1) 0 max(0.001,max(max(P)))])




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
			
			if counter>roll_window_size
				ppc_cols = counter-roll_window_size:counter-1;
				PPC = PPC(:,ppc_cols);
			else
				ppc_cols = 1:size(PPC,2);
			end

			plot(ppc_cols,PPC(1,:),'Color',clrint(1,:));
			for i=2:m
				hold on
				plot(ppc_cols,PPC(i,:),'Color',clrint(i,:));
				hold off
			end
			flabel('Trial','% \Delta',strcat(['Roll ' num2str(roll_window_size) ' % \Delta'...
				'\pi over last \pi']));
		end






		subplot(spr,spc,1);
		scatter(1:num_models,init_p,20.*(1+init_p),pi_colors)
		hold on
		scatter(1:num_models,p,20.*(1+p),pi_colors,'filled')
		hold off
		flabel('j','\pi_j',strcat(['Starting and ending \pi']));

		if isfinite(p_actual)
			y_ax_min_span = 0.5;
			ppdiff=p_actual-p;
			spbi=subplot(spr,spc,5);
			scatter(1:num_models,ppdiff,20.*(1+ppdiff),pi_colors,'filled')
			flabel('j','\Delta',strcat(['\Delta v true']));
			axis(spbi, 'tight')
			currax = axis(spbi);
			if abs(abs(currax(3)) - abs(currax(4))) < y_ax_min_span
				currax(3) = currax(3)-y_ax_min_span;
				currax(4) = currax(4)+y_ax_min_span;
				axis(currax)
			end
			y_ax_min_span = 1;

% 			spbi=subplot(spr,spc,4);
% 			ppdiff=100.*(p_actual-p)./p;
% 			dotweights = 20.*(1+init_p);
% 			sx = 1:num_models;
% 
% 			pct_diff_limit=35;
% 
% 			ndx = ppdiff<pct_diff_limit & ppdiff>-pct_diff_limit;
% 			if sum(ndx)>0
% 				scatter(sx(ndx),ppdiff(ndx),dotweights(ndx),pi_colors(ndx))
% 			end
% 			hold on
% 			ndx = ppdiff>=pct_diff_limit | ppdiff<=-pct_diff_limit;
% 			if sum(ndx)>0
% 				scatter(sx(ndx),pct_diff_limit.*ones(1,sum(ndx)),50.*ones(1,sum(ndx)),pi_colors(ndx),'filled')
% 			end
% 			hold off
% 
% 			axis(spbi, 'tight')
% 			currax = axis(spbi);
% 			if abs(abs(currax(3)) - abs(currax(4))) < y_ax_min_span
% 				currax(3) = currax(3)-y_ax_min_span;
% 				currax(4) = currax(4)+y_ax_min_span;
% 				axis(currax)
% 			end
% 
% 			flabel('j','% \Delta', '% \Delta v true');
		end


		subplot(spr,spc,2);
		plot(ll,'b.-')
		if isfinite(true_loglike)
			hold on
			plot([0 length(ll)], [true_loglike true_loglike], 'r--');
			hold off
		end
		flabel('Trial','l(\pi|x,y)',strcat(['l(\theta)=' num2str(ll(length(ll)))]));
		
% 		if length(ll)>roll_window_size
% 			llchop = ll(length(ll)-roll_window_size:length(ll));
% 			subplot(spr,spc,4);
% 			plot(llchop,'b.-')
% 			if isfinite(true_loglike)
% 				hold on
% 				plot([0 roll_window_size], [true_loglike true_loglike], 'r--');
% 				hold off
% 			end
% 			flabel('Trial','l(\pi|x,y)',strcat([num2str(roll_window_size) ' roll l(\theta)=' num2str(ll(length(ll)))]));
% 		end

	end
