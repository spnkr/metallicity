classdef Observed < handle
	
	properties
		dpath='data/obsdata2.dat';
		name;
		data;
		x,y,lum;
		f_ab;
	end
	
	
	methods(Static)
		
	end
	
	
	methods
		%loading
		function ob = Observed(varargin)
			load_args
			dpath = arg('path',ob.dpath);
			ob.name = arg('name','halo000');
			
			ob.data = load(dpath);
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
					f_ab(i,j) = mo.f_ab(j,ob.x(i),ob.y(i));
				end
			end
			
			ob.f_ab = f_ab;
		end
		
		function save(ob,path)
			save(path,'ob');
			show(['Saved to ' path])
		end
		
		
		
		
		
		
		
		
		
		%em
		function [p,P,norms,plike,clike] = em(ob, varargin)
			load_args
			
			max_iters = arg('max_iters',20);
			min_iters = arg('min_iters',1);
			n = arg('n',length(ob.x));
			min_norm = arg('min_norm',0.01);
			init_str = arg('init','rand(m,1)');
			interactive = arg('interactive',true);
			
			m = size(ob.f_ab,2);

			p = eval(init_str);
			p(p==0) = max(p);
			p = p./sum(p);
			
			P = NaN.*ones(m,max_iters);
			norms = zeros(1,1);
			plike = zeros(1,1);
			clike = zeros(1,1);
			
			w = zeros(n,m);
			
			global im;
			figure(im);
			
			tic;
			P(:,1) = p;
			for counter=1:max_iters
				p0 = p;
				
				if interactive
					ob.plot_progress(norms,plike,clike,P,p,n,m,counter,-counter,init_str,im);
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
				
				plike(counter) = ob.partial_likelihood(p,n,m);
				clike(counter) = ob.complete_likelihood(p,w,n,m);
				
				if norms(counter) < min_norm && counter > min_iters
					%disp('min norm reached; stopping')
					break;
				end
			end
			tmr = toc;
			
			if counter==max_iters
				warning('min norm not reached!')
			end
			
			ob.plot_progress(norms,plike,clike,P,p,n,m,counter,tmr,init_str,im);
			
			im=im+1;
			
			
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
		function l = complete_likelihood(ob,p,w,n,m)
			l = 0;
			for j=1:m
				for i=1:n
					l0 = w(i,j) * log(p(j)*ob.f_ab(i,j));
					if isfinite(l0) && ~isinf(l0)
						l = l+l0;
% 					else
% 						disp(strcat(['error in clike i=' num2str(i) ',j=' num2str(j)]));
					end
				end
			end
		end
		
		function plot_progress(ob,norms,plike,clike,P,p,n,m,counter,tmr,init_str,im)
			figure(im);
			spr=3;spc=2;
			subplot(spr,spc,1);
			plot(norms,'k.-');
			legend(strcat(['norm=' num2str(min(norms))]));
			
			if tmr<0
				flabel('Trial','',strcat(['(In progress) ' num2str(counter) ' runs, n='...
					num2str(n)]));
			else
				flabel('Trial','',strcat([num2str(counter) ' runs in ' ...
					num2str(tmr) 's, n=' num2str(n)]));
			end

			clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
				'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
			
			subplot(spr,spc,2);
			plot(1:size(P,2),P(1,:),strcat([clrs(1) '.-']));
			for i=2:m
				hold on
				plot(1:size(P,2),P(i,:),strcat([clrs(i) '.-']));
				hold off
			end
			flabel('Trial','\pi_j',strcat(['Weights over time']));

			subplot(spr,spc,3);
			scatter(1:25,P(:,1),20.*(1+P(:,1)),10.*(1+P(:,1)))
			flabel('j','\pi_j','Starting weights');
			
			subplot(spr,spc,4);
			scatter(1:25,p,20.*(1+p),10.*(1+p),'filled')
			flabel('j','\pi_j',strcat(['Weights. \pi^{(0)}=' init_str]));
			
			
			subplot(spr,spc,5);
			plot(plike,'b.-')
			flabel('Trial','l(\pi|x,y)',strcat(['Partial Log Likelihood=' num2str(plike(length(plike)))]));
			
			subplot(spr,spc,6);
			plot(clike,'m.-')
			flabel('Trial','l(\pi|x,y,z)',strcat(['Complete Log Likelihood=' num2str(clike(length(clike)))]));
			
			
			
% 			subplot(2,2,4);
% 			scatter(1:25,log(p),20.*(1+p),10.*(1+p),'filled')
% 			flabel('j','log(\pi_j)','Log weights');
		end
		
		function [p, all_p] = em_multi(ob, varargin)
			load_args
			global im;
			xim=im;
			
			count = arg('count',1);
			save_path = arg('save',strcat(['cache/run_auto_' num2str(count) '.mat']));
			
			all_p = zeros(25,count);
			for i=1:count
				im=xim;
				[p,P,norms,plike,clike] = ob.em(cell2mat(varargin));
				all_p(:,i) = p;
				disp(strcat(['Finished run ' num2str(i)]))
			end
			
			save(save_path,'all_p');
			disp(strcat(['Saved all_p to ' save_path]));
			
			ob.plot_weight_changes(struct('all_p',all_p));
		end
		
		
		
		
		
		
		
		
		
		
		%convergence
		function ap = load_em_pis(ob,varargin)
			load_args
			
			pth = arg('path','cache/all_p_50_full_runs.mat');
			
			load(pth);
			ap = all_p;
		end
		
		function plot_weight_changes(ob,varargin)
			load_args
			
			all_p = arg('all_p',NaN);
			
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

			no_std_devs = 2;

			subplot(1,3,1);
			hold on
			for i=1:m
				plot(all_p(i,:),'.','Color',clrs(i,:))

				plot(1:p,mean(all_p(i,:)).*ones(1,p),'-','Color',clrs(i,:))

				plot(1:p,(no_std_devs*std(all_p(1,:))+mean(all_p(i,:))).*ones(1,p),'--','Color',clrs(i,:))
				plot(1:p,(-no_std_devs*std(all_p(1,:))+mean(all_p(i,:))).*ones(1,p),'--','Color',clrs(i,:))
			end


			hold off
			flabel('Trial','\pi_j',[num2str(size(all_p,2)) ' random starting values, +/- 2\sigma, n=all']);


			subplot(1,3,2); 
			vls = zeros(m,5);
			hold on
			for i=1:m
				dv=all_p(i,:);
				mdv = mean(dv);
				sdv = std(dv);
			% 	plot(i,vls(i,:),'.','Color',clrs(i,:))
				errorbar(i,mdv,sdv,'.','Color',clrs(i,:))
			end
			hold off
			flabel('j','\pi_j','Error bars: +/- \sigma');


			subplot(1,3,3); 
			vls = zeros(m,5);
			hold on
			for i=1:m
				dv=all_p(i,:);
				mdv = mean(dv);
				sdv = std(dv);
			% 	plot(i,vls(i,:),'.','Color',clrs(i,:))
				errorbar(i,mdv,2*sdv,'.','Color',clrs(i,:))
			end
			hold off
			flabel('j','\pi_j','Error bars: +/- 2\sigma');
			%plot(1:size(vls,1),vls,'.','Color',clrs(i,:))
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