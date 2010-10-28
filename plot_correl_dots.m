function plot_correl_dots(mi,dot_scale,overlay)
			global im;
			sss=figure(im);
			im=im+1;
			clf(sss);

			hold on
			cv = mi.correl;
			cv=50.*cv./range(range(cv));
			for i=1:mi.num_models
				for j=1:mi.num_models
					[x,y,z] = sphere(20);
					fact=cv(i,j);
					if fact<0
						fc=[223/255 128/255 0/255];
					else
						fc=[39/255 127/255 0/255];
					end
					fact=max(abs(fact),2);
					fact = fact/dot_scale(1);
					
					fact = max(fact,.1);
					fact = fact./dot_scale(2);
					
					if i==j
						fact=0.00001;
					end
					
					ad = mi.pi_est(i)+mi.pi_est(j);
					ad = min(max(ad,0.1),1);
					ead = max(ad-.2,0);
					eclr=min(fc.*1.5,1);
					surf(fact*x+1*i,fact*y+1*j,fact*z,'FaceColor',fc,'EdgeColor',eclr,'FaceAlpha',ad,'EdgeAlpha',ead)
				end
			end

			if 1==11
			[x,y]=meshgrid(1:mi.num_models+1,1:mi.num_models+1);
			z=[];
			m=length(x);
			for i=1:m
				for j=1:m
					if j==m
						ej=mi.pi_est(j-1);
					else
						ej=mi.pi_est(j);
					end
					if i==m
						ei=mi.pi_est(i-1);
					else
						ei=mi.pi_est(i);
					end
					z(i,j) = ei + ej;
				end
			end
			colormap gray
			surf(x-.5,y-.5,-1*ones(length(x),length(x)),'CData',1./max(z,overlay),'EdgeColor','none','FaceColor','flat',...
				'FaceAlpha',0.5)
			end

			hold off
			flabel('\pi_p','\pi_q','Correl of \pi = area of dot; darker background = more \pi_p+\pi_q');
			daspect([1 1 1])
			view(90,90)
			axis([-1 17 -1 17 -2 5])
		end