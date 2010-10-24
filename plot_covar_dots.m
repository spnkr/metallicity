function plot_covar_dots(mi)
			global im;
			sss=figure(im);
			im=im+1;
			clf(sss);

			hold on
			cv = mi.covar;
			cv=50.*cv./range(range(cv));
			for i=1:mi.num_models
				for j=1:mi.num_models
					[x,y,z] = sphere(20);
					fact=cv(i,j);
					if fact<0
						fc=[0.7 0.3 .3];
					else
						fc=[0.6 0.7 .9];
					end
					fact=max(abs(fact),2);
					fact = fact/50;
					ad = mi.pi_est(i)+mi.pi_est(j);
					ad = max(ad,0.1);
					ead = max(ad-.2,0);
					surf(fact*x+1*i,fact*y+1*j,fact*z,'FaceColor',fc,'EdgeColor','none','FaceAlpha',ad,'EdgeAlpha',ead)
				end
			end

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
			surf(x-.5,y-.5,-1*ones(length(x),length(x)),'CData',1./max(z,0.05),'EdgeColor','none','FaceColor','flat',...
				'FaceAlpha',0.5)


			hold off
			flabel('\pi_j','\pi_j','Covar of \pi. Transparency=sum weights, size=covar, blue=+,green=-');
			daspect([1 1 1])
			view(90,90)
			axis([-1 17 -1 17 -2 5])
		end