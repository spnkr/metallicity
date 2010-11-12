function sss=plot_correl_dots(mi,correl,size_is_correl,color_neg)
			global im;
			sss=figure(im);
			im=im+1;
			clf(sss);
			
			adjust_to_baseline=false;
			overlay=false;
			
			cor = abs(correl);
			for i=1:mi.num_models
				for j=1:mi.num_models
					if i==j
						cor(i,j) = 0;
					end
				end
			end

			
			if adjust_to_baseline
			baseline = 1-max(max(cor));
			cor = cor + baseline;
			
			for i=1:mi.num_models
				for j=1:mi.num_models
					if i==j
						cor(i,j) = 0;
					end
				end
			end
			end
			
			hold on
			for i=1:mi.num_models
				for j=1:mi.num_models
					[x,y,z] = sphere(20);
					
					if size_is_correl
						pltitle = 'Correlation of \pi';
						sz = (cor(i,j).^-2).^-1;
						sz = cor(i,j);
						sz = min(max(sz,0.01),.5);
						if i==j
							sz=0.000001;
						end


						if correl(i,j)<0 && color_neg
							clr=[191/255 48/255 0/255];
						else
							clr=[0/255 98/255 223/255];
						end
						
						opac=1;%mi.pi_est(i)+mi.pi_est(j);
					else
						pltitle = 'Correl of \pi. Size of sphere = \pi_p+\pi_q; opacity = correlation';
						sz = mi.pi_est(i)+mi.pi_est(j);
						sz = min(max(sz,0.1),1);
						if i==j
							sz=0.00001;
						end


						if correl(i,j)<0 && color_neg
							clr=[191/255 48/255 0/255];
						else
							clr=[0/255 98/255 223/255];
						end
						
						opac=cor(i,j);
					end
					
					eopac = max(opac./2,0);
					eclr=min(clr.*1.5,1);
					
					surf(sz*x+1*i,sz*y+1*j,sz.*z,'FaceColor',clr,'EdgeColor',eclr,'FaceAlpha',opac,'EdgeAlpha',eopac)
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
			
			m=mi.num_models;
			for i=1:m
				plot([0 m+.5], [i+.5 i+.5], 'Color', [.5 .5 .5])
				plot([i+.5 i+.5], [0 m+.5], 'Color', [.5 .5 .5])
			end

			hold off
			flabel('\pi_p','\pi_q',pltitle);
			daspect([1 1 1])
			axis([.5 m+.5 .5 m+.5])
			view(90,90)
% 			axis([-1 17 -1 17 -2 5])
		end