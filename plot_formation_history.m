function plot_formation_history(mi)
	h=figure(4);
	clf(h);
	hold on;
	m=ceil(sqrt(length(mi.model_skip_ndx)));
	for i=1:m
		plot([0 m+.5], [i+.5 i+.5], 'Color', [.5 .5 .5])
		plot([i+.5 i+.5], [0 m+.5], 'Color', [.5 .5 .5])
	end
	
% 	set(gcf,'Units','normalized');
	
	k=1;
	mm=1;
	clr=[0/255 98/255 223/255];
	eclr=min(clr.*1.5,1);
% 	for i=1:m
% 		for j=1:m
% 			aaaa=k;
% 			aaaa=mm;
% 			aaaa=11;
% 			if length(mi.pi_est)>=mm
% 				if sum(mi.model_skip_ndx==k)==0
% 					[x,y,z] = sphere(20);
% 					sz = mi.pi_est(mm);
% 					sz = min(max(sz,0.01),.8);
% 
% 					surf(sz*x+1*i,sz*y+1*j,sz.*z,'FaceColor',clr,...
% 						'EdgeColor',eclr,'FaceAlpha',1,'EdgeAlpha',.5);
% 
% 
% 					text((i-1)+.6,(j-1)+1.11,1,...
% 						strcat([num2str(round(mi.pi_est(mm)*10000)/100) '%']));
% 					mm=mm+1;
% 				end
% 				k=k+1;
% 			end
% 		end
% 	end

	msk = .9.*ones(m,m);
	X=1;
	Y=1;
	for i=1:length(mi.model_skip_ndx)
		msk(X,Y) = mi.model_skip_ndx(i);
		[i X Y];
		
		Y = Y+1;
		if Y>m
			Y=1;
			X=X+1;
		end
		
	end
	
	xf=0;
	yf=1;
	for i=1:length(mi.model_skip_ndx)
		xf = xf + 1;
		if xf>m
			xf = 1;
			yf = yf + 1;
		end

		if msk(xf,yf)==1
			[x,y,z] = sphere(20);
			sz = mi.pi_est(k);
			sz = min(max(sz,0.01),.8);
			
			
			
			surf(sz*x+1*xf,sz*y+1*yf,sz.*z,'FaceColor',clr,...
				'EdgeColor',eclr,'FaceAlpha',1,'EdgeAlpha',0);


			text(xf+.05,yf-.35,1.1,...
				strcat([num2str(round(mi.pi_est(k)*10000)/100) '%']));
			k=k+1;
		end
	end
	
	hold off
	axis([.5 m+.5 .5 m+.5])
	flabel('Accretion time (Gyr)', 'Mass', ['Formation history - ' mi.filename])
	set(gca,'XTick',0);
	set(gca,'YTick',0);
	view(0,-90)


