function plot_lines_compare(i,sat,sat2,spf)

	subplot(spf,spf,i);
	
	label_plots = true;


	plot(sat(:,1),sat(:,2),'k.','MarkerSize',2)
	hold on
	plot(sat2(:,1),sat2(:,2),'r.','MarkerSize',2)
	hold off
	axis([-3 0 -1 1])

	
	if label_plots
		flabel('Fe/H','\alpha/Fe',['Sat ' num2str(i) ', n=' num2str(size(sat,1))]);
	else
		set(gca,'ytick',[]) 
		set(gca,'xtick',[]) 
	end
		
% 		if 1==11
% 		min_sz=1;
% 		max_sz=100;
% 		wnorm = min(max(.5.*sat(:,4),min_sz),max_sz);
% 		wnorm = max_sz.*wnorm./range(wnorm);
% 		scatter(sat(:,1),sat(:,2),wnorm,[0 0 0],'filled')
% 		flabel('Fe/H','\alpha/Fe',['Sat ' num2str(i) ', n=' num2str(size(sat,1))]);
% 		end
	end