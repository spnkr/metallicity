function plot_lines(i,sat,spf)

	subplot(spf,spf,i);
	
	if 1==1
		

		plot(sat(:,1),sat(:,2),'k.')
		flabel('Fe/H','\alpha/Fe',['Sat ' num2str(i) ', n=' num2str(size(sat,1))]);

		if 1==11
		min_sz=1;
		max_sz=100;
		wnorm = min(max(.5.*sat(:,4),min_sz),max_sz);
		wnorm = max_sz.*wnorm./range(wnorm);
		scatter(sat(:,1),sat(:,2),wnorm,[0 0 0],'filled')
		flabel('Fe/H','\alpha/Fe',['Sat ' num2str(i) ', n=' num2str(size(sat,1))]);
		end
	end