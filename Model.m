classdef Model < handle
	
	properties
		data;
		header='  Fe/H  alpha/Fe  weight  tacc  lsat  nsat  ';
		
		sats;
		n;
	end
	
	
	methods(Static)
		%sat_range=0 to 11
		function mo = generate(sat_range,name)
			%   Fe/H  alpha/Fe  weight  tacc  lsat  nsat  
			mo = Model(NaN,sat_range); 
			save(strcat(['cache/' name '.mat']),'mo')
		end
	end
	
	
	methods
		%--loading
		function mo = Model(n,parts)
			data = zeros(1,6);
			k=1;
			for i=1:length(parts)
				nd = load(strcat(['data/large/datalist_halo2part.' num2str(parts(i)) '.dat']));

				try
					b=size(data,1)+1
					c=b+size(nd,1)-1
					data(b:c,:) = nd;
				catch
					disp('load error')
				end
				if isfinite(n) && size(data,1)>=n
					break;
				end
			end
			
			if isfinite(n) && size(data,1)>=n
				data = data(1:n,:);
			end
			
			data = data(data(:,1)>-3,:);
			
			mo.data = data;
			mo.sats = [min(mo.data(:,6)) max(mo.data(:,6))];
			mo.n = size(data,1);
		end
		
		
		function f=density(mo,nsat,b)
			sat = mo.data(mo.data(:,6)==nsat,:);
			%plot_lines(i,sat,spf);

			%sat = sat(1:50,:);
			data=sat(:,1:3);

			global im;
			fg=figure(im);
			im=im+1;
			clf(fg);
			
			rdata = data(:,1:2);
			
			
			[bandwidth,density,X,Y]=kernel_smooth(rdata,b);%,data(:,3));
			bandwidth
			subplot(1,2,1);
			contour3(X,Y,density,50), hold on
			plot(data(:,1),data(:,2),'r.','MarkerSize',5)
			flabel('Fe/H','\alpha/Fe',['Sat ' num2str(nsat) ', n=' num2str(size(sat,1)) ... 
				', b=[' num2str(bandwidth(1)) ' ' num2str(bandwidth(2)) ']']);

			hold off
			view(-18,44);

			subplot(1,2,2);
			surf(X,Y,density,'LineStyle','none'), view([0,60])
			colormap hot, hold on, alpha(.8)
			set(gca, 'color', 'blue');
			plot(data(:,1),data(:,2),'w.','MarkerSize',5)
			flabel('Fe/H','\alpha/Fe',['Sat ' num2str(nsat) ', n=' num2str(size(sat,1)) ... 
				', b=[' num2str(bandwidth(1)) ' ' num2str(bandwidth(2)) ']']);

			hold off
			
			axis([min(min(X)) max(max(X)) min(min(Y)) ...
			max(max(Y)) min(min(density))-std(std(density)) max(max(density))+std(std(density))])
		
			view(0,90);
		end
		
		function plot_lines(mo,sndx,hide_labels)
			spf=ceil(sqrt(sndx));
			fig;
			
			for i=1:sndx
				sat = mo.data(mo.data(:,6)==i,:);
				plot_lines(i,sat,spf,hide_labels);
			end
		end
	end
	
end