classdef Model < handle
	
	properties
		data;
		header='  Fe/H  alpha/Fe  weight  tacc  lsat  nsat  halo_num';
		
		sats;
		n;
	end
	
	
	methods(Static)
		%sat_range=0 to 11
		function mo = generate(n,sat_range,name)
			%   Fe/H  alpha/Fe  weight  tacc  lsat  nsat  
			mo = Model(n,sat_range); 
			save(strcat(['cache/' name '.mat']),'mo')
			sepr(strcat(['saved to cache/' name '.mat']))
		end
	end
	
	
	methods
		%--loading
		function mo = Model(n,halos)
			data = zeros(1,7);
			parts = 0:22;
			k=1;
			for i=1:length(halos)
				for j=1:length(parts)
					fpl = strcat(['data/halodataall/halodataall/datalist_halo' num2str(halos(i)) 'part.' num2str(parts(j)) '.dat']);
					if exist(fpl) > 0
						disp(strcat(['loading ', fpl]))
						if true
							[h,nd] = hdrload(fpl);

							try
								b=size(data,1)+1
								c=b+size(nd,1)-1
								data(b:c,:) = [nd halos(i)*ones(size(nd,1),1)];
							catch
								disp('load error')
							end
						end
						if isfinite(n) && size(data,1)>=n
							break;
						end
					end
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
		
		function plot_lines(mo,sndx)
			spf=ceil(sqrt(sndx));
			fig;
			
			for i=1:sndx
				sat = mo.data(mo.data(:,6)==i,:);
				plot_lines(i,sat,spf);
			end
		end
		
		function plot_lines_compare(mo,sndx)
			spf=ceil(sqrt(sndx-1));
			fig;
			
			for i=1:sndx-1
				sat = mo.data(mo.data(:,6)==i,:);
				sat2 = mo.data(mo.data(:,6)==i+1,:);
				plot_lines_compare(i,sat,sat2,spf);
			end
		end
		
	end
	
end