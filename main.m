%%
matlabpool open
%% 
matlabpool close

%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;

%1 Fe/H
%2 alpha/Fe
%3 weight
%4 tacc
%5 lsat
%6 nsat  

load('cache/large/mo_40k.mat','mo')


%% 
im=1;
mo.plot_lines(10)



%% 
im=2;
mo.density(11,.9)




















%% 
clc
im=1;

ms=mo.sats(2);
spf=ceil(sqrt(ms));

spf=1;
ms=2;
for i=1:ms
	sat = mo.data(mo.data(:,6)==i,:);
	%plot_lines(i,sat,spf);
	
	%sat = sat(1:50,:);
	data=sat(:,[1,2]);
	
	if 1==1
		fg=figure(1);
		clf(fg);
		b=.9;
		

		[bandwidth,density,X,Y]=kernel_smooth(data,b);
		bandwidth
		subplot(1,2,1);
		contour3(X,Y,density,50), hold on
		plot(data(:,1),data(:,2),'r.','MarkerSize',5)
		flabel('Fe/H','\alpha/Fe',['Sat ' num2str(i) ', n=' num2str(size(sat,1)) ... 
			', b=[' num2str(bandwidth(1)) ' ' num2str(bandwidth(2)) ']']);

		hold off

		subplot(1,2,2);
		surf(X,Y,density,'LineStyle','none'), view([0,60])
		colormap hot, hold on, alpha(.8)
		set(gca, 'color', 'blue');
		plot(data(:,1),data(:,2),'w.','MarkerSize',5)
		flabel('Fe/H','\alpha/Fe',['Sat ' num2str(i) ', n=' num2str(size(sat,1)) ... 
			', b=[' num2str(bandwidth(1)) ' ' num2str(bandwidth(2)) ']']);

		hold off
	end
	
	if 1==11
		fg=figure(2);
		clf(fg);
		[bandwidth,density,X,Y]=kde2d(data);
		bandwidth
		subplot(1,2,1);
		contour3(X,Y,density,50), hold on
		plot(data(:,1),data(:,2),'r.','MarkerSize',5)
		hold off

		subplot(1,2,2);
		surf(X,Y,density,'LineStyle','none'), view([0,60])
		colormap hot, hold on, alpha(.8)
		set(gca, 'color', 'blue');
		plot(data(:,1),data(:,2),'w.','MarkerSize',5)
		hold off
	end
end


%% 













%% generate
%   Fe/H  alpha/Fe  weight  tacc  lsat  nsat  
mo = Model(NaN,0:1); %to 11
size(mo.data)
save('cache/large/mo_40k.mat','mo')







































