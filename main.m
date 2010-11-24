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

motoload='cache/mo_20k.mat';
%motoload='cache/mo_all.mat';
load(motoload,'mo')


%% 
im=1;
mo.plot_lines(113,true);


%% 
im=2;
mo.density(2,.9)




%% 
data = mo.data(mo.data(:,6)==1,:);

%% 
tic
rdata = zeros(43482479,2);
k=1;
m=size(data,1);
for i=1:m
	w = ceil(100*data(i,3));
 	for j=1:w
 		rdata(k,:) = data(i,1:2);
 		k=k+1;
	end
	i
end
toc


%% 

rdata = [];
k=1;
m=size(data,1);
for i=1:m
	w = ceil(100*data(i,3));
 	rdata(k)=w;
 	k=k+1;
end
sum(rdata)

%% 
load('cache/temp/rdata_4m.mat','rdata')

%% 
b=0.9;
nsat=1;
sat=data;

figure(1)
[bandwidth,density,X,Y]=kernel_smooth(data(:,1:2),NaN);%,data(:,3));
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
mo = Model(NaN,0:11); %to 11
size(mo.data)
save('cache/mo_all.mat','mo')


sepr('done')




































