clear
clc
format short g
addpath('etc/');
addpath('data/');

global im;
im=1;
global constant_im;
constant_im=false;

all_halos = [2 5 7 8 9 10 12 14 15 17 20];

mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load(mnames{1});




%% 
[h,d] = hdrload('data/halodata_v1/obsdata_s4t12m77.dat');

%% 
figure(1);
plot(d(:,1),d(:,2),'k.')


%% 
load(strcat(['cache/' 'large/mo_h25']));

%% 
im=1;
mo.plot_lines(1:25)


%% 
im=1;
h=figure(im);
clf(h);
sndx = [2:11];
spf=ceil(sqrt(length(sndx)));
hold on
for i=1:length(sndx)
	sat = mo.data(mo.data(:,6)==sndx(i),:);
% 	sat = sat(1:500,:);
% 	plot_lines(i,sat,spf);
% subplot(spf,spf,i);
	

w=15;


 	%plot(sat(:,1),sat(:,2),'k.','MarkerSize',2)
% 
 colormap summer
 scatter(sat(:,1),sat(:,2),w,ones(1,3).*(i/10),'filled')
	axis([-3 0 -1 1])

	set(gca,'ytick',[]) 
	set(gca,'xtick',[]) 
end
hold off





