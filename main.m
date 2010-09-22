%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;


%% 
obs = load('data/obsdata2.dat');
x = obs(:,1);
y = obs(:,2);
w = obs(:,3);


%% 
models = load('data/modeldata2.dat');
x = models(:,1);
y = models(:,2);
freq = models(:,3);
xbin = models(:,4);
ybin = models(:,5);
lta = models(:,6);
eta = models(:,7);
hma = models(:,8);
lma = models(:,9);

ndx = 1:length(x);

[x(ndx) y(ndx) freq(ndx) xbin(ndx) ybin(ndx) lta(ndx) eta(ndx) hma(ndx) lma(ndx)];




%% 
im=1;
figure(1);
p=1;
for i=min(xbin):max(xbin)
	for j=min(ybin):max(ybin)
		subplot(5,5,p);
		p=p+1;
		ndx = models(:,4)==i & models(:,5)==j & models(:,3)>0;
		mn = min(freq(ndx));
		mx = max(freq(ndx));
		scf = 5 + 10.*(freq(ndx) - mn) ./ (mx-mn);
		scatter(x(ndx),y(ndx),scf,'filled','k');
		title(strcat(['T:' num2str(min(lta(ndx))) '/' num2str(max(lta(ndx))) '; M:' num2str(min(lma(ndx))) '/' num2str(max(hma(ndx)))]))
	end
end



%% 




%% 
scatter_color_by_weight(x,y,w,struct(...
	'title', 'Observed data (blue->green:bright->dim)', 'x','\alpha/Fe', 'y','Fe/H','vary_size',false,'dot_size',3));



