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
clear
clc
format short g
addpath('etc/');
addpath('data/');
load('cache/temp/new2')

data = dat;
grid_size = 0.1;
noise = 0.05;


minx = min(data(:,1));
if minx<-3
	data = data(data(:,1)>=-3,:);
	minx=-3;
end

X = size(data,1);
data(:,1:2) = data(:,1:2) + normrnd(0,noise,X,2);



gridt = [0 1 2 inf];
gridm = [0 1e3 1e6 1e7 1e8 1e9 1e10 inf];
%n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). The last bin
%counts any values of x that match edges(end). Values outside the values in
%edges are not counted.
[nt,tbin] = histc(data(:,4),gridt);
[nm,mbin] = histc(data(:,5),gridm);


data = [data(:,1:6) tbin mbin];
data = data(data(:,7)>0,:);
data = data(data(:,8)>0,:);
X = size(data,1);

tlook = [nt gridt'];
mlook = [nm gridm'];
freq = data(:,3).*data(:,5);
%feh,afe,nfreq,time bin ndx, mass bin ndx, min time, max time, min mass,
%max mass (min/max mass/time are both the min at the moment
data = [data(:,1:2) freq data(:,7) data(:,8)...
	tlook(data(:,7),2) tlook(data(:,7),2) mlook(data(:,8),2) mlook(data(:,8),2)];


%% 
tic
gdata = [];

for i=1:length(gridx)
	for j=1:length(gridy)
		gx = gridx(i)-grid_size/2;
		gy = gridy(j)-grid_size/2;
		mrows=data(data(:,9)==gx & data(:,10)==gy,:);
		if size(mrows,1)>0
			cr = mrows(1,:);
			nfreq = sum(mrows(:,3).*mrows(:,5));
			gdata(size(gdata,1)+1,:) = [gx gy nfreq i j cr(3:6)];
		else
			gdata(size(gdata,1)+1,:) = [gx gy 0 i j 0 0 0 0];
		end
	end
end
toc



%% 



















%% 
clear
clc
format short g
addpath('etc/');
addpath('data/');
load('cache/temp/new2')

data = dat;
grid_size = 0.1;
noise = 0.05;


minx = min(data(:,1));
minx = round(minx*10)/10;
if minx<-3
	data = data(data(:,1)>=-3,:);
	minx=-3;
end
maxx = max(data(:,1));
maxx = round(maxx*10)/10;
miny = min(data(:,2));
miny = round(miny*10)/10;
maxy = max(data(:,2));
maxy = round(maxy*10)/10;
gridx = minx:grid_size:maxx;
gridy = miny:grid_size:maxy;

X = size(data,1);
data(:,1:2) = data(:,1:2) + normrnd(0,noise,X,2);

%n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). The last bin
%counts any values of x that match edges(end). Values outside the values in
%edges are not counted.
[nx,xbin] = histc(data(:,1),gridx);
[ny,ybin] = histc(data(:,2),gridy);


data = [data(:,1:6) xbin ybin];
data = data(data(:,7)>0,:);
data = data(data(:,8)>0,:);

xlook = [nx gridx'];
ylook = [ny gridy'];
data = [data(:,1:8) xlook(data(:,7),2)-grid_size/2 ylook(data(:,8),2)-grid_size/2];

tic
gdata = [];

for i=1:length(gridx)
	for j=1:length(gridy)
		gx = gridx(i)-grid_size/2;
		gy = gridy(j)-grid_size/2;
		mrows=data(data(:,9)==gx & data(:,10)==gy,:);
		if size(mrows,1)>0
			cr = mrows(1,:);
			nfreq = sum(mrows(:,3).*mrows(:,5));
			gdata(size(gdata,1)+1,:) = [gx gy nfreq i j cr(3:6)];
		else
			gdata(size(gdata,1)+1,:) = [gx gy 0 i j 0 0 0 0];
		end
	end
end
toc


%% 


%% 
%weight x lsat
nfreq = data(:,3).*data(:,5);


gdata = [data(:,1:2) nfreq xbin ybin];

%% 





%% 
dat = Halo.read(struct('halos',[5],'xxdata_pont_limit',2000000));




%% 



%% demos
[gdata,noise,grid_size] = Halo.generate_observed_data(struct('n',100,'noise',0.15,'grid_size',0.1,...
	'halos',[2],'data_point_limit',300000));

[dat,no] = Halo.generate_model_data(struct('noise',0.15,'grid_size',0.1,'halos',[2],...
	'data_point_limit',300000));

dat = Halo.read(struct('halos',[2],'data_point_limit',300000));




