%% load data etc
clear
clc
format long g
addpath('etc/');
addpath('data/');
global im;
im=1;
mo = Model.load('cache/models_01.mat');
ob = Observed(struct('name','halo002'));

%% 
mo = Model.generate(0.01,100,struct('save_to','cache/model_temp.mat'));
mo.test_cache(struct('model_no',1,'step_size',0.005));

mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

%% 
ob.plot();



%% 





%% volume checks
vol=0;
k=1;
vols=[];
dx = mo.data{1}(mo.data{1}(:,3)>0,:)
for i=1:size(dx,1)
	vol_of_this_region = (((.05)^2) * dx(i,3));
	vol = vol + vol_of_this_region;
	vols(k,:) = [dx(i,1) dx(i,2) dx(i,3) vol_of_this_region];
	k=k+1;
end
'x,y,density,volume'
vols
vol


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

%binned on a 5x5 grid
XBIN = 0;
YBIN = 0;

%get data for the appropriate bin
data = models(models(:,4)==XBIN & models(:,5)==YBIN,[1 2 3]);

%points that have mass
non_zero_data = data(data(:,3)>0,:);

%find volume
vol=0;
k=1;
vols=[];
for i=1:size(non_zero_data,1)
	vol_of_this_region = (((.05)^2) * non_zero_data(i,3));
	vol = vol + vol_of_this_region;
	vols(k,:) = [non_zero_data(i,1) non_zero_data(i,2) non_zero_data(i,3) vol_of_this_region];
	k=k+1;
end
'x,y,density,volume'
vols
vol

















