%%
matlabpool open
%% 
matlabpool close

%% load data etc
clear
clc
format long g
addpath('etc/');
addpath('data/');
global im;
im=1;
mo = Model.load('cache/models_01.mat');
load('cache/observed.mat');


%% generation
%% 
mo = Model.generate(0.01,100,struct('normalize',false,'save_to','cache/models_temp.mat'));
mo.test_cache(struct('model_no',1,'step_size',0.005));

%% 
ob = Observed(struct('name','halo002'));
ob.load_models(mo);
ob.save('cache/observed_temp.mat');





%% plots
%% 
mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

%% 
ob.plot();



%% em
%% 
[p,P,norms,plike,clike] = ob.em(struct(	'n',100,...
										'min_norm',0.005,...
										'max_iters',25,...
										'min_iters',5,...
										'init','rand(m,1)',...
										'interactive',true));
p

%% 
[p,all_p] = ob.em_multi(struct(	'count',5,...
								'n',10,...
								'min_norm',0.1,...
								'max_iters',75,...
								'interactive',true));
							

							
%% em convergence
%%
ob.plot_weight_changes(struct('path','cache/all_p_50_full_runs.mat'));



%% 
mo = Model.load('cache/models_01_normalized.mat');
load('cache/observed_normalized.mat');

%% 
mo = Model.load('cache/models_01.mat');
load('cache/observed.mat');


try_this_many_times=50;
all_p = zeros(25,try_this_many_times);
for i=1:try_this_many_times
	im=1;
	[p,P,norms,plike,clike] = ob.em(struct(	'Xn',NaN,...
							'min_norm',0.001,...
							'max_iters',75,...
							'Xmin_iters',10,...
							'init','rand(m,1)',...
							'interactive',false));
	all_p(:,i) = p;
	disp(strcat(['Finished run ' num2str(i)]))
end
all_p
sepr('finished all')

%% 
im=1;
[p,P,norms,plike,clike] = ob.em(struct(	'XXn',NaN,...
						'min_norm',0.00001,...
						'max_iters',50,...
						'init','rand(m,1)',...
						'interactive',true));




					
					

%% volume checks
%% 
for i=1:25
	mo.volume(i)
end




%% 
for i=1:25
	x=mo.nonzeros(i);sum(x(:,3))
end



%% 
vol=0;
k=1;
vols=[];
dx = mo.nonzeros(k);
grid_size = .05;

if 1==1
dx = mo.cache{k};
dx = dx(dx(:,1)>0,[2 3 1]);
grid_size = mo.stepsize;
end

for i=1:size(dx,1)
	vol_of_this_region = (((grid_size)^2) * dx(i,3));
	vol = vol + vol_of_this_region;
	vols(k,:) = [dx(i,1) dx(i,2) dx(i,3) vol_of_this_region];
	k=k+1;
end
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

















