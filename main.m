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
mo = Model.load('cache/models_01.mat');
load('cache/observed.mat');





%% simulate
clc
im=1;

m = length(mo.data);

tic

P = [.015 .2 .005 .07 .12 .005 .005 .01 .27 .1 .05 .004 .001 .02 .095 .03];
if sum(P) ~= 1
	error('P must sum to 1')
end

if m ~= length(P)
	error('P and models are different lengths')
end
PC = cumsum(P);

n = 15000;
data = zeros(n,5);
R1 = rand(n,1);
RX = rand(n,1);
RY = rand(n,1);

grid_size = .1;

for i=1:n
	r0 = rand();
	pi_ndx = sum(PC<=r0)+1;
	
	r1 = R1(i);
	nzd = mo.nonzeros(pi_ndx);
	cnzd = cumsum(nzd(:,3));
	nzda = [nzd cnzd];
	
	grid_ndx = sum(cnzd<=r1)+1;
	
	x = nzda(grid_ndx,1);
	y = nzda(grid_ndx,2);
	
	rx = RX(i)*grid_size;
	ry = RY(i)*grid_size;
	
	x = x + rx;
	y = y + ry;
	
	fab = mo.f_ab(pi_ndx,x,y);
	data(i,:) = [x y fab pi_ndx r0];
end

data = data(data(:,3)>0,:);


fig
subplot(1,3,1);
scatter(data(:,1),data(:,2),data(:,3)./sum(data(:,3)),'k','filled');
flabel('Fe/H','\alpha/Fe',[num2str(n) ' generated data points']);

subplot(1,3,2);
scatter(data(:,1),data(:,2),10,'r','filled');
flabel('Fe/H','\alpha/Fe',[num2str(n) ' generated data points']);

subplot(1,3,3);
plot(P,'k.')
flabel('j','\pi_j','True weights');


save(strcat(['cache/generated_' num2str(n) '_auto.mat']),'data');


toc


ob = Observed(struct('name','generated halo','data',data));
ob.load_models(mo);
ob.save('cache/observed_generated.mat');



'doing em'
[p,all_p,norms,clike] = ob.em(struct(	'Xn',100,...
										'min_norm',0.00001,...
										'max_iters',200,...
										'min_iters',5,...
										'init','rand(m,1)',...
										'interactive',true));
p
'finished em'



fig
subplot(1,2,1)
plot(P,'k.')
hold on
plot(p,'g.')
hold off
flabel('j','\pi_j','Actual (black) v predicted (green)');

subplot(1,2,2)
plot(P'-p,'r.')
flabel('j','Real - actual','Differences');


im=im+1;

'trying rand starts'
[p2,all_p] = ob.em_multi(struct(	'count',5,...
									'Xn',10,...
									'min_norm',0.00001,...
									'max_iters',200,...
									'interactive',true));
p2
'finished rand starts'


%% em
%% 
im=3;
[p,P,norms,clike] = ob.em(struct(		'Xn',100,...
										'min_norm',0.003,...
										'max_iters',50,...
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
p










							
%% em convergence
%%
ob.plot_weight_changes(struct('path','cache/all_p_50_full_runs.mat'));

%%  
mo = Model.load('cache/models_01_normalized.mat');
load('cache/observed_normalized.mat');

%% 
mo = Model.load('cache/models_01.mat');
load('cache/observed.mat');


%% 
[p,all_p] = ob.em_multi(struct(	'count',5,...
								'Xn',10,...
								'min_norm',0.0001,...
								'max_iters',75,...
								'interactive',true));
p









%% plots
%% 
mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

%% 
ob.plot();










%% generation
%% 
mo = Model.generate(0.01,100,struct('normalize',true,'save_to','cache/models_temp.mat'));
mo.test_cache(struct('model_no',1,'step_size',0.005));

%% 
ob = Observed(struct('name','halo002'));
ob.load_models(mo);
ob.save('cache/observed_temp.mat');










