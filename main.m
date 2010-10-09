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

[ob, mo] = Observed.load(3,10);





%% 
[p,P,ll,LL] = ob.simulate(mo,struct('sample',10000,...
											'max_seconds',60*30,...
											'interactive',false));



 
multicount = 10;
max_seconds = 60*5;

multicount = 2;
max_seconds = 60*5;

sepr('init with real weights plus or minus 10pct random fluctuations to each')
[best_ll, best_p, all_p, all_ll] = ob.em_multi(struct(	'count',multicount,...
								'Xn',10,...
								'save','cache/p_ll_run_weight_pl_10.mat',...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',ob.p_actual.*0.1.*rand(size(ob.p_actual,1),1)));
print_pi(best_p,mo)


sepr('init with rand weights')
[best_ll, best_p, all_p, all_ll] = ob.em_multi(struct(	'count',multicount,...
								'Xn',10,...
								'save','cache/p_ll_run_weight_rand.mat',...
								'max_seconds', max_seconds, ...
								'interactive',false));
print_pi(best_p,mo)


sepr('init with real weights plus or minus 30pct random fluctuations to each')
[best_ll, best_p, all_p, all_ll] = ob.em_multi(struct(	'count',multicount,...
								'Xn',10,...
								'save','cache/p_ll_run_weight_pl_30.mat',...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',ob.p_actual.*0.3.*rand(size(ob.p_actual,1),1)));
print_pi(best_p,mo)



sepr('init with real weights plus or minus 3pct random fluctuations to each')
[best_ll, best_p, all_p, all_ll] = ob.em_multi(struct(	'count',floor(multicount/4),...
								'Xn',10,...
								'save','cache/p_ll_run_weight_pl_03.mat',...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',ob.p_actual.*0.03.*rand(size(ob.p_actual,1),1)));
print_pi(best_p,mo)
















%% 
for i=1:4
	disp(strcat(['Run ' num2str(i)]))
	[p,P,ll] = ob.em(struct(...
					'max_seconds',15,...
					'init','rand(m,1)',...
					'interactive',true));
	p

	
	im=im+1;
	figure(im);
	im=im+1;
	ll_r = ll(1:2);
	plot(ll_r,'.-')
	legend(strcat(['Final log like: ' num2str(ll_r(length(ll_r)))]), 'Location', 'SouthEast')
	flabel('Trial','Log like', ['Complete log like over 1st ' num2str(length(ll_r)) ' runs (' num2str(i) ')']);

	figure(im);
	im=im+1;
	plot(ll,'.-')
	legend(strcat(['Final log like: ' num2str(ll(length(ll)))]), 'Location', 'SouthEast')
	flabel('Trial','Log like', ['Complete log like over time (' num2str(i) ')']);
	 
end



							
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





%% 
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;

[ob, mo] = Observed.load(3,10);





%% get real log like
p_actual = [	14.467074  ;
				7.4354783  ;
				20.991296  ;
				9.2340355  ;
				33.655754  ;
				8.0191336  ;
				2.4001462  ;
				2.5734037  ;
			   0.11883542  ;
			  0.076990272  ;
			   0.35844400  ;
			   0.25313549  ;
			   0.22529346  ;
			  0.048890239  ;
			  0.024226458  ;
			  0.046833900  ];
p_actual = p_actual./100

n = length(ob.x)
m = size(ob.f_ab,2)


actual_log_like = ob.complete_likelihood(p_actual,n,m)




%% generation
%% 
mo = Model.generate(0.01,100,struct('save_to','cache/models_temp.mat',...
	'path','data/modeldata3.dat','include_blanks',true));
mo.plot();

%% 
mo.test_cache(struct('model_no',1,'step_size',0.005));

%% 
mo = Model.load('cache/models_3.mat');
ob = Observed(struct('name','halo10k','path','data/obsdata2_10000.dat','p_actual',p_actual));
ob.load_models(mo);
ob.save('cache/observed_3_10k.mat');


ob = Observed(struct('name','halo30k','path','data/obsdata2_30000.dat','p_actual',p_actual));
ob.load_models(mo);
ob.save('cache/observed_3_30k.mat');



ob = Observed(struct('name','halo50k','path','data/obsdata2_50000.dat','p_actual',p_actual));
ob.load_models(mo);
ob.save('cache/observed_3_50k.mat');





%% simulate
clc
im=1;

m = length(mo.data);

tic

P = ob.p_actual;%[.015 .2 .005 .07 .12 .005 .005 .01 .27 .1 .05 .004 .001 .02 .095 .03];
if sum(P) ~= 1
	error('P must sum to 1')
end

if m ~= length(P)
	error('P and models are different lengths')
end
PC = cumsum(P);

n = 500;
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
[p,all_p,ll] = ob.em(struct(	'Xn',100,...
										'max_seconds', 60,...
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





