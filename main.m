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
multicount = 2;
max_seconds = 60*60;

[ob, mo] = Observed.load(3,10);




%% templates: em
[ob, mo] = Observed.load(3,10);
[p,P,ll] = ob.em(struct(...
				'max_iters',2,...
				'Xmax_seconds',60*60,...
				'init','rand(m,1)',...
				'interactive',true));
print_pi(p,ll,mo);
ob.plot_differences(p);


%% template: em multi
[ob, mo] = Observed.load(3,10);
multicount = 2;
max_seconds = 5;
max_iters = 100;

sepr('start blah')
ob.em_multi(mo, struct( 'count',multicount,...
						'max_seconds', max_seconds, ...
						'interactive',false));
done('end blah');






%% 
multicount = 10;
max_seconds = 60*20;
[ob, mo] = Observed.load(3,10);
	sepr('start model 3 true')
		init_p_to_use = ob.p_actual;
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',init_p_to_use,...
								'save',strcat(['cache/em_multi_m3_true_' num2str(max_seconds) '_' num2str(multicount) 'x.mat'])));
		done('end model 3 true');


		
		
multicount = 2;
max_seconds = 60*60;
[ob, mo] = Observed.load(5,30);
	sepr('start model 5')
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false));
		done('end model 5');






multicount = 8;
max_seconds = 7*60;
[ob, mo] = Observed.load(3,10);
	sepr('start model 3 true pm 3')
		init_p_to_use = ob.p_actual.*0.03.*rand(size(ob.p_actual,1),1);
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',init_p_to_use,...
								'save',strcat(['cache/em_multi_m3_true_pm3_' num2str(max_seconds) '_' num2str(multicount) 'x.mat'])));
		done('end model 3 true pm 3');
	sepr('start model 3 rand')
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false));
		done('end model 3 rand');







max_seconds = 60*60;
[ob, mo] = Observed.load(3,30);
	sepr('start model 3 30k rand')
		[p,P,ll] = ob.em(struct(...
				'max_seconds',max_seconds,...
				'init','rand(m,1)',...
				'interactive',false));
		print_pi(p,ll,mo);
		ob.plot_differences(p);
		done('end model 3 30k rand');


max_seconds = 60*60*5;
[ob, mo] = Observed.load(3,30);
	sepr('start model 3 30k mega')
		[p,P,ll] = ob.em(struct(...
				'max_seconds',max_seconds,...
				'init','rand(m,1)',...
				'interactive',false));
		print_pi(p,ll,mo);
		ob.plot_differences(p);
		save('cache/mega_30k_run.mat','p','P','ll');
		done('end model 3 30k mega');











%% simulate
[p,P,ll,LL] = ob.simulate(mo,struct('sample',10000,...
											'max_seconds',60*30,...
											'interactive',false));




							
%% em convergence
%%
ob.plot_weight_changes(struct('path','cache/all_p_50_full_runs.mat'));


%% plots
%% 
mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

%% 
ob.plot();



















%% shifts etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;

[ob, mo] = Observed.load(3,10);








%% 
sepr('single dif shift')
load('cache/models_3.mat');
load('cache/observed_3_10k.mat');

ob.x = ob.x + (0.05/2);
ob.y = ob.y + (0.05/2);
ob.load_models(mo);

[p,P,ll] = ob.em(struct(...
				'min_iters',900,...
				'init','rand(m,1)',...
				'interactive',false));
p

[p,P,ll] = ob.em(struct(...
				'min_iters',900,...
				'init','rand(m,1)',...
				'interactive',false));
p



load('cache/models_3_c1.mat');
load('cache/observed_3_10k_c1.mat');


sepr('single very long')
[p,P,ll] = ob.em(struct(...
				'min_iters',5000,...
				'init','rand(m,1)',...
				'interactive',false));
p

sepr('single crazy long')
[p,P,ll] = ob.em(struct(...
				'min_iters',7000,...
				'init','rand(m,1)',...
				'interactive',false));
p

%% 

im=im+1;
multicount = 9;
sepr('multi 1')
max_seconds = 60*45;
[best_ll, best_p, all_p, all_ll] = ob.em_multi(struct(	'count',multicount,...
								'Xn',10,...
								'save','cache/p_ll_run_weight_rand_c.mat',...
								'max_seconds', max_seconds, ...
								'interactive',false));
print_pi(best_p,mo)
ob.plot_differences(best_p)


sepr('multi 2')
im=im+1;
multicount = 9;
sepr('multi')
max_seconds = 60*45;
[best_ll, best_p, all_p, all_ll] = ob.em_multi(struct(	'count',multicount,...
								'Xn',10,...
								'save','cache/p_ll_run_weight_rand_c2.mat',...
								'max_seconds', max_seconds, ...
								'interactive',false));
print_pi(best_p,mo)
ob.plot_differences(best_p)

%% 
im=5;
load('cache/models_3.mat');
load('cache/observed_3_10k.mat');



[p,P,ll] = ob.em(struct(...
				'max_seconds',60*5,...
				'init','rand(m,1)',...
				'interactive',true));
p






