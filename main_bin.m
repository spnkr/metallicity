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


mnames = {'halo3', 'halo5', 'gen1', 'gen2', 'halo3_1600', 'halo5_1600', 'halo3_30k', 'halo3_50k'};
mi = Mixture.load(mnames{1});




%% 
im=1;
mi.plot_correl(true,true);
mi.plot_stdev
mi.plot_zscores



%% 
mi = Mixture(struct(	'save_as','halo3',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_10000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

mi.em(struct('max_iters',100));
mi.update_stats
mi.save





%% 
im=1;
mi = Mixture.load('gen1');
sepr(mi.filename)
mi.em(struct('max_seconds',99999,'p0_eval','rand(num_models,1)','quick_print',10));
mi.update_stats
mi




%% 
[p,ll,P,LL] = ...
	em(mi.x,mi.f,struct(	'num_models',mi.num_models,...
							'p_actual',mi.pi_true,...
							'loglike_true', mi.loglike_true,...
							'interactive',true,...
							'max_seconds',60,...
							'baseline_p',p_bs_target_3,...
							'quick_print',50));

[I,S,V,stdev] = fisher_info(mi.pi_est,mi.f);







%% 
pt = [.48 .32 .1 .1]';
mu = [-1 1 1.5 2];
sigma2 = [1 1 .75 3];

pt = [.68 .32]';
mu = [-1 1];
sigma2 = [1 1];


n=10000;
[x,f] = simulate_n(pt,mu,sigma2, n);
m=size(f,2);

mi = Mixture(struct(	'save_as','gen999',...
						'f',f,...
						'x',x,...
						'pi_true', pt));
mi.models = [mu;sigma2];
mi.save










%% simulate
[p,P,ll,LL] = ob.simulate(mo,struct('sample',500,...
											'max_seconds',60*30,...
											'interactive',false));
p										



































