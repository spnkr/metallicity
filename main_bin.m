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


load('cache/pistar.mat')
load('cache/pistar2.mat')


mnames = {'halo3', 'halo5', 'gen1', 'gen2', 'halo3_1600', 'halo5_1600', 'halo3_30k', 'halo3_50k'};
mi = Mixture.load(mnames{1});



%% 
im=1;
[S,V,correl,stdev] = mi.bootstrap_covariance(pistar);
mi.plot_bootstrap_covar(S,V,correl,stdev,pistar);

%% 
im=6;
mi.bootstrap_plot_bars_both(pistar);

%% 
im=7;
mi.bootstrap_plot_bars(pistar);


%% 
im=8;
mi.plot_info_error_bars(2);

%% 



%% 
MB=1000;
base_mb=0;

'doing para'
size(pistar)

im=-1;
bootstrap_generate(base_mb,mi,MB,10000,struct('Xmax_iters',2,'quick_print',999999,...
											'interactive',false,'do_m_n',true));

'done with para'



'doing NONpara'
size(pistar)

im=-1;
bootstrap_generate(base_mb,mi,MB,10000,struct('Xmax_iters',2,'quick_print',999999,...
											'interactive',false,'do_m_n',false));

'done with NONpara'


%% 
pistar2 = [];
for i=1:2
	mi = Mixture.load(strcat(['temp/halo3_ni_boot_10000_' num2str(i)]));
	pistar2(:,size(pistar2,2)+1) = mi.pi_est;
end

size(pistar2)

save('cache/pistar2.mat','pistar2')

'finished'
clear
clc

load('cache/pistar2.mat')
'final is'
size(pistar2)




%% 
pistarU = pistar2;
im=10;
clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
			'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
fg=figure(im);clf(fg);im=im+1;
hold on;for i=1:size(pistarU,1)
	plot(1:size(pistarU(i,:),2),pistarU(i,:),strcat([clrs(i) '.-']))
end
hold off
flabel('bs trial num','\pi_j','Estimated pi over bootstrap trials')















%% 
clc
mi=Mixture.load('halo3_boot_10000_1');



%% 
clc
im=1;
[bdata] = bootstrap_generate(mi,100);




%% 
im=1;
mi.plot_correl(true,true);
mi.plot_stdev
mi.plot_zscores



%% 
mi2 = Mixture(struct(	'save_as','halo3X',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_10000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

%% 
mi.models_sparse=mi2.models_sparse;












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



































