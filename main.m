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

mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load(mnames{1});


%% mixtures --------------------------------------------------
ii = 3;
hls = [2 5 7 8 9 10 12 14 15 17 20];
mi = Mixture(struct(	'save_as',strcat(['test_halo' num2str(hls(ii)) '']),...
						'model_path','data/mastertemp.dat',...
						'obs_path',strcat(['data/obsdata' num2str(hls(ii)) '.dat']),...
						'pi_true', Mixture.get_pi_true(hls(ii)),...
						'graph',false));

%% 
im=1;
mi.em(struct('max_iters',10,'p0_eval','rand(num_models,1)','quick_print',100,'interactive',true));







%% bootstrapping ---------------------------------------------
bnames = {'bo_halo3_nonpara', 'bo_halo3_para', 'bo_halo3_para_true_pi'};
load(strcat(['cache/' bnames{1}]));

%% 
im=1;
bo.plot_error_bars_both
bo.plot_error_bars
bo.plot_covar
bo.plot_data_coverage('temp/halo3_new_boot_10000_579');

bo = Bootstrap.generate_parametric(	'bo_halo3333',mi,0,5,100,...
									struct(	'max_iters',2,'Ainit_p',mi.pi_true,'Xinteractive',true));

%generate from existant cahced files
mi = Mixture.load('halo3');
bo = Bootstrap('bo_halo3_para_true_pi', mi, 'temp/halo3_boot_10000_', 1:4059);
bo.save













%% models ---------------------------------------------------
%1 Fe/H
%2 alpha/Fe
%3 weight
%4 tacc
%5 lsat
%6 nsat  
monames = {'mo_20k', 'mo_all'};
load(strcat(['cache/' monames{1}]));

%% 
mo.plot_lines(113,true);
mo.density(2,.9)
%generate from raw data files for specified sats
mo = Model.generate(1:2,'mo_something');















%% simulation ---------------------------------
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

mi.em



















