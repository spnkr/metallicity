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
global constant_im;
constant_im=false;

mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load(mnames{1});




%% 
mastertemp = load('data/mastertemp_s4t2m7.dat');

modeldata1 = load('data/modeldata_s4t2m78.dat');
modeldata2 = load('data/modeldata_s4t2m79.dat');



obsdata1 = load('data/obsdata_s4t2m78.dat');
obsdata2 = load('data/obsdata_s4t2m79.dat');


%% 
size(mastertemp)
size(modeldata1)
size(modeldata2)
size(obsdata1)
size(obsdata2)




%% load new data, run em, show/save plots
%make sure the data files don't have headers--must be only numbers
new_halo_name = 'new s4t2m78';%saved as this
model_path = 'data/mastertemp_s4t2m7.dat';
obs_path = 'data/obsdata_s4t2m78.dat';
pi_true = NaN; %replace with nx1 matrix if desired, or leave as-is to skip

tic
im=1;
mi = Mixture(struct(	'save_as',new_halo_name,...
						'model_path',model_path,...
						'obs_path',obs_path,...
						'pi_true', pi_true,...
						'graph',false));


mi.em(struct(	'max_iters',10,...	%max iters to run. change name to something else (e.g. Xmax_iters) to let it run until log like stop changing
				'p0_eval','rand(num_models,1)',... %starting pi value (eval'ed)
				'quick_print',999999999,...		%update UI after this many runs + print current pi state
				'interactive',true));			%show progress (slow)
h = figure(1);
saveas(h,strcat(['media/diag_' mi.filename '.pdf']),'pdf');

mi.update_stats %updates correlation etc. based on em results
mi.save %saves results


mi.plot_info_error_bars(2)
h = figure(2);
saveas(h,strcat(['media/errorbar_' mi.filename '.pdf']),'pdf');


mi.plot_correl(true,true);
h = figure(3);
saveas(h,strcat(['media/correl_' mi.filename '.pdf']),'pdf');


plot_formation_history(mi);
h = figure(4);
saveas(h,strcat(['media/fhist_' mi.filename '.pdf']),'pdf');
toc

sepr('finished');


















%% 
im=1;
hls = [2 3 5 7 8 9 10 12 14 17 20];

for i=1:length(hls)
	im=1;
	ii = hls(i);
	load(strcat(['cache/bo_halo' num2str(ii) '_F_nonpara.mat']))
	
	bo.plot_spread
	h = figure(1);
	saveas(h,strcat(['media_local/spread_f50_nonpara_halo' num2str(ii) '.pdf']),'pdf');
	
	bo.plot_error_bars
	h = figure(2);
	saveas(h,strcat(['media_local/eb_f50_nonpara_halo' num2str(ii) '.pdf']),'pdf');
end







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
bnames = {'bo_halo3_nonpara', 'bo_halo3_para', 'bo_halo3_para_true_pi', 'bo_halo7_mn_init_true', 'bo_halo7_init_true'};
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



















