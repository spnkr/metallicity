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




im=1;
hls = [2 3 5 7 8 9 10 12 14 15 17 20];
B=100;
BN=10000;

tic
parfor i=1:length(hls)
	hndx = hls(i);
	mi = Mixture.load(strcat(['halo' num2str(hndx)]));
	init_p = (mi.pi_true+0.01)./sum(mi.pi_true+0.01);
	bgn_args = struct('Xmax_iters',1,'Xinit_p',init_p,'quick_print',9999999,'interactive',false);
	
	bo = Bootstrap.generate_parametric(strcat(['bo_halo' num2str(hndx) '_para']),mi,0,B,BN,bgn_args);
	bo.pi_true = mi.pi_true;
	bo.save
	%bo.plot_spread
	
	bo = Bootstrap.generate_non_parametric(strcat(['bo_halo' num2str(hndx) '_nonpara']),mi,0,B,BN,bgn_args);
	bo.pi_true = mi.pi_true;
	bo.save
	%bo.plot_spread
	sepr(['Finished halo ' num2str(hndx)])
end
toc

sepr
sepr('Finished all')
sepr



%% 
im=1;
load(strcat(['cache/bo_halo5_nonpara']));
%bo.plot_error_bars;