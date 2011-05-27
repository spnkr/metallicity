clear
clc
format short g
addpath('etc/');
addpath('data/');

global im;
im=1;
global constant_im;
constant_im=false;


%% load saved
mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load('halo2_NEW_05');



%% new way
grid_size = 0.1;
noise = 0.05;


tic
[gridded_model_data,noise,grid_size] = Halo.generate_model_data(struct(...
											'noise',0.15,...
											'grid_size',0.1,...
											'gridm',[0 1e5 1e6 1e7 1e8 1e9 inf],...
											'gridt',[0 2 8 10 12 inf],...
											'halos',[2],...
											'xdata_point_limit',600000));
toc
display('done generate_model_data');

%% 
save('cache/temp/gridded_model_data_15','gridded_model_data')
load('cache/temp/gridded_model_data_15')
load('cache/temp/gridded_model_data_05')

%% 

tic
mi = Mixture(struct(	'save_as',strcat(['halo2_NEW_15']),...
						'model_path','none',...
						'model_obj', gridded_model_data,...
						'obs_path',strcat(['data/obsdata2.dat']),...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

toc
display('done generating density');

tic
mi.em(struct('Xmax_iters',10,'p0_eval','rand(num_models,1)','quick_print',999999999,'interactive',true));
toc
display('done em');
tic
mi.update_stats()
toc
display('done updating stats');

mi.save()


%% 
mi.plot_history();
mi.plot_correl(true,true);
mi.plot_stdev();
mi.plot_zscores();
mi.plot_info_error_bars(1);









%% old data
all_halos = [2 5 7 8 9 10 12 14 15 17 20];

mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load(mnames{1});
















%% ignore below
clear
clc
format short g
addpath('etc/');
addpath('data/');
load('cache/temp/newdat')
im=1;





%% run

[dat,noise,grid_size] = Halo.generate_model_data(struct('noise',0.15,'grid_size',0.1,'halos',[2],...
	'data_point_limit',300000));

[gdata,noise,grid_size] = Halo.generate_observed_data(struct('n',100,'noise',0.15,'grid_size',0.1,...
	'halos',[2],'data_point_limit',300000));

dat = Halo.read(struct('halos',[5],'xxdata_pont_limit',2000000));



%% internal

dat = Halo.read(struct('halos',[5],'xxdata_pont_limit',2000000));
save('cache/temp/newdat','dat');



