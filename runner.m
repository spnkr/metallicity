%%
clear
clc
format short g
addpath('etc/');
addpath('data/');

global im;
im=1;
global constant_im;
constant_im=false;




%% 
mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mnames2 = {'halo2_NEW_05', 'halo2_NEW_15','halo5_NEW_01','halo5_NEW_05','halo2_NEW_15'};
mi = Mixture.load(mnames{1});
mi01 = Mixture.load('halo5_NEW_01');
mi05 = Mixture.load('halo5_NEW_05');
mi15 = Mixture.load('halo5_NEW_15');



%% 
mi01.plot_history();
mi05.plot_history();
mi15.plot_history();




%% 
%load('cache/temp/gridded_model_data_15')
load('cache/temp/gridded_model_data_05')


%% 

tic
[gridded_model_data,noise,grid_size] = Halo.generate_model_data(struct(...
											'noise',0.01,...
											'grid_size',0.1,...
											'gridm',[0 1e5 1e6 1e7 1e8 1e9 inf],...
											'gridt',[0 2 8 10 12 inf],...
											'halos',[5],...
											'xdata_point_limit',600000));
toc
display('done generate_model_data');

%% 
tic
mi = Mixture(struct(	'save_as',strcat(['halo5_NEW_01']),...
						'model_path','none',...
						'model_obj', gridded_model_data,...
						'obs_path',strcat(['data/obsdata5.dat']),...
						'Xpi_true', Mixture.get_pi_true(3),...
						'graph',false));

toc
display('done generating density');
%% 
tic
mi.em(struct('Xmax_iters',10,'p0_eval','rand(num_models,1)','quick_print',999999999,'interactive',true));
toc
display('done em');
tic
mi.update_stats()
toc
display('done updating stats');

%% 
mi.save()

