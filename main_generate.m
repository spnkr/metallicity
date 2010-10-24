
%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;
mnames = {'halo3', 'halo5', 'gen1', 'gen2', 'halo3_1600', 'halo5_1600', 'halo3_30k', 'halo3_50k'};

max_seconds = 99999999999;

mi = Mixture.load('gen1');
sepr(mi.filename)
mi.em(struct('max_seconds',max_seconds,'p0_eval','rand(num_models,1)'));
mi.update_stats
mi
mi.save

mi = Mixture.load('gen2');
sepr(mi.filename)
mi.em(struct('max_seconds',max_seconds,'p0_eval','rand(num_models,1)'));
mi.update_stats
mi
mi.save


mi = Mixture(struct(	'save_as','halo3',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_10000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

mi.em(struct('max_seconds',max_seconds));
mi.update_stats
mi.save
mi.save



mi = Mixture(struct(	'save_as','halo3_1600',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_10000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

mi.em(struct('max_iters',1600,'ll_stop_lookback',1000000));
mi.update_stats
mi.save



mi = Mixture(struct(	'save_as','halo5',...
						'model_path','data/modeldata5.dat',...
						'obs_path','data/obsdata5_30000.dat',...
						'pi_true', Mixture.get_pi_true(5),...
						'graph',false));

mi.em(struct('max_seconds',max_seconds));
mi.update_stats
mi.save

mi = Mixture(struct(	'save_as','halo5_1600',...
						'model_path','data/modeldata5.dat',...
						'obs_path','data/obsdata5_30000.dat',...
						'pi_true', Mixture.get_pi_true(5),...
						'graph',false));

mi.em(struct('max_iters',1600,'ll_stop_lookback',1000000));
mi.update_stats
mi.save




mi = Mixture(struct(	'save_as','halo3_30k',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_30000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

mi.em(struct('max_seconds',max_seconds));
mi.update_stats
mi.save


mi = Mixture(struct(	'save_as','halo3_50k',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_50000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

mi.em(struct('max_seconds',max_seconds));
mi.update_stats
mi.save



sepr
sepr
sepr
sepr
sepr
sepr
sepr
sepr


for i=1:length(mnames)
	mi = Mixture.load(mnames{i})
end
