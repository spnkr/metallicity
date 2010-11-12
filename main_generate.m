%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;

%% 
tic
im=1;
mi = Mixture(struct(	'save_as',strcat(['halo3_r2']),...
						'model_path','data/modeldata3.dat',...
						'obs_path',strcat(['data/obsdata3_10000.dat']),...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));


mi.em(struct('Xmax_iters',2,'p0_eval','rand(num_models,1)','quick_print',999999999,'interactive',false));
h = figure(1);
saveas(h,strcat(['media/diag_' mi.filename '.pdf']),'pdf');

mi.update_stats
mi.save


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




hls = [2 5 7 8 9 10 12 14 15 17 20];
for ii=1:length(hls)
	tic
	im=1;
	try
	mi = Mixture(struct(	'save_as',strcat(['halo' num2str(hls(ii)) '_r2']),...
							'model_path','data/mastertemp.dat',...
							'obs_path',strcat(['data/obsdata' num2str(hls(ii)) '.dat']),...
							'pi_true', Mixture.get_pi_true(hls(ii)),...
							'graph',false));


	mi.em(struct('Xmax_iters',10,'p0_eval','rand(num_models,1)','quick_print',999999999,'interactive',false));
	h = figure(1);
	saveas(h,strcat(['media/diag_' mi.filename '.pdf']),'pdf');

	mi.update_stats
	mi.save


	mi.plot_info_error_bars(2)
	h = figure(2);
	saveas(h,strcat(['media/errorbar_' mi.filename '.pdf']),'pdf');


	mi.plot_correl(true,true);
	h = figure(3);
	saveas(h,strcat(['media/correl_' mi.filename '.pdf']),'pdf');



	plot_formation_history(mi);
	h = figure(4);
	saveas(h,strcat(['media/fhist_' mi.filename '.pdf']),'pdf');
	
	catch
		sepr
		sepr
		sepr
		sepr
		warning('failed!!!!')
		hls(ii)
		sepr
		sepr
		sepr
		sepr
	end
	toc
end





tic
im=1;
mi = Mixture(struct(	'save_as',strcat(['halo5']),...
						'model_path','data/modeldata5.dat',...
						'obs_path',strcat(['data/obsdata5_30000.dat']),...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));


mi.em(struct('Xmax_iters',2,'p0_eval','rand(num_models,1)','quick_print',999999999,'interactive',false));
h = figure(1);
saveas(h,strcat(['media/diag_' mi.filename '.pdf']),'pdf');

mi.update_stats
mi.save


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





tic
im=1;
mi = Mixture(struct(	'save_as',strcat(['halo5_r2']),...
						'model_path','data/modeldata5.dat',...
						'obs_path',strcat(['data/obsdata5_30000.dat']),...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));


mi.em(struct('Xmax_iters',2,'p0_eval','rand(num_models,1)','quick_print',999999999,'interactive',false));
h = figure(1);
saveas(h,strcat(['media/diag_' mi.filename '.pdf']),'pdf');

mi.update_stats
mi.save


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


%% 








































%% 
error('ignore');

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
