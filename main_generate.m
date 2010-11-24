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
ii=2;
mi = Mixture.load(strcat(['halo' num2str(ii)]));
mi.pi_true = [	14.467074       7.4354783       20.988769       0.0000000       0.0000000          ;
      0.0000000       9.2341029       33.655759       8.0197003       0.0000000        ;
      0.0000000       2.4001462       2.5733790       0.0000000      0.11881163        ;
    0.076990272      0.35843001      0.25313093      0.22529881     0.048914533        ;
      0.0000000       0.0000000       0.0000000     0.024224180     0.046834727        ]./100;

pt = [];
k=1;

for x=1:5
	for y=1:5
		
		ndx = (x-1)*5 + y;
		
		if mi.model_skip_ndx(ndx) == 1
			pt(k) = mi.pi_true(x,y);
			k=k+1;
		end
	end
end



mi.pi_true = pt';
mi.save

%% 
ii=20;
mi = Mixture.load(strcat(['halo' num2str(ii)]));

mi.pi_true = [0.0000000       0.0000000       33.588798       0.0000000       ...
       0.0000000       7.6198297       43.265491       2.1882882       6.0608552 ...
      0.18076256      0.87203140      0.84924642       1.0273912       2.2196320 ...
     0.053782485      0.30459020      0.47464266      0.46671873      0.24931823 ...
    0.027163081     0.040342210      0.23066268]./100';

mi.save

%% 



%% 
halos = [2 5 7 8 9 10 12 14 15 17 20];
halos=3;
for i=1:length(halos)
	ii=halos(i);
	mi = Mixture.load(strcat(['halo' num2str(ii)]));
	if 1==11
	im=1;
	figure(1)
	mi.plot_info_error_bars(2)
	h = figure(1);
	saveas(h,strcat(['media_local/errorbar_' mi.filename '.pdf']),'pdf');
	
	im=4;
	figure(4);
	plot_formation_history(mi);
	h = figure(4);
	saveas(h,strcat(['media_local/fhist_' mi.filename '.pdf']),'pdf');
	end
	im=3;
	figure(3);
	
	h = figure(3);
	saveas(h,strcat(['media_local/ed_' mi.filename '.pdf']),'pdf');
end



%% 
clc
h=figure(3);
clf(h);
halos = [2 5 7 8 9 10 12 14 15 17 20];
pir = zeros(22,size(halos,2));
for i=1:length(halos)
	ii=halos(i);
	mi = Mixture.load(strcat(['halo' num2str(ii)]));
	try
	pir(:,i) = (mi.pi_est-mi.pi_true')./mi.pi_true';
	catch
		mi
	end
end
pir=100.*pir;

hold on
for i=1:22
	plot([i i],[min(min(pir(isfinite(pir)))) max(max(pir(isfinite(pir))))],'LineStyle','-','Color',[.8 .8 .8]);
end

plot(pir,'k.','MarkerSize',15)
hold off

axis([0 23 min(min(pir(isfinite(pir)))) max(max(pir(isfinite(pir))))]);

flabel('\pi_j','% \Delta v. true','Percent difference in \pi for all halos');

saveas(h,strcat(['media_local/error_dist.pdf']),'pdf');


%% 
clc
avge=23.778/100;
h=figure(3);
clf(h);
halos = [2 5 7 8 9 10 12 14 15 17 20];
pir = zeros(22,size(halos,2));
for i=1:length(halos)
	ii=halos(i);
	mi = Mixture.load(strcat(['halo' num2str(ii)]));
	try
	pir(:,i) = (mi.pi_est-mi.pi_true')./avge;
	catch
		mi
	end
end
pir=100.*pir;

hold on
for i=1:22
	plot([i i],[min(min(pir(isfinite(pir)))) max(max(pir(isfinite(pir))))],'LineStyle','-','Color',[.8 .8 .8]);
end

plot(pir,'k.','MarkerSize',15,'Color',[.8 .5 .2])
hold off

axis([0 23 min(min(pir(isfinite(pir)))) max(max(pir(isfinite(pir))))]);

flabel('\pi_j','% \Delta v. avg. error',['Percent difference in \pi for all halos over avg error (' num2str(avge*100) '%)']);

saveas(h,strcat(['media_local/error_dist_avg.pdf']),'pdf');


%% 
clc
h=figure(3);
clf(h);
halos = [3];
pir = zeros(16,1);
for i=1:length(halos)
	ii=halos(i);
	mi = Mixture.load(strcat(['halo' num2str(ii)]));
	try
	pir(:,i) = (mi.pi_est-mi.pi_true)./mi.pi_true;
	catch
		mi
	end
end
pir=100.*pir;

hold on
for i=1:16
	plot([i i],[min(min(pir(isfinite(pir)))) max(max(pir(isfinite(pir))))],'LineStyle','-','Color',[.8 .8 .8]);
end


plot(pir,'k.','MarkerSize',15,'Color',[.1 .5 .4])
hold off

axis([0 17 min(min(pir(isfinite(pir)))) max(max(pir(isfinite(pir))))]);

flabel('\pi_j','% \Delta v. true','Percent difference in \pi for halo 3');

saveas(h,strcat(['media_local/error_dist_3.pdf']),'pdf');






%% 
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
