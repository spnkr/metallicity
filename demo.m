%% 
%to run a block of code put the cursor after the %% and press command+enter

%% run this block after you open matlab to load the paths etc.
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;

%these halos are pre-loaded with results
halonames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load(halonames{1});


%% pick a halo here and all code after runs on this halo
halo_ndx = 4;% the index in the list halonames above
mi = Mixture.load(halonames{halo_ndx});



%% 


mxv = 1;

im=1;
h=figure(im);
clf(h)

minw=50;
maxw=100;
hold on
mmm = {'new s4t2m72'};%,'5','7','8','9','10','12','14','15','17','20'};
for i=1:length(mmm)
	abc=strcat('temp/',mmm{i});
	mi = Mixture.load(abc);

	w = (mi.pi_true.*10000).^2;
	w = w./range(w);
	w = minw + (w.*maxw);

	st = (mi.pi_est-mi.pi_true);


	subplot(1,2,1)
	hold on
	plot([0 mi.num_models+1], [0 0],'k--');

	cll = i/length(mmm);
	cll = cll./range(cll);
	cll = 0.3 + (cll.*0.68);
		
	
	
	scatter(1:mi.num_models, st, w, 'filled');
	
	title(mi.filename);
	ylabel('|Est-True|')
	xlabel('j')
 	axis([0 mi.num_models+1 -mxv mxv])
	
	
	
	
	
	subplot(1,2,2)
	hold on
	st = (mi.pi_est-mi.pi_true)./mi.pi_true;
	ndx= st==Inf;
	st(ndx) = 1;

	plot([0 mi.num_models+1], [0 0],'k--');


	scatter(1:mi.num_models, st, w, 'filled');


	title(mi.filename);
	ylabel('|Est-True|/True')
	xlabel('j')
 	axis([0 mi.num_models -mxv mxv])
end
hold off

%% 
%% 
abc='temp/new s4t2m75';
mi = Mixture.load(abc);

mxv = 1;

im=1;
fig;

minw=50;
maxw=100;

w = (mi.pi_true.*10000).^2;
w = w./range(w);
w = minw + (w.*maxw);

clrs = .6.*ones(mi.num_models,3);

subplot(1,2,1)
st = (mi.pi_est-mi.pi_true)./mi.pi_true;
ndx= st==Inf;
st(ndx) = 1;
clrs(ndx,:) = clrs(ndx,:).*rand;

plot([0 mi.num_models+1], [0 0],'k-');

hold on
scatter(1:mi.num_models, st, w, clrs, 'filled');
hold off

title(mi.filename);
ylabel('|Est-True|/True')
xlabel('j')
axis([0 mi.num_models -mxv mxv])



subplot(1,2,2)

st = (mi.pi_est-mi.pi_true);



plot([0 mi.num_models+1], [0 0],'k-');

hold on
scatter(1:mi.num_models, st, w, .6.*ones(mi.num_models,3), 'filled');
hold off
title(mi.filename);
ylabel('|Est-True|')
xlabel('j')
axis([0 mi.num_models+1 -mxv mxv])

%% view info about the selected halo
mi
%pi_est is the best estimate of each pi
%covar, variance, stdev, and correl are information based (as opposed to
%bootstrap based, which is not included here)

% 
% 
% %% temp
% im=15;
% plus_minus_num_std_devs=2;
% mi = Mixture.load('temp/new s4t2m72');
% 
% 
% %% 
% 
% im=15;
% mi.plot_info_error_bars(plus_minus_num_std_devs);
% %% 
% im=15;
% mi.plot_correl(true,true);
% %% 
% im=15;
% plot_formation_history(mi);
%% 
im=15;
mi = Mixture.load('halo2_t');
%mi = Mixture.load('temp/new s4t2m77');

mi.em(struct(	'Xmax_iters',10,...	%max iters to run. change name to something else (e.g. Xmax_iters) to let it run until log like stop changing
				'p0_eval','rand(num_models,1)',... %starting pi value (eval'ed)
				'quick_print',999999999,...		%update UI after this many runs + print current pi state
				'interactive',false));			%show progress (slow)

sepr('finished');

%% plots
im=1; %sets figure window to use

plus_minus_num_std_devs=2;
mi.plot_info_error_bars(plus_minus_num_std_devs);

mi.plot_correl(true,true);

plot_formation_history(mi);






%% run em algo and see interactive results
im=1;
mi.em(struct(	'max_iters',10,...	%max iters to run. change name to something else (e.g. Xmax_iters) to let it run until log like stop changing
				'p0_eval','rand(num_models,1)',... %starting pi value (eval'ed)
				'quick_print',999999999,...		%update UI after this many runs + print current pi state
				'interactive',true));			%show progress (slow)

sepr('finished');
%also see later on how to run em on any data you like








%% load new data, run em, show/save plots
%make sure the data files don't have headers--must be only numbers
new_halo_name = 'testhalo';%saved as this
model_path = 'data/mastertemp.dat';
obs_path = 'data/obsdata2.dat';
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




%% run EM on arbitary data
%in case you don't want to use the Mixture object. 
im=1;
%p=pi_estimate
%ll=best log like
%P = all pi ests over time
%LL = all log likes over time
[p,ll,P,LL] = ...
	em(	mi.x,... nx2 matrix of a/fe,fe/h
		mi.f,... nxm matrix of densities where m=number of non-zero models
		struct(	'num_models',mi.num_models,... %number of non-zero models
							'p_actual',mi.pi_true,... %optional
							'loglike_true', mi.loglike_true,... %optional
							'interactive',true,...
							'max_seconds',10,... %you can use max_seconds instead of max_iters; will stop after this # of seconds
							'quick_print',50));
sepr('finished');
%fisher info, but you need to pass it mi.f--f_ab for each model. mi.f is
%set in Mixture.m's constructor (a method called Mixture)
%[observed_fisher_info,covariance,variance,correl,stdev] = fisher_info(p,mi.f);








