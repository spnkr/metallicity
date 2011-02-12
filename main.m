clear
clc
format short g
addpath('etc/');
addpath('data/');

global im;
im=1;
global constant_im;
constant_im=false;


%mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
%mi = Mixture.load(mnames{1});


%% 
all_halos = [2 5 7 8 9 10 12 14 15 17 20];

Model.generate(NaN,[2 5 7 8 9],'mo_h25789');
load(strcat(['cache/', 'mo_h25789']));

sepr('finished')

%% models ---------------------------------------------------
%1 Fe/H
%2 alpha/Fe
%3 weight
%4 tacc
%5 lsat
%6 nsat  
monames = {'mo_h2','large/mo_h25','large/mo_h257','large/mo_h2578','large/mo_h25789'};
load(strcat(['cache/' monames{3}]));

%% 
Mass = [0 1e5;1e5 1e6;1e6 1e7;1e7 1e8;1e8 1e99];


ndx = mo.curve_ndx_by_mass(0,1e5);

ndx = mo.curve_ndx_by_mass_time(0,1e5,0,10);





%% 
load(strcat(['cache/' monames{5}]));
clc
E=100;
hlosm = 'H:[2,5,7,8,9]';
hlosmfn = 'h25789';

im=1;
msssm = '0_1e5';
[results,mass_time_of_sats,DD,del] = walk_distance(mo,E,1e5,1e6,hlosm);
h = figure(1);
saveas(h,strcat(['media_local/cd_' hlosmfn '_' msssm '.pdf']),'pdf');


im=5;
msssm = '1e5_1e6';
[results,mass_time_of_sats,DD] = walk_distance(mo,E,1e5,1e6,hlosm);
h = figure(5);
saveas(h,strcat(['media_local/cd_' hlosmfn '_' msssm '.pdf']),'pdf');


im=2;
msssm = '1e6_1e7';
[results,mass_time_of_sats,DD] = walk_distance(mo,E,1e6,1e7,hlosm);
h = figure(2);
saveas(h,strcat(['media_local/cd_' hlosmfn '_' msssm '.pdf']),'pdf');


im=3;
msssm = '1e7_1e8';
[results,mass_time_of_sats,DD] = walk_distance(mo,E,1e7,1e8,hlosm);
h = figure(3);
saveas(h,strcat(['media_local/cd_' hlosmfn '_' msssm '.pdf']),'pdf');

im=4;
msssm = '1e8_+';
[results,mass_time_of_sats,DD] = walk_distance(mo,E,1e8,1e99,hlosm);
h = figure(4);
saveas(h,strcat(['media_local/cd_' hlosmfn '_' msssm '.pdf']),'pdf');



%% 
load(strcat(['cache/' monames{2}]));
clc
E=100;
hlosm = 'H:[2,5]';
hlosmfn = 'h25';

im=1;
msssm = '1e5_1e6';
[results,mass_time_of_sats,DD] = walk_distance(mo,E,1e5,1e6,hlosm);

h = figure(1);
saveas(h,strcat(['media_local/cd_' hlosmfn '_' msssm '.pdf']),'pdf');


%% 
im=2;
mo.plot_lines(113);
%mo.density(2,.9)
%generate from raw data files for specified sats
%mo = Model.generate(1:2,'mo_something');


%% 
im=2;
mo.plot_lines_compare(10);












%% 
I=1000;
reps = 10;

model_path = 'data/mastertemp_s4t2m7.dat';
[kld,D_KL] = hellinger(model_path,I,reps,1,2);

model_path = 'data/mastertemp.dat';
[kld,D_KL] = hellinger(model_path,I,reps,1,2);




%% 
%make sure the data files don't have headers--must be only numbers
new_halo_name = 'new s4t2m720';%saved as this
model_path = 'data/mastertemp_s4t2m7.dat';
obs_path = 'data/obsdata_s4t2m720.dat';
%pi_true = [top_left top_right bottom_left bottom_right];
%pi_true = [0 92.8 0.35 6.89]'./100; %replace with nx1 matrix if desired, or leave as-is to skip
pi_true = NaN;





tic
im=1;
mi = Mixture(struct(	'save_as',new_halo_name,...
						'model_path',model_path,...
						'obs_path',obs_path,...
						'pi_true', pi_true,...
						'graph',false));


mi.em(struct(	'X_IGNORE_max_iters',10,...	%max iters to run. change name to something else (e.g. Xmax_iters) to let it run until log like stop changing
				'p0_eval','rand(num_models,1)',... %starting pi value (eval'ed)
				'quick_print',100,...		%update UI after this many runs + print current pi state
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













