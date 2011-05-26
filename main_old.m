clear
clc
format short g
addpath('etc/');
addpath('data/');

global im;
im=1;
global constant_im;
constant_im=false;

all_halos = [2 5 7 8 9 10 12 14 15 17 20];

mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20'};
mi = Mixture.load(mnames{1});



%% 
h=figure(1);clf(h);
ellipse(1,2,pi/8,1,1,'r')



%% 

scale=300;
pval = 0.95;
hi = 10;
mi = Mixture.load(mnames{hi});
hlosmfn = mnames{hi};
% h = figure(1);
% clf(h);

%indice comments by index, not halo #
%1
%indices = [1 2; 1 4; 1 5; 1 6; 2 6; 2 3; 2 4; 2 5; 4 5; 4 6; 8 7;8 9;11 10;12 10;14 9;13 15;11 15];

% 5
%indices = [1 2;1 3;1 4; 1 6;6 5;6 4;8 7; 8 6;12 11;12 10;15 16;16 17;17 18;20 19;19 18;21 22;20 17];


% 10
indices = [1 2;1 6; 1 8;2 3;3 4; 4 5; 5 6;6 7;7 8;8 9;9 10;10 11;11 12;12 13;13 14;14 15;15 16;16 17;17 18;18 19;19 20;20 21;21 22;19 14;...
	14 10;9 6;14 17;18 22];

ROWS = ceil(sqrt(M));
COLS = ceil(sqrt(M));
% ROWS = 3;
% COLS = 3;
M = size(indices,1);
for ndx=1:M
	tic
	ndxi = indices(ndx,1);
	ndxj = indices(ndx,2);
	sig = [mi.covar(ndxi,ndxi) mi.covar(ndxi,ndxj);mi.covar(ndxj,ndxi) mi.covar(ndxj,ndxj)];
	sigi = inv(sig);

	pi_true = [mi.pi_true(ndxi) mi.pi_true(ndxj)];
	pih = [mi.pi_est(ndxi) mi.pi_est(ndxj)];
	pstdev = [mi.stdev(ndxi) mi.stdev(ndxj)];
	pi_eb_ndx = [ndxi ndxj];

	rnk = rank(sigi);
	if rnk ~= 2;warning('rank of inv(sigma) < 2');end
	cval = chi2inv(pval,rnk);


	map = zeros(scale.^2,2);
	for k = 1:scale
		map((k-1)*scale+1:k*scale,1) = k;
		map((k-1)*scale+1:k*scale,2) = linspace(1,scale,scale)';
	end



	map=map./scale;
	omap = map;
	map(:,1) = map(:,1)-pih(1);
	map(:,2) = map(:,2)-pih(2);

	filled = zeros(1,2);
	z = 1;
	pp=length(map);
	for i = 1:pp
		a = map(i,:)';
		if (a' * sigi * a) < cval
			filled(z,:) = omap(i,:);
			z = z+1;
		end
	end



% 	subplot(ROWS,COLS,ndx);
	h = figure(1);
	clf(h);
	plot(filled(:,1),filled(:,2),'Color', [0.8 0.8 0.8], 'Marker', '.', 'LineStyle', 'none');
	axis([0 1 0 1]);

	hold on

	h1=plot(pi_true(1),pi_true(2),'ko','MarkerSize',10,'LineWidth',2);
	h2=plot(pih(1),pih(2),'gx','MarkerSize',10,'LineWidth',2,'Color',[10 198 28]./255);
	legend([h1 h2],strcat(['True \pi: ' num2str(round(pi_true(1)*10000)/100) '%, ' num2str(round(pi_true(2)*10000)/100) '%']),...
		strcat(['Est \pi: ' num2str(round(pih(1)*10000)/100) '%, ' num2str(round(pih(2)*10000)/100) '%']));

	title(strcat(['\pi_' num2str(ndxi) ' and \pi_' num2str(ndxj) '(' num2str(pval*100) '%)']));
	hold off
	saveas(h,strcat(['media_local/el_h' hlosmfn '_' num2str(ndxi) 'v' num2str(ndxj) '.pdf']),'pdf');

	toc
end

'DONE'
%% 
im=10;
mi.plot_info_error_bars(1)


im=20;
mi.plot_correl(true,true);

%% 
clc
tic
I=10000;
reps = 10;


mpaths = {'modeldata_s4t10m7','modeldata_s4t12m7','modeldata_s4t8m7','modeldata_s4t2m7'};


for i=1:length(all_halos)
	im=1;
	plot_hellinger(I,reps,all_halos(i),mpaths);
	h = figure(1);
	saveas(h,strcat(['media_local/hd_' num2str(all_halos(i)) '_' num2str(I) '.pdf']),'pdf');
end

mpaths = {'mastertemp_s4t10m7','mastertemp_s4t12m7','mastertemp_s4t8m7','mastertemp_s4t2m7'};
im=1;
plot_hellinger(I,reps,-1,mpaths);
h = figure(1);
saveas(h,strcat(['media_local/hd_master_' num2str(I) '.pdf']),'pdf');

toc
%% 




%% 


Model.generate(NaN,[2 5 7 8 9 10 12 14 15 17 20],'mo_h25789101214151720');
load(strcat(['cache/', 'mo_h25789101214151720']));

clc
E=10000;
hlosm = 'H:[2,5,7,8,9,10,12,14,15,17,20]';
hlosmfn = 'h25789101214151720';


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




%% models ---------------------------------------------------
%1 Fe/H
%2 alpha/Fe
%3 weight
%4 tacc
%5 lsat
%6 nsat  
monames = {'mo_h25789101214151720','mo_h2','large/mo_h25','large/mo_h257','large/mo_h2578','large/mo_h25789'};
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













