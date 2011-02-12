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
Model.generate(NaN,[2 5 7 8],'mo_h2578');
load(strcat(['cache/', 'mo_h2578']));


%% models ---------------------------------------------------
%1 Fe/H
%2 alpha/Fe
%3 weight
%4 tacc
%5 lsat
%6 nsat  
monames = {'mo_h2','large/mo_h25','large/mo_h257','large/mo_h2578'};
load(strcat(['cache/' monames{2}]));



%% 
l=1000
y=betarnd(2,5,1,l);
y=sort(y);
plot(1:l,y);



%% 
im=1;
clc
results = {};
tic
M = max(mo.sats)-1;
E = 1000;
for i=1:M
	j = i+1;
	
	x = mo.data(mo.data(:,6)==i,:);
	x = sortrows(x,1);
	x(:,3) = x(:,3)./sum(x(:,3));
	y = mo.data(mo.data(:,6)==j,:);
	y(:,3) = y(:,3)./sum(y(:,3));
	y = sortrows(y,1);
	
	xc = cumsum(x(:,3));
	yc = cumsum(y(:,3));
	
	D = zeros(1,E);
	
	for e=1:(E-1)
		xndx = sum(xc<=(e/E))+1;
		yndx = sum(yc<=(e/E))+1;
		
		D(e) = sqrt((x(xndx,1)-y(yndx,1))^2+(x(xndx,2)-y(yndx,2))^2);
	end
	
	results{i} = struct('i',i,'j',j,'sum',sum(D),'normalized_sum',sum(D)/E,'D',D);
	strcat(['finished for ',num2str(i),'/',num2str(j)])
end

toc

del = zeros(1,M);
for k=1:M
	del(k) = results{k}.normalized_sum;
end

%% 
figure(1)
subplot(1,2,1)
plot(1:M,del,'k.-');
title(strcat('Distance between curves, E=', num2str(E)));
xlabel('Sat index');
ylabel('D(i,i-1)');



subplot(1,2,2)
plot(1:M,log(del),'r.-');
title(strcat('Log distance between curves, E=', num2str(E)));
xlabel('Sat index');
ylabel('log(D(i,i-1))');




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
[kld,D_KL] = kl_divergence(model_path,I,reps,1,2);

model_path = 'data/mastertemp.dat';
[kld,D_KL] = kl_divergence(model_path,I,reps,1,2);




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













