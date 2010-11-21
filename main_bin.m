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

load('cache/pistar_orig.mat')
load('cache/pistar.mat')
load('cache/pistar_mn.mat')


mnames = {'halo3', 'halo2','halo5','halo7','halo8','halo9','halo10','halo12','halo14','halo15','halo17','halo20',...
	'gen1', 'gen2', 'halo3_1600', 'halo5_1600', 'halo3_30k', 'halo3_50k'};
mi = Mixture.load(mnames{1});



%% 












%% 

clc
x=[];
x(:,1) = mi.pi_true;
x(:,2) = mi.pi_est;
for i=1:3
	mii = Mixture.load('temp/halo3_new_boot_10000_', num2str(i));
	x(:,i+2) = mii.pi_est;
end
rounder(x.*100)


%% 
im=1;
mii = Mixture.load('temp/halo3_new_boot_10000_', num2str(1));
mii.em(struct('init_p',NaN,'quick_print',10,'interactive',true));





%% 
clc
h=figure(2);
clf(h);
subplot(1,2,1);
fi=1;
ndxa=1:size(mi.x,1);
ndx = mi.f(ndxa,fi)>0;
F=mi.f(ndx,fi);
scatter(mi.x(ndx,1), mi.x(ndx,2), max(1,F), [.1 .1 .2], 'filled');
flabel('x','y',['halo3'])
halo3_F = [min(F) max(F)]
axis([-3 0 -.3 .6])


subplot(1,2,2);
fi=1;
ndx = mii.f(ndxa,fi)>0;
F=mii.f(ndx,fi);
scatter(mii.x(ndx,1), mii.x(ndx,2), max(1,F), [.5 .4 .1], 'filled');
flabel('x','y',['GEN'])
gen_F = [min(F) max(F)]
axis([-3 0 -.3 .6])



%% 
im=1;
[S,V,correl,stdev] = mi.bootstrap_covariance(pistar_mn);
mi.plot_bootstrap_covar(S,V,correl,stdev,pistar_mn);

%% 
im=8;
mi.bootstrap_plot_bars_both(pistar);

%% 
im=9;
mi.bootstrap_plot_bars(pistar);
title('pistar')

im=8;
mi.bootstrap_plot_bars(pistar_mn);
title('mn')


%% 
im=8;
mi.plot_info_error_bars(2);



%% 
round([mean(pistar,2) mean(pistar_mn,2) mi.pi_true].*100.*1000)./1000


%% 
MB=1;
base_mb=593;

tic
im=1;
pmt = bootstrap_generate(base_mb,mi,MB,5000,struct('max_iters',200,'quick_print',999999,...
											'interactive',true,'do_m_n',false,'init_p',mi.pi_true));

difference_in_pct=rounder(100.*([abs(pmt-mi.pi_est)]),1000)
toc



%% 
mi_10k = Mixture.load(strcat(['temp/halo3_new_boot_10000_' num2str(587)]));
mi_100 = Mixture.load(strcat(['temp/halo3_new_boot_100_' num2str(588)]));
mi_50 = Mixture.load(strcat(['temp/halo3_new_boot_50_' num2str(589)]));
mi_1k = Mixture.load(strcat(['temp/halo3_new_boot_1000_' num2str(590)]));
mi_2k = Mixture.load(strcat(['temp/halo3_new_boot_2000_' num2str(591)]));


%% 
pistar = [];
for i=1:586
	mi = Mixture.load(strcat(['temp/halo3_new_boot_10000_' num2str(i)]));
	pistar(:,size(pistar,2)+1) = mi.pi_est;
end

size(pistar)

save('cache/pistar.mat','pistar')

'finished'
clear
clc

load('cache/pistar.mat')
'final is'
size(pistar)


%% 
pistar_mn = [];
for i=1:1000
	mi = Mixture.load(strcat(['temp/halo3_new_mn_boot_10000_' num2str(i)]));
	pistar_mn(:,size(pistar_mn,2)+1) = mi.pi_est;
end

size(pistar2)

save('cache/pistar_mn.mat','pistar_mn')

'finished'
clear
clc

load('cache/pistar_mn.mat')
'final is'
size(pistar_mn)




%% 
pistarU = pistar_mn;
im=10;
clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
			'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
fg=figure(im);clf(fg);im=im+1;
hold on;for i=1:size(pistarU,1)
	plot(1:size(pistarU(i,:),2),pistarU(i,:),strcat([clrs(i) '.-']))
end
hold off
flabel('bs trial num','\pi_j','Estimated pi over bootstrap trials')







%% 
im=1;
[S,V,correl,stdev] = mi.bootstrap_covariance(pistar_mn);
mi.correl = correl;
mi.variance = V;
mi.stdev=stdev;
mi.covar=S;







%% 
clc
mi=Mixture.load('halo3_boot_10000_1');



%% 
clc
im=1;
[bdata] = bootstrap_generate(mi,100);




%% 
im=1;
mi.plot_correl(true,true);
mi.plot_stdev
mi.plot_zscores






%%
figure(10)
plot(1,2)

%% 
im=1;
mi = Mixture.load('gen1');
sepr(mi.filename)
mi.em(struct('max_seconds',99999,'p0_eval','rand(num_models,1)','quick_print',10));
mi.update_stats
mi




%% 
[p,ll,P,LL] = ...
	em(mi.x,mi.f,struct(	'num_models',mi.num_models,...
							'p_actual',mi.pi_true,...
							'loglike_true', mi.loglike_true,...
							'interactive',true,...
							'max_seconds',10,...
							'quick_print',50));

[I,S,V,stdev] = fisher_info(mi.pi_est,mi.f);







%% 
pt = [.48 .32 .1 .1]';
mu = [-1 1 1.5 2];
sigma2 = [1 1 .75 3];

pt = [.68 .32]';
mu = [-1 1];
sigma2 = [1 1];


n=10000;
[x,f] = simulate_n(pt,mu,sigma2, n);
m=size(f,2);

mi = Mixture(struct(	'save_as','gen999',...
						'f',f,...
						'x',x,...
						'pi_true', pt));
mi.models = [mu;sigma2];
mi.save










%% simulate
[p,P,ll,LL] = ob.simulate(mo,struct('sample',500,...
											'max_seconds',60*30,...
											'interactive',false));
p										



































