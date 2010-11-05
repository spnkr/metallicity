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


load('cache/pistar.mat')

mnames = {'halo3', 'halo5', 'gen1', 'gen2', 'halo3_1600', 'halo5_1600', 'halo3_30k', 'halo3_50k'};
mi = Mixture.load(mnames{1});




%% 
clc
[S,V,correl,stdev] = bootstrap_covar(mi,pistar);


%% 
w = max(mi.pi_est.*2000,10);

fg=figure(1)
clf(fg)
hold on
scatter(1:mi.num_models,mi.variance,w,[0 0 0],'filled')
scatter(1:mi.num_models,V,w,[.5 .1 .3])
hold off
flabel('j','Variance',['Filled: Info variance. Open: Boostrap variance (n=' num2str(size(pistar,2)) ')'])
axis([0 size(pistar,1) min(min(min([mi.variance;V])),-.01) max(max(max([mi.variance;V])),.01)])


fg=figure(2)
clf(fg)
hold on
scatter(1:mi.num_models,mi.stdev,w,[0 0 0],'filled')
scatter(1:mi.num_models,stdev,w,[.5 .1 .3])
hold off
flabel('j','Standard deviation',['Filled: Info \sigma. Open: Boostrap \sigma (n=' num2str(size(pistar,2)) ')'])
axis([0 size(pistar,1) min(min(min([mi.stdev;stdev])),-.01) max(max(max([mi.stdev;stdev])),.01)])


fg=figure(3)
clf(fg)

stddelta=mi.stdev-stdev;
scatter(1:mi.num_models,stddelta,w,[0 0 0],'filled')

flabel('\pi','Standard deviation delta',['Info v. Boostrap \sigma (n=' num2str(size(pistar,2)) ')'])
axis([0 size(pistar,1) min(min(stddelta),-.1) max(max(stddelta),.1)])



fg=figure(4)
clf(fg)

w = zeros(mi.num_models,mi.num_models);
for i=1:mi.num_models
	for j=1:mi.num_models
		w(i,j) = mi.pi_est(i) + mi.pi_est(j);
	end
end
w=max(200.*reshape(w,1,size(w,1)^2),10);

hold on
stddelta=reshape(mi.correl,1,size(mi.correl,1)^2)-reshape(correl,1,size(correl,1)^2);
scatter(1:mi.num_models^2,stddelta,w,[0 0 0],'filled')
hold off
flabel('j','Correlation delta',['Info v. Boostrap Correlation (n=' num2str(size(pistar,2)) ')'])
axis([0 length(stddelta) min(min(stddelta),-.1) max(max(stddelta),.1)])





fg=figure(5)
clf(fg)
rc = ceil(sqrt(mi.num_models));
for i=1:mi.num_models
	subplot(rc,rc,i);
	hist(pistar(i,:))
end



fg=figure(6)
clf(fg)
conf = zeros(mi.num_models,3);
B = size(pistar,2);
clev=0.95;
lndx = ceil(B*((1-clev)/2))
undx = ceil(B - lndx)
for i=1:mi.num_models
	sv=sort(pistar(i,:));
	mu=mi.pi_est(i);
	conf(i,1) = mu;
	conf(i,2) = sv(lndx);
	conf(i,3) = sv(undx);
end

% errorbar(1:mi.num_models,conf(:,1),conf(:,2),conf(:,3),'k.');
hold on

for i=1:mi.num_models
	plot(i.*ones(1,2), conf(i,[2 3]), 'b.-');
end

plot(mi.pi_est,'k.')
plot(mi.pi_true,'rx')
hold off
flabel('j','\pi',[num2str(clev*100) '% Boostrap error bars (n=' num2str(size(pistar,2)) '). x=true, .=est'])



fg=figure(7)
clf(fg)

errorbar(1:mi.num_models,mi.pi_est,3.*mi.stdev,'c.');
hold on
errorbar(1:mi.num_models,mi.pi_est,2.*mi.stdev,'g.');
errorbar(1:mi.num_models,mi.pi_est,mi.stdev,'k.');
plot(mi.pi_true,'rx')
hold off
flabel('j','\pi',['Info error bars. black=1x, green=2x, cyan=3x. x=true, .=est'])





%% 
MB=5000;
base_mb=4059;

'starting with size'
size(pistar)

im=-1;
bootstrap_generate(base_mb,mi,MB,10000,struct('XXmax_iters',2,'quick_print',999999,...
											'interactive',false));

'done with bootstrap_covar'


%% 
pistar = [];
for i=1:4059
	mi = Mixture.load(strcat(['temp/halo3_boot_10000_' num2str(i)]));
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
im=10;
clrs = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k' ...
			'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r' 'g' 'b' 'c' 'm' 'y' 'k'];
fg=figure(im);clf(fg);im=im+1;
hold on;for i=1:size(pistar,1)
	plot(1:size(pistar(i,:),2),pistar(i,:),strcat([clrs(i) '.-']))
end
hold off
flabel('bs trial num','\pi_j','Estimated pi over bootstrap trials')


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
mi2 = Mixture(struct(	'save_as','halo3X',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_10000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

%% 
mi.models_sparse=mi2.models_sparse;












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
							'max_seconds',60,...
							'baseline_p',p_bs_target_3,...
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



































