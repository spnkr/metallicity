
%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;


mnames = {'halo3', 'halo5', 'gen1', 'gen2'};
mi = Mixture.load(mnames{1});



im=1;
ns=[250 500 1000 1500 2000 3000 5000 7000 9000 9969];
pval=[];
ll=[];
llt=[];
for k=1:length(ns)
	i = ns(k);
	mi = Mixture.load('halo3');
	mi.save_as(strcat(['halo3_n' num2str(i)]));
	sepr(mi.filename);
	mi.em(struct('max_iters',2,'n',i,'p0_eval','rand(num_models,1)','quick_print',NaN));
	mi.update_stats
	mi.save
	pval(k) = mi.lrt_chi2_sig;
	ll(k) = mi.est_loglike;
	llt(k) = mi.true_loglike;
end


im=im+1;
fig
subplot(1,3,1);
plot(ns,lrt,'.-')
flabel('n','LRT p-value','p-value LRT with different sample sizes');
subplot(1,3,2);
plot(ns,ll,'.-')
flabel('n','l(\pi)','l(\pi) with different sample sizes');
subplot(1,3,3);
plot(ns,llt,'.-')
flabel('n','l(\pi true)','l(\pi true) with different sample sizes');
