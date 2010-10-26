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


mnames = {'halo3', 'halo5', 'gen1', 'gen2', 'halo3_1600', 'halo5_1600', 'halo3_30k', 'halo3_50k'};
mi = Mixture.load(mnames{1});

%% 
clc
im=1;
mi = Mixture.load('halo3')
mi.plot_covar
mi.plot_stdev


%% 
clc
zscore = (mi.pi_est-mi.pi_true)./mi.stdev
vx=[mi.pi_est mi.pi_true mi.stdev].*100
[(1:16)' vx zscore]



%% 
im=1;
mi.plot_zscores




%% 
figure(3)
sc=mi.covar./range(range(mi.covar));
sc=abs(log(sc));
cl=sc;
x=1:mi.num_models;
y=x;
plot(x,y,5,cl,'filled')
flabel('j','j','Information based Covar \pi')

%% 






figure(9)
[x,y]=meshgrid(1:mi.num_models,1:mi.num_models);
tri = delaunay(x,y);
z=[];
for i=1:length(x)
	for j=1:length(y)
		z(i,j) = mi.pi_est(i) + mi.pi_est(j);
	end
end
colormap hot
surf(x,y,-1*ones(length(x),length(x)),'CData',1./max(z,0.05),'EdgeColor','none','FaceColor','flat')




%% 
figure(7)
[x,y]=meshgrid(1:mi.num_models,1:mi.num_models);
tri = delaunay(x,y);
z = mi.covar;
trimesh(tri,x,y,z)


%% 
figure(5)
waterfall(mi.covar)


%% 
mi = Mixture(struct(	'save_as','halo3',...
						'model_path','data/modeldata3.dat',...
						'obs_path','data/obsdata3_10000.dat',...
						'pi_true', Mixture.get_pi_true(3),...
						'graph',false));

mi.em(struct('max_iters',100));
mi.update_stats
mi.save





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
























%% 
clc
pi = [1 2 3]'
f = [10 11 12 13;21 22 23 24;31 32 33 34]'
n=4;
m=3;

[I,S,V] = fisher_info(pi,f);
I
S
V

%% 

I = observed_info(pi,f);
S = observed_covar(I)

s_3_1 = S(1,1)+S(2,1)
s_3_2 = S(1,2)+S(2,2)

s_1_3 = S(1,1)+S(1,2)
s_2_3 = S(2,1)+S(2,2)

s_3_3 = S(1,1)+S(1,2)+S(2,1)+S(2,2)

SS=S;
SS(3,1)=s_3_1;
SS(3,2)=s_3_2;
SS(3,3)=s_3_3;
SS(1,3)=s_1_3;
SS(2,3)=s_2_3;

SS

SS-S
norm(SS-S)
%% 

sepr

denom = (f*pi).^2

d1_1 = sum(1./denom)

d1_1=0;
for i=1:n
	d1_1 = d1_1 + (f(i,1)/(denom(i)*f(i,1)));
end
d1_1

d1_1 = sum(f(:,1)./(denom.*f(:,1)))


d1_2 = sum(f(:,1)./(denom.*f(:,2)))

d1_2=0;
for i=1:n
	d1_2 = d1_2 + (f(i,1)/(denom(i)*f(i,2)));
end
d1_2



%% 
denom = zeros(n,1);

for i=1:n
	for k=1:m
		denom(i) = denom(i) + pi(k)*f(i,k);
	end
	denom(i) = denom(i)^-2;
end
denom


for i=1:n
	denom(i) = f(i,:)*pi;
end
denom


denom = (f*pi).^-2;

denom




