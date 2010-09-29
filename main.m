%%
matlabpool open
%% 
matlabpool close

%% load data etc
clear
clc
format long g
addpath('etc/');
addpath('data/');
global im;
im=1;
mo = Model.load('cache/models_01.mat');
load('cache/observed.mat');



%% 
mo = Model.generate(0.01,100,struct('normalize',false,'save_to','cache/models_temp.mat'));
mo.test_cache(struct('model_no',1,'step_size',0.005));

%% 
mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

mo.volume(1)

%% 
ob = Observed(struct('name','halo002'));
ob.load_models(mo);
ob.save('cache/observed_temp.mat');

%% 
ob.plot();
[p,P,norms,plike,clike] = ob.em(struct(	'n',100,...
										'min_norm',0.005,...
										'max_iters',25,...
										'min_iters',5,...
										'init','rand(m,1)',...
										'interactive',true));
p

%% 


%% 
mo = Model.load('cache/models_01_normalized.mat');
load('cache/observed_normalized.mat');

%% 
mo = Model.load('cache/models_01.mat');
load('cache/observed.mat');


try_this_many_times=50;
all_p = zeros(25,try_this_many_times);
for i=1:try_this_many_times
	im=1;
	[p,P,norms,plike,clike] = ob.em(struct(	'Xn',NaN,...
							'min_norm',0.001,...
							'max_iters',75,...
							'Xmin_iters',10,...
							'init','rand(m,1)',...
							'interactive',false));
	all_p(:,i) = p;
	disp(strcat(['Finished run ' num2str(i)]))
end
all_p
sepr('finished all')

%% 
im=1;
[p,P,norms,plike,clike] = ob.em(struct(	'XXn',NaN,...
						'min_norm',0.00001,...
						'max_iters',50,...
						'init','rand(m,1)',...
						'interactive',true));




					
					

%% 
figure(1)
scatter(1:25,p)
for i=1:25
	text(i,p(i),strcat([' ' num2str(mo.props{i}(1)) '-' num2str(mo.props{i}(2)) 'Gyr']),'FontSize',8)
	disp(sprintf(strcat([num2str(p(i)*100) '\t' num2str(mo.props{i}(1)) '-' num2str(mo.props{i}(2)) 'Gyr\t' ...
		num2str(mo.props{i}(3)) '-' num2str(mo.props{i}(4)) 'M'])))
end

%% 
fg=figure(2)
clf(fg)
m = size(all_p,1);
p = size(all_p,2);

clrs = rand(m,3);

no_std_devs = 2;

subplot(1,2,1);
hold on
for i=1:m
	plot(all_p(i,:),'.','Color',clrs(i,:))
	
	plot(1:p,mean(all_p(i,:)).*ones(1,p),'-','Color',clrs(i,:))
	
	plot(1:p,(no_std_devs*std(all_p(1,:))+mean(all_p(i,:))).*ones(1,p),'--','Color',clrs(i,:))
	plot(1:p,(-no_std_devs*std(all_p(1,:))+mean(all_p(i,:))).*ones(1,p),'--','Color',clrs(i,:))
end


hold off
flabel('Trial','\pi_j','Random starting values, +/- 2\sigma, n=all');


subplot(1,2,2); 
vls = zeros(m,5);
hold on
for i=1:m
	dv=all_p(i,:);
	mdv = mean(dv);
	sdv = std(dv);
	vls(i,:) = [min(dv) -sdv+mdv mdv mdv+sdv max(dv)];
% 	plot(i,vls(i,:),'.','Color',clrs(i,:))
	errorbar(i,mdv,sdv,'.','Color',clrs(i,:))
end
hold off
flabel('j','\pi_j','Error bars for final weights over 50 random starts, n=all');
%plot(1:size(vls,1),vls,'.','Color',clrs(i,:))


%% 
figure(3)

plot([1 2 3 4],[13 20 38 51],'k.-')


%% 
max_ob_ndx=10;
ob.x = ob.x(1:max_ob_ndx);
ob.y = ob.y(1:max_ob_ndx);
ob.lum = ob.lum(1:max_ob_ndx);

%% 
max_iters = 20;
m = length(mo.data);
n = length(ob.x);

p = (1/m).*ones(m,1);
P = NaN.*ones(m+1,max_iters);%first row is norm

w0 = zeros(n,m); %col 1 is for p1, col2 for p2, etc.
w = zeros(n,m);

min_norm = 0.01;

tic
for counter=1:max_iters
	p0 = p;
	
	for j=1:m
		for i=1:n
			p_f_ak_bk = zeros(m,1);
			
			x = ob.x(i);
			y = ob.y(i);

			for k=1:m
				p_f_ak_bk(k) = p0(k).*mo.f_ab(k,x,y);
			end

			w(i,j) = (p0(j).*mo.f_ab(j,x,y)) ./ sum(p_f_ak_bk);
		end
		
		p(j) = sum(w(:,j))/n;
	end
	
	P(:,counter) = [norm(p0-p);p]
	
	if abs(norm(p0-p)) < min_norm
		'min norm reached; stopping'
		break;
	end
end
tmr = toc;

norms = P(1,isfinite(P(1,:)))';
figure(1)

plot(norms,'k.-');
legend(strcat(['norm=' num2str(min(norms))]));

flabel('trial','',strcat([num2str(counter) ' runs in ' num2str(tmr) 's']));
%% 







































%% volume checks
for i=1:25
	x=mo.nonzeros(i);sum(x(:,3))
end



%% 
vol=0;
k=1;
vols=[];
dx = mo.nonzeros(k);
grid_size = .05;

if 1==1
dx = mo.cache{k};
dx = dx(dx(:,1)>0,[2 3 1]);
grid_size = mo.stepsize;
end

for i=1:size(dx,1)
	vol_of_this_region = (((grid_size)^2) * dx(i,3));
	vol = vol + vol_of_this_region;
	vols(k,:) = [dx(i,1) dx(i,2) dx(i,3) vol_of_this_region];
	k=k+1;
end
vol


%% 
models = load('data/modeldata2.dat');
x = models(:,1);
y = models(:,2);
freq = models(:,3);
xbin = models(:,4);
ybin = models(:,5);
lta = models(:,6);
eta = models(:,7);
hma = models(:,8);
lma = models(:,9);

%binned on a 5x5 grid
XBIN = 0;
YBIN = 0;

%get data for the appropriate bin
data = models(models(:,4)==XBIN & models(:,5)==YBIN,[1 2 3]);

%points that have mass
non_zero_data = data(data(:,3)>0,:);

%find volume
vol=0;
k=1;
vols=[];
for i=1:size(non_zero_data,1)
	vol_of_this_region = (((.05)^2) * non_zero_data(i,3));
	vol = vol + vol_of_this_region;
	vols(k,:) = [non_zero_data(i,1) non_zero_data(i,2) non_zero_data(i,3) vol_of_this_region];
	k=k+1;
end
'x,y,density,volume'
vols
vol

















