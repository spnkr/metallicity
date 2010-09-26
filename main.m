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
mo = Model.generate(0.01,100,struct('save_to','cache/model_temp.mat'));
mo.test_cache(struct('model_no',1,'step_size',0.005));

mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

%% 
ob = Observed(struct('name','halo002'));
ob.load_models(mo);
ob.save('cache/observed.mat');

ob.plot();


%% 
try_this_many_times=10;
all_p = zeros(25,try_this_many_times);
for i=1:try_this_many_times
im=i+1;
[p,P,norms,plike,clike] = ob.em(struct(	'n',2000,...
						'min_norm',0.005,...
						'max_iters',50,...
						'init','rand(m,1)',...
						'interactive',false));
all_p(:,i) = p;
end
all_p


%% 
im=1;
[p,P,norms,plike,clike] = ob.em(struct(	'XXn',NaN,...
						'min_norm',0.00001,...
						'max_iters',50,...
						'init','rand(m,1)',...
						'interactive',true));




					
					

%% 
figure(2)
scatter(1:25,p)
for i=1:25
	text(i,p(i),strcat([' ' num2str(mo.props{i}(1)) '-' num2str(mo.props{i}(2)) 'Gyr']),'FontSize',8)
	disp(sprintf(strcat([num2str(p(i)*100) '\t' num2str(mo.props{i}(1)) '-' num2str(mo.props{i}(2)) 'Gyr\t' ...
		num2str(mo.props{i}(3)) '-' num2str(mo.props{i}(4)) 'M'])))
end

%% 
figure(2)
plot(all_p(1,:),'k.-')
hold on
for i=2:size(all_p,2)
	plot(all_p(i,:),'k.-')
end
hold off
flabel('Trial','\pi_j','Random starting values, n=2000');


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
vol=0;
k=1;
vols=[];
dx = mo.data{1}(mo.data{1}(:,3)>0,:)
for i=1:size(dx,1)
	vol_of_this_region = (((.05)^2) * dx(i,3));
	vol = vol + vol_of_this_region;
	vols(k,:) = [dx(i,1) dx(i,2) dx(i,3) vol_of_this_region];
	k=k+1;
end
'x,y,density,volume'
vols
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

















