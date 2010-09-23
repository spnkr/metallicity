%% load data etc
clear
clc
format long g
addpath('etc/');
addpath('data/');
global im;
im=1;
load('cache/models_01.mat');


%% 
stepsize = 0.001;
[X,Y] = meshgrid(-3:stepsize:0, -.7:stepsize:.9);
Z = zeros(size(X,1),size(X,2));
p = size(Z,1);
q = size(Z,2);

profile on


tic
for i=1:p
	for j=1:q
		xs = X(i,j);
		ys = Y(i,j);

		a = mo.f_ab(1,xs,ys);

		if isfinite(a)
			Z(i,j) = a;
		end
	end
end
cache=toc



tic
 for i=1:p
 	for j=1:q
 		xs = X(i,j);
 		ys = Y(i,j);
 
 		a = mo.f_ab_live(1,xs,ys);
 
 		if isfinite(a)
 			Z(i,j) = a;
 		end
 	end
 end
live=toc


profile off
profile report


%% 
'comparing all'
sigma = 0;
stepsize = 0.5;
[X,Y] = meshgrid(-3:stepsize:0, -.7:stepsize:.9);
Z = zeros(size(X,1),size(X,2));
p = size(Z,1);
q = size(Z,2);

for i=1:p
	for j=1:q
		xs = X(i,j);
		ys = Y(i,j);

		a = mo.f_ab(1,xs,ys);
		b = mo.f_ab_live(1,xs,ys);
		v = norm(a-b);
		sigma = sigma + v;
		
		if v ~= 0
			err=[a b xs ys]
		end
	end
end

sigma
'done comparing all'


%% 
'generating 0.01 models'
tic
mo = Model(struct('step',0.01,'precision',100));
mo.save('cache/models_01.mat');
toc
'saved 01'

%% 
'generating 0.001 models'
tic
mo = Model(struct('step',0.001,'precision',1000));
mo.save('cache/models_001.mat');
toc
'done generating 0.001 model'



%% 
xs = -0.36
ys = 0.04
a = mo.f_ab(1,xs,ys)
b = mo.f_ab_live(1,xs,ys)
v = norm(a-b)

%% 

%% 
stepsize = 0.01;
[X,Y] = meshgrid(-3:stepsize:0, -.7:stepsize:.9);
Z = zeros(size(X,1),size(X,2));
p = size(Z,1);
q = size(Z,2);

profile on


for i=1:p
	for j=1:q
		xs = X(i,j);
		ys = Y(i,j);

		a = mo.f_ab(1,xs,ys);

		if isfinite(a)
			Z(i,j) = a;
		end
	end
end

% for i=1:p
% 	for j=1:q
% 		xs = X(i,j);
% 		ys = Y(i,j);
% 
% 		a = mo.f_ab_live(1,xs,ys);
% 
% 		if isfinite(a)
% 			Z(i,j) = a;
% 		end
% 	end
% end


profile off
profile report

%% 
obs = load('data/obsdata2.dat');
x = obs(:,1);
y = obs(:,2);
w = obs(:,3);
ndx = x>-3;
x=x(ndx);
y=y(ndx);
w=w(ndx);

%% 
scatter_color_by_weight(x,y,log(w),struct(...
	'title', 'halo002', 'x','Fe/H', 'y','\alpha/Fe','colormap','winter','dot_size',11));
axis([-3 0 -.5 1])



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



%% 















%% 
for i=1:25
subplot(5,5,i);
view(0,90);
axis([-3 0 -.5 1]);
end


%% 
grid_results = {};
stepsize = 0.01;
im=1;
fg = figure(im);
clf(fg,'reset')

spndx=1;
for n=min(ybin):max(ybin)
	for m=min(xbin):max(xbin)
		subplot(5,5,spndx);
		
		%--------------------------------------------
		
		[X,Y] = meshgrid(-3:stepsize:0, -.7:stepsize:.9);
		Z = zeros(size(X,1),size(X,2));


		bin_data = models(models(:,4)==m & models(:,5)==n,:);

		%[min(bin_data(:,9)) max(bin_data(:,9)) min(bin_data(:,8)) max(bin_data(:,8));...
		%min(bin_data(:,6)) max(bin_data(:,6)) min(bin_data(:,7)) max(bin_data(:,7))]

		x = bin_data(:,1);
		y = bin_data(:,2);
		freq = bin_data(:,3);

		p = size(Z,1);
		q = size(Z,2);



		for i=1:p
			for j=1:q
				xs = X(i,j);
				ys = Y(i,j);

				a = bin_data(	bin_data(:,1)==max(x(x<=xs))...
								& ...
								bin_data(:,2)==max(y(y<=ys)) ...
							,3);

				if isfinite(a)
					Z(i,j) = a;
				end
			end
		end
		
		AlphaData = Z;
		AlphaData(AlphaData>0) = 1;
		
		colormap winter
		surf(X,Y,Z,'AlphaData',AlphaData,'AlphaDataMapping','scaled','FaceAlpha','flat','EdgeAlpha',0,...
			'CData',log(Z));
		
% 		hold on
% 		plot3(bin_data(:,1),bin_data(:,2),bin_data(:,3),'r.')
% 		hold off
		flabel('Fe/H','\alpha/Fe','F_{a,b}')
		view(0,90);
		axis([-3 0 -.5 1]);
		
		obj = struct('time',strcat([num2str(min(bin_data(:,6))) '-' num2str(min(bin_data(:,7))) ' Gyr']),...
			'mass',strcat([num2str(min(bin_data(:,9))) '-' num2str(min(bin_data(:,8))) ' M']),...
			'X',X,'Y',Y,'Z',Z,'x',x,'y',y,'freq',freq,'data',bin_data);
		grid_results{spndx} = obj;
		spndx=spndx+1;
		
		%--------------------------------------------
		title(strcat(['' num2str(min(bin_data(:,6))) '-' num2str(min(bin_data(:,7))) ...
			' Gyr; M=' num2str(min(bin_data(:,9)))  '-' num2str(min(bin_data(:,8))) ]))
	end
end

save('cache/grid_results_0.001.mat','grid_results');


%% 
figure(3)
load('cache/grid_results_0.01.mat');
gx = grid_results{1};
AlphaData = gx.Z;
		AlphaData(AlphaData>0) = 1;
		
		colormap winter
		surf(gx.X,gx.Y,gx.Z,'AlphaData',AlphaData,'AlphaDataMapping','scaled','FaceAlpha','flat','EdgeAlpha',0,...
			'CData',log(gx.Z));
		
% 		hold on
% 		plot3(bin_data(:,1),bin_data(:,2),bin_data(:,3),'r.')
% 		hold off
		flabel('Fe/H','\alpha/Fe','F_{a,b}')
		view(0,90);
		axis([-3 0 -.5 1]);

%% 
figure(1);
[X,Y] = meshgrid(-3:.2:0, -.7:.2:.9);
%round because otherwise the grid can't use ==
X = round(X.*10)./10;
Y = round(Y.*10)./10;

bin_data = models(models(:,4)==0 & models(:,5)==0,:);
bin_data = bin_data(bin_data(:,3)>0,:);
%round because otherwise the grid can't use ==
bin_data(:,1) = round(bin_data(:,1).*10)./10;
bin_data(:,2) = round(bin_data(:,2).*10)./10;
%plot3(bin_data(:,1),bin_data(:,2),bin_data(:,3),'r.')

%Z = X .* exp(-X.^2 - Y.^2)
m = size(X,1);
n = size(Y,2);
Z = zeros(size(X,1),size(X,2));
for i=1:m
	for j=1:n
		a = bin_data(bin_data(:,1)==X(i,j) & bin_data(:,2)==Y(i,j),3);
		if isfinite(a)
			Z(i,j) = a;
		end
	end
end


surf(X,Y,Z)
hold on
plot3(bin_data(:,1),bin_data(:,2),bin_data(:,3),'r.')
hold off
flabel('Fe/H','\alpha/Fe','F_{a,b}')





%% 
clc
figure(1);
[X,Y] = meshgrid(1:1:2,1:1:5)
Z = X .* exp(-X.^2 - Y.^2)
sepr


%% 
im=1;
figure(1);
p=1;
for i=min(xbin):max(xbin)
	for j=min(ybin):max(ybin)
		subplot(5,5,p);
		p=p+1;
		ndx = models(:,4)==i & models(:,5)==j & models(:,3)>0;
		mn = min(freq(ndx));
		mx = max(freq(ndx));
		scf = 5 + 10.*(freq(ndx) - mn) ./ (mx-mn);
		scatter(x(ndx),y(ndx),scf,'filled','k');
		title(strcat(['T=' num2str(min(eta(ndx))) ',' num2str(min(lta(ndx))) ...
			'; M=' num2str(min(lma(ndx)))  ',' num2str(min(hma(ndx))) ]))
	end
end



%% 







