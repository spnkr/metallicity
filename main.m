%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;


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

ndx = 1:length(x);

[x(ndx) y(ndx) freq(ndx) xbin(ndx) ybin(ndx) lta(ndx) eta(ndx) hma(ndx) lma(ndx)];





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
figure(2);
[X,Y] = meshgrid(-3:.05:0, -.7:.05:.9);
Z = zeros(size(X,1),size(X,2));


bin_data = models(models(:,4)==0 & models(:,5)==0,:);
%bin_data = bin_data(bin_data(:,3)>0,:);
x = bin_data(:,1);
y = bin_data(:,2);
freq = bin_data(:,3);

xb = -3:.2:0;
yb = -.7:.2:.9;

m = length(xb);
n = length(yb);

p = size(Z,1);
q = size(Z,2);

x_lower = min(xb);
y_lower = min(yb);

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


surf(X,Y,Z)
hold on
plot3(bin_data(:,1),bin_data(:,2),bin_data(:,3),'r.')
hold off
flabel('Fe/H','\alpha/Fe','F_{a,b}')



%% 
old
for i=1:m
	if i>1
		x_lower = xb(i-1);
	end
	x_upper = xb(i);
	for j=1:n
		if j>1
			y_lower = yb(j-1);
		end
		y_upper = yb(j);
		
		a = bin_data(	bin_data(:,1)>=x_lower &  bin_data(:,1)<x_upper...
						& ...
						bin_data(:,2)>=y_lower &  bin_data(:,2)<y_upper ...
					,3);
		
		if isfinite(a)
			Z(i,j) = a;
		end
	end
end



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







