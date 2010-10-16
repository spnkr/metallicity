
clear all
fid = fopen('data/obsdata2_10000.dat','r'); x = fscanf(fid,'%f',[3,inf]); x= x';
figure(1); plot(x(:,1),x(:,2),'.'); xlabel('[\alpha/Fe]'); ylabel('[Fe/H]'); 

xmin = -3; xmax = 0; ymin = -0.5; ymax = 1;
idx = ((x(:,1) >= xmin) & (x(:,1) <= xmax) & (x(:,2) >= ymin) & (x(:,2) <= ymax));

x  = x(idx,:);
figure(2); plot(x(:,1),x(:,2),'.'); xlabel('[\alpha/Fe]'); ylabel('[Fe/H]'); 



fid = fopen('data/modeldata2.dat','r'); data = fscanf(fid,'%f',[9,inf]); data = data';
figure(3); plot(data(:,1),data(:,2),'.'); xlabel('[\alpha/Fe]'); ylabel('[Fe/H]'); 



% return;
r = sqrt(2*0.05*0.05);
t = (1/8:1/4:1)'*2*pi; x1 = r*sin(t); x2 = r*cos(t);

k = 0; s = 0; 
for i=0:4,
    for j =0:4,
        k = k + 1;
        idx = ((data(:,4)==j) & (data(:,5)==i) & (data(:,3)~=0));
        distx = data(idx==1,:); n = length(distx);
        subplot(5,5,k); axis([-3 0 -0.5 1]);
        if (n ~= 0) s = s + 1; end
        for l = 1:n,
            subplot(5,5,k); fill(distx(l,1)+x1,distx(l,2)+x2,[distx(l,3) 0 0]); hold on;
        end
        xlabel('[\alpha/Fe]'); ylabel('[Fe/H]'); axis([-3 0 -0.5 1]);
        % return;
    end
end
T = s;
xn = 30; yn = 17;
xgrid = -2.95:0.1:-0.05;
ygrid = -0.65:0.1:0.95;

% return
k=0; sq = zeros(xn,yn);
f = zeros(T,xgrid,ygrid);
indic= zeros(25,1); m = 0;
for i=0:4,
    for j =0:4,
        m = m + 1;
        idx = ((data(:,4)==j) & (data(:,5)==i));
        distx = data(idx==1,:); n = sum(distx(:,3)~=0);
        if (n==0) continue; end
        l = 0; k = k + 1;
        indic(m) = 1;
        for ix = 1:xn,
            for iy = 1:yn,
                l = l + 1;
                sq(ix,iy) = distx(l,3);
                f(k,ix,iy) = distx(l,3);
            end
        end
        figure(5); subplot(5,5,m); mesh(xgrid,ygrid,sq'); hold on;
    end
end

f = 100*f;
fval = EvalDens(x(:,1),x(:,2),f,xgrid,ygrid,T,xn,yn);
id = (sum(fval')~=0); fval = fval(id,:); 

% return
%The EM algorithm
u = gamrnd(1/T,1,T,1);
pi_0 = u/sum(u); 
pi_1 = ones(T,1)/T; pi_0 = pi_1 + 1;
eps = 0.0001;
k = 0;
N = length(fval);
w = zeros(N,T);

figure(10);
subplot(1,2,1);
plot(1,1);
hold on;

while (norm(pi_0 - pi_1) >= eps)
    k = k + 1;
    s = 0;
    for i = 1:N,
        w(i,:) = pi_1.*fval(i,:)'; 
        tot = sum(w(i,:));
        w(i,:) = w(i,:)/tot;
        s = s + log(tot);
    end
    pi_0 = pi_1; 
    pi_1 = sum(w)'/N
    subplot(1,2,1); plot(k,pi_1,'.');
	
    subplot(1,2,2); plot(k,s,'.');
	
    if(k > 500) break; end
end
hold off;
