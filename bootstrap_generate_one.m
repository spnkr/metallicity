function [data] = bootstrap_generate_one(mi, n, varargin)
	load_args
	
	if length(varargin)>0
		doplot=cell2mat(varargin(1));
	else
		doplot=false;
	end
	
	models = mi.models_sparse;
	m = size(models,1);
% 	P = mi.pi_true;
	P = mi.pi_est;
	M=size(P,1);
	grid_size = mi.bin_step;

	if ~isfinite(P);error('no true weights found');end


	PC = cumsum(P);

	data = zeros(n,5);
	R1 = rand(n,1);

	blur_factor=grid_size/2;

	RX = rand(n,1);
	ndx=randperm(length(RX));
	ndx=ndx(1:length(ndx)/2);
	RX(ndx)=RX(ndx)-1;
	RX = RX.*blur_factor;

	RY = rand(n,1);
	ndx=randperm(length(RY));
	ndx=ndx(1:length(ndx)/2);
	RY(ndx)=RY(ndx)-1;
	RY = RY.*blur_factor;

	for i=1:n
		try
		r0 = rand();
		pi_ndx = sum(PC<=r0)+1;
		if pi_ndx>M
			pi_ndx=M;
			disp(strcat(['bad rand ' num2str(i)]));
		end

		r1 = R1(i);
		nzd = models(models(:,4)==pi_ndx & models(:,3)>0 & models(:,1)>-3,1:3);

		cnzd = cumsum(nzd(:,3));
		nzda = [nzd cnzd];

		grid_ndx = sum(cnzd<=r1)+1;

		x = nzda(grid_ndx,1);
		y = nzda(grid_ndx,2);

		x = x + RX(i);
		y = y + RY(i);

		%fab = mo.f_ab(pi_ndx,x,y);
		fab = ones(size(x,1),size(x,2));
		data(i,:) = [x y fab pi_ndx r0];
		catch
			i
		end
	end

	data = data(data(:,3)>0,:);

	global im;
	if doplot
		fig;
		scatter(data(:,1),data(:,2),data(:,3)./sum(data(:,3)),'k','filled');
		flabel('Fe/H','\alpha/Fe',[num2str(n) ' generated data points']);
	end

end