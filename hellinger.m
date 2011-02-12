function [d_kl,D_KL] = hellinger(model_path,I,reps,pi_1,pi_2)
	new_halo_name = 'KLD';%saved as this
	obs_path = 'data/obsdata_s4t2m720.dat';
	pi_true = NaN;


	
	im=1;
	mi0 = Mixture(struct(	'save_as',new_halo_name,...
							'model_path',model_path,...
							'obs_path',obs_path,...
							'pi_true', pi_true,...
							'graph',false));


	xr = range(mi0.xrange)+mi0.bin_step;
	yr = range(mi0.yrange)+mi0.bin_step;
	minx = min(mi0.xrange)-(mi0.bin_step/2);
	miny = min(mi0.yrange)-(mi0.bin_step/2);
	
	

	
	D_KL = zeros(1,reps);
	
	for i=1:reps
		
		xyc = repmat([xr yr],I,1);
		xyc = rand(I,2).*xyc + repmat([minx miny],I,1);

		mi = Mixture(struct(	'save_as',new_halo_name,...
								'model_path',model_path,...
								'x',xyc,...
								'pi_true', pi_true,...
								'graph',false));


		d_kl = 0;
		K = size(mi.x,1); %num valid (non-zero) points
		for k=1:K
			c = (sqrt(mi.f(k,2+pi_1))-sqrt(mi.f(k,2+pi_2))).^2;
			if isfinite(c)
				d_kl = d_kl + c;
			end
		end


		d_kl = d_kl/(2*I);
		D_KL(i) = d_kl;
	end
	
	D_KL'
	
	d_kl = mean(D_KL);
	
	sepr(['H(p,q) = ',num2str(d_kl),' for I=',num2str(I),', repeated ',num2str(reps),' times for ',num2str(pi_1),' v. ',num2str(pi_2),' models.  [',model_path,']'])