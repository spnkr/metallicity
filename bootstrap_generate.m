function ps0 = bootstrap_generate(base_mb, mi, B, ns, varargin)
	load_args
	do_m_n = arg('do_m_n',false);
	
	if do_m_n
		file_base = 'temp/halo3_new_mn_boot_';
		sepr('m-out-of-n');
	else
		file_base = 'temp/halo3_new_boot_';
		sepr('generated');
	end
	
	global pistarGL;

	models = mi.models;

	pistar = zeros(mi.num_models,B);

	for b=1:B
		tic
		if do_m_n
			bdata = bootstrap_generate_m_n(mi,ns);
		else
			bdata = bootstrap_generate_one(mi,ns);
			bdata = sortrows(bdata,4);
		end
		bmi = Mixture(struct(...
				'save_as',strcat([file_base num2str(ns) '_' num2str(b+base_mb)]),...
				'model_path',mi.model_path,...
				'x',bdata,...
				'pi_true', mi.pi_est));

		fig
		if size(bmi.x,2)>2
			scatter(bmi.x(:,1),bmi.x(:,2),bmi.x(:,3)./sum(bmi.x(:,3)),'k','filled');
		else
			scatter(bmi.x(:,1),bmi.x(:,2),10,'k','filled');
		end
		flabel('Fe/H','\alpha/Fe',['Bootstrap ' num2str(b) ' generated data points']);
		axis([-3 1 -.7 1]);

		ps0 = bmi.em(cell2mat(varargin));
		pistar(:,b) = ps0;
		pistarGL = pistar;
		bmi.save;
		b
		toc
	end

