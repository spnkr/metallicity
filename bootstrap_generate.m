function bootstrap_generate(base_mb, mi, B, ns, varargin)
	load_args

	global pistarGL;

	models = mi.models;

	pistar = zeros(mi.num_models,B);

	for b=1:B
		tic
		bdata = bootstrap_generate_one(mi,ns);
		bmi = Mixture(struct(...
				'save_as',strcat(['temp/halo3_ni_boot_' num2str(ns) '_' num2str(b+base_mb)]),...
				'model_path',mi.model_path,...
				'x',bdata,...
				'pi_true', mi.pi_true));

		fig
		scatter(bmi.x(:,1),bmi.x(:,2),bmi.x(:,3)./sum(bmi.x(:,3)),'k','filled');
		flabel('Fe/H','\alpha/Fe',['Bootstrap ' num2str(b) ' generated data points']);
		axis([-3 1 -.7 1]);

		ps0 = bmi.em(cell2mat(varargin));
		pistar(:,b) = ps0;
		pistarGL = pistar;
		bmi.save;
		b
		toc
	end

