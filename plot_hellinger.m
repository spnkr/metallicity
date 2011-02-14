function plot_hellinger(I,reps,halonum,mpaths)
	results = {};

	for i=1:length(mpaths)
		if halonum==-1
			fnn = mpaths{i};
			fnn2=fnn;
		else
			fnn = strcat(mpaths{i},num2str(halonum));
			fnn2 = strcat(mpaths{i}, ' Halo', num2str(halonum));
		end
		
		model_path = strcat(['data/halodata_v1/',fnn,'.dat']);
		[hd1,HD1] = hellinger(model_path,I,reps,1,2);
		[hd2,HD2] = hellinger(model_path,I,reps,2,3);
		[hd3,HD3] = hellinger(model_path,I,reps,3,4);

		res = struct('name',fnn2,'model_path',model_path,...
			'hd1',hd1,'HD1',HD1,'hd2',hd2,'HD2',HD2,'hd3',hd3,'HD3',HD3);

		results{length(results)+1} = res;
	end


	fig
	MM=length(results);
	for i=1:MM
		subplot(ceil(sqrt(MM)),ceil(sqrt(MM)),i)
	% 	plot(1:3,[results{i}.hd1 results{i}.hd2 results{i}.hd3],'r.')
		plot(repmat(1:3,1,reps),[results{i}.HD1 results{i}.HD2 results{i}.HD3],'bo')
		title(results{i}.name)
		ylabel('Hellinger distance of i,i+1')
		xlabel(strcat(['Model index. I=', num2str(I), '; ', num2str(reps), 'x']))
		axis([0 4 0 1])
	end