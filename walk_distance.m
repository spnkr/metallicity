function [results,mass_time_of_sats,DD,del] = walk_distance(mo,E,min_mass,max_mass,hlosm)

	data = mo.curve_by_mass_sorted_by_time(min_mass,max_mass);
	M = size(unique(data(:,6)),1)-1;
	N = size(data,1);
	
	results = {};
	
	sndx_ordered = [];
	mass_time_of_sats = [];
	
	%this must be slow...
	last_sndx = -1;
	tic
	for i=1:N
		if data(i,6) ~= last_sndx
			sndx_ordered(length(sndx_ordered)+1) = data(i,6);
			mass_time_of_sats(size(mass_time_of_sats,1)+1,:) = data(i,4:6);
			last_sndx = data(i,6);
		end
	end
	toc
	
	tic
	for i=1:M
		j = i+1;

		c1 = mo.data(mo.data(:,6)==sndx_ordered(i),:);
		c1 = sortrows(c1,1);
		c1(:,3) = c1(:,3)./sum(c1(:,3));

		c2 = mo.data(mo.data(:,6)==sndx_ordered(j),:);
		c2(:,3) = c2(:,3)./sum(c2(:,3));
		c2 = sortrows(c2,1);

		c1c = cumsum(c1(:,3));
		c2c = cumsum(c2(:,3));

		D = zeros(1,E);

		for e=1:(E-1)
			c1ndx = sum(c1c<=(e/E))+1;
			c2ndx = sum(c2c<=(e/E))+1;

			D(e) = sqrt((c1(c1ndx,1)-c2(c2ndx,1))^2+(c1(c1ndx,2)-c2(c2ndx,2))^2);
		end
		
		timei = mass_time_of_sats(mass_time_of_sats(:,3)==sndx_ordered(i),1);
		timej = mass_time_of_sats(mass_time_of_sats(:,3)==sndx_ordered(j),1);
		
		results{i} = struct('i',i,'j',j,'sat_i',sndx_ordered(i),'sat_j',sndx_ordered(j),...
			'sum',sum(D),'normalized_sum',sum(D)/E,'D',D,...
			'time_i',timei,...
			'mass_i',mass_time_of_sats(mass_time_of_sats(:,3)==sndx_ordered(i),2),...
			'time_j',timej,...
			'mass_j',mass_time_of_sats(mass_time_of_sats(:,3)==sndx_ordered(j),2),...
			'avg_time',(timei+timej)/2);
		strcat(['finished for ',num2str(i),'/',num2str(j)])
	end

	toc

	del = zeros(M,3);
	for k=1:M
		del(k,:) = [results{k}.normalized_sum results{k}.avg_time results{k}.j];
	end

	DD = sum(del(:,1));
	
	fig;
	subplot(1,2,1)
	plot(del(:,3),del(:,1),'k.-')
	title(strcat(hlosm, '. Distance (\Sigma=',num2str(DD),'), E=', num2str(E)));
	xlabel('Sat index (ordered by time; equally spaced)');
	ylabel('D(i,i-1)');



	subplot(1,2,2)
	plot(del(:,2),del(:,1),'b.-')
	title(strcat('Difference by avg time. Mass from ',...
		num2str(min(mass_time_of_sats(:,2))), ' to ', num2str(max(mass_time_of_sats(:,2)))));
	
	
	xlabel(strcat(['Time. From ' num2str(min(mass_time_of_sats(:,1))) ' to ' num2str(max(mass_time_of_sats(:,1)))]));
	ylabel('D(i,i-1)');
	
	%plot(1:M,log(del),'r.-');
	%title(strcat('Log distance, E=', num2str(E)));
	
	
	
	
	
	
	
	
	
	
	
	
	
	