function show(inputs)
	
	if length(inputs) > 0
		vs = inputs;
		vals = {};
		
		all_chars = true;
		for i=1:length(vs)
			v = vs(i);
			if ~ischar(v)
				all_chars = false;
			end
		end
		
		
		if all_chars || strcmp(class(inputs),'char')
			disp(strcat(inputs));
		else
			for i=1:length(vs)
				v = vs{i};
				if ischar(v)
					vals{i} = v;
				else
					vals{i} = num2str(v);
				end
			end
			disp(strcat(vals));
		end
		
		
	end
	