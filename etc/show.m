function show(varargin)
	if length(varargin) > 0
		vs = varargin;
		vals = {};
		
		all_chars = true;
		for i=1:length(vs)
			v = vs(i);
			if ~ischar(v)
				all_chars = false;
			end
		end
		


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
	