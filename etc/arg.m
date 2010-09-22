%loads the paramter of the given name from the variable called args in
%scope. returns NaN if it didn't exist, or the second param (default val)
function [v,key] = arg(key, varargin)
	if evalin('caller',['isfield(args, ''', key, ''')'])
		v = evalin('caller',['args.', key]);
	else
		if(length(varargin)>0)
			v = cell2mat(varargin(1));
		else
			v = NaN;
		end
	end
	