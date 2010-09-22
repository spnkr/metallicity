%if there were extra parameters passed put them into an object called args
%so that the function arg(key,[default value]) works.
function load_args()
	if evalin('caller', 'exist(''varargin'')')
		assignin('caller', 'args', evalin('caller', 'cell2mat(varargin)'));
	else
		assignin('caller', 'args', struct());
	end