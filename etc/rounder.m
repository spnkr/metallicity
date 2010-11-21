function y=rounder(x,varargin)
	if length(varargin)>0
		rb=cell2mat(varargin(1));
	else
		rb=100;
	end
	y = round(x.*rb)./rb;