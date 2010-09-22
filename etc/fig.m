function [out]= fig()
	global im;
	
	if im > 0
		figure(im);
		im=im+1;
		out=true;
	else
		out=false;
	end
	