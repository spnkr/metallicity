function [out]= fig()
	global im;
	global constant_im;
	
	if im > 0
		figure(im);
		if ~constant_im
			im=im+1;
		end
		out=true;
	else
		out=false;
	end
	