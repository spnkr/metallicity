function print_pi(p,mo)
	for i=1:length(p)
		disp(sprintf(strcat([num2str(i) '\t' num2str(100*p(i)) '\tT=' num2str(mo.props{i}(1)) ...
			'-' num2str(mo.props{i}(2)) '\tM=' num2str(mo.props{i}(3)) '-' num2str(mo.props{i}(4))   ])))
	end