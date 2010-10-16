%% regenerate all
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;

p_actual = [	14.467074  ;
				7.4354783  ;
				20.991296  ;
				9.2340355  ;
				33.655754  ;
				8.0191336  ;
				2.4001462  ;
				2.5734037  ;
			   0.11883542  ;
			  0.076990272  ;
			   0.35844400  ;
			   0.25313549  ;
			   0.22529346  ;
			  0.048890239  ;
			  0.024226458  ;
			  0.046833900  ];
p_actual = p_actual./100;



mo = Model.generate(0.01,100,struct('path','data/modeldata3.dat','include_blanks',false,...
	'normalize',false,'shift_data',false));

mo = Model.load('cache/models_3.mat');


ob = Observed(struct('name','halo_3_10k','path','data/obsdata2_10000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 10;
ob.x = ob.x + (0.1/2);
ob.y = ob.y + (0.1/2);
ob.load_models(mo);
Observed.save(ob);

ob = Observed(struct('name','halo_3_30k','path','data/obsdata2_30000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 30;
ob.x = ob.x + (0.1/2);
ob.y = ob.y + (0.1/2);
ob.load_models(mo);
Observed.save(ob);


ob = Observed(struct('name','halo_3_50k','path','data/obsdata2_50000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 50;
ob.x = ob.x + (0.1/2);
ob.y = ob.y + (0.1/2);
ob.load_models(mo);
Observed.save(ob);



mo = Model.generate(0.01,100,struct('path','data/modeldata5.dat','include_blanks',false,...
	'normalize',false,'shift_data',false));

mo = Model.load('cache/models_5.mat');

ob = Observed(struct('name','halo_5_30k','path','data/obsdata5_30000.dat','p_actual',NaN));
ob.model_number = 5;
ob.point_count = 30;
ob.x = ob.x + (0.1/2);
ob.y = ob.y + (0.1/2);
ob.load_models(mo);
Observed.save(ob);