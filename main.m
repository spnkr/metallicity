%%
matlabpool open
%% 
matlabpool close

%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('data/');
global im;
im=1;
multicount = 2;
max_seconds = 60*60;


[ob3, mo3] = Observed.load(3,10);
[ob5, mo5] = Observed.load(5,30);

init_p_val = ones(16,1)./16;
p_my_target_3 = [    0.13093;   0.0011212;     0.23474;    0.067527;     0.32525;    0.083172;    0.016957;    0.031472;  1.0601e-36;   0.0007638;   0.0027027;  2.9913e-10;   0.0036901;  1.1769e-09; 4.5421e-115;   0.0022481];
p_bs_target_3 = [    0.15993;    0.056266;     0.20149;     0.09207;     0.34038;    0.089855;    0.025727;    0.021953;  1.9384e-37; 6.7372e-119;   0.0042462;   0.0018148;   0.0054376;  6.0639e-09;  0.00016809;  0.00066224];


p_5_my_best = [    53.849;      7.423;     22.389;      8.459;     1.5405;    0.50683;    0.31515;     3.2273;    0.64169;    0.31798;    0.40411;    0.51391;    0.37232;   0.039342; 0.00085221];


%% 



%% 
clc
constl=9222;
unconstl=9228.5;

t = -2*(uc-c);

pval = 1-chi2cdf(t,15)


%%
im=1;
p_b=ob2.em(struct(	'max_seconds',60,...
							'p',init_p_val,...
							'interactive',true,...
							'interactive_print_interval',1));
p_b





%% 
im=2;
p_b=ob5.em(struct(	'max_seconds',60,...
							'p',init_p_val,...
							'interactive',true,...
							'interactive_print_interval',1));
p_b


%% 
im=2;
tic
p_me=ob.em(struct(	'max_iters',2156,...
					'interactive',true,...
					'p', init_p_val,...
					'baseline_p',p_bs_target_3));
p_me
toc

%% 
im=3;
tic
p_me2=ob.em(struct(	'max_iters',5000,...
					'interactive',true,...
					'p', init_p_val,...
					'baseline_p',p_bs_target_3));
p_me2
toc


%% templates: em
[ob, mo] = Observed.load(3,10);
[p,P,ll] = ob.em(struct(...
				'max_iters',2,...
				'Xmax_seconds',60*60,...
				'init','rand(m,1)',...
				'interactive',true));
print_pi(p,ll,mo);
ob.plot_differences(p);


%% template: em multi
[ob, mo] = Observed.load(3,10);
multicount = 2;
max_seconds = 5;
max_iters = 100;

sepr('start blah')
ob.em_multi(mo, struct( 'count',multicount,...
						'max_seconds', max_seconds, ...
						'interactive',false));
done('end blah');






%% 
multicount = 10;
max_seconds = 60*20;
[ob, mo] = Observed.load(3,10);
	sepr('start model 3 true')
		init_p_to_use = ob.p_actual;
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',init_p_to_use,...
								'save',strcat(['cache/em_multi_m3_true_' num2str(max_seconds) '_' num2str(multicount) 'x.mat'])));
		done('end model 3 true');


		

multicount = 2;
max_seconds = 60*60;
[ob, mo] = Observed.load(5,30);
	sepr('start model 5')
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false));
		done('end model 5');






multicount = 8;
max_seconds = 7*60;
[ob, mo] = Observed.load(3,10);
	sepr('start model 3 true pm 3')
		init_p_to_use = ob.p_actual.*0.03.*rand(size(ob.p_actual,1),1);
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false,...
								'p',init_p_to_use,...
								'save',strcat(['cache/em_multi_m3_true_pm3_' num2str(max_seconds) '_' num2str(multicount) 'x.mat'])));
		done('end model 3 true pm 3');
	sepr('start model 3 rand')
		ob.em_multi(mo, struct( 'count',multicount,...
								'max_seconds', max_seconds, ...
								'interactive',false));
		done('end model 3 rand');







max_seconds = 60*60;
[ob, mo] = Observed.load(3,30);
	sepr('start model 3 30k rand')
		[p,P,ll] = ob.em(struct(...
				'max_seconds',max_seconds,...
				'init','rand(m,1)',...
				'interactive',false));
		print_pi(p,ll,mo);
		ob.plot_differences(p);
		done('end model 3 30k rand');


max_seconds = 60*60*5;
[ob, mo] = Observed.load(3,30);
	sepr('start model 3 30k mega')
		[p,P,ll] = ob.em(struct(...
				'max_seconds',max_seconds,...
				'init','rand(m,1)',...
				'interactive',false));
		print_pi(p,ll,mo);
		ob.plot_differences(p);
		save('cache/mega_30k_run.mat','p','P','ll');
		done('end model 3 30k mega');











%% simulate
[p,P,ll,LL] = ob.simulate(mo,struct('sample',10000,...
											'max_seconds',60*30,...
											'interactive',false));




							
%% em convergence
%%
ob.plot_weight_changes(struct('path','cache/all_p_50_full_runs.mat'));


%% plots
%% 
mo.plot_single(struct('model',1,'step_size',0.01,'precision',100))
mo.plot();
mo.plot(struct('overlay',true,'observed',ob));
mo.plot_single(struct('model',1));

%% 
ob.plot();



















%% shifts etc
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

grid_size_fc = 0.1;

mo = Model.generate(0.01,100,struct('path','data/modeldata3.dat','include_blanks',false,...
	'normalize',false,'shift_data',false,'model_number',3));

mo = Model.load('cache/models_3.mat');


ob = Observed(struct('name','halo_3_10k','path','data/obsdata2_10000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 10;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);

ob = Observed(struct('name','halo_3_30k','path','data/obsdata2_30000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 30;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);


ob = Observed(struct('name','halo_3_50k','path','data/obsdata2_50000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 50;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);




%% 
im=1;
tic
p_b=ob.em_bodhi(struct(	'max_seconds',60*60,...
						'baseline_p',p_bs_target_3,...
						'p',init_p_val,...
						'interactive',true,...
						'model_path','data/modeldata2.dat',...
						'obs_path','data/obsdata2_10000.dat'));
p_b
save('cache/p_b_1h.mat','p_b')
toc



mo = Model.generate(0.01,100,struct('path','data/modeldata3.dat','include_blanks',false,...
	'normalize',false,'shift_data',false,'model_number',3));

mo = Model.load('cache/models_3.mat');


ob = Observed(struct('name','halo_3_10k','path','data/obsdata2_10000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 10;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);


init_p_val = ones(16,1)./16;
p_my_target_3 = [    0.13093;   0.0011212;     0.23474;    0.067527;     0.32525;    0.083172;    0.016957;    0.031472;  1.0601e-36;   0.0007638;   0.0027027;  2.9913e-10;   0.0036901;  1.1769e-09; 4.5421e-115;   0.0022481];
p_bs_target_3 = [    0.15993;    0.056266;     0.20149;     0.09207;     0.34038;    0.089855;    0.025727;    0.021953;  1.9384e-37; 6.7372e-119;   0.0042462;   0.0018148;   0.0054376;  6.0639e-09;  0.00016809;  0.00066224];
[ob, mo] = Observed.load(3,10);

im=2;
tic
p_me=ob.em(struct(	'max_seconds',60*60,...
					'interactive',true,...
					'p', init_p_val,...
					'baseline_p',p_bs_target_3));
p_me
toc
save('cache/p_3_me_1h_eql_start.mat','p_me')
sepr('done em')






im=3;
tic
p_me_rnd=ob.em(struct(	'max_seconds',60*60,...
					'interactive',true,...
					'baseline_p',p_bs_target_3));
p_me
toc
save('cache/p_me_1h_rnd.mat','p_me_rnd')
sepr('done em w rand')









mo = Model.generate(0.01,100,struct('path','data/modeldata3.dat','include_blanks',false,...
	'normalize',true,'shift_data',false,'model_number',3));
ob = Observed(struct('name','halo_3_10k','path','data/obsdata2_10000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 10;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);


im=4;
tic
p_me_nrm=ob.em(struct(	'max_seconds',60*60,...
						'interactive',true,...
						'p', init_p_val,...
						'baseline_p',p_bs_target_3));
p_me_nrm
toc
save('cache/p_3_me_1h_normalized_fixed.mat','p_me_nrm')
sepr('done normalized em')



im=5;
tic
p_me_nrm_rnd=ob.em(struct(	'max_seconds',60*60,...
							'interactive',true,...
							'baseline_p',p_bs_target_3));
p_me_nrm
toc
save('cache/p_3_me_1h_normalized_rand.mat','p_me_nrm_rnd')
sepr('done normalized em rand')




sepr('start model 5')
mo = Model.generate(0.01,100,struct('path','data/modeldata5.dat','include_blanks',false,...
	'normalize',false,'shift_data',false,'model_number',5));

mo = Model.load('cache/models_5.mat');

ob = Observed(struct('name','halo_5_30k','path','data/obsdata5_30000.dat','p_actual',NaN));
ob.model_number = 5;
ob.point_count = 30;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);



%% 


[ob, mo] = Observed.load(5,30);

ob.complete_likelihood(ob.p_actual,n,m);



%% 
im=6;
tic
p_me_5=ob.em(struct(	'max_seconds',60*60,...
					'interactive',true));
p_me_5
toc
save('cache/p_5_me_1h_eql_start.mat','p_me_5')
sepr('done em for 5')






sepr('done all')














%% 
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

grid_size_fc = 0.1;

% 
% mo = Model.generate(0.01,100,struct('path','data/modeldata3.dat','include_blanks',false,...
% 	'normalize',true,'shift_data',false,'model_number',3));
ob = Observed(struct('name','halo_3_10k','path','data/obsdata2_10000.dat','p_actual',p_actual));
ob.model_number = 3;
ob.point_count = 10;
ob.x = ob.x + (grid_size_fc/2);
ob.y = ob.y + (grid_size_fc/2);
ob.load_models(mo);
Observed.save(ob);


%% 
init_p_val = ones(16,1)./16;
[ob, mo] = Observed.load(3,10);


im=3;
tic
p_me=ob.em(struct(	'max_iters',100,...
					'interactive',true,...
					'p', init_p_val,...
					'model_path','data/modeldata2.dat',...
					'obs_path','data/obsdata2_10000.dat'));
p_me
toc





im=2;
tic
p_me=ob.em(struct(	'max_iters',2156,...
					'interactive',true,...
					'p', init_p_val,...
					'model_path','data/modeldata2.dat',...
					'obs_path','data/obsdata2_10000.dat'));
p_me
toc






%% 
[ob, mo] = Observed.load(5,30);
im=3;
tic
p_me=ob.em(struct(	'max_iters',2156,...
					'interactive',true,...
					'model_path','data/modeldata5.dat',...
					'obs_path','data/obsdata5_30000.dat'));
p_me
toc






