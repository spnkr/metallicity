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

[ob2, mo2] = Observed.load(2,10);
[ob5, mo5] = Observed.load(5,30);

init_p_val = ones(16,1)./16;

%% model 2 realization 2
im=1;
p_b=ob2.em(struct(	'max_seconds',60,...
					'p',init_p_val,...
					'interactive',true,...
					'interactive_print_interval',1));
p_b


%% model 5
im=2;
p_b=ob5.em(struct(	'max_seconds',60,...
					'p',init_p_val,...
					'interactive',true,...
					'interactive_print_interval',1));
p_b




























%% 
clc
constl=9222;
unconstl=9228.5;

t = -2*(uc-c);

pval = 1-chi2cdf(t,15)




