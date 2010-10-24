%% load data etc
clear
clc
format short g
addpath('etc/');
addpath('tools/');
addpath('data/');
global im;
im=1;



%% 
im=1;
mu1=0;
var1=1;
mu2=1;
var2=2;

x = generate_mixture([mu1 var1],[mu2 var2],10000,0.7);
plot_mixture(x);
em_mixture(x, mu1, var1, mu2, var2);




