classdef Galaxy < handle
	
	properties
		data; %[vlos r theta x1 x2 err pm tub sel_ndx];
		
		rotated_by=0;
		r_core=0;
		include_unbound=true;
		
		vel_disp;
		name,color,n;
		
		
		rotation_candidates;
		
		major_axis_pa;	%loaded data is already corrected for this; informational only
		
		
		omega=0; % from bisector test
		r0;%bisector
		
		
		bisector_w_str = 'pm./(err.^2+sigma_0^2)';
		
		%modelling ones
		nadwat_show_data_points=true;
		nadwat_2d_b,nadwat_1d_b,nadwat_step;
		nadwat_2d_num_reductions;
	end
	
	
	methods(Static)
		function [leo, sculptor, fornax, draco, sextans, carina, galaxies] = load_all()
			leo = load_leo_data();
			sculptor = load_sculptor_data();
			fornax = load_fornax_data();
			draco = load_draco_data();
			sextans = load_sextans_data();
			carina = load_carina_data();
			galaxies = [leo;sextans;sculptor;fornax;draco;carina];
		end
		
		function [sim1, sim2, sim3, galaxies] = load_all_sim(varargin)
			load_args
			r_core = arg('r_core',0);
			
			va = struct('file','');
			va.plane = arg('plane', 'vel');
			va.file = '100';
			
			
			[x1, x2, th, ra, vlos, tub, pm] = load_sim(va);
			sim1 = Galaxy('Sim 1',0,vlos,ra,th,zeros(length(vlos),1),x1,x2,pm,tub);
			
			
			
			va.file = '110';
			[x1, x2, th, ra, vlos, tub, pm] = load_sim(va);
			sim2 = Galaxy('Sim 2',0,vlos,ra,th,zeros(length(vlos),1),x1,x2,pm,tub);
			
			
			
			va.file = '120';
			[x1, x2, th, ra, vlos, tub, pm] = load_sim(va);
			sim3 = Galaxy('Sim 3',0,vlos,ra,th,zeros(length(vlos),1),x1,x2,pm,tub);
			
			galaxies = [sim1;sim2;sim3];
		end
		
		
		function [data] = sort_by_vals(sb,data)
			if strcmp(sb,'r')
				sndx = 2;
			elseif strcmp(sb,'x1')
				sndx = 4;
			end
			
			data = sortrows(data,sndx);
		end
	end
	
	
	methods
		function g = Galaxy(name,vel_disp,vlos,r,theta,err,x1,x2,pm,tub)
			g.name=name;
			g.vel_disp=vel_disp;
			
			%x1,x2,vlos,r,err,theta,pm;
			
			if ~isfinite(tub)
				tub = zeros(length(r),1);
			end
			
			ndx = ones(length(r),1);
			
			g.data = [vlos r theta x1 x2 err pm tub ndx];
		end
		

	
		function a = vlos(g);a = g.data(logical(g.data(:,9)),1);end
		function a = r(g);a = g.data(logical(g.data(:,9)),2);end
		function a = theta(g);a = g.data(logical(g.data(:,9)),3);end
		function a = x1(g);a = g.data(logical(g.data(:,9)),4);end
		function a = x2(g);a = g.data(logical(g.data(:,9)),5);end
		function a = err(g);a = g.data(logical(g.data(:,9)),6);end
		function a = pm(g);a = g.data(logical(g.data(:,9)),7);end
		function a = tub(g);a = g.data(logical(g.data(:,9)),8);end

		
		%
		%
		% bisector
		%
		%
		function [nu, sigma] = null_nu_sigma_0(g)
			[nu, sigma] = null_nu_sigma_0(g.vlos, g.err, g.pm, g.vel_disp, g.bisector_w_str);
		end
		
		function [omega, b_obs, r0, sig_level, all, best] = bisector(g, rho_range, varargin)
			
			%if you change this then you also need to change the
			%probability memberships to be all 1's so that the nu, sigma
			%calculation will be 1/n, not a weighted sum
			%w_str = 'pm./(err.^2+sigma_0^2)';
			
			va = cell2mat(varargin);
			va.name = g.name;
			
			[b_obs, omega, r0, sig_level, all, best] = bisector(rho_range, g.bisector_w_str,...
				g.vlos, g.theta, g.r, g.err, g.pm, g.vel_disp, va);
			g.omega = omega;
			g.r0 = r0;
			
		end
		
		
		
		
		
		
		
		
		%
		%
		% nadwat
		%
		%		
		%Nadaraya-Watson kernel regression.
		%pass a number second to use a different bandwidth
		%set self.nadwat_show_data_points to show scaled data points
		function out = nadwat(g, varargin)
			load_args
			
			b = arg('b', g.nadwat_2d_b);
			shrink = arg('shrink', g.nadwat_2d_num_reductions);
			cropped_alpha = arg('cropped_alpha', 0);
			step = arg('step', g.nadwat_step); %10 is ok
			show_data_points = arg('show_data_points', true);
			
			ndx = g.pm>0;
			
			px1 = g.x1;
			px1 = px1(ndx);
			px2 = g.x2;
			px2 = px2(ndx);
			pvlos = g.vlos;
			pvlos = pvlos(ndx);
			
			out = nadwat(px1, px2, pvlos, b,...
				step, g.name, show_data_points, cropped_alpha, shrink);
		end
		
		%first opt is bandwidth, set to -1 to use default. second opt is
		%the step size, which will control accuracy of residuals. default
		%is 10.
		function [xh,yh] = nadwat1d(g, varargin)
			%step size of 10 is fast and accurate enough
			if length(varargin)>0 & cell2mat(varargin(1)) ~= -1
				b=cell2mat(varargin(1));
			else
				b=g.nadwat_1d_b;
			end
			
			if length(varargin)>1
				nws = cell2mat(varargin(2));
			else
				nws=g.nadwat_step;
			end
			
			[xh,yh] = nadwat1d(g.x1, g.vlos, b, nws, g.name, g.nadwat_show_data_points);
		end
		
		
		
		
		
		function [sse_r, sse_rcos, sse_x1_cos, sse_xi_gamma] = find_sses(g,om_step,varargin)
			load_args
			omegas = arg('omegas',0:om_step:2*pi);
			
			global im;
			%xim=im;
			%im=-1;
			
			sse_r = g.isotonic_r(struct('truncate',[10 10]));
			[sse_rcos, lam, nu_rcos, omega_rcos] = g.isotonic_rtheta(struct('truncate',[10 10],'omegas',omegas));
			[sse_x1_cos, lam, vd2, omega_x1_cos] = g.isotonic_x1(struct('truncate',[10 10],'omegas',omegas));
			[sse_xi_gamma, yhat, nu_xi_gamma, xi, gamma, omega_xi_gamma] = g.isotonic_xi_gamma(struct('truncate',[10 10],'omegas',omegas));
			
			%im=xim;
			
			strcat(['SSEs for ' g.name ': r:' num2str(round(sse_r)) ', rcos:' num2str(round(sse_rcos)) ...
				', x1cos:' num2str(round(sse_x1_cos)) ', xi_gamma:' num2str(round(sse_xi_gamma)) ''])
			
			strcat(['nus for ' g.name ': rcos:' num2str(round2(nu_rcos)) ', xi_gamma:' num2str(round2(nu_xi_gamma))])
			
			strcat(['omegas for ' g.name ': rcos:' num2str(round2(omega_rcos)) ...
				', x1cos:' num2str(round2(omega_x1_cos)), ', xi_gamma:' num2str(round2(omega_xi_gamma))])
			
			%[sse_r, sse_rcos, sse_x1_cos, sse_xi_gamma]
		end
		
		
		
		
		
		%
		%
		% testing
		%
		%
		function [pvalue, sse_obs, sses] = test_permute_rtheta(g, varargin)
			global im;
			load_args
			
			omegas = arg('omegas',0);
			truncate = arg('truncate', [0 10]);
			vd2 = arg('vel_disp',g.vel_disp)^2;
			num_bins = arg('num_bins', 10);
			use_nu = true;
			
			prefix = arg('prefix','');
			silent = arg('silent',false);
			fast = arg('fast', false);
			
			omega_free = arg('omega_free', true);
			omega_per_bin_free = arg('omega_per_bin_free', false);
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			
			if silent
				doplot=false;
			elseif gui
				doplot=true;
			else
				if fig
					doplot=true;
				else
					doplot=false;
				end
			end
			
			r = g.r;
			%--
			g.sort_by('r');
			iterations = arg('iterations', 100);
			
			vx = cell2mat(varargin);
			sse_obs = g.isotonic_rtheta(vx);
			
			sses = zeros(iterations,1);
			
			
			r=g.r;
			r = r(g.pm>0);
			
			bin_size = floor(length(r)/num_bins);
			
			if omega_per_bin_free
				for i=1:iterations
					sse_binned = 0;
					for k=1:num_bins
						st = ((k-1)*bin_size)+1;
						en = k*bin_size;

						if k==num_bins
							en = length(r);
						end

						bin_ndx = st:en;
						m = length(bin_ndx);

						perm_ndx = randperm(m);

						bin_theta = g.theta;
						bin_theta = bin_theta(bin_ndx);
						bin_theta = bin_theta(perm_ndx);
						
						pvlos = g.vlos;
						pvlos = pvlos(bin_ndx);
						
						pr = g.r;
						pr = pr(bin_ndx);
						
						perr = g.err;
						perr = perr(bin_ndx);
						
						ppm = g.pm;
						ppm = ppm(bin_ndx);
						
						if omega_free
							sse_addition = isotonic_rtheta_omega_free(g.name, omegas,...
								pvlos, pr, bin_theta, perr, ...
								ppm, vd2, truncate, use_nu,...
								doplot, gui, gui_graph0, fast);
						else
							sse_addition = isotonic_rtheta(pvlos, pr, bin_theta, perr, ...
								ppm, vd2, omegas(1), truncate, use_nu);
						end

						sse_binned = sse_binned + sse_addition;
					end
					sses(i) = sse_binned;
				end
			else
				for i=1:iterations
					m = length(r);

					bin_vlos = zeros(m,1);
					bin_r = zeros(m,1);
					bin_theta = zeros(m,1);
					bin_err = zeros(m,1);
					bin_pm = zeros(m,1);

					for k=1:num_bins
						st = ((k-1)*bin_size)+1;
						en = k*bin_size;

						if k==num_bins
							en = length(r);
						end

						bin_ndx = st:en;
						m = length(bin_ndx);

						perm_ndx = randperm(m);

						bin_theta_int = g.theta;
						bin_theta_int = bin_theta_int(bin_ndx);
						bin_theta_int = bin_theta_int(perm_ndx);
						
						pvlos = g.vlos;
						
						pr = g.r;
						
						perr = g.err;
						
						ppm = g.pm;
						
						
						bin_vlos(bin_ndx) = pvlos(bin_ndx);
						bin_r(bin_ndx) = pr(bin_ndx);
						bin_theta(bin_ndx) = bin_theta_int;
						bin_err(bin_ndx) = perr(bin_ndx);
						bin_pm(bin_ndx) = ppm(bin_ndx);
					end

					if omega_free
						sses(i) = isotonic_rtheta_omega_free(g.name, omegas,...
							bin_vlos, bin_r, bin_theta, bin_err, bin_pm, vd2, truncate, use_nu,...
							doplot, gui, gui_graph0, fast);
					else
						sses(i) = isotonic_rtheta(bin_vlos, bin_r, bin_theta, bin_err, ...
							bin_pm, vd2, omegas(1), truncate, use_nu);
					end

				end
			end
			
			
		
			
			pvalue = sum(sses<sse_obs)/length(sses);
			
			if fig
				hist(sses);
				hold on;
				plot([sse_obs],[0],'r.', 'MarkerSize', 50);
				hold off;
				flabel(['SSE'],['Freq'],[prefix g.name ' (' num2str(iterations) '/' num2str(num_bins) ') ' ...
					'p-value=' num2str(pvalue) ', sse obs=' ...
					num2str(sse_obs)]);
			end
			
			%--
		end
		
		function [pvalue, sse_obs, sses] = test_permute_rtheta_2(g, varargin)
			global im;
			load_args
			
			omega = arg('omega',0);
			truncate = arg('truncate', [0 10]);
			vd2 = arg('vel_disp',g.vel_disp)^2;
			use_nu = true;
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			
			
			r = g.r;
			
			num_bins = arg('num_bins', 10);
			
			g.sort_by('r');
			iterations = arg('iterations', 100);
			
			sse_obs = isotonic_rtheta(g.vlos, g.r, g.theta, g.err, g.pm, vd2, omega, truncate, use_nu);
			
			sses = zeros(iterations,1);
			
			
			r=g.r;
			r = r(g.pm>0);
			
			bin_size = floor(length(r)/num_bins);
					
					
			for i=1:iterations
				m = length(r);
				
				bin_vlos = zeros(m,1);
				bin_r = zeros(m,1);
				bin_theta = zeros(m,1);
				bin_err = zeros(m,1);
				bin_pm = zeros(m,1);
				
				for k=1:num_bins
					st = ((k-1)*bin_size)+1;
					en = k*bin_size;

					if k==num_bins
						en = length(r);
					end

					bin_ndx = st:en;
					m = length(bin_ndx);
					
					perm_ndx = randperm(m);
					
					bin_theta_int = g.theta(bin_ndx);
					bin_theta_int = bin_theta_int(perm_ndx);
					
					bin_vlos(bin_ndx) = g.vlos(bin_ndx);
					bin_r(bin_ndx) = g.r(bin_ndx);
					bin_theta(bin_ndx) = bin_theta_int;
					bin_err(bin_ndx) = g.err(bin_ndx);
					bin_pm(bin_ndx) = g.pm(bin_ndx);
				end
				
				sses(i) = isotonic_rtheta(bin_vlos, bin_r, bin_theta, bin_err, ...
						bin_pm, vd2, omega, truncate, use_nu);
			end
			
			pvalue = sum(sses<sse_obs)/length(sses);
			
			if fig
				hist(sses);
				hold on;
				plot([sse_obs],[0],'r.', 'MarkerSize', 50);
				hold off;
				flabel(['SSE'],['Freq'],['\lambda test i=' num2str(iterations) '. ' g.name ...
					', p-value=' num2str(pvalue) ', sse obs=' ...
					num2str(sse_obs)]);
			end
			
		end

		
		
		
		function [pvalue, sse_obs, sses] = test_permute_r(g,varargin)
			global im;
			load_args
			
			truncate = arg('truncate', [10 10]);
			
			vd2 = arg('vel_disp',g.vel_disp)^2;
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			
			r = g.r;
			
			iterations = arg('iterations', 100);
			
			
			
			
			
			
			sse_obs = isotonic_r(g.vlos, g.r, g.err, g.pm, vd2, truncate);
			
			sses = zeros(iterations,1);
			
			m = length(g.r);
			for i=1:iterations
				ndx = randperm(m);
				pvlos = g.vlos;
				pvlos = pvlos(ndx);
				sses(i) = isotonic_r(pvlos, g.r, g.err, g.pm, vd2, truncate);
			end
			
			pvalue = sum(sses<sse_obs)/length(sses);
			
			
			if fig
				hist(sses);
				hold on;
				plot([sse_obs],[0],'r.', 'MarkerSize', 50);
				hold off;
				flabel(['SSE'],['Freq'],['R test i=' num2str(iterations) '. ' g.name ...
					', p-value=' num2str(pvalue) ', sse obs=' ...
					num2str(sse_obs)]);
			end
			
		end
		
		
		
		
		%
		%
		% isotonic
		%
		%	
		function [sse, lam, vd2, r] = isotonic_r(g, varargin)
			global im;
			load_args
			
			truncate = arg('truncate', [0 0]);
			
			vd2 = arg('vel_disp',g.vel_disp)^2;
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			
			r = g.r;
			
			if gui
				doplot=true;
			else
				if fig
					doplot=true;
				else
					doplot=false;
				end
			end
			
			[sse, lam, nu, r] = isotonic_r(g.vlos, g.r, g.err, g.pm, vd2, truncate);
			
			
			if doplot
				if gui
					axes(gui_graph0);
				else
					%subplot(1,2,1)
				end
				plot(r,lam,'b');
				flabel(['r'],['\mu'],['R best. ' g.name ', nu=' num2str(nu) ', sse=' num2str(sse) ...
					]);
			end
			
			
		end
		
				
		% tries every possible omega to find the best
		function [sse, lam, vd2, omega, x1] = isotonic_x1(g, varargin)
			global im;
			load_args
			
			omegas = arg('omegas',0:0.1:2*pi);
			truncate = arg('truncate', [0 0]);
			
			vd2 = arg('vel_disp',g.vel_disp)^2;
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			
			x1 = g.x1;
			
			if gui
				doplot=true;
			else
				if fig
					doplot=true;
				else
					doplot=false;
				end
			end
			
			sse = NaN;
			all_sse = zeros(length(omegas),1);
			
			for k=1:length(omegas)
				om = omegas(k);
				[sse0, lam0, omega0, x10] = isotonic_x1(g.vlos, g.r, g.theta, g.err, g.pm, vd2, om,...
					truncate);
				
				all_sse(k) = sse0;
				
				%ignore flat estimators where flat is rise of lte 1e-12
				if abs(lam0(length(lam0)) - lam0(1)) > 1e-12 && (~isfinite(sse) || sse0 < sse)
					sse=sse0;
					lam=lam0;
					omega=omega0;
					x1=x10;
				end
				
				xim=im;
				if length(omegas)<200 && doplot
					if gui
						axes(gui_graph0);
					else
						subplot(1,2,1)
					end
					plot(x10,lam0,'k');
					flabel(['x_1'],['\lambda'],['X_1 (all). ' g.name ' \omega=' num2str(omega0) ', sse=' num2str(sse0) ...
						', \sigma_0=' num2str(sqrt(vd2))]);
					pause(.03);
				end
				im=xim;
			end
			
			
			
			if doplot
				if gui
					axes(gui_graph0);
				else
					subplot(1,2,1)
				end
				plot(x1,lam,'b');
				flabel(['x_1'],['\lambda'],['X_1 best. ' g.name ' \omega=' num2str(omega) ', sse=' num2str(sse) ...
					', \sigma_0=' num2str(sqrt(vd2))]);
			end
			
			
			if length(omegas)>1 && doplot
				if gui
					axes(gui_graph1);
				else
					subplot(1,2,2)
				end
				plot(omegas, all_sse,'b');
				flabel(['\omega'],['SSE'],['X_1. ' g.name, ' all sses']);
			end
			
		end
		
		% tries every possible omega to find the best
		function [sse, lam, nu, omega, theta, r] = isotonic_rtheta(g, varargin)
			global im;
			load_args
			
			omegas = arg('omegas',0:1:2*pi);
			truncate = arg('truncate', [0 0]);
			vd2 = arg('vel_disp',g.vel_disp)^2;
			use_nu = arg('use_nu', true);
			
			fast = arg('fast',false);
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			silent = arg('silent',false);
			
			
			if silent
				doplot=false;
			elseif gui
				doplot=true;
			else
				if fig
					doplot=true;
				else
					doplot=false;
				end
			end
			
			[sse, lam, nu, omega, theta, r, all_sse] = isotonic_rtheta_omega_free(g.name, omegas,...
				g.vlos, g.r, g.theta, g.err, g.pm, vd2, truncate, use_nu,...
				doplot, gui, gui_graph0,fast);
			
			
			
			if doplot
				if gui
					axes(gui_graph0);
				else
					subplot(1,2,1)
				end
				plot(r,lam,'r');
				flabel(['r'],['\lambda'],['Cosine best. ' g.name ' \omega= ' num2str(omega) ', sse=' num2str(sse)...
					', \nu= ' num2str(nu) ', \sigma_0=' num2str(sqrt(vd2))]);
			end
			
			
			if length(omegas)>1 && doplot
				if gui
					axes(gui_graph1);
				else
					subplot(1,2,2)
				end
				plot(omegas, all_sse,'r');
				flabel(['\omega'],['SSE'],['Cosine. ' g.name, ' all sses']);
			end
		end
		
		
		
		
		

		
		
		function [sse, yhat, nu, xi, gamma, omega, theta, r] = isotonic_xi_gamma(g, varargin)
			global im;
			imo = im + 1;
			im = im + 1;
			load_args
			
			
			omegas = arg('omegas',0:0.1:2*pi);
			truncate = arg('truncate', [0 0]);
			vd2 = arg('vel_disp',g.vel_disp)^2;
			
			gui = arg('gui',false);
			gui_graph0 = arg('gui_graph_0',NaN);
			gui_graph1 = arg('gui_graph_1',NaN);
			gui_graph2 = arg('gui_graph_2',NaN);
			gui_graph3 = arg('gui_graph_3',NaN);
			color = arg('color','k');
			
			g.constrict_by(struct('include_unbound',false));
			g.sort_by('r');
			
			vx = cell2mat(varargin);
			vx.omegas = [0];
			vx.use_nu = true;
			vx.silent = true;
			
			[sse_xi, xi, nu, omega, theta, r] = g.isotonic_rtheta(vx);
			
			vlos_subtracted = g.vlos - nu - xi.*cos(g.theta);
			
			g.plot_iso_est(r, 'r', xi, '\xi', '\xi', gui, gui_graph0, [3 1], imo, color);
			legend(strcat(['\nu=' num2str(round(nu*10)/10) '; \omega=' num2str(round(omega*10)/10)]));
			
			
			
			
			%now do gamma
			if gui
				doplot=true;
			else
				if fig
					doplot=true;
				else
					doplot=false;
				end
			end
			
			sse_gamma = NaN;
			all_sse_gamma = zeros(length(omegas),1);
			for k=1:length(omegas)
				om = omegas(k);
				[sse0, gamma0, nu0_not_used, omega0, theta0, r0, vd2] = ...
					isotonic_rtheta(vlos_subtracted, g.r, g.theta, g.err, g.pm, vd2, om, truncate, false);
				
				all_sse_gamma(k) = sse0;
				%ignore flat estimators where flat = rise of 1e-12
				if abs(gamma0(length(gamma0)) - gamma0(1)) > 1e-12 && (~isfinite(sse_gamma) || sse0 < sse_gamma)
					sse_gamma=sse0;
					gamma=gamma0;
					omega=omega0;
					theta=theta0;
					r=r0;
				end
				
				xim=im;
				if length(omegas)<100 && doplot
					g.plot_iso_est(r0, 'r', gamma0, '\gamma', '\gamma', gui, gui_graph0, [3 2], imo, color);
					legend(strcat(['\omega=' num2str(round(omega0*10)/10)]));
					
					pause(.03);
				end
				im=xim;
			end
			
			if im>0 
				im=im+1;
			end
			
			g.plot_iso_est(r, 'r', gamma, '\gamma', 'Final \gamma', gui, gui_graph0, [3 2], imo, color);
			legend(strcat(['\omega=' num2str(round(omega0*10)/10)]));
			
			
			g.plot_iso_est(omegas, '\omega', all_sse_gamma, 'SSE', '\gamma SSE', gui, gui_graph1, [3 3], imo, color);
			legend(strcat(['SSE=' num2str(round(sse_gamma*10)/10)]));
			
			
			yhat = nu + xi.*cos(g.theta) + gamma.*cos(g.theta);
			sse = sum((g.vlos - yhat).^2);
			
			
			if doplot || fig
				if gui
					axes(gui_graph0)
				end
				plot(r,xi,'b')
				hold on
				plot(r,gamma,'k')
				hold off
				title(strcat([g.name '. SSE=' num2str(sse) '; nu=' num2str(nu)]))
				legend('\xi: rotation', '\gamma: streaming')
				xlabel('r')
				ylabel('\xi/\gamma')
			end
			
			if ~gui
				if fig
					subplot(2,2,1);
					plot(r,yhat,'b.');
					flabel('r','\nu + \xi(r)cos(\theta) + \gamma(r)cos(\theta)',g.name);

					subplot(2,2,2);
					plot(cos(theta),yhat,'b.')
					flabel('cos(\theta)','\nu + \xi(r)cos(\theta) + \gamma(r)cos(\theta)',g.name)

					subplot(2,2,3);
					plot(r.*cos(theta),yhat,'b.')
					flabel('r*cos(\theta)','\nu + \xi(r)cos(\theta) + \gamma(r)cos(\theta)',g.name)

					subplot(2,2,4);
					plot(theta,yhat,'b.')
					flabel('\theta','\nu + \xi(r)cos(\theta) + \gamma(r)cos(\theta)',g.name)
				end
			end
			
		end
		

				
		function plot_iso_est(g,x,xlabel,f,funclabel,title,gui,gui_graph,subplots,imo,color)
			global im;
			
			if im < 1
				return
			end
			
			if gui
				doplot=true;
			else
				im=imo;
				if fig
					doplot=true;
				else
					doplot=false;
				end
			end
			
			if doplot
				if gui
					axes(gui_graph);
				elseif length(subplots) > 1
					subplot(1,subplots(1),subplots(2));
				end
				plot(x,f,color);
				flabel([xlabel],[funclabel],[g.name '. ' title]);
			end
		end
		
				
		%
		%
		% plotting
		%
		%	
		%[vlos r theta x1 x2 err pm tub sel_ndx];
		function plot_x1_x2(g,varargin)
			global im;
			load_args
			
			color_by = arg('color_by','vlos');
			color_scale = arg('color_scale', 'all');
			
			curcolmap = colormap;
			curmapsize = size(curcolmap,1);
				
			if strcmp(color_by,'vlos')
				if strcmp(color_scale,'all')
					minz = min(g.data(:,1)); maxz = max(g.data(:,1));
				else
					minz = min(g.vlos); maxz = max(g.vlos);
				end
				
				mapz = round(1 + (g.vlos - minz) ./ (maxz-minz) .* (curmapsize-1));
			else
				if strcmp(color_scale,'all')
					minz = max(g.data(:,7)); maxz = min(g.data(:,7));
				else
					minz = max(g.pm); maxz = min(g.pm);
				end
				
				mapz = round(1 + (g.pm - minz) ./ (maxz-minz) .* (curmapsize-1));
			end
			
			mapz(~isfinite(mapz))=1;
			
			figure(im)
			im=im+1;
			scatter3(g.x1, g.x2, g.vlos, 20, curcolmap(mapz,:), 'filled')
			xlabel('x_1')
			ylabel('y_1')
			zlabel('Velocity')
			title(strcat([g.name ' member stars']))
		end
		
		function plot_r_theta(g,varargin)
			global im;
			load_args
			
			color_by = arg('color_by','vlos');
			color_scale = arg('color_scale', 'all');
			
			curcolmap = colormap;
			curmapsize = size(curcolmap,1);
				
			if strcmp(color_by,'vlos')
				if strcmp(color_scale,'all')
					minz = min(g.data(:,1)); maxz = max(g.data(:,1));
				else
					minz = min(g.vlos); maxz = max(g.vlos);
				end
				
				mapz = round(1 + (g.vlos - minz) ./ (maxz-minz) .* (curmapsize-1));
			else
				if strcmp(color_scale,'all')
					minz = max(g.data(:,7)); maxz = min(g.data(:,7));
				else
					minz = max(g.pm); maxz = min(g.pm);
				end
				
				mapz = round(1 + (g.pm - minz) ./ (maxz-minz) .* (curmapsize-1));
			end
			
			mapz(~isfinite(mapz))=1;
			
			figure(im)
			im=im+1;
			scatter3(g.r, g.theta, g.vlos, 20, curcolmap(mapz,:), 'filled')
			xlabel('r')
			ylabel('\theta')
			zlabel('Velocity')
			title(strcat([g.name ' member stars']))
		end
		
		%plots rotation curve using binning, not sure about +/- stuff
		function plot_rotation_curve(g,varargin)
			load_args
			
			nbins = arg('nbins',[12]);
			global im;
			
			colors = ['b' 'r' 'c' 'y' 'k'];
			if im>0
				figure(im);
				im=im+1;
				xlabel('r bin number (evenly distributed by count)');
				ylabel('v');
				title(strcat(['Rotation curve for ' g.name ' with ' num2str(nbins) ' bins.']));
				
				hold on
				
				for j=1:length(nbins)
					ndx=g.pm>.9;
					mu = [];
					s = [];
					r=g.r(ndx);
					v=g.vlos(ndx);
					sgm=g.err(ndx);
					for i=1:nbins(j)
						mx = floor((i*length(r))/nbins(j));
						mn = floor(((i-1)*length(r))/nbins(j))+1;
						
						v_bin = v(mn:mx);
						
						mu(i) = mean(v_bin);
						s(i) = sqrt(mean((v_bin-mu(i)).^2));%vel disp
						%s(i) = (mean((sgm(mn:mx))));%average measure error
					end
					errorbar(linspace(min(r),max(r),nbins(j)), mu, s, strcat([colors(j) '.']))
				end
				
				plot([min(r) max(r)], [mean(v) mean(v)], 'k-')
				
				hold off
			end
		end
		
		%plot data for simulation analysis
		function plot_sim_data(g,varargin)
			load_args
			show_all_rotations = arg('show_all_rotations',false);
			gui = arg('gui',false);
			
			%gui
			%global im;
			if gui || fig
				subplot(2,3,1);plot(g.x1,g.x2,'.');
				
				if ~gui
					title(strcat(['n=' num2str(length(g.x1))]));
				end
				
				xlabel('X_1'); ylabel('X_2');hold off;

				subplot(2,3,2);plot3(g.x1,g.x2,g.vlos,'.');
				xlabel('X_1'); ylabel('X_2'); zlabel('Vlos'); hold off;

				subplot(2,3,3);plot(g.x1,g.vlos,'.'); xlabel('X_1'); ylabel('Vlos'); hold off;
				subplot(2,3,4);plot(g.x2,g.vlos,'.');  xlabel('X_2'); ylabel('Vlos'); hold off;
				subplot(2,3,5);plot(g.r,g.vlos,'.');  xlabel('R'); ylabel('Vlos'); hold off;

				subplot(2,3,6);ksdensity(g.vlos);

				if show_all_rotations
					for om = 0:0.05:2*pi,
						x1_n = g.r.*cos(g.theta - om);
						x2_n = g.r.*sin(g.theta - om);
						figure(im); axis([-25 25 -40 100]); plot(x1_n,g.vlos,'.');title(num2str(om));
						xlabel('ra*cos(th - om)'); ylabel('Vlos'); hold off;
						pause(0.5);
					end
					im=im+1;
				end
			end
		end
		
		
		
		%
		%
		% utilities
		%
		%
		
		%pass 0 to remove core cutoff
		function limit_to_core(g,r_core)
			g.constrict_by(struct('r_core',r_core));
		end
		
		
		%rotate the galaxy as it was originally positioned to the passed number of radians
		function rotate_to(g, omega_radians)
			g.data(:,3) = g.data(:,3) - g.rotated_by;%reset to zero
			g.data(:,3) = g.data(:,3) + omega_radians;%rotate
			
			g.rotated_by = omega_radians;
			
			g.data(:,4) = g.data(:,2).*cos(g.data(:,3));
			g.data(:,5) = g.data(:,2).*sin(g.data(:,3));
		end
		
		
		function constrict_by(g,varargin)
			
			load_args
			
			r_core = arg('r_core', g.r_core);
			include_unbound = arg('include_unbound', g.include_unbound);
			
			g.r_core=r_core;
			g.include_unbound=include_unbound;
			
			
			if r_core ~= 0
				if include_unbound
					idx = g.data(:,2) < r_core;
				else
					idx = g.data(:,2) < r_core & g.data(:,8) == 0;
				end
			else
				if include_unbound
					idx = isfinite(g.data(:,2));
				else
					idx = isfinite(g.data(:,2)) & g.data(:,8) == 0;
				end
			end
			
			g.data(:,9) = idx;
		end
		
		
		
		function sort_by(g,sb)
			g.data = Galaxy.sort_by_vals(sb,g.data);
		end

		
		
		
	end
	
end