classdef Halo < handle
	
	properties
		halos;
		header;
		noise;
		grid_size;
	end
	
	methods(Static)
		
		function [gdata,noise,grid_size] = generate_model_data(varargin)
			data = Halo.read(cell2mat(varargin));
			load_args
			grid_size = arg('grid_size',0.1);
			noise = arg('noise');
			plots = arg('plots',false);
			gridt = arg('gridt',NaN);
			gridm = arg('gridm',NaN);
			
			if ~isfinite(gridt)
				error('you need to provide both gridt and gridm')
			end
			if ~isfinite(gridm)
				error('you need to provide both gridt and gridm')
			end
			
			minx = min(data(:,1));
			if minx<-3
				data = data(data(:,1)>=-3,:);
				minx=-3;
			end

			X = size(data,1);
			data(:,1:2) = data(:,1:2) + normrnd(0,noise,X,2);



			%n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). The last bin
			%counts any values of x that match edges(end). Values outside the values in
			%edges are not counted.
			[nt,tbin] = histc(data(:,4),gridt);
			[nm,mbin] = histc(data(:,5),gridm);


			data = [data(:,1:6) tbin mbin];
			data = data(data(:,7)>0,:);
			data = data(data(:,8)>0,:);
			X = size(data,1);

			tlook = [nt gridt'];
			mlook = [nm gridm'];
			freq = data(:,3).*data(:,5);
			%feh,afe,nfreq,time bin ndx, mass bin ndx, min time, max time, min mass,
			%max mass (min/max mass/time are both the min at the moment
			data = [data(:,1:2) freq data(:,7) data(:,8)...
				tlook(data(:,7),2) tlook(data(:,7),2) mlook(data(:,8),2) mlook(data(:,8),2)];





			tic
			gdata = [];
			for ii=1:length(gridt)
				for jj=1:length(gridm)
					ndx = data(:,4)==ii & data(:,5)==jj;
					sd = data(ndx,:);
					minx = min(sd(:,1));
					minx = round(minx*10)/10;
					if minx<-3
						sd = sd(sd(:,1)>=-3,:);
						minx=-3;
					end
					maxx = max(sd(:,1));
					maxx = round(maxx*10)/10;
					miny = min(sd(:,2));
					miny = round(miny*10)/10;
					maxy = max(sd(:,2));
					maxy = round(maxy*10)/10;
					gridx = minx:grid_size:maxx;
					gridy = miny:grid_size:maxy;

					X = size(sd,1);

					%n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). The last bin
					%counts any values of x that match edges(end). Values outside the values in
					%edges are not counted.
					[nx,xbin] = histc(sd(:,1),gridx);
					[ny,ybin] = histc(sd(:,2),gridy);


					sd = [sd(:,1:9) xbin ybin];
					sd = sd(sd(:,10)>0,:);
					sd = sd(sd(:,11)>0,:);

					xlook = [nx gridx'];
					ylook = [ny gridy'];
					sd = [sd(:,1:11) xlook(sd(:,10),2)-grid_size/2 ylook(sd(:,11),2)-grid_size/2];



					for i=1:length(gridx)
						for j=1:length(gridy)
							gx = gridx(i)-grid_size/2;
							gy = gridy(j)-grid_size/2;
							mrows=sd(sd(:,12)==gx & sd(:,13)==gy,:);
							if size(mrows,1)>0
								cr = mrows(1,:);
								nfreq = sum(mrows(:,3));
								gdata(size(gdata,1)+1,:) = [gx gy nfreq cr(4:9)];
							else
								gdata(size(gdata,1)+1,:) = [gx gy 0 0 0 0 0 0 0];
							end
						end
					end

				end
			end
			
			
		end
		function [gdata,noise,grid_size] = generate_observed_data(varargin)
			data = Halo.read(cell2mat(varargin));
			load_args
			grid_size = arg('grid_size',0.1);
			noise = arg('noise');
			plots = arg('plots',false);
			N = arg('n',1);
			
			minx = min(data(:,1));
			minx = round(minx*10)/10;
			if minx<-3
				data = data(data(:,1)>=-3,:);
				minx=-3;
			end
			maxx = max(data(:,1));
			maxx = round(maxx*10)/10;
			miny = min(data(:,2));
			miny = round(miny*10)/10;
			maxy = max(data(:,2));
			maxy = round(maxy*10)/10;
			gridx = minx:grid_size:maxx;
			gridy = miny:grid_size:maxy;

			X = size(data,1);
			data(:,1:2) = data(:,1:2) + normrnd(0,noise,X,2);
			
			%todo:cdf or some other way to select other than first N
			gdata = data(1:N,1:2);
		end
		function data = read(varargin)
			load_args
			halos = arg('halos',[]);
			data_point_limit = arg('data_point_limit',NaN);
			
			data = zeros(1,6);
			for i=1:length(halos)
				h=halos(i);
				for k=0:50
					fn = strcat(['data/large/halodatalist/datalist_halo'...
						num2str(h) 'part.' num2str(k) '.dat']);
					if exist(fn,'file')
						show(strcat(['loading ' fn]));
						[hdr,d]=hdrload(fn);
						if ~isfinite(data)
							data = d;
						else
							data = [data;d];
						end
					else
						show(strcat(['can''t find' fn 'stopping']));
						break
					end
					if isfinite(data_point_limit) && size(data,1)>data_point_limit
						show('truncating to',data_point_limit,'data points');
						data = data(1:data_point_limit,:);
						break;
					end
				end
			end
		end
	end
	
	methods
		
	end
end