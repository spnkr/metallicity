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
			
			if plots;fig;plot(data(:,1),data(:,2),'b.');end
			
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
			
			if plots;fig;plot(data(:,1),data(:,2),'b.');end
			
			
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