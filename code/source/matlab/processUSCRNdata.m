function processUSCRNdata(years, months, orbits)
% function processUSCRNdata
%   
%   Create in-situ dataset (.nc) with USCRN soil moisture data and MODIS 
%   NDVI data. The raw USCRN data is filetered by quality control flags.
%   It downloads USCRN data if not exist. Interpolated NDVI data is used
%   and converted to VWC. A dynamic (time-varying) vegetation structure 
%   for the discrete scatterer model is modeled based on VWC, vegetation
%   parameters, and scaling parameters. If the dataset exists, it only
%   plots the data.
%
%   See also main.
%
%   Copyright (C) 2023 Seho Kim
%     
%   This file is part of SOCRATES-Retrieval.
%     
%   SOCRATES-Retrieval is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%     
%   SOCRATES-Retrieval is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%     
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%   Version 1.0

startMonth = months(1);
endMonth = months(2);

FILLVALUE = -9999;

dataColor = {'r.', 'g.', 'b.', 'c.', 'm.'};

% Directories
dir_in = './inputs/';
dir_geometry = './data/geometry/';
dir_USCRN = './data/USCRN/';
dir_processed = './data/processed/';
dir_NDVI = './data/MODIS-NDVI/';

% Read the list of selected stations
temp = readlines(strcat(dir_in, 'station_in.txt'),"EmptyLineRule","skip");
nStation = str2double(temp(1));
station_in = temp(2:end);

% Read data of the selected stations from the station table
station_info = readcell(strcat(dir_in, 'station_list.xlsx'), 'NumHeaderLines', 1);
station_name = station_info(:,1);
[~, idx_station_in] = ismember(station_in, station_name);
station_vegType = lower(station_info(idx_station_in,2));
station_lat = cell2mat(station_info(idx_station_in,3));
station_lon = cell2mat(station_info(idx_station_in,4));
station_soilDepth = str2num(cell2mat(station_info(idx_station_in,6)));
station_soilType = 0.01*cell2mat(station_info(idx_station_in, 7:9));

nYear = length(years);
nOrb = length(orbits);
nLayer = size(station_soilDepth, 2);
nChar = 20; % Dimension of varaible 'vegType' of NetCDF

% Stem factor
stem_factor = [15.96, 19.15, 7.98, 12.77, 12.77, 3.0, 1.5, 4.0, 3.0, 1.5, 4.0, 3.5, 6.49, 3.25, 0, 0];

% Nominal vegetation parameters for grass modeled as an elliptical disk
nVar.grass = 6; % Number of variables for vegetation model
dim1.grass = 0.2; % Disk semi-major axis (m)
dim2.grass = 0.01; % Disk semi-minor axis (m)
dim3.grass = 0.0002; % Disk thickness (m)
m_v.grass = 0.5; % Volumetric Moisture (g/cm3)
d.grass = 0.425; % Layer thickness (m)
density.grass = 577; % Number density
vf.grass = density.grass * pi*dim1.grass*dim2.grass*dim3.grass; % Volume fraction
W.grass = vf.grass * m_v.grass*1e3 * d.grass; % Vegetation water content

nVeg = 1;
nVegVar = nVar.grass;

%% Compute VWC based on dynamic vegetation model
for iStation = 1 : nStation
    % Read pre-processed MODIS NDVI data for each station
    filename_NDVI = strcat('MODIS_NDVI_', station_in(iStation), '_interp.txt');
    fullpath_NDVI = strcat(dir_NDVI, filename_NDVI);
    temp = readmatrix(fullpath_NDVI);
    yyyy{iStation} = temp(:,1);
    doy{iStation} = temp(:,2);
    ANC.NDVIdate{iStation} = datetime(num2str(yyyy{iStation}), 'InputFormat', 'yyyy') + caldays(doy{iStation}-1);
    ANC.NDVI{iStation} = temp(:,3);
    % Masking fillvalue
    ANC.NDVI{iStation}(ANC.NDVI{iStation}==FILLVALUE) = NaN;
    % Compute VWC from NDVI
    ANC.VWC.grass{iStation} = (1.9134*ANC.NDVI{iStation}.^2 - 0.3215*ANC.NDVI{iStation}) ...
        + stem_factor(10) * (ANC.NDVI{iStation} - 0.1) / (1 - 0.1);

    % Compute dynamic vegetation parameters using the scaling method
    scalingParams = [0.1, 0.3, 0.1, 0.1, 0.1, 0.3];
    alpha = ANC.VWC.grass{iStation} / W.grass;
    ANC.density.grass{iStation} = alpha .^scalingParams(1) * density.grass;
    ANC.dim1.grass{iStation} = alpha .^scalingParams(2)* dim1.grass;
    ANC.dim2.grass{iStation} = alpha .^scalingParams(3) * dim2.grass;
    ANC.dim3.grass{iStation} = alpha .^scalingParams(4) * dim3.grass;
    ANC.m_v.grass{iStation} = alpha .^scalingParams(5) * m_v.grass;
    ANC.d.grass{iStation} = alpha .^scalingParams(6) * d.grass;
    
    % Compute VWC from modeled time-varying vegetation
    ANC.VWC.grass_new{iStation} = ANC.density.grass{iStation}*pi.*ANC.dim1.grass{iStation}.*ANC.dim2.grass{iStation}.*ANC.dim3.grass{iStation}.*ANC.m_v.grass{iStation}*1e3.*ANC.d.grass{iStation};

    % Assign variables
    VWCref{iStation} = ANC.VWC.grass_new{iStation};
    vegVar{iStation} = [ANC.density.grass{iStation}, ANC.dim1.grass{iStation}, ANC.dim2.grass{iStation}, ANC.dim3.grass{iStation}, ANC.m_v.grass{iStation}, ANC.d.grass{iStation}];
end

%% WRITE NETCDF
for iYear = 1 : nYear
    for iOrb = 1 : nOrb
        prev_year_str = num2str(years(iYear) - 1);
        curr_year_str = num2str(years(iYear));
        dir_in1 = strcat(dir_USCRN, prev_year_str, '/');
        dir_in2 = strcat(dir_USCRN, curr_year_str, '/');
        dir_out = strcat(dir_processed, curr_year_str, '/');
        % Create directories if not exist
        if(~exist(dir_out, 'dir'))
            mkdir(dir_out);
        end
        if(~exist(dir_in1, 'dir'))
            mkdir(dir_in1);
        end
        if(~exist(dir_in2, 'dir'))
            mkdir(dir_in2);
        end
        % Incidence angle data
        filename_inc = strcat('inc_', station_in, '_', orbits{iOrb}, '.mat');
        for iStation = 1 : nStation
            filename_in1 = strcat('CRNH0203-', prev_year_str, '-', station_in(iStation), '.txt');
            filename_in2 = strcat('CRNH0203-', curr_year_str, '-', station_in(iStation), '.txt');
            filename_out = sprintf('USCRN_hourly_%sM%dM%d_%s_%s', curr_year_str, startMonth, endMonth, orbits{iOrb}, station_in(iStation));
            
            fullpath_out = strcat(dir_out, filename_out, '.nc');
            
            fprintf('Creating %s ... ', fullpath_out);
            % Create USCRN NetCDF file if not exist
            if(~exist(fullpath_out, 'file'))
                % Read USCRN data
                % Data from the last hour of the previous year would be needed if starting month is January
                if ~exist(strcat(dir_in1, filename_in1), 'file') && startMonth == 1
                    % Download USCRN hourly data if not exist
                    file_temp = fopen(strcat(dir_in1, filename_in1), 'w+');
                    options = weboptions('Timeout', 150);
                    data = webread(strcat('https://www1.ncdc.noaa.gov/pub/data/uscrn/products/hourly02/', prev_year_str, '/', filename_in1), options);
                    fprintf(file_temp, '%s', data);
                    fclose(file_temp);                
                end
                if ~exist(strcat(dir_in2, filename_in2), 'file')
                    % Read USCRN data on website
                    file_temp = fopen(strcat(dir_in2, filename_in2), 'w+');
                    options = weboptions('Timeout', 150);
                    data = webread(strcat('https://www1.ncdc.noaa.gov/pub/data/uscrn/products/hourly02/', curr_year_str, '/', filename_in2), options);
                    fprintf(file_temp, '%s', data);
                    fclose(file_temp);                
                end
                if startMonth == 1
                    temp1 = readmatrix(strcat(dir_in1, filename_in1));
                    temp2 = readmatrix(strcat(dir_in2, filename_in2));
                    temp = [temp1; temp2];
                else
                    temp = readmatrix(strcat(dir_in2, filename_in2));
                end
                
                utcDate = temp(:,2); % YYYYMMDD
                utcTime = temp(:,3); % HHmm
                utcDatetime{iStation} = datetime(utcDate,'ConvertFrom','yyyymmdd') + timeofday(datetime(num2str(utcTime,'%04d'),'Format','HHmm'));
                
                % Index for data within months of interest
                idx_valid = utcDatetime{iStation} >= datetime(years(iYear), startMonth, 1, 0, 0, 0) & utcDatetime{iStation} <= datetime(years(iYear), endMonth, eomday(years(iYear), endMonth), 24, 00, 00);
    
                % Headers of USCNR data - refer to https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/readme.txt
                utcDate = utcDate(idx_valid);
                utcTime = utcTime(idx_valid);
                utcDatetime{iStation} = utcDatetime{iStation}(idx_valid);
                jd = juliandate(utcDatetime{iStation});
    
                lst_date = temp(idx_valid,4); % YYYYMMDD
                lst_time = temp(idx_valid,5); % HHmm
                T_air_calc = temp(idx_valid,9); % Celsius
                T_air_hr_avg = temp(idx_valid,10);
                T_air_max = temp(idx_valid,11);
                T_air_min = temp(idx_valid,12);
                P_calc = temp(idx_valid,13); % mm
                T_sur = temp(idx_valid,21);
                T_sur_max = temp(idx_valid,23);
                T_sur_min = temp(idx_valid,25);
                VSM = temp(idx_valid, 29:33); % m^3/m^3
                T_soil = temp(idx_valid, 34:38); % Celsius
    
                nSample = size(utcDatetime{iStation}, 1);
                
                % Load incidence angle data
                load(strcat(dir_geometry,filename_inc(iStation)));
    
                % Read vegetation data
                idxYear = year(ANC.NDVIdate{iStation}) == years(iYear);
                VWCdaily = VWCref{iStation}(idxYear);
                vegVarDaily = vegVar{iStation}(idxYear,:);
                VWChourly{iStation} = [repelem(VWCdaily, 24); VWCdaily(end)];
                vegVarHourly{iStation} = [repelem(vegVarDaily, 24, 1); vegVarDaily(end,:)];
    
                % Set fill value
                VSM(VSM == -99) = FILLVALUE;
                T_air_min(T_air_min == -9999) = FILLVALUE;
                T_soil(T_soil == -9999) = FILLVALUE;
                P_calc(P_calc == -9999) = FILLVALUE;
                VWChourly{iStation}(isnan(VWChourly{iStation})) = FILLVALUE;
                vegVarHourly{iStation}(isnan(vegVarHourly{iStation})) = FILLVALUE;
    
                nStep = 1;
                idx_sample = 1:nStep:nSample;
                nSample = length(idx_sample);
                inc_selected = inc_selected(idx_sample,:);
                utcDatetime{iStation} = utcDatetime{iStation}(idx_sample);
                jd = jd(idx_sample);
                utcDate = utcDate(idx_sample);
                utcTime = utcTime(idx_sample);
                lst_date = lst_date(idx_sample);
                lst_time = lst_time(idx_sample);
                VSM = VSM(idx_sample,:);
                T_soil = T_soil(idx_sample,:);
                T_air_min = T_air_min(idx_sample);
                P_calc = P_calc(idx_sample);
                VWChourly{iStation} = VWChourly{iStation}(idx_sample,:);
                vegVarHourly{iStation} = vegVarHourly{iStation}(idx_sample,:);
                
                % Write NetCDF
                nccreate(fullpath_out, 'lat', 'FillValue', FILLVALUE, 'Datatype', 'single', 'Format', 'netcdf4');    
                nccreate(fullpath_out, 'lon', 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'soilDepth', 'Dimensions', {'layers', nLayer}, 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'soilType', 'Dimensions', {'soilSeparates', 3}, 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'juliandate', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'utcDate', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'utcTime', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'lstDate', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'lstTime', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'soilMoist', 'Dimensions', {'samples', nSample, 'layers', nLayer}, 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'soilTemp', 'Dimensions', {'samples', nSample, 'layers', nLayer}, 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'airTempMin', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'prec', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'single');
                nccreate(fullpath_out, 'vegType', 'Dimensions', {'chars', nChar}, 'Datatype', 'char');
                nccreate(fullpath_out, 'incMUOS', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'incORBCOMM', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'incGPS', 'Dimensions', {'samples', nSample}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'VWCref', 'Dimensions', {'samples', nSample, 'vegs', nVeg}, 'FillValue', FILLVALUE, 'Datatype', 'double');
                nccreate(fullpath_out, 'vegVar', 'Dimensions', {'samples', nSample, 'vegVars', nVegVar}, 'FillValue', FILLVALUE, 'Datatype', 'double');
     
                if(nSample > 0)
                    ncwrite(fullpath_out, 'lat', station_lat(iStation));
                    ncwrite(fullpath_out, 'lon', station_lon(iStation));
                    ncwrite(fullpath_out, 'soilDepth', station_soilDepth(iStation,:));
                    ncwrite(fullpath_out, 'soilType', station_soilType(iStation,:));
                    ncwrite(fullpath_out, 'juliandate', jd); 
                    ncwrite(fullpath_out, 'utcDate', utcDate); 
                    ncwrite(fullpath_out, 'utcTime', utcTime); 
                    ncwrite(fullpath_out, 'lstDate', lst_date); 
                    ncwrite(fullpath_out, 'lstTime', lst_time);
                    ncwrite(fullpath_out, 'soilMoist', VSM);
                    ncwrite(fullpath_out, 'soilTemp', T_soil);
                    ncwrite(fullpath_out, 'airTempMin', T_air_min);
                    ncwrite(fullpath_out, 'prec', P_calc);
                    ncwrite(fullpath_out, 'vegType', cell2mat(station_vegType(iStation,:)));
                    ncwrite(fullpath_out, 'incORBCOMM', inc_selected(:,1));
                    ncwrite(fullpath_out, 'incMUOS', inc_selected(:,2));
                    ncwrite(fullpath_out, 'incGPS', inc_selected(:,3));
                    ncwrite(fullpath_out, 'VWCref', VWChourly{iStation});
                    ncwrite(fullpath_out, 'vegVar', vegVarHourly{iStation});
                end
                fprintf('complete.\n');
            else
                % Read USCRN NetCDF
                utcDate = ncread(fullpath_out, 'utcDate'); 
                utcTime = ncread(fullpath_out, 'utcTime'); 
                utcDatetime{iStation} = datetime(utcDate,'ConvertFrom','yyyymmdd') + timeofday(datetime(num2str(utcTime,'%04d'),'Format','HHmm'));
                VSM = ncread(fullpath_out, 'soilMoist');
                T_soil = ncread(fullpath_out, 'soilTemp');
                T_air_min = ncread(fullpath_out, 'airTempMin');
                P_calc = ncread(fullpath_out, 'prec');
                clear inc_selected
                inc_selected(:,1) = ncread(fullpath_out, 'incORBCOMM');
                inc_selected(:,2) = ncread(fullpath_out, 'incMUOS');
                inc_selected(:,3) = ncread(fullpath_out, 'incGPS');
                VWChourly{iStation} = ncread(fullpath_out, 'VWCref');
                vegVarHourly{iStation} = ncread(fullpath_out, 'vegVar');
                fprintf('exists.\n');
            end
            
            % Plot VSM, Temperature, Precipitation, and VWC
            VSM(VSM == FILLVALUE) = NaN;
            T_air_min(T_air_min == FILLVALUE) = NaN;
            T_soil(T_soil == FILLVALUE) = NaN;
            P_calc(P_calc == FILLVALUE) = NaN;
            VWChourly{iStation}(VWChourly{iStation} == FILLVALUE) = NaN;
            vegVarHourly{iStation}(vegVarHourly{iStation} == FILLVALUE) = NaN;
    
            fullpath_plot = strcat(dir_out, filename_out, '.png');
            fprintf('Ploting %s ... ', fullpath_plot);
            if(~exist(fullpath_plot,'file'))
                hf = figure('visible', 'off', 'Position', [0 0 1600 700]);
                subplot(4,1,1); hold on; grid on;
                yyaxis left
                for iLayer = 1 : 5
                    plot(utcDatetime{iStation}, VSM(:,iLayer), dataColor{iLayer}, 'LineStyle', 'None');
                end
                ylabel('VSM [m^3m^{-3}]');
                yyaxis right
                scatter(utcDatetime{iStation}, P_calc, 'k.');
                ylabel('Precipitation [mm]');
                xlabel('Date');
                legend('5cm', '10cm', '20cm', '50cm', '100cm', 'Preci', 'Location', 'northeastoutside');
                subplot(4,1,2); hold on; grid on;
                for iLayer = 1 : 5
                    plot(utcDatetime{iStation}, T_soil(:,iLayer), dataColor{iLayer}, 'LineStyle', 'None');
                end
                scatter(utcDatetime{iStation}, T_air_min, 'k.');
                ylabel('Temp. [C]');
                xlabel('Date');
                legend('5cm', '10cm', '20cm', '50cm', '100cm', 'min air', 'Location', 'northeastoutside');
                subplot(4,1,3); hold on; grid on;
                plot(utcDatetime{iStation}, VWChourly{iStation}, 'LineWidth', 1);
                legend('VWC', 'Location', 'northeastoutside')
                ylabel('VWC [kg/m^2]');
                xlabel('Date');
                subplot(4,1,4); hold on; grid on;
                for iTx = 1 : 3
                    plot(utcDatetime{iStation}, inc_selected(:,iTx), '.', 'LineStyle', 'none');
                end
                legend('Orbcomm', 'MUOS', 'GPS', 'Location', 'northeastoutside')
                ylabel('Incidence Angle [deg]');
                xlabel('Date');
                print(hf, fullpath_plot, '-dpng');
                fprintf('complete.\n');
            else
                fprintf('exists.\n');
            end
        end
    end
end

end