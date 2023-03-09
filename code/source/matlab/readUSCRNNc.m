function data = readUSCRNNc(path, station_name, year, months, orbit)
% function readUSCRNNc
%   
%   Reads USCRN dateset created by processUSCRNdata.
%
%   See also estimate_b_parameter, plotSynObs, calcStatAndPlot.
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

    data.station_name_space = strrep(station_name,'_', ' ');
    for iStation = 1:length(station_name)
        fullpath_USCRN = sprintf('%s%d/USCRN_hourly_%dM%dM%d_%s_%s.nc', ...
            path, year, year, months(1), months(2), orbit, station_name{iStation});
        if(exist(fullpath_USCRN, 'file')>0)
            % Read USCRN NetCDF
            info_USCRN{iStation} = ncinfo(fullpath_USCRN);
            station_lat{iStation} = ncread(fullpath_USCRN, 'lat');
            station_lon{iStation} = ncread(fullpath_USCRN, 'lon');
            station_soilDepth{iStation} = ncread(fullpath_USCRN, 'soilDepth');
            station_soilFractions{iStation} = ncread(fullpath_USCRN, 'soilType');
            utcDate{iStation} = ncread(fullpath_USCRN, 'utcDate'); 
            utcTime{iStation} = ncread(fullpath_USCRN, 'utcTime'); 
            soilMoist{iStation} = ncread(fullpath_USCRN, 'soilMoist');
            soilTemp{iStation} = ncread(fullpath_USCRN, 'soilTemp');
            airTempMin{iStation} = ncread(fullpath_USCRN, 'airTempMin');
            prec{iStation} = ncread(fullpath_USCRN, 'prec');
            vegType{iStation} = ncread(fullpath_USCRN, 'vegType');
            incMUOS{iStation} = ncread(fullpath_USCRN, 'incMUOS');
            incORBCOMM{iStation} = ncread(fullpath_USCRN, 'incORBCOMM');
            incGPS{iStation} = ncread(fullpath_USCRN, 'incGPS');

            nLayer{iStation} = length(station_soilDepth{iStation});
            nSample{iStation} = size(soilMoist{iStation}, 1);
            utcDatetime{iStation} = datetime(utcDate{iStation},'ConvertFrom','yyyymmdd') + timeofday(datetime(num2str(utcTime{iStation},'%04d'),'Format','HHmm'));

            vegVar{iStation} = ncread(fullpath_USCRN, 'vegVar');
            VWCref{iStation} = ncread(fullpath_USCRN, 'VWCref');
        end
    end
    data.info_USCRN = info_USCRN;
    data.station_name = station_name;
    data.station_lat = station_lat;
    data.station_lon = station_lon;
    data.station_soilDepth = station_soilDepth;
    data.station_soilFractions = station_soilFractions;
    data.utcDate = utcDate;
    data.utcTime = utcTime;
    data.utcDatetime = utcDatetime;
    data.soilMoist = soilMoist;
    data.soilTemp = soilTemp;
    data.airTempMin = airTempMin;
    data.prec = prec;
    data.vegType = vegType;
    data.incMUOS = incMUOS;
    data.incORBCOMM = incORBCOMM;
    data.incGPS = incGPS;
    data.nLayer = nLayer;
    data.nSample = nSample;
    data.vegVar = vegVar;
    data.VWCref = VWCref;
end