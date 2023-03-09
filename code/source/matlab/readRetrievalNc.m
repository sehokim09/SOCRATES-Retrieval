function data = readRetrievalNc(path, station_name, year, months, orbit, period, window, stddev, polRx, flagVeg)
% function readRetrievalNc
%   
%   Reads retrieved data created by SOCRATES-Retrieval.
%
%   See also calcStatAndPlot.
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

    nStation = length(station_name);
    nPol = length(polRx);
 
    for iStation = 1 : nStation
        filename_ret = sprintf('USCRN_%dM%dM%d_p%dw%dstd%d_%s_%s_ret.nc', year, months(1), months(2), period, window, stddev*100, orbit, station_name{iStation});
        fullpath_ret = strcat(path, filename_ret);
        % Read Retrieval NetCDF
        if(exist(fullpath_ret, 'file')>0)
            info_ret{iStation} = ncinfo(fullpath_ret);
            flagQcRet{iStation} = ncread(fullpath_ret, 'flagQcRet');
            ind_productTime{iStation} = ncread(fullpath_ret, 'productTimeIndex');
            ind_QC_GOOD{iStation} = flagQcRet{iStation} == 0;
            nRet(iStation) = sum(ind_QC_GOOD{iStation});
        
            ind_productTime{iStation} = ind_productTime{iStation} + 1; % Indexes start from 1 in MATLAB
            ind_obsTime{iStation} = ncread(fullpath_ret, 'obsTimeIndex')';
            ind_obsTime{iStation} = ind_obsTime{iStation}(ind_QC_GOOD{iStation}, :) + 1; % Indexes start from 1 in MATLAB
            soilMoistRet{iStation} = ncread(fullpath_ret, 'soilMoistRet')';
            VWCRet{iStation} = ncread(fullpath_ret, 'VWCRet');

            if flagVeg
                for iPol = 1 : nPol
                    switch(polRx(iPol))
                        case 'R'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_R_V_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_R_V_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_R_V_err')';
                        case 'L'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_R_V_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_R_V_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_R_V_err')';
                        case 'X'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_X_V_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_X_V_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_X_V_err')';
                        case 'Y'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_X_V_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_X_V_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_X_V_err')';
                    end                    
                end
            else
                for iPol = 1 : nPol
                    switch(polRx(iPol))
                        case 'R'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_R_B_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_R_B_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_R_B_err')';
                        case 'L'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_R_B_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_R_B_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_R_B_err')';
                        case 'X'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_X_B_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_X_B_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_CO_X_B_err')';
                        case 'Y'
                            refl_ref{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_X_B_ref')';
                            refl_est{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_X_B_est')';
                            refl_err{iStation, iPol} = ncread(fullpath_ret, 'reflect_X_X_B_err')';
                    end                   
                end
            end
        
            penDepthHRet{iStation} = ncread(fullpath_ret, 'penDepthHRet')';
            penDepthVRet{iStation} = ncread(fullpath_ret, 'penDepthVRet')';
            reflErr{iStation} = ncread(fullpath_ret, 'reflErr')';
            
            typePOMERet{iStation} = ncread(fullpath_ret, 'typePOME');
            saTime{iStation} = ncread(fullpath_ret, 'sa_time');
            saIter{iStation} = ncread(fullpath_ret, 'sa_iter');
            saJump{iStation} = ncread(fullpath_ret, 'sa_jump');
            saCostMin{iStation} = ncread(fullpath_ret, 'sa_costMin');
            runtime{iStation} = ncread(fullpath_ret, 'runtime');      
                    
            saTimeMean(iStation) = mean(saTime{iStation}, 'omitnan');
        end
    end

    data.info_ret = info_ret;
    data.nRet = nRet;
    data.flagQcRet = flagQcRet;
    data.ind_productTime = ind_productTime;
    data.ind_QC_GOOD = ind_QC_GOOD;
    data.ind_obsTime = ind_obsTime;
    data.soilMoistRet = soilMoistRet;
    data.VWCRet = VWCRet;
    data.refl_ref = refl_ref;
    data.refl_est = refl_est;
    data.refl_err = refl_err;
    data.penDepthHRet = penDepthHRet;
    data.penDepthVRet = penDepthVRet;
    data.reflErr = reflErr;
    data.typePOMERet = typePOMERet;
    data.saTime = saTime;
    data.saIter = saIter;
    data.saJump = saJump;
    data.saCostMin = saCostMin;
    data.runtime = runtime;
    data.saTimeMean = saTimeMean;
end