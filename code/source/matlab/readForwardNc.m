function data = readForwardNc(path, station_name, year, months, orbit)
% function readForwardNc
%   
%   Reads synthetic observations created by SOCRATES-Retrieval.
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

    for iStation = 1 : length(station_name)
        filename_forward = sprintf('USCRN_hourly_%dM%dM%d_%s_%s_forward.nc', ...
            year, months(1), months(2), orbit, station_name{iStation});
        fullpath_forward = strcat(path, filename_forward);

        if(exist(fullpath_forward, 'file')>0)
            % Read Forward NetCDF
            info_forward{iStation} = ncinfo(fullpath_forward);
            flagQcForward{iStation} = ncread(fullpath_forward, 'flagQcForward');
            typePOME{iStation} = ncread(fullpath_forward, 'typePOME');
            surPOME{iStation} = ncread(fullpath_forward, 'surPOME');
            bottPOME{iStation} = ncread(fullpath_forward, 'bottPOME');
            avgPOME{iStation} = ncread(fullpath_forward, 'avgPOME');
            soilDepthPOME{iStation} = ncread(fullpath_forward, 'soilDepthPOME');
            soilMoistPOME{iStation} = ncread(fullpath_forward, 'soilMoistPOME');

            reflCoefR_V_R{iStation} = ncread(fullpath_forward, 'reflCoefCoR_R_V');
            reflCoefI_V_R{iStation} = ncread(fullpath_forward, 'reflCoefCoI_R_V');
            reflCoefR_V_X{iStation} = ncread(fullpath_forward, 'reflCoefCoR_X_V');
            reflCoefI_V_X{iStation} = ncread(fullpath_forward, 'reflCoefCoI_X_V');
            reflCoefR_V_L{iStation} = ncread(fullpath_forward, 'reflCoefXR_R_V');
            reflCoefI_V_L{iStation} = ncread(fullpath_forward, 'reflCoefXI_R_V');
            reflCoefR_V_Y{iStation} = ncread(fullpath_forward, 'reflCoefXR_X_V');
            reflCoefI_V_Y{iStation} = ncread(fullpath_forward, 'reflCoefXI_X_V');

            reflCoefR_B_R{iStation} = ncread(fullpath_forward, 'reflCoefCoR_R_B');
            reflCoefI_B_R{iStation} = ncread(fullpath_forward, 'reflCoefCoI_R_B');
            reflCoefR_B_X{iStation} = ncread(fullpath_forward, 'reflCoefCoR_X_B');
            reflCoefI_B_X{iStation} = ncread(fullpath_forward, 'reflCoefCoI_X_B');
            reflCoefR_B_L{iStation} = ncread(fullpath_forward, 'reflCoefXR_R_B');
            reflCoefI_B_L{iStation} = ncread(fullpath_forward, 'reflCoefXI_R_B');
            reflCoefR_B_Y{iStation} = ncread(fullpath_forward, 'reflCoefXR_X_B');
            reflCoefI_B_Y{iStation} = ncread(fullpath_forward, 'reflCoefXI_X_B');

            refl_V_R{iStation} = ncread(fullpath_forward, 'reflect_CO_R_V');
            refl_V_L{iStation} = ncread(fullpath_forward, 'reflect_X_R_V');
            refl_V_X{iStation} = ncread(fullpath_forward, 'reflect_CO_X_V');
            refl_V_Y{iStation} = ncread(fullpath_forward, 'reflect_X_X_V');

            refl_B_R{iStation} = ncread(fullpath_forward, 'reflect_CO_R_B');
            refl_B_L{iStation} = ncread(fullpath_forward, 'reflect_X_R_B');
            refl_B_X{iStation} = ncread(fullpath_forward, 'reflect_CO_X_B');
            refl_B_Y{iStation} = ncread(fullpath_forward, 'reflect_X_X_B');

            penDepthH{iStation} = ncread(fullpath_forward, 'penDepthH');
            penDepthV{iStation} = ncread(fullpath_forward, 'penDepthV');
            
            VWCrefmodel{iStation} = ncread(fullpath_forward, 'VWCrefmodel');
            argH{iStation} = [ncread(fullpath_forward, 'argHr') + 1i* ncread(fullpath_forward, 'argHi')];
            argV{iStation} = [ncread(fullpath_forward, 'argVr') + 1i* ncread(fullpath_forward, 'argVi')];

            VODH{iStation} = ncread(fullpath_forward, 'VODH');
            VODV{iStation} = ncread(fullpath_forward, 'VODV');
        end
    end
    data.info_forward = info_forward;
    data.flagQcForward = flagQcForward;
    data.typePOME = typePOME;
    data.surPOME = surPOME;
    data.bottPOME = bottPOME;
    data.avgPOME = avgPOME;
    data.soilDepthPOME = soilDepthPOME;
    data.soilMoistPOME = soilMoistPOME;
    data.reflCoefR_V_R = reflCoefR_V_R;
    data.reflCoefI_V_R = reflCoefI_V_R;
    data.reflCoefR_V_L = reflCoefR_V_L;
    data.reflCoefI_V_L = reflCoefI_V_L;
    data.reflCoefR_V_X = reflCoefR_V_X;
    data.reflCoefI_V_X = reflCoefI_V_X;
    data.reflCoefR_V_Y = reflCoefR_V_Y;
    data.reflCoefI_V_Y = reflCoefI_V_Y;
    data.reflCoefR_B_R = reflCoefR_B_R;
    data.reflCoefI_B_R = reflCoefI_B_R;
    data.reflCoefR_B_L = reflCoefR_B_L;
    data.reflCoefI_B_L = reflCoefI_B_L;
    data.reflCoefR_B_X = reflCoefR_B_X;
    data.reflCoefI_B_X = reflCoefI_B_X;
    data.reflCoefR_B_Y = reflCoefR_B_Y;
    data.reflCoefI_B_Y = reflCoefI_B_Y;
    data.refl_V_R = refl_V_R;
    data.refl_V_L = refl_V_L;
    data.refl_V_X = refl_V_X;
    data.refl_V_Y = refl_V_Y;
    data.refl_B_R = refl_B_R;
    data.refl_B_L = refl_B_L;
    data.refl_B_X = refl_B_X;
    data.refl_B_Y = refl_B_Y;
    data.penDepthH = penDepthH;
    data.penDepthV = penDepthV;
    data.VWCrefmodel = VWCrefmodel;
    data.argH = argH;
    data.argV = argV;
    data.VODH = VODH;
    data.VODV = VODV;
end