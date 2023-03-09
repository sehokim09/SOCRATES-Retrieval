function plotSynObs(years, months, orbits)
% function plotSynObs
%   
%   Reads and plots USCRN data and synthetic observations.
%
%   See also main, readUSCRNNc, readForwardNc.
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

% Directories
dir_in = './inputs/';
dir_processed = './data/processed/';
dir_forward = './results/forward/';
dir_figures = './figures/forward/';

% Read the list of selected stations
temp = readlines(strcat(dir_in, 'station_in.txt'),"EmptyLineRule","skip");
nStation = str2double(temp(1));
station_in = temp(2:end);
station_name_space = strrep(station_in,'_', ' ');

% Read data of the selected stations from the station table
station_info = readcell(strcat(dir_in, 'station_list.xlsx'), 'NumHeaderLines', 1);
station_name_list = station_info(:,1);
[~, idx_station_in] = ismember(station_in, station_name_list);
station_soilDepth = str2num(cell2mat(station_info(idx_station_in,6)));
nLayer = length(station_soilDepth(1,:));
idx_layer = floor((station_soilDepth(1,:)) / 0.01) + 1;

freq_MHz = [137; 255; 370; 1575.42];
freqOfInterest = [1,2,3,4]; % Change this to pick frequencies of interest
freqStr = strcat(num2str(freq_MHz(freqOfInterest)), 'MHz');
nSoOp = length(freqOfInterest);
nYear = length(years);
nOrb = length(orbits);

% Read data
pol = {'RHCP','LHCP','X-pol','Y-pol'};
dataColor = {'r.', 'g.', 'b.', 'c.', 'm.'};  
fprintf("Saving plots in %s ... ", dir_figures);
for iYear = 1 : nYear
    for iOrb = 1 : nOrb
        true = readUSCRNNc(dir_processed, station_in, years(iYear), months, orbits{iOrb});
        forward = readForwardNc(dir_forward, station_in, years(iYear), months, orbits{iOrb});
        for iStation = 1 : nStation
            soilMoistPOMEref_insituDepth{iStation} = forward.soilMoistPOME{iStation}(:,idx_layer);
            refl_ref{iStation,1} = forward.refl_V_R{iStation};
            refl_ref{iStation,2} = forward.refl_V_L{iStation};
            refl_ref{iStation,3} = forward.refl_V_X{iStation};
            refl_ref{iStation,4} = forward.refl_V_Y{iStation};
            VWCrefmodelSum{iStation} = sum(forward.VWCrefmodel{iStation}, 2);
            
            % Plot Soil moisture and reflectivity
            hf = figure('visible', 'off', 'Position', [0 0 1600 1200]);
            ht = tiledlayout(7, 1, 'TileSpacing', 'loose', 'Padding','loose');
        
            nexttile; hold on; grid on;
            for iLayer = 1 : nLayer
                hp_sm(iLayer) = plot(true.utcDatetime{iStation}, soilMoistPOMEref_insituDepth{iStation}(:,iLayer), dataColor{iLayer}, 'LineStyle', 'None');
            end
            ylabel('VSM (m^3m^{-3})');
            legend(hp_sm, '5cm', '10cm', '20cm', '50cm', '100cm', 'Location', 'northeastoutside');
            title('Reference Soil Moisture')
            for iPol = 1 : 4
                nexttile; hold on; grid on;
                for iSoOp = 1 : nSoOp
                    hp_refl(iSoOp) = plot(true.utcDatetime{iStation}, refl_ref{iStation,iPol}(:,iSoOp), dataColor{iSoOp}, 'LineStyle', 'None');
                end
                ylabel('\Gamma_{ref} (linear)');
                legend(hp_refl, freqStr, 'Location', 'northeastoutside');
                title(strcat('Reflectivity Reference - ', pol{iPol}))
            end
            nexttile; hold on; grid on;
            hp_VWC(1) = plot(true.utcDatetime{iStation}, VWCrefmodelSum{iStation}, dataColor{1}, 'LineStyle', 'None');
            ylabel('VWC');
            legend(hp_VWC, 'ref', 'Location', 'northeastoutside');
            title('Vegetation Water Content')
        
            nexttile; hold on; grid on;
            hp_inc(1) = plot(true.utcDatetime{iStation}, true.incORBCOMM{iStation}, dataColor{1}, 'LineStyle', 'None');
            hp_inc(2) = plot(true.utcDatetime{iStation}, true.incMUOS{iStation}, dataColor{2}, 'LineStyle', 'None');
            hp_inc(3) = plot(true.utcDatetime{iStation}, true.incGPS{iStation}, dataColor{3}, 'LineStyle', 'None');
            ylabel('\theta_s');
            legend(hp_inc, 'ORBCOMM','MUOS','GPS', 'Location', 'northeastoutside');
            title('Incidence Angles')
            
            sgtitle(sprintf('%s, %d', station_name_space{iStation}, years(iYear)));
            print(hf, sprintf('%sUSCRN_%dM%dM%d_%s_%s_TS', dir_figures, years(iYear), months(1), months(2), orbits{iOrb}, station_in{iStation}), '-dpng');
        end
    end
end
fprintf("Complete.\n");

end