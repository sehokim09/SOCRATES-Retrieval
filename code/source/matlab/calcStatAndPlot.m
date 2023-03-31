function calcStatAndPlot(year, months, orbit, polRx, freq, period, window, stddev, flagVeg, flagPlot)
% function calcStatAndPlot
%
%   Calculates, plots, and saves error statistics of retrieved data for a
%   given system parameters.
%
%   See also main, analyzeRetrieval, readUSCRNNc, readForwardNc,
%   readRetrievalNc.
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
dir_inverse = './results/inverse/';
dir_figures = './figures/inverse/';
dir_ret = sprintf('p%dw%df%s%sstd%d%s/', period, window, strrep(num2str(freq), ' ', ''), polRx, stddev*100, orbit);
dir_out_figure = strcat(dir_figures, dir_ret);
dir_out_stat = strcat(dir_inverse, dir_ret);
if(~exist(dir_out_figure, 'dir'))
    mkdir(dir_out_figure);
end
if(~exist(dir_out_stat, 'dir'))
    mkdir(dir_out_stat);
end

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
str_soilDepth = strcat(num2str(station_soilDepth(1,:)'*100),'cm');
idx_layer = floor((station_soilDepth(1,:)) / 0.01) + 1;

nPol = length(polRx);
nSoOp = length(freq);

freq_MHz = [137; 255; 370; 1575.42];
str_freq = strcat(num2str(freq_MHz(freq)), 'MHz');
    
% Read NetCDF
fprintf("Computing statistics ... ");
true = readUSCRNNc(dir_processed, station_in, year, months, orbit);
forward = readForwardNc(dir_forward, station_in, year, months, orbit);
ret = readRetrievalNc([dir_inverse, dir_ret], station_in, year, months, orbit, period, window, stddev, polRx, flagVeg);

% Post-processing 
for iStation = 1 : nStation
    idx.flagQcForward{iStation} = forward.flagQcForward{iStation} == 0;
    forward.VWCrefmodel{iStation}(~idx.flagQcForward{iStation},:) = NaN;
    VWCrefmodelSum{iStation} = sum(forward.VWCrefmodel{iStation}, 2);
    nSublayer = length(forward.soilDepthPOME{iStation});
    % SM POME reference
    soilMoistPOMEref{iStation} = forward.soilMoistPOME{iStation}(ret.ind_productTime{iStation}, :);
    soilMoistPOMEref{iStation} = soilMoistPOMEref{iStation}(ret.ind_QC_GOOD{iStation},1:nSublayer);
    idx.validref{iStation} = ~any(isnan(soilMoistPOMEref{iStation}),2);
    soilMoistPOMEref{iStation} = soilMoistPOMEref{iStation}(idx.validref{iStation},:);
    soilMoistPOMEref_insituDepth{iStation} = soilMoistPOMEref{iStation}(:,idx_layer);
    % SM POME estimate
    soilMoistPOMEest{iStation} = ret.soilMoistRet{iStation}(ret.ind_QC_GOOD{iStation},1:nSublayer);
    soilMoistPOMEest{iStation} = soilMoistPOMEest{iStation}(idx.validref{iStation},:);
    soilMoistPOMEest_insituDepth{iStation} = soilMoistPOMEest{iStation}(:, idx_layer);
    % VWC reference
    VWCrefmodelSum{iStation} = VWCrefmodelSum{iStation}(ret.ind_productTime{iStation});
    VWCrefmodelSum{iStation} = VWCrefmodelSum{iStation}(ret.ind_QC_GOOD{iStation});
    VWCrefmodelSum{iStation} = VWCrefmodelSum{iStation}(idx.validref{iStation});
    % VWC estimate
    VWCest{iStation} = ret.VWCRet{iStation}(ret.ind_QC_GOOD{iStation});
    VWCest{iStation} = VWCest{iStation}(idx.validref{iStation});
    % Reflectivity refernece and estimate
    for iPol = 1 : nPol
        refl_ref{iStation,iPol} = ret.refl_ref{iStation,iPol}(ret.ind_QC_GOOD{iStation},:);
        refl_ref{iStation,iPol} = refl_ref{iStation,iPol}(idx.validref{iStation},:);
        refl_est{iStation,iPol} = ret.refl_est{iStation,iPol}(ret.ind_QC_GOOD{iStation},:);
        refl_est{iStation,iPol} = refl_est{iStation,iPol}(idx.validref{iStation},:);
        % RMSE Reflectivity
        RMSE_refl(iStation,iPol,:) = sqrt(mean((refl_ref{iStation,iPol} - refl_est{iStation,iPol}).^2,1));
    end
    % Datetime of retrieval
    utcDatetimeRet{iStation} = true.utcDatetime{iStation}(ret.ind_productTime{iStation});
    utcDatetimeRet{iStation} = utcDatetimeRet{iStation}(ret.ind_QC_GOOD{iStation});
    utcDatetimeRet{iStation} = utcDatetimeRet{iStation}(idx.validref{iStation});
    % Incidence angle
    inc.ret{iStation} = [true.incORBCOMM{iStation}(ret.ind_productTime{iStation}), true.incMUOS{iStation}(ret.ind_productTime{iStation}),...
                         true.incMUOS{iStation}(ret.ind_productTime{iStation}), true.incGPS{iStation}(ret.ind_productTime{iStation})];
    inc.ret{iStation} = inc.ret{iStation}(ret.ind_QC_GOOD{iStation},:);
    inc.ret{iStation} = inc.ret{iStation}(idx.validref{iStation},:);
    % RMSE of Soil Moisture - all depths
    nRet(iStation,1) = sum(idx.validref{iStation});
    bias.SM(iStation,:) = mean(soilMoistPOMEest{iStation} - soilMoistPOMEref{iStation}, 1);
    RMSE(iStation,:) = sqrt(mean((soilMoistPOMEest{iStation} - soilMoistPOMEref{iStation}).^2, 1));
    % RMSE from surface to specific depth
    for iLayer = 1 : nSublayer
        RMSE_top(iStation,iLayer) = sqrt(mean(RMSE(iStation,1:iLayer).^2,2));
    end
    % RMSE of Vegetation Water Content
    bias.VWC(iStation) = mean(VWCest{iStation} - VWCrefmodelSum{iStation}, 1);
    RMSE_VWC(iStation,:) = sqrt(mean((VWCest{iStation} - VWCrefmodelSum{iStation}).^2, 1));
end

% RMSE of Soil Moisture - in-situ depths
bias.SM_insituDepth = bias.SM(:,idx_layer);
RMSE_insituDepth = RMSE(:,idx_layer);  
% RMSE top depth averaged over stations
RMSE_top_mean = sqrt(mean(RMSE_top.^2,1));
% RMSE top in-situ depths
RMSE_top_mean_insituDepth = RMSE_top_mean(idx_layer);

fprintf('done.\n');
%% OUTPUTS
dateStr = sprintf('%dM%dM%d', year, months(1), months(2));
if nStation < 3
    colDivisor = nStation/3;
else
    colDivisor = 1;
end
rowDivisor = ceil(nStation/3)/3;
if flagPlot
    dataColor = {'r.', 'g.', 'b.', 'c.', 'm.'};  
    fprintf('Saving plots in %s ... ', dir_out_figure);
    % One-to-one Matchups - Soil Moisture
    soilMoist_max = 0.6;
    for iStation = 1 : nStation
        soilMoist_max = max([soilMoist_max, max(max([soilMoistPOMEest_insituDepth{iStation}, soilMoistPOMEref_insituDepth{iStation}]))]);
    end
    xVSM = 0:0.001:soilMoist_max;    
    hf = figure('visible', 'off', 'Position', [0 0 2500*colDivisor 1250*rowDivisor]);
    ht = tiledlayout('flow', 'TileSpacing', 'Compact');
    for iStation = 1 : nStation
        nexttile; hold on;
        for iLayer = 1 : nLayer
            scatter(soilMoistPOMEref_insituDepth{iStation}(:,iLayer), soilMoistPOMEest_insituDepth{iStation}(:,iLayer), dataColor{iLayer});
        end
        plot(xVSM, xVSM, 'k--', 'LineWidth', 1);
        axis([0 soilMoist_max 0 soilMoist_max])
        grid on;
        legend(strcat(str_soilDepth, ' RMSE=',num2str(squeeze(RMSE_insituDepth(iStation,:)'), '%.3f'), ', Bias=', num2str(squeeze(bias.SM_insituDepth(iStation,:)'), '%.3f')), 'Location', 'northeastoutside');
        title(station_name_space{iStation})
    end
    xlabel(ht, 'Modeled Soil Moisture [m^3m^{-3}]', 'FontSize', 11); ylabel(ht, 'Retrieved Soil Moisture [m^3m^{-3}]', 'FontSize', 11);
    str_title = sprintf('RMSE: (5cm) %.3f, (10cm) %.3f, (20cm) %.3f, (50cm) %.3f, (100cm) %.3f', ...
        RMSE_top_mean_insituDepth(1), RMSE_top_mean_insituDepth(2), RMSE_top_mean_insituDepth(3), RMSE_top_mean_insituDepth(4), RMSE_top_mean_insituDepth(5));
    title(ht, str_title, 'FontSize', 11);
    print(hf, sprintf('%sUSCRN_%s_%s_SM_matchups', dir_out_figure, dateStr, orbit), '-dpng');
    
    % One-to-one Matchups - Vegetation Water Content
    VWC_max = 4;
    for iStation = 1 : nStation
        VWC_max = max([VWC_max, max(max([VWCrefmodelSum{iStation}, VWCest{iStation}]))]);
    end
    xVWC = 0:0.01:VWC_max;    
    hf = figure('visible', 'off', 'Position', [0 0 2500*colDivisor 1250*rowDivisor]);
    ht = tiledlayout('flow', 'TileSpacing', 'Compact');
    for iStation = 1 : nStation
        nexttile; hold on;
        scatter(VWCrefmodelSum{iStation}, VWCest{iStation}, '.k');
        plot(xVWC, xVWC, 'k--', 'LineWidth', 1);
        axis([0 VWC_max 0 VWC_max])
        grid on;
        legend(['RMSE = ', num2str(RMSE_VWC(iStation), '%.3f'), ', Bias=', num2str(bias.VWC(iStation), '%.3f')], 'Location', 'northeastoutside');
        title(station_name_space{iStation})
    end
    xlabel(ht, 'Modeled Vegetation Water Content [kg/m^{3}]', 'FontSize', 11); ylabel(ht, 'Retrieved Vegetation Water Content [kg/m^{3}]', 'FontSize', 11);
    print(hf, sprintf('%sUSCRN_%s_%s_VWC_matchups', dir_out_figure, dateStr, orbit), '-dpng');
    
    % One-to-one Matchups - Reflectivity
%     refl_max = 0.3;
%     for iPol = 1 : nPol
%         for iStation = 1 : nStation
%             refl_max = max([refl_max, max(max([refl_ref{iStation,iPol}, refl_est{iStation,iPol}]))]);
%         end
%         xRefl = 0:0.01:refl_max;    
%         hf = figure('visible', 'off', 'Position', [0 0 800 700]);
%         ht = tiledlayout('flow', 'TileSpacing', 'Compact');
%         for iStation = 1 : nStation
%             nexttile; hold on;
%             for iSoOp = 1 : nSoOp 
%                 scatter(refl_ref{iStation,iPol}(:,iSoOp), refl_est{iStation,iPol}(:,iSoOp), dataColor{iSoOp});    
%             end
%             plot(xRefl, xRefl, 'k--', 'LineWidth', 1);
%             axis([0 refl_max 0 refl_max])
%             grid on;
%             % legend(num2str(squeeze(RMSE_refl(iStation,iPol,:)), '%.4f'));
%             title(station_name_space{iStation})
%         end
%         xlabel(ht, 'Reference Reflectivity', 'FontSize', 11); ylabel(ht, 'Estimated Reflectivity', 'FontSize', 11);
%         legend(str_freq, 'Location', 'southoutside', 'orientation', 'horizontal');
%         print(hf, sprintf('%sUSCRN_%s_%s_%s_reflectivity_matchups', dir_out_figure, dateStr, orbit, polRx(iPol)), '-dpng');
%     end
    
    % Time series of soil moisture and VWC
    for iStation = 1 : nStation
        hf = figure('visible', 'off', 'Position', [0 0 1600 1200]);
        ht = tiledlayout(3, 1, 'TileSpacing', 'loose', 'Padding','loose');
    
        nexttile; hold on; grid on;
        for iLayer = 1 : nLayer
            hp_sm(iLayer) = plot(utcDatetimeRet{iStation}, soilMoistPOMEref_insituDepth{iStation}(:,iLayer), dataColor{iLayer}, 'LineStyle', 'None');
        end
        ylabel('VSM (m^3m^{-3})');
        legend(hp_sm, '5cm', '10cm', '20cm', '50cm', '100cm', 'Location', 'northeastoutside');
        title('Reference Soil Moisture')
    
        nexttile; hold on; grid on;
        for iLayer = 1 : nLayer
            hp_sm(iLayer) = plot(utcDatetimeRet{iStation}, soilMoistPOMEest_insituDepth{iStation}(:,iLayer), dataColor{iLayer}, 'LineStyle', 'None');
        end
        ylabel('VSM (m^3m^{-3})');
        legend(hp_sm, '5cm', '10cm', '20cm', '50cm', '100cm', 'Location', 'northeastoutside');
        title('Retrieved Soil Moisture')
    
        nexttile; hold on; grid on;
        hp_VWC(1) = plot(utcDatetimeRet{iStation}, VWCrefmodelSum{iStation}, dataColor{1}, 'LineStyle', 'None');
        hp_VWC(2) = plot(utcDatetimeRet{iStation}, VWCest{iStation}, dataColor{2}, 'LineStyle', 'None');
        ylabel('VWC');
        legend(hp_VWC, 'ref', 'est', 'Location', 'northeastoutside');
        title('Vegetation Water Content')
        
        sgtitle(strcat(station_name_space{iStation}, ', ', year));
        print(hf, sprintf('%sUSCRN_%s_%s_%s_TS_ret', dir_out_figure, dateStr, orbit, station_in{iStation}), '-dpng');
    end
    
    % POME vs Retrieval (entire column)
    hf = figure('visible', 'off', 'Position', [0 0 2500*colDivisor 1250*rowDivisor]);
    ht = tiledlayout('flow', 'TileSpacing', 'Compact');
    for iStation = 1 : nStation
        nexttile; hold on;
        plot(squeeze(RMSE(iStation, :)), forward.soilDepthPOME{iStation}*100, 'k-', 'LineWidth', 1)
        xlim([0 0.35])
        set(gca, 'ydir', 'reverse')
        grid on;
        legend(['Bias=', num2str(mean(bias.SM(iStation,:)),'%.3f')], 'Location', 'northeastoutside')
        title(station_name_space{iStation})
    end
    xlabel(ht, 'RMSE (m^3/m^{-3})', 'FontSize', 11)
    ylabel(ht, 'depth (cm)', 'FontSize', 11)
    str_title = sprintf('RMSE: (5cm) %.3f, (10cm) %.3f, (20cm) %.3f, (50cm) %.3f, (100cm) %.3f',...
        RMSE_top_mean_insituDepth(1), RMSE_top_mean_insituDepth(2), RMSE_top_mean_insituDepth(3), RMSE_top_mean_insituDepth(4), RMSE_top_mean_insituDepth(5));
    title(ht, str_title, 'FontSize', 11);
    print(hf, sprintf('%sUSCRN_%s_%s_SMP_RMSE', dir_out_figure, dateStr, orbit), '-dpng');

    % Simulated Annealing Performance
    for iStation = 1 : nStation
        hf = figure('visible', 'off');
        subplot(4,1,1); hold on; grid on;
        bar(utcDatetimeRet{iStation}, ret.saIter{iStation}(ret.ind_QC_GOOD{iStation}));
        ylabel('# of iterations');
        
        subplot(4,1,2); hold on; grid on;
        bar(utcDatetimeRet{iStation}, ret.saJump{iStation}(ret.ind_QC_GOOD{iStation}));
        ylabel('# of Jumps');
        
        subplot(4,1,3); hold on; grid on;
        bar(utcDatetimeRet{iStation}, ret.saTime{iStation}(ret.ind_QC_GOOD{iStation})/60);
        ylabel('CPU time (min)');
        title(sprintf("Elapsed real time: %.2f min", ret.runtime{iStation}/60)); 
        
        subplot(4,1,4); hold on; grid on;
        bar(utcDatetimeRet{iStation}, log10(ret.saCostMin{iStation}(ret.ind_QC_GOOD{iStation})));
        ylabel('log(cost)');
        xlabel('Date');
        print(hf, sprintf('%sUSCRN_%s_%s_%s_SA', dir_out_figure, dateStr, orbit, station_in{iStation}), '-dpng');
    end
    fprintf('done.\n');
end

% Write statistics
fprintf('Saving statistics in %s ... ', dir_out_stat);
table_out = table(RMSE,'RowNames', station_in, 'DimensionNames', {'Soil_Depth', 'Data'});
filename_table = sprintf('%sUSCRN_%s_p%dw%dstd%d_%s_RMSE_SMP_alldepth.txt', dir_out_stat, dateStr, period, window, stddev*100, orbit);
writetable(table_out, filename_table, 'Delimiter', '\t', 'WriteRowNames', 1);

RMSE = sqrt(mean(RMSE_insituDepth.^2,1))';
table_out = table(RMSE, 'RowNames', cellstr(str_soilDepth), 'DimensionNames', {'Soil_Depth', 'Data'});
filename_table = sprintf('%sUSCRN_%s_p%dw%dstd%d_%s_RMSE_SMP.txt', dir_out_stat, dateStr, period, window, stddev*100, orbit);
writetable(table_out, filename_table, 'Delimiter', '\t', 'WriteRowNames', 1);

table_out = table(RMSE_VWC,'RowNames', station_in, 'DimensionNames', {'Station', 'Data'});
filename_table = sprintf('%sUSCRN_%s_p%dw%dstd%d_%s_RMSE_VWC.txt', dir_out_stat, dateStr, period, window, stddev*100, orbit);
writetable(table_out, filename_table, 'Delimiter', '\t', 'WriteRowNames', 1);

nProfiles = [nRet; sum(nRet)];
table_out = table(nProfiles,'RowNames', [station_in; {'Total'}], 'DimensionNames', {'Station', 'Data'});
filename_table = sprintf('%sUSCRN_%s_p%dw%dstd%d_%s_nProfiles.txt', dir_out_stat, dateStr, period, window, stddev*100, orbit);
writetable(table_out, filename_table, 'Delimiter', '\t', 'WriteRowNames', 1);
fprintf('done.\n');
end