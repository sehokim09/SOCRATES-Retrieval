function estimate_b_parameter(years, months, orbit)
% function estimate_b_parameter
%   
%   Reads USCRN data and synthetic observations to estimate the b parameter
%   relating VOD and VWC for each frequency and polarizatoin. It plots and
%   saves the b parameters.
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

%% Directories
dir_in = './inputs/';
dir_forward = './results/forward/';
dir_processed = './data/processed/';
dir_out = './results/';
dir_figures = './figures/';
if ~exist(dir_figures, 'dir')
    mkdir(dir_figures);
end

%% PARAMETERS
trainDataRatio = 0.1;

freq_MHz = [137; 255; 370; 1575.42];
freqOfInterest = [1,2,3,4];
freqStr = strcat(num2str(freq_MHz(freqOfInterest)), 'MHz');

nYear = length(years);
nOrb = length(orbit);
nFreq = length(freqOfInterest);

% Read the list of selected stations
temp = readlines(strcat(dir_in, 'station_in.txt'),"EmptyLineRule","skip");
nStation = str2double(temp(1));
station_in = temp(2:end);

%% READ DATA
VWC.concat = [];
VOD.concat.H = [];
VOD.concat.V = [];
VOD.concat.R = [];
VOD.concat.L = [];
incAng.concat = [];
true = cell(nYear,1);
forward = cell(nYear,1);

surPOME = [];
bottPOME = [];
avgPOME = [];
for iYear = 1 : nYear
    for iOrb = 1 : nOrb
        true{iYear,iOrb} = readUSCRNNc(dir_processed, station_in, years(iYear), months, orbit{iOrb});
        forward{iYear,iOrb} = readForwardNc(dir_forward, station_in, years(iYear), months, orbit{iOrb});
        for iStation = 1 : nStation
            % Quality flag
            idx.flagQcForward{iYear,iOrb,iStation} = forward{iYear,iOrb}.flagQcForward{iStation} == 0;
            forward{iYear,iOrb}.VWCrefmodel{iStation}(~idx.flagQcForward{iYear,iOrb,iStation},:) = NaN;
            % Incidence angle
            incAng.raw{iYear,iOrb,iStation} = [true{iYear,iOrb}.incORBCOMM{iStation}, true{iYear,iOrb}.incMUOS{iStation},...
                                          true{iYear,iOrb}.incMUOS{iStation}, true{iYear,iOrb}.incGPS{iStation}];
            incAng.qc{iYear,iOrb,iStation} = incAng.raw{iYear,iOrb,iStation}(idx.flagQcForward{iYear,iOrb,iStation});
            % Ratio of vegetation reflectivity and bare soil reflectivity
            GammaRatio.L{iYear,iOrb,iStation} = forward{iYear,iOrb}.refl_V_L{iStation} ./ forward{iYear,iOrb}.refl_B_L{iStation};
            GammaRatio.R{iYear,iOrb,iStation} = forward{iYear,iOrb}.refl_V_R{iStation} ./ forward{iYear,iOrb}.refl_B_R{iStation};
            GammaRatio.X{iYear,iOrb,iStation} = forward{iYear,iOrb}.refl_V_X{iStation} ./ forward{iYear,iOrb}.refl_B_X{iStation};
            GammaRatio.Y{iYear,iOrb,iStation} = forward{iYear,iOrb}.refl_V_Y{iStation} ./ forward{iYear,iOrb}.refl_B_Y{iStation};
            % Vegetation optical depth
            VOD.H{iYear,iOrb,iStation} = -0.5 * cosd(incAng.raw{iYear,iOrb,iStation}) .* log(GammaRatio.Y{iYear,iOrb,iStation});        
            VOD.V{iYear,iOrb,iStation} = -0.5 * cosd(incAng.raw{iYear,iOrb,iStation}) .* log(GammaRatio.X{iYear,iOrb,iStation});        
            VOD.R{iYear,iOrb,iStation} = -0.5 * cosd(incAng.raw{iYear,iOrb,iStation}) .* log(GammaRatio.R{iYear,iOrb,iStation});
            VOD.L{iYear,iOrb,iStation} = -0.5 * cosd(incAng.raw{iYear,iOrb,iStation}) .* log(GammaRatio.L{iYear,iOrb,iStation});
            % Filtered VOD
            VOD.Hqc{iYear,iOrb,iStation} = VOD.H{iYear,iOrb,iStation}(idx.flagQcForward{iYear,iOrb,iStation},:);
            VOD.Vqc{iYear,iOrb,iStation} = VOD.V{iYear,iOrb,iStation}(idx.flagQcForward{iYear,iOrb,iStation},:);
            VOD.Rqc{iYear,iOrb,iStation} = VOD.R{iYear,iOrb,iStation}(idx.flagQcForward{iYear,iOrb,iStation},:);
            VOD.Lqc{iYear,iOrb,iStation} = VOD.L{iYear,iOrb,iStation}(idx.flagQcForward{iYear,iOrb,iStation},:);
            % Vegetation water content
            VWC.raw{iYear,iOrb,iStation} = sum(forward{iYear,iOrb}.VWCrefmodel{iStation}, 2);
            % Filtered VWC
            VWC.qc{iYear,iOrb,iStation} = VWC.raw{iYear,iOrb,iStation}(idx.flagQcForward{iYear,iOrb,iStation});
            % Data concatenation
            VWC.concat = [VWC.concat; VWC.qc{iYear,iOrb,iStation}];
            VOD.concat.H = [VOD.concat.H; VOD.Hqc{iYear,iOrb,iStation}];
            VOD.concat.V = [VOD.concat.V; VOD.Vqc{iYear,iOrb,iStation}];
            VOD.concat.R = [VOD.concat.R; VOD.Rqc{iYear,iOrb,iStation}];
            VOD.concat.L = [VOD.concat.L; VOD.Lqc{iYear,iOrb,iStation}];
            incAng.concat = [incAng.concat; incAng.qc{iYear,iOrb,iStation}];
            
            surPOME = [surPOME; forward{iYear,iOrb}.surPOME{iStation}(idx.flagQcForward{iYear,iOrb,iStation})];
            bottPOME = [bottPOME; forward{iYear,iOrb}.bottPOME{iStation}(idx.flagQcForward{iYear,iOrb,iStation})];
            avgPOME = [avgPOME; forward{iYear,iOrb}.avgPOME{iStation}(idx.flagQcForward{iYear,iOrb,iStation})];
        end
    end
end

%% RANDOM SELECTION OF TRAINING DATA
% Total number of data
nData = length(VWC.concat);
% Sampling training data uniformly at random
[VWC.train, idx.train] = datasample(VWC.concat, round(trainDataRatio*nData));
VOD.train.H = VOD.concat.H(idx.train,:);
VOD.train.V = VOD.concat.V(idx.train,:);
VOD.train.R = VOD.concat.R(idx.train,:);
VOD.train.L = VOD.concat.L(idx.train,:);

%% GENERATE THE b PARAMETER
% Linear regression
b.linReg.H = VWC.train \ VOD.train.H;
b.linReg.V = VWC.train \ VOD.train.V;
b.linReg.R = VWC.train \ VOD.train.R;
b.linReg.L = VWC.train \ VOD.train.L;
% Calculated VOD for training data
VOD.calc.H = b.linReg.H .* VWC.train;
VOD.calc.V = b.linReg.V .* VWC.train;
VOD.calc.R = b.linReg.R .* VWC.train;
VOD.calc.L = b.linReg.L .* VWC.train;

%% ERROR STATISTICS
% Root mean square error
b.RMSE.H = sqrt(mean((VOD.calc.H - VOD.train.H).^2, 1));
b.RMSE.V = sqrt(mean((VOD.calc.V - VOD.train.V).^2, 1));
b.RMSE.R = sqrt(mean((VOD.calc.R - VOD.train.R).^2, 1));
b.RMSE.L = sqrt(mean((VOD.calc.L - VOD.train.L).^2, 1));

% R-squared
b.Rsq.H = 1 - sum((VOD.calc.H - VOD.train.H).^2,1)./sum((VOD.calc.H - mean(VOD.train.H)).^2,1);
b.Rsq.V = 1 - sum((VOD.calc.V - VOD.train.V).^2,1)./sum((VOD.calc.V - mean(VOD.train.V)).^2,1);
b.Rsq.R = 1 - sum((VOD.calc.R - VOD.train.R).^2,1)./sum((VOD.calc.R - mean(VOD.train.R)).^2,1);
b.Rsq.L = 1 - sum((VOD.calc.L - VOD.train.L).^2,1)./sum((VOD.calc.L - mean(VOD.train.L)).^2,1);

%% WRITE DATA
filename_table = strcat(dir_in, 'b_parameters.txt');
if ~exist(filename_table, 'file')
    b_linReg_H = b.linReg.H(freqOfInterest)';
    b_linReg_V = b.linReg.V(freqOfInterest)';
    b_linReg_R = b.linReg.R(freqOfInterest)';
    b_linReg_L = b.linReg.L(freqOfInterest)';
    
    RMSE_H = b.RMSE.H(freqOfInterest)';
    RMSE_V = b.RMSE.V(freqOfInterest)';
    RMSE_R = b.RMSE.R(freqOfInterest)';
    RMSE_L = b.RMSE.L(freqOfInterest)';
    
    table_out = table(b_linReg_R, b_linReg_L, b_linReg_V, b_linReg_H,...
                       RMSE_R, RMSE_L, RMSE_V, RMSE_H,...
                      'RowNames', cellstr(freqStr));
    writetable(table_out, filename_table, 'Delimiter', '\t', 'WriteRowNames', 1);
end
if ~exist(strcat(dir_out, 'b_parameters.mat'), "file")
    save(strcat(dir_out, 'b_parameters.mat'), "b", "VWC", "VOD", '-mat')
end

%% PLOT
if ~exist(strcat(dir_figures, 'b_linearRegression.png'), "file")
    load(strcat(dir_out, 'b_parameters.mat'));
    dataColor = {'r.', 'g.', 'b.', 'c.', 'm.'};  
    lineColor = {'r', 'g', 'b', 'c', 'm'};
    
    % One-to-one Matchups - parameter b
    maxVal = max(max([VOD.calc.H, VWC.train]));
    xCenter = 0:0.01:maxVal;
    hf = figure('visible', 'off', 'Position', [0 0 1400 700]);
    ht = tiledlayout(2, 2, 'TileSpacing', 'Compact');
    % H-pol
    nexttile; hold on; grid on;
    for iFreq = 1 : nFreq
        hp(iFreq) = scatter(VWC.train, VOD.train.H(:,iFreq), dataColor{iFreq});  
    end
    tau = b.linReg.H .* xCenter';
    plot(xCenter, xCenter, 'k--', 'LineWidth', 1);
    for iFreq = 1 : nFreq
        plot(xCenter, tau(:,iFreq), lineColor{iFreq}, 'LineWidth', 1);
    end
    axis([0 maxVal 0 maxVal])
    legend(hp, strcat(freqStr, ', b=', num2str(b.linReg.H(freqOfInterest)', '%.2f'), ', RMSE=',  num2str(b.RMSE.H(freqOfInterest)', '%.4f')), 'Location', 'northwest');
    title('H-pol');
    % V-pol
    nexttile; hold on; grid on;
    for iFreq = 1 : nFreq
        hp(iFreq) = scatter(VWC.train, VOD.train.V(:,iFreq), dataColor{iFreq});  
    end
    tau = b.linReg.V .* xCenter';
    plot(xCenter, xCenter, 'k--', 'LineWidth', 1);
    for iFreq = 1 : nFreq
        plot(xCenter, tau(:,iFreq), lineColor{iFreq}, 'LineWidth', 1);
    end
    axis([0 maxVal 0 maxVal])
    legend(hp, strcat(freqStr, ', b=', num2str(b.linReg.V(freqOfInterest)', '%.2f'), ', RMSE=',  num2str(b.RMSE.V(freqOfInterest)', '%.4f')), 'Location', 'northwest');
    title('V-pol');
    % R-pol
    nexttile; hold on; grid on;
    for iFreq = 1 : nFreq
        hp(iFreq) = scatter(VWC.train, VOD.train.R(:,iFreq), dataColor{iFreq});  
    end
    tau = b.linReg.R .* xCenter';
    plot(xCenter, xCenter, 'k--', 'LineWidth', 1);
    for iFreq = 1 : nFreq
        plot(xCenter, tau(:,iFreq), lineColor{iFreq}, 'LineWidth', 1);
    end
    axis([0 maxVal 0 maxVal])
    legend(hp, strcat(freqStr, ', b=', num2str(b.linReg.R(freqOfInterest)', '%.2f'), ', RMSE=',  num2str(b.RMSE.R(freqOfInterest)', '%.4f')), 'Location', 'northwest');
    title('R-pol');
    % L-pol
    nexttile; hold on; grid on;
    for iFreq = 1 : nFreq
        hp(iFreq) = scatter(VWC.train, VOD.train.L(:,iFreq), dataColor{iFreq});  
    end
    tau = b.linReg.L .* xCenter';
    plot(xCenter, xCenter, 'k--', 'LineWidth', 1);
    for iFreq = 1 : nFreq
        plot(xCenter, tau(:,iFreq), lineColor{iFreq}, 'LineWidth', 1);
    end
    axis([0 maxVal 0 maxVal])
    legend(hp, strcat(freqStr, ', b=', num2str(b.linReg.L(freqOfInterest)', '%.2f'), ', RMSE=',  num2str(b.RMSE.L(freqOfInterest)', '%.4f')), 'Location', 'northwest');
    title('L-pol');
    
    xlabel(ht, 'Vegetation Water Content [kg/m^{3}]', 'FontSize', 11); ylabel(ht, 'Vegetation Optical Depth', 'FontSize', 11);
    title(ht, 'Linear Regression Relation Between VOD & VWC');
    print(hf, strcat(dir_figures, 'b_linearRegression'), '-dpng');
end

end