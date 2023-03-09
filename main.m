% script main
%   
%   Guides work flow of retrieval simulation using SOCRATES-Retrieval. It
%   is recommended to run one section at a time. Some sections may take
%   several hours. Read carefully and understand the functionality before
%   running each section.
%
%   See also processUSCRNdata, runForwardModelMC, estimate_b_parameter, 
%   plotSynObs, runInverseModelMC, plotRetrieval.
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

clear all; close all;
addpath code/source/matlab/

%% Process source data
years = [2016, 2017]; % Set years [2016, 2017, 2018, 2019, 2020]
months = [1, 12]; % Set starting month and ending month
orbits = {'ISS', 'SSO'}; % Set orbit 'ISS', 'SSO'

processUSCRNdata(years, months, orbits);

%% Compile SMAT-retrieval
system('make clean; make');

%% Generate synthetic observations for all freq/pol
% It will generate hourly synthetic reflectivities on bare soil &
% vegetation at all four frequencies and co/cross-pol. No need to run
% multiple times once the data is generated.
% NOTE: This may take several hours and will overwrite any existing
% files. Use with caution.
runForwardModel(years, months, orbits);

%% Estimate the b parameter for all freq/pol
estimate_b_parameter(years, months, orbits);

%% Plot synthetic observations
plotSynObs(years, months, orbits);

%% Run retrieval simulation
% Set SoOp-R system parameters for Monte Carlo simulations using the
% synthetic observations generated above. Refer to Kim et al. 2023
% NOTE: This may take several hours and will overwrite any existing
% files. Use with caution.
freq = {[1,2,3,4]}; % Frequency combinations, 1: 137MHz, 2: 255MHz, 3:370MHz, 4:1575.42MHz
% freq = {[1,2,3,4], [1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4]};
polRx = {'XY'}; % Receive polarization, 'L', 'XY'
stddev = 0.04; % Measurement error parameter eta
period = 168; % Period P in hours
window = 12; % Observation time window T in hours
np = 24; % Number of processors

% Run Monte Carlo simulations for different SoOp-R parameters
runInverseModel(np, years, months, orbits, freq, polRx, stddev, period, window)

%% Plot retrieval results
plotRetrieval(years, months, orbits, freq, polRx, stddev, period, window)
