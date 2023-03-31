% script main
%   
%   Guides work flow of retrieval simulation using SOCRATES-Retrieval. It
%   is recommended to run one section at a time. Some sections may take
%   several hours. Read instructions carefully and understand the
%   functionality before running each section.
%
%   See also processUSCRNdata, runForwardModelMC, estimate_b_parameter, 
%   plotSynObs, runInverseModelMC, analyzeRetrieval.
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

% Create necessary directories
dir_obj = './code/object';
dir_forward = './results/forward';
dir_inverse = './results/inverse';
if ~exist(dir_obj, 'dir')
    mkdir(dir_obj);
end
if ~exist(dir_forward, 'dir')
    mkdir(dir_forward);
end
if ~exist(dir_inverse, 'dir')
    mkdir(dir_inverse);
end

% Compile SOCRATES-Retrieval
system('make clean; make');

%% Step 1: Process source data
% Select stations of interest from ./inputs/station_list.xlsx and put them
% into ./inputs/station_in.txt. Note that the first line of the text file 
% must be the total number of the stations of interest listed in the file.
% The name of each station must be provided from the second line separated
% by a line.

% Set years, [2016, 2017, 2018, 2019, 2020]
years = [2016];

% Set start month and end month. The start month must be less than the end
% month
months = [1, 12];

% Set orbit, 'ISS', 'SSO'
orbits = {'ISS'};

% Start processing
processUSCRNdata(years, months, orbits);

%% Step 2-1: Generate synthetic observations for all freq/pol
% It will generate hourly synthetic reflectivities on bare soil &
% vegetation at all four frequencies and all polarizations.
% NOTE: This may take several hours and will overwrite any existing
% files. Use with caution.

runForwardModel(years, months, orbits);

%% Step 2-2: Plot synthetic observations
plotSynObs(years, months, orbits);

%% Step 3: Estimate the b parameter for all freq/pol
estimate_b_parameter(years, months, orbits);

%% Step 4-1: Run retrieval simulation
% Set SoOp-R system parameters for Monte Carlo simulations using the
% synthetic observations generated above. Refer to Kim et al. 2023
% NOTE: This may take several hours and will overwrite any existing
% files. Use with caution. You can cancel the simulation by ctrl+c in
% Command Window

% Frequency combinations, 1: 137MHz, 2: 255MHz, 3:370MHz, 4:1575.42MHz
% syntax: {[1,2,3,4], [1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4]};
freq = {[1,2,3,4]};
% Receive polarization, 'L', 'XY'
polRx = {'XY'};
% Measurement error parameter eta
stddev = 0.04;
% Period P in hours
period = 168;
% Observation time window T in hours
window = 12;
% Number of processors. Adjust this value according to your system.
np = 24;

% Run Monte Carlo simulations for different SoOp-R parameters
runInverseModel(np, years, months, orbits, freq, polRx, stddev, period, window)

%% Step 4-2: Analyze retrieval results
analyzeRetrieval(years, months, orbits, freq, polRx, stddev, period, window)
