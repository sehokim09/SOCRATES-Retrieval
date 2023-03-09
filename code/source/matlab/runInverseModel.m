function runInverseModel(np, years, months, orbits, freq, polRx, stddev, period, window)
% function runInverseModel
%   
%   Updates input files for SOCRATES-Retrieval and run it as a inverse mode
%   for each system configuration and station to retrieve SM and VWC.
%   SOCRATES-retrieval requires forward data created by runForwardModel.
%   This may take several hours to complete the process.
%
%   See also main, overwriteLineInFile.
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

%% PARAMETERS
% Transmitter parameters
f_MHz = {'137', '255', '370', '1575.42'};
EIRP = {'1', '27', '37', '32'};
polRxForward = {'R', 'R', 'R', 'R'};
alt_km = {'750', '35786', '35786', '20180'};
bw_kHz = {'25', '25', '5000', '4092'};

nYear = length(years);
nOrb = length(orbits);
nPol = length(polRx);
nDev = length(stddev);
nFreq = length(freq);
nWindow = length(window);

%% GET DIRECTORIES
dir_in = './inputs/';
dir_out = './results/inverse/';
dir_inp_main = strcat(dir_in, 'Inp_Main.txt');
dir_inp_fixed = strcat(dir_in, 'Inp_Fixed.txt');
dir_inp_retrieval = strcat(dir_in, 'Inp_Retrieval.txt');
dir_inp_station = strcat(dir_in, 'station_in.txt');
dir_inp_station_temp = strcat(dir_in, 'station_in_temp.txt');

% Read the list of selected stations
temp = readlines(dir_inp_station,"EmptyLineRule","skip");
nStation = str2double(temp(1));
station_in = temp(2:end);

% Create temporary station input file
fileID = fopen(dir_inp_station_temp, 'w');
fprintf(fileID, '1\n');

line_inp_main = [9, 10];
line_inp_retrieval = [2, 3, 4];
line_inp_fixed = [19, 25, 31, 37];
for iFreq = 1 :nFreq
    for iYear = 1 : nYear
        % Update Inp_Main.txt
        str = sprintf('%d  %d                         ! Starting Year, Ending Year\n', years(iYear), years(end));
        overwriteLineInFile(dir_inp_main,line_inp_main(1),str);
        str = sprintf('%d     %d                           ! Starting Month, Ending Month\n', months(1), months(2));
        overwriteLineInFile(dir_inp_main,line_inp_main(2),str);
        for iOrb = 1 : nOrb
            for iPol = 1 : nPol
                for iDev = 1 : nDev
                    for iWindow = 1 : nWindow
                        folder_out = sprintf('p%dw%df%s%sstd%d%s', period, window(iWindow), strrep(num2str(freq{iFreq}), ' ', ''), polRx{iPol}, stddev(iDev)*100, orbits{iOrb});
                        if ~exist(strcat(dir_out, folder_out), "dir")
                            % Make output directory
                            system(['mkdir ', dir_out, folder_out])
                        end
                        % Update Inp_Retrieval.txt
                        str = sprintf('%s                             ! Receive Polarization(R,L,X,Y,RL,RX,RY,LX,LY,XY)\n', polRx{iPol});
                        overwriteLineInFile(dir_inp_retrieval,line_inp_retrieval(1),str);
                        str = sprintf('%.2f                          ! Standard Deviation of Reflectivity Error Ratio\n', stddev(iDev));
                        overwriteLineInFile(dir_inp_retrieval,line_inp_retrieval(2),str);
                        str = sprintf('%d  %d                         ! Measurement Period [hr], Time window [hr]\n', period, window(iWindow));
                        overwriteLineInFile(dir_inp_retrieval,line_inp_retrieval(3),str);
                        % Update Inp_Fixed.txt
                        str = sprintf('%s                           ! (For retrieval only) Orbit type (ISS, SSO)\n', orbits{iOrb});
                        overwriteLineInFile(dir_inp_fixed,11,str);
                        str = sprintf('%d                             ! Number of Transmitters\n', length(freq{iFreq}));
                        overwriteLineInFile(dir_inp_fixed,17,str);
                        for iLine = 1 : length(freq{iFreq})
                            str = sprintf('%s                           ! Transmitter Frequency [MHz]\n', f_MHz{freq{iFreq}(iLine)});
                            overwriteLineInFile(dir_inp_fixed,line_inp_fixed(iLine),str);
                            str = sprintf('%s                           ! Transmitter Altitude [km]\n', alt_km{freq{iFreq}(iLine)});
                            overwriteLineInFile(dir_inp_fixed,line_inp_fixed(iLine)+1,str);
                            str = sprintf('%s    %s                        ! Transmitter EIRP [dB], Polarization (R,L,X,Y)\n', EIRP{freq{iFreq}(iLine)}, polRxForward{freq{iFreq}(iLine)});
                            overwriteLineInFile(dir_inp_fixed,line_inp_fixed(iLine)+2,str);
                            str = sprintf('%s                            ! Bandwidth of channel [kHz]\n', bw_kHz{freq{iFreq}(iLine)});
                            overwriteLineInFile(dir_inp_fixed,line_inp_fixed(iLine)+4,str);
                        end
                        
                        for iStation = 1 : nStation
                            filename_ret = sprintf('USCRN_%dM%dM%d_p%dw%d_%s_std%d_%s_ret.nc', years(iYear), months(1), months(2), period, window(iWindow), station_in{iStation}, stddev(iDev)*100, orbits{iOrb});
                            fullpath = strcat(dir_out, folder_out, '/', filename_ret);
                            fprintf('Starting simulated retrieval %s ... \n', fullpath);
                            if exist(fullpath, "file")
                                fprintf('exist.\n');
                            else
                                % Update station_list.txt
                                str = sprintf('%s\n', station_in{iStation});
                                overwriteLineInFile(dir_inp_station_temp,2,str);

                                % Run simulation
                                fprintf('run SMAT-retrieval.\n');
                                status = system(['mpirun -np ', num2str(np), ' SMAT-retrieval i ', folder_out]);
                                if status == 0
                                    fprintf('%s complete.\n', fullpath);
                                else
                                    fprintf('Error: SMAT-retrieval %s\n', fullpath);
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Delete the temporary file
delete ./inputs/station_in_temp.txt

end