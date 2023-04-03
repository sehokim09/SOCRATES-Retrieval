function runForwardModel(dirs, years, months, orbits)
% function runForwardModelMC
%   
%   Updates input files for SOCRATES-Retrieval and run it as a forward mode
%   for each orbit to generate synthetic observations.
%   SOCRATES-retrieval requires USCRN dataset created by processUSCRNdata.
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

% Directories
dir_inp_main = strcat(dirs.in, 'Inp_Main.txt');
dir_inp_fixed = strcat(dirs.in, 'Inp_Fixed.txt');

for iOrb = 1 : length(orbits)
    fprintf('Starting synthetic observation generation ... @ %s orbit\n', orbits{iOrb});
    % Update Inp_Main.txt
    str = sprintf('%d  %d                         ! Starting Year, Ending Year\n', years(1), years(end));
    overwriteLineInFile(dir_inp_main,9,str);
    str = sprintf('%d     %d                           ! Starting Month, Ending Month\n', months(1), months(2));
    overwriteLineInFile(dir_inp_main,10,str);
    % Update Inp_Fixed.txt
    str = sprintf('%s                           ! (For retrieval only) Orbit type (ISS, SSO)\n', orbits{iOrb});
    overwriteLineInFile(dir_inp_fixed,11,str);
    % Run simulation
    system([dirs.exe, 'SOCRATES-Retrieval f']);
end

end