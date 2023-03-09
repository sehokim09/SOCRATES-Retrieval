function plotRetrieval(years, months, orbits, freq, polRx, stddev, period, window)
% function plotRetrieval
%
%   Calculates, plots, and saves error statistics of retrieved data by
%   calling calcStatAndPlot in loops.
%
%   See also main, calcStatAndPlot.
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

    flagVeg = true;
    flagPlot = true;
    nYear = length(years);
    nOrb = length(orbits);
    nPol = length(polRx);
    nStd = length(stddev);
    nWindow = length(window);
    nFreq = length(freq);
    
    for iYear = 1 : nYear
        for iOrb = 1 : nOrb
            for iPol = 1 : nPol
                for iWindow = 1 : nWindow
                    for iStd = 1 : nStd
                        for iFreq = 1 : nFreq
                            calcStatAndPlot(years(iYear), months, orbits{iOrb}, polRx{iPol}, freq{iFreq}, period, window(iWindow), stddev(iStd), flagVeg, flagPlot);
                        end
                    end
                end
            end
        end
    end
    disp('Complete.')
end
