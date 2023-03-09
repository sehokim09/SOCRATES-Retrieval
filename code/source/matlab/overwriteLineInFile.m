function overwriteLineInFile(Filename,LineNum,Text)
% function overwriteLineInFile
%   
%   Overwrites the line LineNum of Filename with Text. 
%
%   See also runForwardModel, runInverseModel.
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

    infile = fopen(Filename);
    outfile = fopen('MatlabTemp.txt','w');
    
    % Copy infile to outfile, up to line to be replaced
    for i=1:(LineNum-1)
        line = fgets(infile,512);
        fprintf(outfile,'%s',line);
    end
    
    % Read and discard line from infile
    line = fgets(infile,512);

    % Write replacement to outfile
    fprintf(outfile,'%s',Text);
    
    % Copy remainder of infile to outfile
    while(~feof(infile))
        line = fgets(infile,512);
        fprintf(outfile,'%s',line);
    end
    
    fclose(infile);
    fclose(outfile);         
    
    % Overwrite old file with modified file
    line = sprintf('mv MatlabTemp.txt %s',Filename);
    system(line);
end