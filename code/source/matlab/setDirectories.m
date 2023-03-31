function dirs = setDirectories(dir_top)

% object
dirs.obj = strcat(dir_top, 'code/object/');

% data
dirs.geometry = strcat(dir_top, 'data/geometry/');
dirs.processed = strcat(dir_top, 'data/processed/');
dirs.USCRN = strcat(dir_top, 'data/USCRN/');
dirs.NDVI = strcat(dir_top, 'data/MODIS-NDVI/');

% input
dirs.in = strcat(dir_top, 'inputs/');

% output
dirs.out = strcat(dir_top, 'results/');
dirs.forward = strcat(dir_top, 'results/forward/');
dirs.inverse = strcat(dir_top, 'results/inverse/');
dirs.figures = strcat(dir_top, 'figures/');

% executable file
dirs.exe = dir_top;

if ~exist(dirs.obj, 'dir')
    mkdir(dirs.obj);
end
if ~exist(dirs.forward, 'dir')
    mkdir(dirs.forward);
end
if ~exist(dirs.inverse, 'dir')
    mkdir(dirs.inverse);
end

end