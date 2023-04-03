function dirs = setDirectories(dir_top)

% object
dirs.obj = strcat(dir_top, 'code/object/');

% data
dirs.geometry = strcat(dir_top, 'data/geometry/');
dirs.NDVI = strcat(dir_top, 'data/MODIS-NDVI/');

% input
dirs.in = strcat(dir_top, 'code/inputs/');

% output
dirs.out = strcat(dir_top, 'results/');
dirs.USCRN = strcat(dir_top, 'results/USCRN/');
dirs.processed = strcat(dir_top, 'results/data-processed/');
dirs.forward = strcat(dir_top, 'results/forward/');
dirs.inverse = strcat(dir_top, 'results/inverse/');
dirs.figures = strcat(dir_top, 'results/figures/');

% executable file
dirs.exe = strcat(dir_top, 'results/');

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