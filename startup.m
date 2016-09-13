function startup

if ispc
    addpath('M:\fieldtrip');
    addpath('M:\fieldtrip\qsub');
    addpath('M:\dyncon');
    addpath('M:\MATLAB');
    cd('M:\dyncon\erfosc');
else
    addpath('/home/electromag/matves/fieldtrip');
    addpath('/home/electromag/matves/fieldtrip/qsub');
    addpath('/home/electromag/matves/MATLAB');
    addpath('home/electromag/matves/dyncon');
    addpath('/home/electromag/matves/Data/bbcm/');
    
    cd('home/electromag/matves/dyncon/erfosc');
end
ft_defaults

end