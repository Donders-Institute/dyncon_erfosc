function startup

if ispc
    addpath('H:\common\matlab\fieldtrip');
    addpath('H:\common\matlab\fieldtrip/qsub');
    addpath('M:\MATLAB');
    cd 'M:\MATLAB'
else
    addpath('/home/electromag/matves/fieldtrip');
    addpath('/home/electromag/matves/fieldtrip/qsub');
    addpath('/home/electromag/matves/MATLAB');
    
    addpath /home/electromag/matves/Data/bbcm/
    
    cd '/home/electromag/matves/MATLAB'
end
ft_defaults

end