function startup

if ispc
    addpath('M:\fieldtrip');
    addpath('M:\fieldtrip\qsub');
    addpath('M:\dyncon');
    addpath('M:\MATLAB');
    addpath('M:\MATLAB\cellfunction');
    cd('M:\dyncon\erfosc');
else
    addpath('/home/electromag/matves/fieldtrip');
    addpath('/home/electromag/matves/fieldtrip/qsub');
    addpath('/home/electromag/matves/MATLAB');
    addpath('/home/electromag/matves/dyncon');
    addpath('/home/electromag/matves/Data/bbcm/');
    addpath('/home/electromag//matves/MATLAB/cellfunction/');
    % Dss and ASEO
    addpath('/home/electromag/matves/MATLAB/dss_aseo/peakfit');
    addpath('/home/electromag/matves/MATLAB/dss_aseo/ASEO');
    addpath('/home/electromag/matves/MATLAB/dss_aseo/dss2_1-0');
    
    cd('/home/electromag/matves/dyncon/erfosc');
end
ft_defaults

end