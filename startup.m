function startup

if ispc
    addpath('P:\3011085.02\scripts\fieldtrip');
    addpath('P:\3011085.02\scripts\fieldtrip\qsub');
    addpath('P:\3011085.02\scripts\dyncon');
    addpath('M:\MATLAB');
    addpath('M:\MATLAB\cellfunction');
    cd('P:\3011085.02\scripts\erfosc');
else
    addpath('/project/3011085.02/scripts/fieldtrip');
    addpath('/project/3011085.02/scripts/fieldtrip/qsub');
    
    addpath('/project/3011085.02/scripts/erfosc');
    addpath('/home/electromag/matves/MATLAB');
    addpath('/home/electromag/matves/MATLAB/cellfunction/');
    % Dss and ASEO
    addpath('/home/electromag/matves/MATLAB/dss_aseo/peakfit');
    addpath('/home/electromag/matves/MATLAB/dss_aseo/ASEO');
    addpath('/home/electromag/matves/MATLAB/dss_aseo/dss2_1-0');
   
    cd('/project/3011085.02/scripts/erfosc');
end
ft_defaults

end