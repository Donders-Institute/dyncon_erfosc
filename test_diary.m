function test_diary

%% start diary
workSpace = whos;
diaryname = sprintf('/project/3011085.02/scripts/erfosc/tmpDiary_%s.txt', datestr(now, 'dd.mm.yyyy_HH:MM:SS.FFF')); % working directory
diary(diaryname) % save command window output
fname = mfilename('fullpath')
datetime

fid = fopen(fullfile([fname '.m']));
tline = fgets(fid); % returns first line of fid
while ischar(tline) % at the end of the script tline=-1
    disp(tline) % display tline
    tline = fgets(fid); % returns the next line of fid
end
fclose(fid);

for i = 1:numel(workSpace) % list all workspace variables
    workSpace(i).name % list the variable name
    printstruct(eval(workSpace(i).name)) % show its value(s)
end

%% start body
variable1 = rand(1,10);
variable2 = variable1+1;

%% end diary

filename = '/project/3011085.02/scripts/erfosc/test_diary';
save(fullfile([filename '.mat']), 'variable1', 'variable2') % example saving
diary off
movefile(diaryname, fullfile([filename '.txt']));