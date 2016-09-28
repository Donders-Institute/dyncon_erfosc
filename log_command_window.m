%% start function
workSpace = whos;
diary('tmpDiary') % save command window output
fname = mfilename('fullpath') % display function path
datetime % display date and time

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



%% end function
save(fullfile([filename '.mat']), 'variable1', 'variable2') % example saving
diary off % toggly diary
movefile('tmpDiary', fullfile([filename '.txt'])); % give diary file same name as saved file