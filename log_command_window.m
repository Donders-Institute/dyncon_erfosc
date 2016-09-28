%% start function
workSpace = whos;
diary('tmpDiary') % save command window output
mfilename('fullpath')
datetime
for i = 1:numel(workSpace) % list all workspace variables
    workSpace(i).name % list the variable name
    printstruct(eval(workSpace(i).name)) % show its value(s)
end

%% start body



%% end function
save(fullfile([filename '.mat']), 'variable1', 'variable2') % example saving
diary off % toggly diary
movefile('tmpDiary', fullfile([filename '.txt'])); % give diary file same name as saved file