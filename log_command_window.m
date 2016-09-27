%% start function
workSpace = whos;
diary('tmpDiary') % save command window output
for i = 1:numel(workSpace) % list all workspace variables
    eval(workSpace(i).name)
end

%% start body



%% end function
save(fullfile([filename '.mat']), 'variable1', 'variable2') % example saving
diary off % toggly diary
movefile('tmpDiary', fullfile([filename '.txt'])); % give diary file same name as saved file