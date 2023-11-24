function [fpaths, fname, allinfo]=path2fldrfiles(searchFolder, searchPattern)
%% Return full paths to folders or files within a given search Folder whose filenames 
% match the given string pattern. 
% filetypes example = '*.m' if unspecified then there are . and .. that
% show up in the fldr_info file

if exist("searchPattern","var")
    allinfo=dir(fullfile(searchFolder, strcat('**/',searchPattern)));
else
    allinfo = dir(searchFolder);
    allinfo(ismember( {allinfo.name}, {'.', '..'})) = [];
    searchPattern='';
end
fpaths=join([{allinfo.folder}', {allinfo.name}'],filesep,2);
fname={allinfo.name}';
fprintf("Found %d files matching the pattern %s\n", length(fpaths), searchPattern )

end
