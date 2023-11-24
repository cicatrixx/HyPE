function fpathlist=path2fldrfiles(fldrpath,filetypes)
% filetypes example = '*.m' if unspecified then there are . and .. that
% show up in the fldr_info file
fldr_info  = dir(fullfile(fldrpath, filetypes));
fpathlist = fullfile(fldrpath, {fldr_info.name} )';
end
