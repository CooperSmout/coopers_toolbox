function [subjects,files,newFolder,dataFolder,conditions] = prepDataFiles(filetype,dataFolder,newFolderName)


% function to identify subject data files within one folder that are all named
% according to same subj_xxx.xxx format, and create new folder if specified

% INPUTS
% filetype: 'xxx'
% dataFolder: string containing target directory (if left out will manually locate folder)
% newFolderName: string containing new folder name

% define sol (slash)
if ispc % windows
    sol = '\';
else % mac
    sol = '/';
end

% identify target folder
if nargin<2 || isempty(dataFolder)
    
    dataFolder = uigetdir(cd, ['Find the folder containing your ' filetype ' files...']);
    dataFolder = [dataFolder sol];
    if isequal(dataFolder,0)  % user clicked on cancel
        return;
    end
elseif ~exist(dataFolder,'dir')
    error('folder does not exist')
end

% Create new folder?
if nargin>2
    newFolder = [dataFolder newFolderName sol];
    mkdir(newFolder)
else
    newFolder=[];
end

% directory listing
all_files = dir(dataFolder);
[sorted_files,sorted_index] = sortrows({all_files.name}');

% remove files that aren't target filetype (and include underscore)
files = { };
for current_file = 1:length(sorted_files)
    extension_locn = strfind(sorted_files{current_file}, ['.' filetype]);
    if ~isempty(extension_locn) && strcmp(sorted_files{current_file}(extension_locn+1:end),filetype) && ~strcmp(sorted_files{current_file},'settings.mat')
        files = [files sorted_files(current_file)];
    end
end

% only continue if target files present
if isempty(files)
    warning(['No ' filetype ' files found in specified folder'])
    fclose('all');
    subjects = {};
    return;
end


% remove .xxx extensions
files = strrep(files, ['.' filetype], '');


% get unique list of participants and conditions
for current_file = 1:length(files)
    filename = files{current_file};
    delimeter = strfind(filename, '_');
    if delimeter
        participant{current_file} = filename(1:delimeter(1)-1);
        condition{current_file} = filename(delimeter(1)+1:end);
    else
        participant{current_file} = filename;
        condition{current_file} = 'Conditions not specified in filenames';
    end
end
subjects = unique(participant);
conditions = unique(condition);





