% bulk rename files

function renameFiles(folder)

if ~strcmp(folder(end),'/')
    folder = [folder '/'];
end
files = dir(folder);
for F = 3:length(files)
    filename = files(F).name;
    delimeter = strfind(filename, '_');
    newname = [filename(1:delimeter-1) filename(end-3:end)];
    movefile([folder filename],[folder newname]);
end
    