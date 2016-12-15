function []=crafs_multiple_refinements(filename)

file=sprintf('ROUTINES/%s.txt',filename);
fid = fopen(file);

command = fgets(fid);

while ischar(command)
    eval(command)
    command = fgets(fid);
end

fclose(fid);

end