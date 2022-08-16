function clear_files(file_name)
    ff=fopen(file_name,'w');
    fclose(ff);
end