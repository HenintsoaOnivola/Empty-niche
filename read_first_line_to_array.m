function the_arr=read_first_line_to_array(fil_name)
    fil=fopen(fil_name);
    line=fgetl(fil);
    fclose(fil);
    if ~isempty(line)
        the_cell_char=regexp(line,'\t','split');
        the_cell_char=the_cell_char(~ismember(the_cell_char,''));
        the_arr=cellfun(@str2num,the_cell_char);
    else
        the_arr=[];
    end
end
    