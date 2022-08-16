function the_arr=read_last_line_to_array(fil)
    ff=fopen(fil);
    fil_cont=textscan(ff,'%s','Delimiter','\n');
    fclose(ff);
    fil_cont=fil_cont{1};
    if ~isempty(fil_cont)
        the_char=fil_cont{end};
        the_cell_char=regexp(the_char,'\t','split');
        the_cell_char=the_cell_char(~ismember(the_cell_char,''));
        the_arr=cellfun(@str2num,the_cell_char);
    else
        the_arr=[];
    end
end