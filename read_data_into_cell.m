function the_cell_char=read_data_into_cell(fil_name)
    fil=fopen(fil_name);
    fil_cont=textscan(fil,'%s','Delimiter','\n');
    fclose(fil);
    fil_cont=fil_cont{1};
    the_cell_char=cell(length(fil_cont),1);
    for el=1:length(fil_cont)
        the_char=fil_cont{el};
        the_cell_char_el=regexp(the_char,'\t','split');
        the_cell_char_el=the_cell_char_el(~ismember(the_cell_char_el,''));
        the_cell_char{el}=cellfun(@str2num,the_cell_char_el);
    end
    the_cell_char=the_cell_char(~cellfun('isempty',the_cell_char));
    
end