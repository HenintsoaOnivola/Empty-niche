function the_arr=read_file_to_array(fil_name)
    fil=fopen(fil_name);
    fil_cont=textscan(fil,'%s','Delimiter','\n');
    fclose(fil);
    fil_cont=fil_cont{1};
    if ~isempty(fil_cont)
        the_char=fil_cont{1};
        the_cell_char=regexp(the_char,'\t','split');
        the_cell_char=the_cell_char(~ismember(the_cell_char,''));
        the_first_arr=cellfun(@str2num,the_cell_char);
        the_arr=zeros(length(fil_cont),length(the_first_arr));
        for each_l=1:length(fil_cont)
            the_char=fil_cont{each_l};
            the_cell_char=regexp(the_char,'\t','split');
            the_cell_char=the_cell_char(~ismember(the_cell_char,''));
            the_arr(each_l,:)=cellfun(@str2num,the_cell_char);
        end
    else
        the_arr=[];
    end
end