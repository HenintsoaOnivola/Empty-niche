function the_mat=cell_to_mat_fill_nan(the_cell)
    the_rows = cellfun(@numel,the_cell);
    the_mat=nan(length(the_cell),max(the_rows));
    for k=1:length(the_cell) 
        the_mat(k,1:the_rows(k))=the_cell{k};
    end
end