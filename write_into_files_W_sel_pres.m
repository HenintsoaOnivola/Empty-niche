function write_into_files_W_sel_pres(path_file_name,the_data)
    if (isa(the_data,'char'))
        ffil=fopen(path_file_name,'Wb');
        fprintf(ffil,the_data);
        fclose(ffil);
    end
        
    if (isa(the_data,'cell'))
        ffil=fopen(path_file_name,'Wb');
        for ind_row=1:length(the_data)
            element_row=the_data{ind_row};
            for ind_cl=1:length(element_row)
                %fprintf(ffil,'%s',num2str(element_row(ind_cl)));
                fprintf(ffil,'%0.20f',element_row(ind_cl));
                fprintf (ffil,'\t');
            end
            fprintf(ffil,'\n');
        end
        fclose(ffil);
    end
    if (isa(the_data,'double'))
        ffil=fopen(path_file_name,'Wb');
        for ind_row=1:size(the_data,1)
            element_row=the_data(ind_row,:);
            for ind_cl=1:length(element_row)
                %fprintf(ffil,'%s',num2str(element_row(ind_cl)));
                fprintf(ffil,'%0.20f',element_row(ind_cl));
                fprintf (ffil,'\t');
            end
            fprintf(ffil,'\n');
        end
        fclose(ffil);
    end
    
end
