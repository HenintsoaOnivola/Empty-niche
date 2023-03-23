function generate_traits_phase_portrait(rA,rP,const,h,x0,y0,kx,ky,log_sig_c_an,log_sig_m,log_sig_an)
    
    sig_c_an=exp(log_sig_c_an);
    sig_c_pl=sig_c_an;
    sig_m=exp(log_sig_m);
    
    sig_an=log_sig_an;
    sig_pl=sig_an;
    
    %%start with defining all possible ancestral pairs over the trait range
    anc_an=2:0.1:4;
    anc_an=[anc_an 4.05];
    anc_pl=1:0.1:3;

    all_pairs=zeros(5,2);
    count_pair=1;
    for i=1:length(anc_an)
        for j= anc_pl
            if rem(i,2)==1
                all_pairs(count_pair,:)=[anc_an(i) j];
            else
                all_pairs(count_pair,:)=[anc_an(i)-0.05 j];
            end
            count_pair=count_pair+1;
        end
    end

    for i=2:0.05:3
        for j=2:0.05:3
            if i>=j
                all_pairs(count_pair,:)=[i j];
                count_pair=count_pair+1;
            end
        end
    end

    all_pairs=unique(all_pairs,'rows');

    %start solving the eco-evolutionary dynamics for each ancestral pairs
    %unti ecological stability
    file_content=cell(8*400,1);
    file_content{1}=[2 4 1 3];
    count_file_line=2;
    %start the integration over initial values
    for pp=1:size(all_pairs,1)
        ini_an=all_pairs(pp,1);
        ini_pl=all_pairs(pp,2);
        disp ([ini_an ini_pl]);
        %check that it is inside the viable range
        if carrying_capacity(kx,x0,sig_an,ini_an)>0 && carrying_capacity(ky,y0,sig_pl,ini_pl)>0
            disp ('inside the range');
            infinite=200;
            negative_pop=false;
            convergence=false;
            count_integrate=0;
            time_plot=zeros(7e005,1);
            time_plot(:)=NaN;
            sol_plot=zeros(infinite*5000,4);
            sol_plot(:)=NaN;
            time_index=0;
            previous_t_end=0;
            ini_cond=[1 1 ini_an ini_pl];
            while ( (convergence==false) && count_integrate<=infinite )
                [time,sol,TE,YE,IE]=solve_eco_evo_phase_portrait(ini_cond,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m);
                if count_integrate==infinite
                    disp ('not convergent');
                    IE=1;
                end
                current_time_range=time+previous_t_end;
                time_plot(time_index+1:time_index+numel(time))=current_time_range;
                sol_plot(time_index+1:time_index+numel(time),:)=sol;
                previous_t_end=time_plot(time_index+numel(time));
                time_index=time_index+numel(time);

                sel_grad=selection_gradient_function([sol(end,1) sol(end,2) sol(end,3) sol(end,4)],rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,1,1);
                sel_grad_an=sel_grad{1};
                sel_grad_pl=sel_grad{2};

                the_threshold=10^(-8);

                if abs(sel_grad_an)<the_threshold && abs(sel_grad_pl)<the_threshold
                    convergence=true;
                end

                ini_cond=sol(end,:);
                count_integrate=count_integrate+1;
            end

            time_plot=time_plot(~isnan(time_plot));
            sol_plot=sol_plot(1:numel(time_plot),:);

            %%%%%%%%%%%%%%%
            %write into the files
            if negative_pop==false

                ini_conds=[ini_an ini_pl];
                an_pop=sol_plot(:,1);
                pl_pop=sol_plot(:,2);
                an_way=sol_plot(:,3);
                pl_way=sol_plot(:,4);

                file_content{count_file_line+0}=ini_conds;
                file_content{count_file_line+1}=time_plot;
                file_content{count_file_line+2}=an_pop;
                file_content{count_file_line+3}=pl_pop;
                file_content{count_file_line+4}=an_way;
                file_content{count_file_line+5}=pl_way;
                file_content{count_file_line+6}='not extinct';
                count_file_line=count_file_line+7;
            else
                disp ('negative');
                ini_conds=[ini_an ini_pl];
                an_pop=sol_plot(:,1);
                pl_pop=sol_plot(:,2);
                an_way=sol_plot(:,3);
                pl_way=sol_plot(:,4);

                file_content{count_file_line+0}=ini_conds;
                file_content{count_file_line+1}=time_plot;
                file_content{count_file_line+2}=an_pop;
                file_content{count_file_line+3}=pl_pop;
                file_content{count_file_line+4}=an_way;
                file_content{count_file_line+5}=pl_way;
                file_content{count_file_line+6}='extinct';
                count_file_line=count_file_line+7;
            end


        else
            disp ('outside the range');
            ini_conds=[ini_an ini_pl];
            file_content{count_file_line+0}=ini_conds;
            file_content{count_file_line+1}='outside range';
            file_content{count_file_line+2}='outside range';
            file_content{count_file_line+3}='outside range';
            file_content{count_file_line+4}='outside range';
            file_content{count_file_line+5}='outside range';
            file_content{count_file_line+6}='outside range';
            count_file_line=count_file_line+7;
        end
    end
    file_name=strcat('logSig_an=',num2str(log_sig_an),'logSig_m=',num2str(log_sig_m),'logSig_C=',num2str(log_sig_c_an),'.dat');
    ffil=fopen(strcat('./Results_phase_portrait/',file_name),'w');
    fclose(ffil);

    new_file_content=file_content(1:count_file_line-1);
    file_name='data_for_FigB1.dat';
    ffil=fopen(file_name,'a');
    for el_file=1:length(new_file_content)
        if isa(new_file_content{el_file},'double')==1
            the_line=new_file_content{el_file};
            for ind=1:length(the_line)
                fprintf(ffil,'%s',num2str(the_line(ind)));
                fprintf (ffil,'\t');
            end
            fprintf(ffil,'\n');
        end
        if isa(new_file_content{el_file},'char')==1
            fprintf(ffil,'%s',new_file_content{el_file});
            fprintf(ffil,'\n');
        end
    end
    fclose(ffil);

end