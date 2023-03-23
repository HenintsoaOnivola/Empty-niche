function res_array=define_native_invader_function(an_tr,pl_tr,sep)

    sep=unique(sep);
    an_id=cell(length(an_tr),1);
    pl_id=cell(length(an_tr),1);

    previous=1;
    for s=sep
        for i=previous:s
            an=an_tr{i};
            an_id_ar=zeros(1,length(an));
            pl=pl_tr{i};
            pl_id_ar=zeros(1,length(pl));
            if previous~=1
                if i==previous
                    if length(an)>length(an_tr{i-1})
                        [consid_el, consid_ind]=setdiff(an,an_tr{i-1});
                        consid_ind=transpose(consid_ind);
                        [same_el,same_ind1,same_ind2]=intersect(an_tr{i-1},an,'stable');
                        an_id_ar(transpose(same_ind2))=an_id{i-1}(same_ind1);

                        for j=consid_ind
                            logical1=find(abs(an_tr{i-1}-an(j)+0.005)<1e-9);
                            if ~isempty(logical1)
                                an_id_ar(j)=an_id{i-1}(logical1);
                            end
                            logical2=find(abs(an_tr{i-1}-an(j)-0.005)<1e-9);
                            if ~isempty(logical2)
                                an_id_ar(j)=an_id{i-1}(logical2);
                            end
                            logical3=find(abs(an_tr{i-1}-an(j))<1e-9);
                            if ~isempty(logical3)
                                an_id_ar(j)=an_id{i-1}(logical3);
                            end

                            if (isempty(logical1) && isempty(logical2) && isempty(logical3))
                                an_id_ar(j)=1;
                            end
                        end

                    else
                        for j=1:length(an)
                            logical=find(abs(an_tr{i-1}-an(j))<1e-9);
                            an_id_ar(j)=an_id{i-1}(logical(1));
                        end

                    end

                    if length(pl)>length(pl_tr{i-1})
                        [consid_el, consid_ind]=setdiff(pl,pl_tr{i-1});
                        consid_ind=transpose(consid_ind);
                        [same_el,same_ind1,same_ind2]=intersect(pl_tr{i-1},pl,'stable');
                        pl_id_ar(transpose(same_ind2))=pl_id{i-1}(same_ind1);
                        for j=consid_ind
                            logical1=find(abs(pl_tr{i-1}-pl(j)+0.005)<1e-10);
                            if ~isempty(logical1)
                                pl_id_ar(j)=pl_id{i-1}(logical1);
                            end
                            logical2=find(abs(pl_tr{i-1}-pl(j)-0.005)<1e-10);
                            if ~isempty(logical2)
                                pl_id_ar(j)=pl_id{i-1}(logical2);
                            end
                            if (isempty(logical1) && isempty(logical2))
                                pl_id_ar(j)=1;
                            end
                        end

                    else
                        for j=1:length(pl)
                            logical=find(abs(pl_tr{i-1}-pl(j))<1e-9);
                            pl_id_ar(j)=pl_id{i-1}(logical(1));
                        end

                    end
                else
                    an_id_ar=an_id{i-1};
                    pl_id_ar=pl_id{i-1};
                end
            end
            an_id{i}=an_id_ar;
            pl_id{i}=pl_id_ar;

        end
  
        previous=s+1;

    end
    res_array={an_id,pl_id};

end


