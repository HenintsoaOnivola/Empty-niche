function merged_value=merge(the_trait,the_pop)
    morph_dictionary=containers.Map('KeyType','double','ValueType','double');
    morph_dictionary(the_trait(1))=the_pop(1);
    for each_el=2:length(the_trait)
        each_tr=the_trait(each_el);
        current_keys=morph_dictionary.keys();
        current_keys=cell2mat(current_keys(:));
      
        %equal_index=abs(each_tr-current_keys)<=10^(-4);
        equal_index=find(abs(current_keys-each_tr)<=10^(-2));
        if sum(equal_index)~=0
            equal_index=equal_index(1);
            the_key=current_keys(equal_index);
            morph_dictionary(the_key)=morph_dictionary(the_key)+the_pop(each_el);
            
        else
            morph_dictionary(each_tr)=the_pop(each_el);
        end
    end
    
    keys=morph_dictionary.keys();
    new_tr=zeros(1,length(keys));
    new_ab=zeros(1,length(keys));
          
    for ii=1:length(keys)
        new_tr(ii)=keys{ii};
        new_ab(ii)=morph_dictionary(keys{ii});
    end
    merged_value={new_tr;new_ab};
end