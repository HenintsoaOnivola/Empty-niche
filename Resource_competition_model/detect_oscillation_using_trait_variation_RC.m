function variation = detect_oscillation_using_trait_variation_RC(all_solutions,animal)
    aa=all_solutions(:,animal+1:2*animal);
    for each_aa=1:size(aa,2)
        variation=false;
        num_max=0;
        num_min=0;
        the_trait=aa(:,each_aa);
        for t=2:length(the_trait)-1
            if the_trait(t)>the_trait(t-1) && the_trait(t)>the_trait(t+1)
                num_max=num_max+1;
            end
            if the_trait(t)<the_trait(t-1) && the_trait(t)<the_trait(t+1)
                num_min=num_min+1;
            end
           
        end
        
        if num_max>3 && num_min>3
            variation=true;
        end
    end
    
end