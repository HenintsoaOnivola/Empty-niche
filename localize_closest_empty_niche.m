function res=localize_closest_empty_niche(the_trait,res_pop,other_res_pop,res_tr,other_res_tr,r,h,const,sig_c,sig,x,kx,sig_m)

    left_inter=the_trait-0.5:0.001:the_trait;
    right_inter=the_trait:0.001:the_trait+0.5;
    
    fit_left=fitness_land_function(left_inter,res_pop,other_res_pop,res_tr,other_res_tr,r,h,const,sig_c,sig,x,kx,sig_m);
    fit_right=fitness_land_function(right_inter,res_pop,other_res_pop,res_tr,other_res_tr,r,h,const,sig_c,sig,x,kx,sig_m);
    
    pos_left=left_inter(fit_left>0.01);
    pos_right=right_inter(fit_right>0.01);
    
    dist_left=the_trait-max(pos_left);
    dist_right=min(pos_right)-the_trait;
    
    if isempty(dist_left)
        dist_left=Inf;
    end
    if isempty(dist_right)
        dist_right=Inf;
    end
    
    if dist_left<=dist_right
        res='left';
    else
        res='right';
    end
        
end