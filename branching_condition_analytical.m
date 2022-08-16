function info=branching_condition_analytical(INI,rA,rP,h,const,sig_c_an,sig_c_pl,sig_an,sig_pl,x0,y0,kx,ky,sig_m,animal,plant)

    ben=@(sig_m,x1,y1) exp(-(1/(2*sig_m^2))*((x1-y1).^2));   
  
    ben_prime=@(sig_m,x1,y1) -(1/sig_m^2)*(x1-y1).*ben(sig_m,x1,y1);
    
    ben_second=@(sig_m,x1,y1) -(1/sig_m^2)*(ben(sig_m,x1,y1)+(ben_prime(sig_m,x1,y1).*(x1-y1)));

    % competition kernal 
    alphaA=@(sig_c_an,x1,x2) exp(-(1/(2*sig_c_an^2))*((x1-x2).^2));
    alphaP=@(sig_c_pl,x1,x2) exp(-(1/(2*sig_c_pl^2))*((x1-x2).^2));
    
    alphaA_prime=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(x1-x2).*alphaA(sig_c_an,x1,x2);
    alphaP_prime=@(sig_c_pl,x1,x2) -(1/sig_c_pl^2)*(x1-x2).*alphaP(sig_c_pl,x1,x2);
    
    alphaA_second=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(alphaA(sig_c_an,x1,x2)+(alphaA_prime(sig_c_an,x1,x2).*(x1-x2)));
    alphaP_second=@(sig_c_pl,x1,x2) -(1/sig_c_pl^2)*(alphaP(sig_c_pl,x1,x2)+(alphaP_prime(sig_c_pl,x1,x2).*(x1-x2)));
    
   
    % carrying capacity function for animal
    function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_prime(k0,x0,a,x)
        K=-2*k0*(1/a^2)*(x-x0);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_second(k0,x0,a,x)
        K=-2*k0*(1/a^2);
        K=repmat(K,size(x));
        K(abs(x-x0)>a)=0;
    end
  
   
    function comp =competition_term_animal(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_plant(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaP(sig_c_pl,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_animal_prime(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA_prime(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_plant_prime(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaP_prime(sig_c_pl,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    
    function comp =competition_term_animal_second(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA_second(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function comp =competition_term_plant_second(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaP_second(sig_c_pl,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function benef=benefit(the_traits, the_other_traits)
        %return the matrix Bij
        RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
        SS=repmat(the_other_traits,[length(the_traits),1]);
        benef=ben(sig_m,RR,SS);
    end


    function benef_pr=benefit_prime(the_traits, the_other_traits)
        %return the matrix Bij
        RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
        PP=repmat(the_other_traits,[length(the_traits),1]);
        benef_pr=ben_prime(sig_m,RR,PP);
    end

  

    function benef_second=benefit_second(the_traits, the_other_traits)
        %return the matrix Bij
        RR=repmat(transpose(the_traits),[1,length(the_other_traits)]);
        PP=repmat(the_other_traits,[length(the_traits),1]);
        benef_second=ben_second(sig_m,RR,PP);
    end

   
    function disr_full=disruptive_selection(r,carry_cap,carry_cap_pr,carry_cap_sc,competition_term,competition_term_pr,competition_term_sc,benef_mat,benef_pr,benef_sc,the_other_pop)
        
        disr_comp_1=(competition_term_sc.*carry_cap - carry_cap_sc.*competition_term);
        disr_comp_2=(carry_cap).^(2);
        disr_comp_3=competition_term_pr.*carry_cap-carry_cap_pr.*competition_term;
        disr_comp_4=2*carry_cap_pr.*carry_cap;
        disr_comp=-r*( ((disr_comp_1.*disr_comp_2)-(disr_comp_3.*disr_comp_4))./((carry_cap).^(4)));
        
        num_mutu=the_other_pop*transpose(benef_mat);
        num_mutu_pr=the_other_pop*transpose(benef_pr);
        num_mutu_sc=the_other_pop*transpose(benef_sc);
        denum_mutu=1+h*the_other_pop*transpose(benef_mat);
        denum_mutu_pr=h*the_other_pop*transpose(benef_pr);
        denum_mutu_sc=h*the_other_pop*transpose(benef_sc);
        disr_mutu_num=((num_mutu_sc.*denum_mutu-denum_mutu_sc.*num_mutu).*(denum_mutu.^2))  - (2*denum_mutu.*denum_mutu_pr.*(num_mutu_pr.*denum_mutu-denum_mutu_pr.*num_mutu));
        disr_mutu_denum=denum_mutu.^(4);
        disr_mutu=disr_mutu_num./disr_mutu_denum;
        disr=disr_comp+const*disr_mutu;
        disr_full=[transpose(disr) transpose(disr_comp) transpose(const*disr_mutu)];

    end


    
    A=INI(1:animal);
    P=INI(animal+1:animal+plant);
    x=INI(animal+plant+1:2*animal+plant);
    y=INI(2*animal+plant+1:2*animal+2*plant);
  
    carry_cap_an=k_x(kx,x0,sig_an,x);
    carry_cap_an_pr=k_x_prime(kx,x0,sig_an,x);
    %disp (carry_cap_an_pr);
    carry_cap_an_sc=k_x_second(kx,x0,sig_an,x);
    %disp (carry_cap_an_sc);
    competition_term_an=competition_term_animal(x,A);
    competition_term_an_pr=competition_term_animal_prime(x,A);
    competition_term_an_sc=competition_term_animal_second(x,A);
    
    carry_cap_pl=k_x(ky,y0,sig_pl,y);
    carry_cap_pl_pr=k_x_prime(ky,y0,sig_pl,y);
    carry_cap_pl_sc=k_x_second(ky,y0,sig_pl,y);
    competition_term_pl=competition_term_plant(y,P);
    competition_term_pl_pr=competition_term_plant_prime(y,P);
    competition_term_pl_sc=competition_term_plant_second(y,P);
    
   

  
    benef_mat_an=benefit(x,y);
    benef_mat_pl=benefit(y,x);
    benef_mat_an_pr=benefit_prime(x, y);
    benef_mat_pl_pr=benefit_prime(y, x);
    benef_mat_an_sc=benefit_second(x, y);
    benef_mat_pl_sc=benefit_second(y, x);
   
  
    disr_an=disruptive_selection(rA,carry_cap_an,carry_cap_an_pr,carry_cap_an_sc,competition_term_an,competition_term_an_pr,competition_term_an_sc,benef_mat_an,benef_mat_an_pr,benef_mat_an_sc,P);
    disr_pl=disruptive_selection(rP,carry_cap_pl,carry_cap_pl_pr,carry_cap_pl_sc,competition_term_pl,competition_term_pl_pr,competition_term_pl_sc,benef_mat_pl,benef_mat_pl_pr,benef_mat_pl_sc,A);
    
    
    
    info={disr_an disr_pl};
end