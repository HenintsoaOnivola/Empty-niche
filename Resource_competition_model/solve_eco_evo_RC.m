function [time, solution,TE,YE,IE]=solve_eco_evo_RC(INI,rA,sig_c_an,sig_an,x0,kx,animal,past_time,the_count_integrate,count_sys)

  

    % competition kernal 
    alphaA=@(sig_c_an,x1,x2) exp(-(1/(2*sig_c_an^2))*((x1-x2).^2));
    alphaA_prime=@(sig_c_an,x1,x2) -(1/sig_c_an^2)*(x1-x2).*alphaA(sig_c_an,x1,x2);

    function K=k_x(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
    end

    function K=k_x_prime(k0,x0,a,x)
        K=-2*k0*(1/a^2)*(x-x0);
        K(abs(x-x0)>a)=0;
    end
    
    time_scale_and_mut_rate=10^(-6);
    extinction_threshold=10^(-8);
  
    function comp =competition_term(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end


    function comp =competition_term_prime(the_traits, the_pops)
        %traits and the_pops should be row vectors, and comp as well
        mat_comp=repmat(the_traits,[length(the_traits),1]);
        comp=alphaA_prime(sig_c_an,mat_comp,transpose(mat_comp));
        comp=the_pops*comp;
    end

    function gr_r=growth_rate(the_traits,the_pops)
     
        gr_r=the_pops.*(rA  - ( rA*competition_term(the_traits, the_pops)./k_x(kx,x0,sig_an,the_traits)) ) ;
    end

   
    function sel_gr=selection_gradient(the_traits,the_pops)
        
        first_term_num=-rA*((competition_term_prime(the_traits, the_pops).*k_x(kx,x0,sig_an,the_traits)) - (k_x_prime(kx,x0,sig_an,the_traits).*competition_term(the_traits, the_pops)));
        first_term_denum=k_x(kx,x0,sig_an,the_traits).^2;
        first_term=first_term_num./first_term_denum;
        
        sel_gr=the_pops.*(first_term);
    end
  

    options=odeset('Events', @events);
    
    function [value,isterminal,direction]=events(t,y)
        value=ones(3,1);
        y=transpose(y);
        
        %checking for non extinction
        density=y(1:animal);
        traits=y(animal+1:length(y));
      
        %extinction event
        %%%%%%%%%%%%%%%%%%%%%%
        if sum(density<=extinction_threshold)>=1
            value(1)=0;
        else
            value(1)=1;
        end
        
        
        %equilibrium event
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sel_gr=selection_gradient(traits,density);
        count_conv=0;
        for i=1:animal
            if abs(sel_gr(i))<10^(-8)
                count_conv=count_conv+1;
            end
        end

        if count_sys==1
            if count_conv==(animal)  && past_time+t>100
                value(2)=0;
            else
                value(2)=1;
            end
        else
        
            if count_conv==animal 
                value(2)=0;
            else
                value(2)=1;
            end
        end
    
        
        % not integrated or when oscillations are detected
        %%%%%%%%%%%%%%%%%%%%%%%%
        if the_count_integrate==0
            maximum_time=240;
        else
            maximum_time=240;
        end
        
        if maximum_time-toc<0
            value(3)=0;
        else
            value(3)=1;
        end
        
        isterminal=ones(size(value));
        direction=zeros(size(value));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
 
    function dxdt=ode_function_model(t,IN)
        
        A=IN(1:animal);
        x=IN(animal+1:end);
        
        A=transpose(A);
        x=transpose(x);
    
        dA=growth_rate(x,A);
      
        dx=time_scale_and_mut_rate*selection_gradient(x,A);
        
        dxdt=[dA dx];
        dxdt=transpose(dxdt);
        
    end
    [time,solution,TE,YE,IE]=ode15s(@ode_function_model,[0,5000],INI,options);      
end