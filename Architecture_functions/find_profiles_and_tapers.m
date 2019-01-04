function [ radius_taper volume_taper volume_coeff Sail_cylinders, Sail_radius,Volume_profile,Tot_Sail_area] = find_profiles_and_tapers(radius,cyl_length,comps,centres,height,top_of_trunk,PLOT)
    % Sail area and profile
    %Total sail area is the area projected in the xy plane.
    %Area = 2r*lcos(z/abs(vector))

    x_comp=comps(:,1); y_comp=comps(:,2); z_comp=comps(:,3);
    Sail_area=2*radius.*cyl_length.*cos(z_comp./(sqrt(x_comp.^2+y_comp.^2+z_comp.^2)));
    cyl_volume=cyl_length.*pi.*(radius.^2);
    Tot_Sail_area=sum(Sail_area);
    
    bad_n=100;
    for n=5:100;
        if n>=bad_n; continue; end
        h_inc=1:n;
        %height_incs=h_inc.*height./n-0.5*height./n;    
        for i=h_inc
            scaled_height=centres(:,3)/height;
            temp_rows=find(scaled_height>(i-1)/n & scaled_height<=i/n);
            Volume_profile_temp(i)=sum(cyl_volume(temp_rows));
        end
        
        %  BREAK HERE IF N volume profile has zeros
        if length(find(Volume_profile_temp==0))>=1;    bad_n=n;   continue;   end
        Volume_profile=Volume_profile_temp;
        %clearvars Sail_radius_x Sail_radius_y Sail_radius Sail_cylinders Volume_profile model2 v
        % =======================================
        for i=h_inc  % These don't get affected by the too high n
            scaled_height=centres(:,3)/height;
            temp_rows=find(scaled_height>(i-1)/n & scaled_height<=i/n);
            if length(temp_rows)>0
                Sail_radius_x(i)=max(centres(temp_rows,1))-min(centres(temp_rows,1));
                Sail_radius_y(i)=max(centres(temp_rows,2))-min(centres(temp_rows,2));
                Sail_radius(i)=mean([Sail_radius_x(i) Sail_radius_y(i)]);
            end
            Sail_cylinders(i)=sum(Sail_area(temp_rows));
        end
        % Mass taper
        h_inc=1:n;  
        height_incs=h_inc.*height./n-0.5*height./n; % mid points in the windows up the tree
        L2=(height-height_incs)./height;  % redefine downwards vector
        fit2=fitlm(log(L2),log(Volume_profile),'exclude',find(Volume_profile==0));
        model2=exp(fit2.Coefficients.Estimate(1)).*L2.^fit2.Coefficients.Estimate(2);
        volume_coeff=exp(fit2.Coefficients.Estimate(1));
        volume_taper=fit2.Coefficients.Estimate(2);
    end % end of loop over n's

    % =======================================
    % Diameter taper
    trunk=1:top_of_trunk; %temp to make fit
    L1=(height-centres(trunk,3))./height;
    fit1=fitlm(log(L1),log(radius(trunk))); 
    radius_taper=fit1.Coefficients.Estimate(2);
    model1=exp(fit1.Coefficients.Estimate(1)).*L1.^fit1.Coefficients.Estimate(2);
        
        
    if PLOT==1
        close all
        subplot(1,2,1)
        plot(Sail_cylinders,height_incs,'-o');   legend boxoff;  ylabel('Height (m)'); hold on;
        plot(Volume_profile,height_incs,'--')
        legend('Projected cylinder area','Volume profile','Location','SouthEast'); legend boxoff;   ylabel('Height (m)')

        if 1==1
            subplot(1,2,2)
            plot(Sail_radius_x,height_incs,'-o');   hold on;
            plot(Sail_radius_y,height_incs,'-x')
            plot(Sail_radius,height_incs,'--','LineWidth',1.2)
            legend('Sail radius x','Sail radius y','Sail radius mean','Location','SouthEast'); legend boxoff;   ylabel('Height (m)')
        end
    end
     if PLOT==2
         %%
        close all
        color1=brewermap(8,'blues');
        color2=brewermap(8,'dark2');
        subplot(1,2,1)
        plot(smooth(Volume_profile,4),height_incs,'--')
        ylabel('Height (m)'); hold on;
        
        subplot(1,2,2)
        rows0=find(Sail_radius==0);
        Sail_radius_x(rows0)=Sail_cylinders(rows0);
        Sail_radius_y(rows0)=Sail_cylinders(rows0);
        Sail_radius(rows0)=Sail_cylinders(rows0);
        plot(smooth(Sail_radius_x,4),height_incs,'Color',color1(4,:),'LineWidth',1.2);   hold on;
        plot(smooth(Sail_radius_y,4),height_incs,'Color',color1(4,:),'LineWidth',1.2)
        plot(smooth(Sail_radius,4),height_incs,'--','Color',color2(8,:),'LineWidth',1.2)
        box off; axis off; 
        %legend('Sail radius x','Sail radius y','Sail radius mean','Location','SouthEast'); legend boxoff;   ylabel('Height (m)')
        %%
     end
     if PLOT==3
         %%
        close all
        color1=brewermap(8,'oranges');
        color2=brewermap(8,'dark2');
        fig=figure;
        set(fig,'defaultAxesColorOrder',[color2(8,:); color2(2,:)]);
        Sail_radius(find(Sail_radius==0))=Sail_cylinders(find(Sail_radius==0));
        plot(height_incs,smooth(Sail_radius,4),'Color',color2(8,:),'LineWidth',1.2)
        h=ylabel('Sail area (m^2)','rotation',270); 
        xticks ''; ylim([0 25])
        yticks([0 25]); ytickangle(270); 
        set(gca, 'FontName', 'Helvetica','FontSize', 12) ;
        axis ij
        
        yyaxis right
        plot(height_incs,(smooth(Volume_profile,4)),'Color',color2(2,:),'LineWidth',1.2)
        ylabel('Volume (m^3)','rotation',270);       
        box off; xticks ''; ylim([0 4])
        yticks([0 4]); ytickangle(270); 
        axis ij
        %legend('Sail radius x','Sail radius y','Sail radius mean','Location','SouthEast'); legend boxoff;   ylabel('Height (m)')
        %%
     end
     if PLOT==5
         %%
        close all
        color1=brewermap(8,'blues');
        color2=brewermap(8,'dark2');
        h1=plot(smooth(Volume_profile,3),height_incs,'Color',color2(2,:),'LineWidth',1.6); hold on
        Sail_radius(find(Sail_radius==0))=Sail_cylinders(find(Sail_radius==0));
        h2=plot(smooth(Sail_radius,3)./5,height_incs,'Color',color2(8,:),'LineWidth',1.6);
        %area(cat(2,Sail_radius_x',Sail_radius_y') ,height_incs')
        box off; axis off
        legend([h2 h1],'','','Location','East'); legend boxoff;   
          set(gca, 'FontName', 'Helvetica','FontSize', 12) ;
        %%
     end
end

