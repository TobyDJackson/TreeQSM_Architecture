function [CoV, CoV_stem, CoV_crown] = find_CoVs(centres,cyl_volume,branch_order,PLOT)
    % Centre of volume calculations
    CoV=sum(centres(:,[1 2 3]).*cyl_volume./sum(cyl_volume));
    CoV_stem=sum(centres(branch_order==0,[1 2 3]).*cyl_volume(branch_order==0)./sum(cyl_volume(branch_order==0)));
    CoV_crown=sum(centres(branch_order~=0,[1 2 3]).*cyl_volume(branch_order~=0)./sum(cyl_volume(branch_order~=0)));
    
    if PLOT==1
        scatter3(centres(:,1),centres(:,2),centres(:,3),cyl_volume*100);  hold on
        h1=scatter3(CoV(:,1),CoV(:,2),CoV(:,3),100,'x','red','LineWidth',2);
        h2=scatter3(CoV_stem(:,1),CoV_stem(:,2),CoV_stem(:,3),100,'x','black','LineWidth',2);
        h3=scatter3(CoV_crown(:,1),CoV_crown(:,2),CoV_crown(:,3),100,'x','green','LineWidth',2);
        legend([h1 h2 h3],'Whole tree','Stem','Crown','Location','SouthEast'); legend boxoff;
    end
    %% OLD calculation below here
    if 1==2;
        n_bins=100;
        scaled_height=centres(:,3)/height;
        for h_inc=1:n_bins
            rows=find(scaled_height>(h_inc-1)/n_bins & scaled_height<=h_inc/n_bins);
            Vol(h_inc)=sum(cyl_volume(rows));
        end
        Vol(1)=cyl_volume(1);
        cumVol=cumsum(Vol);
        temp=find(cumVol>=Tot_volume/2);
        h_ramp=linspace(0,height,n_bins);
        CoV_height=h_ramp(temp(1));
    end
end

