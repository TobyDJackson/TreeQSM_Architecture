function [max_crown_width,crown_height,crown_area, vol_asym, rel_rad_asym] = calculate_CrownAsymmetry(num_cyls,h,branch_order,x,y,z,cyl_volume,PLOT)
    % Crown Asymmetry measure-------------------------------
    theta=zeros(num_cyls,1);
    rr=zeros(num_cyls,8);
    vol=zeros(num_cyls,8);
    crown_height=max(h)-min(h(branch_order==1));
   % cen_x=x-x(1);  cen_y=y-y(1);
    cen_x=x-mean(x);  cen_y=y-mean(y);
    cyl_radius=sqrt(cen_x.^2+cen_y.^2);   
    
    %This part JUST calculates the angle theta (radians) of the cylinder w.r.t the 1st beam coordinates
    for i=1:num_cyls 
         if cen_x(i)>0  
            if cen_y(i)>0 %Top right quadrant, theta=0 to pi/2
                theta(i)=atan(abs(cen_x(i)/cen_y(i)));
            else          %Bottom right quadrant, theta=pi/2 to pi
                theta(i)=pi/2+atan(abs(cen_y(i)/cen_x(i)));
            end
         else
            if cen_y(i)<0 %Bottom left quadrant, theta=pi to 3pi/2
                theta(i)=pi+atan(abs(cen_x(i)/cen_y(i)));
            else          %Top left quadrant, theta=3pi/2 to 2pi
                theta(i)=3*pi/2+atan(abs(cen_y(i)/cen_x(i)));
            end
         end
    end %End loop over num_cyls i - we have now calculated theta(i) (radians)

    
    for sec=1:8 %loop over the sections (number of sections can change)
        temp=find(theta>(sec-1)*pi/4 & theta<sec*pi/4); 
        rr(temp,sec)=cyl_radius(temp); %This reads in the radius of each cylinder into the column corresponding to its segment 
        vol(temp,sec)=cyl_volume(temp);
        [crown_radius(1,sec) rows(1,sec)]=max(rr(:,sec)); %finds the max of each column - the radius in each of the 8 directions
        crown_volume(1,sec)=sum(vol(:,sec)); %sums the total branch volume in each section / column
    end %end loop over sections sec
       
    if PLOT==1
         scatter(cen_x,cen_y,cyl_volume*10000)
         title('Crown dimensions')
         hold on
         plot(0,0,'r+','MarkerSize',20,'LineWidth',3)
         plot(cen_x(rows),cen_y(rows),'b+','MarkerSize',10,'LineWidth',3)
    end
    for sec=2:8
        mean_x=mean([cen_x(rows(sec-1)),cen_x(rows(sec))]);
        mean_y=mean([cen_y(rows(sec-1)),cen_y(rows(sec))]);
        len_x=cen_x(rows(sec-1))-cen_x(rows(sec));
        len_y=cen_y(rows(sec-1))-cen_y(rows(sec));
        areas(1,sec-1)=0.5*sqrt(len_x.^2+len_y.^2)*sqrt(mean_x.^2+mean_y.^2);
        if PLOT==1;
            plot([cen_x(rows(sec-1)),cen_x(rows(sec))],[cen_y(rows(sec-1)),cen_y(rows(sec))])
            plot([0,mean_x],[0,mean_y])
        end
    end
    mean_x=mean([cen_x(rows(1)),cen_x(rows(8))]);
    mean_y=mean([cen_y(rows(1)),cen_y(rows(8))]);
    len_x=cen_x(rows(1))-cen_x(rows(8));
    len_y=cen_y(rows(1))-cen_y(rows(8));
    areas(1,8)=0.5*sqrt(len_x.^2+len_y.^2)*sqrt(mean_x.^2+mean_y.^2);
    if PLOT==1
        plot([cen_x(rows(1)),cen_x(rows(8))],[cen_y(rows(1)),cen_y(rows(8))])
        plot([0,mean_x],[0,mean_y])
        pause
        hold off
    end
    
    crown_area=sum(areas);
    max_crown_radius=max(crown_radius);
    rad_asym=max(crown_radius)-min(crown_radius); %This one isnt scaled by the tree size
    rel_rad_asym=(max(crown_radius)-min(crown_radius))/mean(crown_radius); %This one is scaled by crown size
    
    %Max crown radius
    crown_width(1)=crown_radius(1)+crown_radius(5);
    crown_width(2)=crown_radius(2)+crown_radius(6);
    crown_width(3)=crown_radius(3)+crown_radius(7);
    crown_width(4)=crown_radius(4)+crown_radius(8);
    max_crown_width=max(crown_width);
    
    vol_asym=max(crown_volume)-min(crown_volume); %This one isnt scaled by the tree size
    rel_vol_asym=(max(crown_volume)-min(crown_volume))/mean(crown_volume); %This one is scaled by crown size
    %--------------------end of Asymmetry-----------------------------
end

