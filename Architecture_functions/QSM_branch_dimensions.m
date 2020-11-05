function branch_dimensions = QSM_branch_dimensions(CylData)
% This function splits the QSM into 10 bins vertically and 10 bins
% horizontally and calculates branch diameters for each
% Could also save an 2D image of the raster of branch diameters
% Could also repeat for angle and length

    [num_cyls,radius,cyl_length,x,y1,z,x_comp,y_comp,z_comp,comps,centres,h,parent_index,ext_index,...
    branch_id,branch_order,running_no,TreeHeight,cyl_volume,tree_volume,canopy_volume,...
    canopy_vol_ratio,h_vol,Tot_volume,dbh] = architecture_definitions_TABLE(CylData);
    
    height=max(h);
    x_radius=max(centres(:,1))-min(centres(:,1))/2;
    scaled_height=centres(:,3)/height;
    scaled_x=abs(centres(:,1))/x_radius;
    %scaled_x=scaled_x-min(scaled_x); 
    
    n=10;
    vbr_min=nan(1,n);
    vbr_max=nan(1,n);
    vbr_mean=nan(1,n);
    vbr_median=nan(1,n);
    hbr_min=nan(1,n);
    hbr_max=nan(1,n);
    hbr_mean=nan(1,n);
    hbr_median=nan(1,n);
    
    for i=1:n  % height
        temp_rows=find(scaled_height>(i-1)/n & scaled_height<=i/n);
        if length(temp_rows)>0
            vbr_min(i)=min(radius(temp_rows));
            vbr_median(i)=median(radius(temp_rows));  
            vbr_mean(i)=mean(radius(temp_rows));
            vbr_max(i)=max(radius(temp_rows));
        end
    end
    
    for i=1:n      
        temp_rows=find(scaled_x>(i-1)/n & scaled_x<=i/n);
        if length(temp_rows)>0
            hbr_min(i)=min(radius(temp_rows));
            hbr_median(i)=median(radius(temp_rows));  
            hbr_mean(i)=mean(radius(temp_rows));
            hbr_max(i)=max(radius(temp_rows));
        end   
    end
    
    branch_dimensions=struct('height',height,'crown_radius',x_radius,'hbr_min',hbr_min,'hbr_max',hbr_max,...
        'hbr_median',hbr_median,'hbr_mean',hbr_mean,'vbr_min',vbr_min,'vbr_max',vbr_max,...
        'vbr_median',vbr_median,'vbr_mean',vbr_mean);
        
    %for i=1:n  
    %    for j=1:n % loop over x slices
    %        temp_rows=find(scaled_height>(i-1)/n & scaled_height<=i/n & scaled_x>(j-1)/n & scaled_x<=j/n);
    %        if length(temp_rows)>0
    %            min_branch_radii(i,j)=min(radius(temp_rows));
    %            median_branch_radii(i,j)=median(radius(temp_rows));  
    %           mean_branch_radii(i,j)=mean(radius(temp_rows));
    %            max_branch_radii(i,j)=max(radius(temp_rows));
    %        end
    %    end
    %end
    
    
    %image(flip(max_branch_radii*100)); colorbar
    
end

