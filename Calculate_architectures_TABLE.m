function [ architectures ] = Calculate_architectures_TABLE( QSM,PLOT )
    
   
    %%  Definitions and simple measures
    [num_cyls,radius,cyl_length,x,y,z,x_comp,y_comp,z_comp,comps,centres,h,parent_index,ext_index,...
    branch_id,branch_order,running_no,height,cyl_volume,tree_volume,canopy_volume,...
    canopy_vol_ratio,h_vol,Tot_volume,dbh] = architecture_definitions_TABLE(QSM);

    
    %% Detailed measures 
    %n=find_vertical_segments(centres,height);
    %n=10; 
    %PLOT=0;  %top_of_trunk=8;
    %[ rTaper, vTaper, vCoeff, Sail_cylinders, Sail_radius,Volume_profile,Tot_Sail_area]=...
    %find_profiles_and_tapers(radius,cyl_length,comps,centres,height,top_of_trunk,n,PLOT);
    %rTaper=real(rTaper);
    
    %%
    %[CoV, CoV_stem, CoV_crown] = find_CoVs(centres,cyl_volume,branch_order,PLOT); CoV_height=CoV(:,3);
    
    %%
    [max_crown_width,crown_height,crown_area, vol_asym, rel_rad_asym] = calculate_CrownAsymmetry(num_cyls,h,branch_order,x,y,z,cyl_volume,PLOT)
    %AspectRatio=max_crown_width/crown_height;
    %rel_CrownDepth=crown_height./height;
    CrownDensity=canopy_volume/crown_area;
    
    %%
    [mean_angle] = calculate_MeanAngle(running_no,parent_index,x_comp,y_comp,z_comp,PLOT);

    %% This will fail if the extension cylinder definitions have been changed
    ON=1;
    if ON==1
       pf=calculate_PathFraction(ext_index,cyl_length,parent_index,PLOT);
    else     
        pf=nan;   
    end

    %% Create output table
    architectures=table(height, dbh, tree_volume, pf,  max_crown_width, crown_area, rel_rad_asym, mean_angle, ...
          canopy_vol_ratio, CrownDensity, num_cyls,...
        'variablenames',{'TreeHeight','DBH','TotalVolume','PathFraction','MaxCrownWidth','CrownArea','CrownAsym','MeanAngle', ...
       'CanopyVolumeRatio','CrownDensity','num_cyls'});

      
end

