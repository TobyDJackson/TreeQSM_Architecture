function [ architectures Sail_area_profile ] = Calculate_architectures_STRUCT( QSM,PLOT )
    

    
    %%  Definitions and simple measures
    [num_cyls,radius,cyl_length,x,y,z,x_comp,y_comp,z_comp,comps,centres,h,parent_index,ext_index,...
    branch_id,branch_order,running_no,height,cyl_volume,tree_volume,canopy_volume,...
    canopy_vol_ratio,h_vol,Tot_volume,dbh] = architecture_definitions(QSM);

    %% Detailed measures    
     n=10; PLOT=1;  top_of_trunk=8;
    [ rTaper vTaper vCoeff Sail_cylinders, Sail_radius,Volume_profile,Tot_Sail_area] = ...
    find_profiles_and_tapers(radius,cyl_length,comps,centres,height,top_of_trunk,0);

     %%
    [CoV, CoV_stem, CoV_crown] = find_CoVs(centres,cyl_volume,branch_order,PLOT);
    %%
     [max_crown_width,crown_height,crown_area, vol_asym, rel_rad_asym] = calculate_CrownAsymmetry(num_cyls,h,branch_order,x,y,z,cyl_volume,PLOT);
     %%
     [mean_angle] = calculate_MeanAngle(running_no,parent_index,x_comp,y_comp,z_comp,PLOT);

    %% This will fail if the extension cylinder definitions have been changed
     ON=1;
     if ON==1
        [pf] = calculate_PathFraction(ext_index,cyl_length,parent_index,PLOT);
     else;     pf=nan;   end


  
   %% Create output table
   architectures=table(height, dbh, tree_volume, pf, crown_height, max_crown_width,crown_area, rel_rad_asym, mean_angle, ...
                                    Tot_Sail_area,  canopy_vol_ratio, CoV_height,num_cyls,'variablenames',...
                                 {'TreeHeight', 'DBH', 'TotalVolume', 'PathFraction',  'CrownHeight', 'MaxCrownWidth','CrownArea','RelRadAsym', ...
                                 'MeanAngle',  'TotalSailArea' ,'CanopyVolumeRatio' ,'CoV_height','num_cyls'});

       

end

