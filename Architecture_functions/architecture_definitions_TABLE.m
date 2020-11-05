function [num_cyls,radius,cyl_length,x,y,z,x_comp,y_comp,z_comp,comps,centres,h,parent_index,ext_index,...
    branch_id,branch_order,running_no,height,cyl_volume,tree_volume,canopy_volume,...
    canopy_vol_ratio,h_vol,Tot_volume,dbh] = architecture_definitions_TABLE(QSM)


    num_cyls=length(QSM(:,1)); 
    radius=QSM(:,1);
    cyl_length=QSM(:,1);
    x=QSM(:,3)-QSM(1,3);
    y=QSM(:,4)-QSM(1,4);
    z=QSM(:,5)-QSM(1,5);
    x_comp=QSM(:,6);
    y_comp=QSM(:,7);
    z_comp=QSM(:,8);
    comps=QSM(:,6:8);
    parent_index=QSM(:,9);
    ext_index=QSM(:,10);
    branch_id=QSM(:,11);
    branch_order=QSM(:,12);
    running_no=QSM(:,13);
    height=max(z)-min(z);         
    centres=cat(2,x+0.5*cyl_length.*x_comp,y+0.5*cyl_length.*y_comp,z+0.5*cyl_length.*z_comp);
    cyl_volume=cyl_length.*pi.*(radius.^2);
    h=centres(:,3);
    tree_volume=sum(cyl_volume);
    canopy_volume=sum(cyl_volume(find(branch_order>=1)));
    stem_volume=sum(cyl_volume(find(branch_order==0)));
    %canopy_vol_ratio=tree_volume/canopy_volume;
    canopy_vol_ratio=canopy_volume./stem_volume; % this is the same as K
    h_vol=height/canopy_volume; %Height to canopy volume ratio
    Tot_volume=sum(cyl_volume);
    [dbh] = find_dbh(radius,branch_order,h,cyl_length,z_comp);


end

