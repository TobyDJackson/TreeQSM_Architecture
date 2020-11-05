function [num_cyls,radius,cyl_length,x,y,z,x_comp,y_comp,z_comp,comps,centres,h,parent_index,ext_index,...
    branch_id,branch_order,running_no,height,cyl_volume,tree_volume,canopy_volume,...
    canopy_vol_ratio,h_vol,Tot_volume,dbh] = architecture_definitions_STRUCT(QSM)


    num_cyls=length(QSM.cylinder.radius); 
    radius=QSM.cylinder.radius;
    cyl_length=QSM.cylinder.length;
    x=QSM.cylinder.start(:,1)-QSM.cylinder.start(1,1);
    y=QSM.cylinder.start(:,2)-QSM.cylinder.start(1,2);
    z=QSM.cylinder.start(:,3)-QSM.cylinder.start(1,3);
    x_comp=QSM.cylinder.axis(:,1);
    y_comp=QSM.cylinder.axis(:,2);
    z_comp=QSM.cylinder.axis(:,3);
    comps=QSM.cylinder.axis;
    parent_index=QSM.cylinder.parent;
    ext_index=QSM.cylinder.extension;
    branch_id=QSM.cylinder.branch;
    branch_order=QSM.cylinder.BranchOrder;
    running_no=QSM.cylinder.PositionInBranch;
    height=max(z)-min(z);         
    centres=cat(2,x+0.5*cyl_length.*x_comp,y+0.5*cyl_length.*y_comp,z+0.5*cyl_length.*z_comp);
    cyl_volume=cyl_length.*pi.*(radius.^2);
    h=centres(:,3);
    tree_volume=sum(cyl_volume);
    canopy_volume=sum(cyl_volume(find(branch_order>=1)));
    stem_volume=sum(cyl_volume(find(branch_order==0)));
    %canopy_vol_ratio=tree_volume/canopy_volume;
    canopy_vol_ratio=canopy_volume./stem_volume;
    h_vol=height/canopy_volume; %Height to canopy volume ratio
    Tot_volume=sum(cyl_volume);
    [dbh] = find_dbh(radius,branch_order,h,cyl_length,z_comp);


end

