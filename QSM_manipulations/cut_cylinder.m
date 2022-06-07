function [cylinder1] = cut_cylinder(cylinder,cyls_out)
    cylinder1.radius=cylinder.radius(cyls_out,:);
    cylinder1.length=cylinder.length(cyls_out,:);
    cylinder1.start=cylinder.start(cyls_out,:);
    cylinder1.axis=cylinder.axis(cyls_out,:);
    cylinder1.parent=cylinder.parent(cyls_out,:);
    cylinder1.extension=cylinder.extension(cyls_out,:);
    cylinder1.added=cylinder.added(cyls_out,:);
    %cylinder1.UnmodRadius=cylinder.UnmodRadius(cyls_out,:);
    cylinder1.branch=cylinder.branch(cyls_out,:);
    cylinder1.BranchOrder=cylinder.BranchOrder(cyls_out,:);
    cylinder1.PositionInBranch=cylinder.PositionInBranch(cyls_out,:);
end

