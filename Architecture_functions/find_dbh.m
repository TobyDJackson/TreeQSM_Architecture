function [dbh] = find_dbh(radius,branch_order,h,cyl_length,z_comp)
%Finds the diameter at approx 1.3m - dbh.
    r=radius(find(branch_order==0 & h<=1.3 & h+cyl_length.*z_comp>=1.3)); %This finds the radius at 1.3m
    if length(r)>=2 
        r=r(1); %Just choose the first - the lowest
    end
    if length(r)==0 %This defaults to the first beam if the above case was empty.
        r=radius(1);
    end
    dbh=2*r; 
end

