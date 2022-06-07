function [CylData_new] = add_buttress_cylinder(CylData,h,radius,PLOT)

    %1.       radius
    %2.       length
    %3.       x-coordinate of the starting point
    %4.       y-coordinate of the starting point
    %5.       z-coordinate of the starting point
    %6.       x-component of the cylinder axis
    %7.       y-component of the cylinder axis
    %8.       z-component of the cylinder axis

    CylData_new=cat(1,nan(1,size(CylData,2)),CylData);
    CylData_new(1,1)=radius;
    CylData_new(1,2)=h;
    CylData_new(2:end,5)=CylData_new(2:end,5)+h;
    CylData_new(1,3:5)=[0 0 0];
    CylData_new(1,6:8)=CylData(1,6:8);

    %9.      index (row number in this file) of the parent cylinder
    %10.    index (row number in this file) of the extension cylinder
    %11.    branch (row number in the branch data-file) of the cylinder
    %12.    branch order of the branch the cylinder belongs
    %14.    Indication if the cylinder is added after normal cylinder fitting (=1 if added)
    CylData_new(1,9)=0;
    CylData_new(2:end,9)=CylData_new(2:end,9)+1;
    CylData_new(1,10)=0;
    CylData_new(2:end,10)=CylData_new(2:end,10)+1;
    CylData_new(1,11)=CylData_new(2,11);
    CylData_new(1,12)=CylData_new(2,12);
    CylData_new(1,14)=CylData_new(2,14);

    %13.    running number of the cylinder in the branch it belongs
    CylData_new(1,13)=0;
    branch1_rows=find(CylData_new(:,11)==1)
    CylData_new(branch1_rows,13)=CylData_new(branch1_rows,13)+1;
    if PLOT==1
        subplot(1,2,1)
        plot_cylinder_model(CylData,1,20,1)
        subplot(1,2,2)
        plot_cylinder_model(CylData_new,1,20,1)
    end   
end

