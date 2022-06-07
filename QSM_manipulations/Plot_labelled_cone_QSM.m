function Plot_labelled_cone_QSM(cylinder,label_these_beams,fig_no)

    len=length(cylinder.radius);
    radius=cylinder.radius;
    length_cyl=cylinder.length;
    cylinder.start=cylinder.start-cylinder.start(1,:);
    xyz_pos=cylinder.start;
    xyz_vec=cylinder.axis;
    %Vector for half way along beam
    location=double(xyz_pos+0.5*([length_cyl.*xyz_vec(:,1) length_cyl.*xyz_vec(:,2) length_cyl.*xyz_vec(:,3)]) );
    x_offset=double((ones(len,1)-xyz_vec(:,1)).*2.*radius);
    y_offset=double((ones(len,1)-xyz_vec(:,2)).*2.*radius);
    z_offset=double((ones(len,1)-xyz_vec(:,3)).*2.*radius);

    %selection=zeros(len,1); 
    %for i=1:length(color_these)
    %    selection(color_these(i))=1; 
    %end
    %plot_cylinder_model_OLD(cylinder,fig_no,20,0.2,selection)   
    plot_cone_model(cylinder,fig_no,20,1)   
    
     %these labels are by row number, so it starts at 1 as in original qsm,
     %not at 2 as in abaqus. I need to add one
     for i=1:length(label_these_beams)
        label=num2str((label_these_beams(i))');
        xpos=location(label_these_beams(i),1)+x_offset(label_these_beams(i));
        ypos=location(label_these_beams(i),2);
        zpos=location(label_these_beams(i),3)+z_offset(label_these_beams(i));
        
        text(xpos,ypos,zpos,label)
     end
            
end



