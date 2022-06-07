function Plot_labelled_QSM(CylData,label_these_beams,color_these,fig_no)

    len=length(CylData);
    xyz_pos=CylData(:,3:5)-repmat(CylData(1,3:5),len,1);
    CylData(:,3:5)=xyz_pos; %This centres the tree on zero.
    xyz_vec=CylData(:,6:8);
    %Vector for half way along beam
    location=double(xyz_pos+0.5*([CylData(:,2).*xyz_vec(:,1) CylData(:,2).*xyz_vec(:,2) CylData(:,2).*xyz_vec(:,3)]) );
    x_offset=double((ones(len,1)-xyz_vec(:,1)).*2.*CylData(:,1));
    y_offset=double((ones(len,1)-xyz_vec(:,2)).*2.*CylData(:,1));
    z_offset=double((ones(len,1)-xyz_vec(:,3)).*2.*CylData(:,1));

    selection=zeros(len,1); 
    for i=1:length(color_these)
        selection(color_these(i))=1; 
    end
    plot_cylinder_model_OLD(CylData,fig_no,20,0.2,selection)   

    
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



