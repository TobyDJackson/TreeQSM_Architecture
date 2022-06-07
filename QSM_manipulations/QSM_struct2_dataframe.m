function [CylData] = QSM_struct2_dataframe(QSM)


   C = QSM.cylinder;
   format long
       
    part1=table(C.radius ,C.length,C.start(:,1),C.start(:,2),C.start(:,3), C.axis(:,1),C.axis(:,2),C.axis(:,3),...
    'VariableNames', {'Radius','Length','x_coord','y_coord','z_coord','x_vector','y_vector','z_vector'});
    
    part2=table(single(C.parent), single(C.extension), single(C.branch) ,single(C.BranchOrder), single(C.PositionInBranch), single(C.UnmodRadius), C.added,...
       'VariableNames',{'Parent_cyl','Extension_cyl','Branch_id','Branch_order','Position_in_branch','Unmodified_radius','Added'});
       
     CylData=cat(2,part1,part2);
       
end

