function [mean_angle] = calculate_MeanAngle(running_no,parent_index,x_comp,y_comp,z_comp,PLOT)
% Mean insertion angle ---------------
    rows1=find(running_no==1);
    parents=parent_index(rows1);
    angle=zeros(length(rows1),1);
    for i=1:length(rows1)
        if parents(i)==0 continue
        end
        xyz(i,[1:3])=[x_comp(rows1(i)) y_comp(rows1(i)) z_comp(rows1(i))] ;
        xyz_par(i,[1:3])=[x_comp(parents(i)) y_comp(parents(i)) z_comp(parents(i))] ;
        angle(i)=asind(dot(xyz(i,:),xyz_par(i,:))./(norm(xyz(i,:))*norm(xyz_par(i,:))));
    end
    if PLOT==1
        hist(angle)
        title('Branching Angles')
        pause
    end
    mean_angle=mean(abs(angle));
end

