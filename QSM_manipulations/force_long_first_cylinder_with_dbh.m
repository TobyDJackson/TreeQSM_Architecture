function [CylData1] = force_long_first_cylinder_with_dbh(cylinder_data,h,DBH,PLOT)

% I need to add an if !
limit=5;
temp_5=cylinder_data(1:limit,:);
temp_CylData=cylinder_data;

%% Find point h
al=sqrt(cylinder_data(1:limit,3).^2+cylinder_data(1:limit,4).^2+cylinder_data(1:limit,5).^2); % absolute length
above_h=find(al>h);
temp=cat(2,linspace(cylinder_data(above_h(1)-1,3),cylinder_data(above_h(1),3),100)',...
    linspace(cylinder_data(above_h(1)-1,4),cylinder_data(above_h(1),4),100)',...
    linspace(cylinder_data(above_h(1)-1,5),cylinder_data(above_h(1),5),100)');
al_temp=sqrt(temp(:,1).^2+temp(:,2).^2+temp(:,3).^2);
[val pos]=min(abs(al_temp-h));
al_exact=al_temp(pos);
coord_h=temp(pos,:); %coordinates of the point at which length = h

%% Replace first beam
new_beam_vector=(temp(pos,:)-[0 0 0])./al_temp(pos);
new_beam=[DBH al_temp(pos) 0 0 0 new_beam_vector 0 2 1 0 1 0];
skipped_nodes=above_h(1)-2;

%% Decide if we need to modify beam above and do it (I have to define a beam either way)
next_len=al(above_h(1))-al_temp(pos);

if above_h(1) ~= limit & next_len<2 %then I want to modify - skip the next node
    disp('Updated again')
    next_node=above_h(2);
    skipped_nodes=skipped_nodes+1;
    next_beam_radius=mean([cylinder_data(above_h(1),1) cylinder_data(above_h(2),1)]);
    next_beam_length=sqrt(sum((cylinder_data(next_node,3:5)-temp(pos,1:3)).^2));
    next_beam_vector=(cylinder_data(next_node,3:5)-temp(pos,1:3))./next_beam_length;
    next_beam1=[next_beam_radius  next_beam_length temp(pos,1:3) next_beam_vector...
        cylinder_data(next_node,9)-skipped_nodes cylinder_data(next_node,10)-skipped_nodes cylinder_data(next_node,11:14)]; 
    temp_CylData(:,9)=temp_CylData(:,9)-skipped_nodes+1;
    temp_CylData(:,10)=temp_CylData(:,10)-skipped_nodes+1;
    CylData1=cat(1,new_beam,next_beam1,temp_CylData(next_node:end,:));
    %plot_cylinder_model(CylData1,1,20,1)
else %otherswise just stick with the next node as is - so I need to define a beam from h to next node
    disp('Im going in the else')
    next_node=above_h(1);
    next_beam_radius=cylinder_data(above_h(1)-1,1);
    next_beam_length=sqrt(sum((cylinder_data(next_node,3:5)-temp(pos,1:3)).^2));
    next_beam_vector=(cylinder_data(next_node,3:5)-temp(pos,1:3))./next_beam_length;
    next_beam1=[next_beam_radius  next_beam_length temp(pos,1:3) next_beam_vector...
        cylinder_data(next_node,9)-skipped_nodes cylinder_data(next_node,10)-skipped_nodes cylinder_data(next_node,11:14)];    
    temp_CylData(:,9)=temp_CylData(:,9)-(skipped_nodes-1);
    temp_CylData(:,10)=temp_CylData(:,10)-(skipped_nodes-1);
    CylData1=cat(1,new_beam,next_beam1,temp_CylData(next_node:end,:));
end
if PLOT==1
    subplot(1,2,1)
    plot_cylinder_model(temp_5,1,20,1)
    subplot(1,2,2)
    plot_cylinder_model(CylData1(1:5,:),1,20,1)
    %plot_cylinder_model(CylData1,1,20,1)
    %pause
end

    
end

