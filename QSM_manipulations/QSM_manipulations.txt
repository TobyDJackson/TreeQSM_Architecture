% This is a collections of functions which manipulate QSMs
% They work mostly with old versions of the QSM (or only the CylData) but some of them transform between versions

function [CylData,BranchData] = QSM_struct2array(QSM)

           C = QSM.cylinder;
           format long
           CylData = cat(2,C.radius ,C.length, C.start, C.axis, single(C.parent), single(C.extension), single(C.branch) ,...
               single(C.BranchOrder), single(C.PositionInBranch), C.added);
           B = QSM.branch;
           BranchData = cat(2,B.order, B.parent, B.volume, B.length, B.angle);

      %CylData = single([Rad Len Sta Axe CPar CExt BoC Added]);
       % BranchData = single([BOrd BPar BVol BLen BAng]);

end

function [CylData] = QSM_struct2_dataframe(QSM)


   C = QSM.cylinder;
   format long
       
    part1=table(C.radius ,C.length,C.start(:,1),C.start(:,2),C.start(:,3), C.axis(:,1),C.axis(:,2),C.axis(:,3),...
    'VariableNames', {'Radius','Length','x_coord','y_coord','z_coord','x_vector','y_vector','z_vector'});
    
    part2=table(single(C.parent), single(C.extension), single(C.branch) ,single(C.BranchOrder), single(C.PositionInBranch), single(C.UnmodRadius), C.added,...
       'VariableNames',{'Parent_cyl','Extension_cyl','Branch_id','Branch_order','Position_in_branch','Unmodified_radius','Added'});
       
     CylData=cat(2,part1,part2);
       
end

function [QSM] = QSM_cell2struct(CylData,BranchData)
    %This function takes the CylData and BranchData that come out of the QSM sometimes and
    %saves them as a struct.
   
    
    Rad=CylData(:,1); Len=CylData(:,2); Sta=CylData(:,3:5); Axe=CylData(:,6:8); CPar=CylData(:,9); CExt=CylData(:,10); BoC=CylData(:,11:13); Added=CylData(:,14); 
    BOrd=BranchData(:,1); BPar=BranchData(:,2); BVol=BranchData(:,3); BLen=BranchData(:,4); BAng=BranchData(:,5); 
    
    cylinder=struct('radius',Rad,'length',Len,'start',Sta,'axis',Axe,'parent',CPar,'extension',CExt,...
         'BranchOrder',BoC(:,2),'branch',BoC(:,1),'PositionInBranch',BoC(:,3),'added',Added);
     
     bases=find(cylinder.PositionInBranch==1);
     BDiam=cylinder.radius(bases)*2;
     
     branch=struct('order',BOrd,'parent',BPar,'volume',BVol,'length',BLen,'angle',BAng,'diameter',BDiam);
     
     QSM=struct('cylinder',cylinder,'branch',branch);
        
end

function [QSM] = QSM_bits2struct(Rad,Len,Sta,Axe,CPar,CExt,BoC,Added,BOrd,BPar,BVol,BLen,BAng)
%This function takes the pieces that come out of the QSM sometimes and
%saves them as a struct.

     cylinder=struct('radius',Rad,'length',Len,'start',Sta,'axis',Axe,'parent',CPar,'extension',CExt,...
         'BranchOrder',BoC(:,2),'branch',BoC(:,1),'PositionInBranch',BoC(:,3),'added',Added);
     
     bases=find(cylinder.PositionInBranch==1);
     BDiam=cylinder.radius(bases)*2;
        branch=struct('order',BOrd,'parent',BPar,'volume',BVol,'length',BLen,'angle',BAng,'diameter',BDiam);
        QSM=struct('cylinder',cylinder,'branch',branch);
        
end



function [CylData, BranchData] = QSM_bits2cell(Rad,Len,Sta,Axe,CPar,CExt,BoC,Added,BOrd,BPar,BVol,BLen,BAng)
%This function takes the pieces that come out of the QSM sometimes and
%saves them as a cell
        CylData = single([Rad Len Sta Axe CPar CExt BoC Added]);  
        BranchData = single([BOrd BPar BVol BLen BAng]);
        %ModelData={CylData, BranchData}; %repackages data in cell array for simplify code   

end


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

function [ CylData_simplified_added ] = add_simplified_area( CylData, CylData_updated )
	% This functions is intended to simulate the effect of leaves by adding area to the outer cylinders

	[ Arch_Values_original Sail Arch_Names ] = Calculate_architectures_2017( CylData );
	Sail_Original=Arch_Values_original(10);
	Tot_Original=sum(CylData(:,1).*CylData(:,2));

	[ Arch_Values_simplified Sail Arch_Names ] = Calculate_architectures_2017( CylData_updated  );
	Sail_simplified=Arch_Values_simplified(10)
	Tot_simplified=sum(CylData_updated(:,1).*CylData_updated(:,2));

	[length(CylData) Sail_Original Tot_Original length(CylData_updated) Sail_simplified Tot_simplified]
	CylData_simplified_added=[length(CylData) Sail_Original Tot_Original length(CylData_updated) Sail_simplified Tot_simplified];

end


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

% This file is part of TREEQSM.
% 
% TREEQSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TREEQSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

function QSM = simplify_qsm_STRUCT(QSM,MaxBranchOrder,SmallRadii,ReplaceIterations,PLOT)

% ---------------------------------------------------------------------
% SIMPLIFY_QSM.M   Simplifies cylinder QSMs by restricting the maximum
%                       branching order, by removing thin branches, and by 
%                       replacing two concecutive cylinders with a longer cylinder
%
% Version 1.10
% Latest update     16 Aug 2017
%
% Copyright (C) 2015-2017 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Inputs:
% QSM               Model file, output of treeqsm.m, must contains only one model
% MaxBranchOrder    Maximum branching order, higher order branches removed
% SmallRadii        Minimum acceptable radius for a branch at its base
% ReplaceIterations Number of iterations for replacing two concecutive
%                       cylinders inside one branch with one longer cylinder  
%
% Output:
% Modified QSM      NOTICE: Only cylinder data is modified.
if PLOT==1
plot_cone_model(QSM.cylinder,1,20,1)
end
%pause(0.1)

%% Maximum branching order
c = QSM.cylinder;
nc = size(c.radius,1);
disp([num2str(nc),' cylinders originally'])
% Cylinders with branching order up to MaxBranchOrder
SmallOrder = c.BranchOrder <= MaxBranchOrder; 
N = fieldnames(c);
n = max(size(N));
for i = 1:n
    c.(N{i}) = c.(N{i})(SmallOrder,:);
end

%nc = nnz(SmallOrder);
%disp([num2str(nc),' cylinders after branching order simplification'])

%% Small branches
if nargin >= 3 && SmallRadii > 0
    % Determine child cylinders
    nc = size(c.radius,1);
%     CChi = cell(nc,1);
%     for i = 1:nc
%         P = QSM.cylinder.parent(i);
%         if P > 0
%             PE = c.extension(P);
%             if PE ~= i
%                 CChi{P} = [CChi{P}; i];
%             end
%         end
%     end
    % Determine child branches
    BPar = QSM.branch.parent;
    nb = size(BPar,1);
    BChi = cell(nb,1);
    for i = 1:nb
        P = BPar(i);
        if P > 0
            BChi{P} = [BChi{P}; i];
        end
    end
    
    % Remove branches whose radii is too small compared to its parent
    Keep = true(nc,1);
    Pass = true(nb,1);
    for i = 1:nb
        if Pass(i)
            if QSM.branch.diameter(i) < SmallRadii
                B = i;
                BC = BChi{B};
                while ~isempty(BC)
                    B = [B; BC];
                    BC = vertcat(BChi{BC});
                end
                Pass(B) = false;
                m = length(B);
                for k = 1:m
                    Keep(c.branch == B(k)) = false;
                end
            end
        end
    end

    % Modify topology information
    Ind = (1:1:nc)';
    m = nnz(Keep);
    Ind(Keep) = (1:1:m)';
    I = c.parent > 0;
    c.parent(I) = Ind(c.parent(I));
    I = c.extension > 0;
    c.extension(I) = Ind(c.extension(I));
        
    % Update/reduce cylinders
    for i = 1:n
        c.(N{i}) = c.(N{i})(Keep,:);
    end

    %nc = nnz(Keep);
    %disp([num2str(nc),' cylinders after small branch simplification'])
end
    
%% Cylinder replacing
if nargin == 4
    
    % Determine child cylinders
    nc = size(c.radius,1);
    CChi = cell(nc,1);
    for i = 1:nc
        P = c.parent(i);
        if P > 0
            PE = c.extension(P);
            if PE ~= i
                CChi{P} = [CChi{P}; i];
            end
        end
    end
    
    % Replace cylinders
    for j = 1:ReplaceIterations
        
        nc = size(c.radius,1);
        Ind = (1:1:nc)';
        Keep = false(nc,1);
        i = 1;
        while i <= nc
            t = 1;
            while i+t <= nc && c.branch(i+t) == c.branch(i)
                t = t+1;
            end
            Cyls = (i:1:i+t-1)';
            S = c.start(Cyls,:);
            A = c.axis(Cyls,:);
            L = c.length(Cyls);
            if t == 1 % one cylinder in the branch
                Keep(i) = true;
            elseif ceil(t/2) == floor(t/2) % even number of cylinders in the branch
                I = (1:2:t)'; % select 1., 3., 5., ...
                % Correct radii, axes and lengths
                E = S(end,:)+L(end)*A(end,:);
                S = S(I,:);
                m = length(I);
                if m > 1
                    A = [S(2:end,:); E]-S(1:end,:);
                else
                    A = E-S(1,:);
                end
                L = sqrt(sum(A.*A,2));
                A = [A(:,1)./L A(:,2)./L A(:,3)./L];
                cyls = Cyls(I);
                Keep(cyls) = true;
                V = pi*c.radius(Cyls).^2.*c.length(Cyls);
                J = (2:2:t)';
                V = V(I)+V(J);
                R = sqrt(V./L/pi);
                c.radius(cyls) = R;
                
            else % odd number of cylinders
                I = [1 2:2:t]'; % select 1., 2., 4., 6., ...
                % Correct radii, axes and lengths
                E = S(end,:)+L(end)*A(end,:);
                S = S(I,:);
                l = L(1);
                a = A(I,:);
                m = length(I);
                if m > 2
                    a(2:end,:) = [S(3:end,:); E]-S(2:end,:);
                else
                    a(2,:) = E-S(2,:);
                end
                A = a;
                L = sqrt(sum(A.*A,2));
                L(1) = l;
                A(2:end,:) = [A(2:end,1)./L(2:end) A(2:end,2)./L(2:end) A(2:end,3)./L(2:end)];
                cyls = Cyls(I);
                Keep(cyls) = true;
                V = pi*c.radius(Cyls).^2.*c.length(Cyls);
                J = (3:2:t)';
                V = V(I(2:end))+V(J);
                R = sqrt(V./L(2:end)/pi);
                c.radius(cyls(2:end)) = R;
            end
            
            if t > 1
                % Modify cylinders
                c.length(cyls) = L;
                c.axis(cyls,:) = A;
                % Correct branching/topology information
                c.PositionInBranch(cyls) = (1:1:m)';
                c.extension(cyls) = [cyls(2:end); 0];
                c.parent(cyls(2:end)) = cyls(1:end-1);
                par = c.parent(cyls(1));
                if par > 0 && ~Keep(par)
                    par0 = c.parent(par);
                    if Keep(par0) && c.extension(par0) == par
                        c.parent(cyls(1)) = par0;
                    end
                end
                
                % Correct child branches
                chi = vertcat(CChi{Cyls});
                if ~isempty(chi)
                    par = c.parent(chi);
                    J = Keep(par);
                    par = par(~J)-1;
                    c.parent(chi(~J)) = par;
                    
                    par = c.parent(chi);
                    rp = c.radius(par);
                    sp = c.start(par,:);
                    ap = c.axis(par,:);
                    lc = c.length(chi);
                    sc = c.start(chi,:);
                    ac = c.axis(chi,:);
                    ec = sc+[lc.*ac(:,1) lc.*ac(:,2) lc.*ac(:,3)];
                    m = length(chi);
                    for k = 1:m
                        [d,V,h,B] = distances_to_line(sc(k,:),ap(k,:),sp(k,:));
                        V = V/d;
                        sc(k,:) = sp(k,:)+rp(k)*V+B;
                    end
                    ac = ec-sc;
                    [ac,lc] = normalize(ac);
                    c.length(chi) = lc;
                    c.start(chi,:) = sc;
                    c.axis(chi,:) = ac;
                end 
            end
            
            i = i+t;
        end
        % Change topology (parent, extension) indexes
        m = nnz(Keep);
        Ind(Keep) = (1:1:m)';
        I = c.parent > 0;
        c.parent(I) = Ind(c.parent(I));
        I = c.extension > 0;
        c.extension(I) = Ind(c.extension(I));
        
        % Update/reduce cylinders
        for i = 1:n
            c.(N{i}) = c.(N{i})(Keep,:);
        end
        
        if j < ReplaceIterations
            % Determine child cylinders
            nc = size(c.radius,1);
            CChi = cell(nc,1);
            for i = 1:nc
                P = c.parent(i);
                if P > 0
                    PE = c.extension(P);
                    if PE ~= i
                        CChi{P} = [CChi{P}; i];
                    end
                end
            end
        end
    end
end

QSM.cylinder = c;

nc = size(c.radius,1);
disp([num2str(nc),' cylinders after simplification'])
if PLOT==1
    plot_cone_model(QSM.cylinder,2,20,1)
end
%pause

function ModelData = simplify_qsm_CELL(ModelData,MaxBranchOrder,SmallRadii,ReplaceIterations)
disp('Im in here now')
% ---------------------------------------------------------------------
% SIMPLIFY_MODEL.M   Simplifies cylinder model by destricting the maximum
%                       branching order, removing branches that have small
%                       radii compared to the radii of the parent branch,
%                       and by replacing two concecutive cylinders with 
%                       one longer cylinder
%
% Version 1.1
% Author        Pasi Raumonen
% Created       12 Nov 2015
%
% The distribution and commercial use of this software is prohibited
% without the permission of Pasi Raumonen. 
%
% Copyright:    Pasi Raumonen
% ---------------------------------------------------------------------
%
% Inputs:
% ModelData         ModelData file, output of qsm_tree.m
% MaxBranchOrder    Maximum branching order, higher order branches removed
%% TOBY HAS UPDATED THIS to remove radii under a given radius limit, not ratio to parent.
% SmallRadii        Relative radius for removal: 
%                       if R_child < SmallRadii*R_parent, then remove
%                       branch and all its child branches
% ReplaceIterations Number of iterations for replacing two concecutive
%                       cylinders inside one branch with one longer cylinder  
%
% Output:
% Modified ModelData-file. NOTICE: Only cylinder data CylData is modified



%% Maximum branching order
CylData = ModelData{1}; % Cylinder data from ModelData
nc = size(CylData,1);
disp([num2str(nc),' cylinders originally'])

%plot_cylinder_model(CylData,1,20,1)
%pause

SmallOrder = CylData(:,12) <= MaxBranchOrder; % Select cylinders from branches with branching order upto MaxBranchOrder
CylData = CylData(SmallOrder,:);
nc = size(CylData,1);
disp([num2str(nc),' cylinders after branching order simplification'])


%% Small branches
if nargin >= 3 && SmallRadii > 0
    % Determine child cylinders
    %toby_count=1
    nc = size(CylData,1);
    CChi = cell(nc,1);
    for i = 1:nc
        P = CylData(i,9);
        if P > 0
            PE = CylData(P,10);
            if PE ~= i
                CChi{P} = [CChi{P}; i];
            end
        end
    end
    % Determine child branches
    BranchData = ModelData{2};
    nb = size(BranchData,1);
    BChi = cell(nb,1);
    for i = 1:nb
        P = BranchData(i,2);
        if P > 0
            BChi{P} = [BChi{P}; i];
        end
    end
    
    % Remove branches whose radii is too small compared to its parent
    Keep = true(nc,1);
    for i = 1:nc
        if Keep(i)
            C = CChi{i};
            I = CylData(C,1) < SmallRadii;
            if any(I)
                C = C(I);
                m = length(C);
                for j = 1:m
                    B = CylData(C(j),11);
                    BC = BChi{B};
                    while ~isempty(BC)
                        B = [B; BC];
                        BC = vertcat(BChi{BC});
                    end
                    nb = length(B);
                    for k = 1:nb
                        J = CylData(:,11) == B(k);
                        Keep(J) = false;
                    end
                end
            end
        end
    end
    % Modify topology information
    Ind = (1:1:nc)';
    m = nnz(Keep); %Total number to keep. Keep is a logical deleting vector
    Ind(Keep) = (1:1:m)'; %Change index so the new ones only count up to new total.
    I = CylData(:,9) > 0; %Logical for all that don't come out of thin air! All 1's for a good tree.
    CylData(I,9) = Ind(CylData(I,9));
    I = CylData(:,10) > 0; %Logical for all that have an extension cylinder 
    %% Toby added a bit here
    %if max(CylData(I,10))>length(Ind) %Then we are going to get an error message
    %    CylData(find(CylData(:,10)==max(CylData(:,10))),10)=length(Ind);
    %end
    CylData(I,10) = Ind(CylData(I,10));
        
    CylData = CylData(Keep,:);
    nc = size(CylData,1);
    disp([num2str(nc),' cylinders after small branch simplification'])
    %toby=10
end
    
%% Cylinder replacing
if nargin >= 4
    % Determine child cylinders
    nc = size(CylData,1);
    CChi = cell(nc,1);
    for i = 1:nc
        P = CylData(i,9);
        if P > 0
            PE = CylData(P,10);
            if PE ~= i
                CChi{P} = [CChi{P}; i];
            end
        end
    end

    
    
    %% Replace cylinders
    for j = 1:ReplaceIterations
        nc = size(CylData,1);
        Ind = (1:1:nc)';
        Keep = false(nc,1);
        i = 1;
        while i <= nc
            t = 1;
            while i+t <= nc && CylData(i+t,11) == CylData(i,11)
                t = t+1;
            end
            Cyls = (i:1:i+t-1)';
            S = CylData(Cyls,3:5);
            A = CylData(Cyls,6:8);
            L = CylData(Cyls,2);
            if t == 1 % one cylinder in the branch
                Keep(i) = true;
            elseif ceil(t/2) == floor(t/2) % even number of cylinders in the branch
                I = (1:2:t)'; % select 1., 3., 5., ...
                % Correct radii, axes and lengths
                E = S(end,:)+L(end)*A(end,:);
                S = S(I,:);
                if length(I) > 1
                    A = [S(2:end,:); E]-S(1:end,:);
                else
                    A = E-S(1,:);
                end
                L = sqrt(sum(A.*A,2));
                A = [A(:,1)./L A(:,2)./L A(:,3)./L];
                cyls = Cyls(I);
                m = length(I);
                Keep(cyls) = true;
                V = pi*CylData(Cyls,1).^2.*CylData(Cyls,2);
                J = (2:2:t)';
                V = V(I)+V(J);
                R = sqrt(V./L/pi);
                CylData(cyls,1) = R;
                CylData(cyls,2) = L;
                CylData(cyls,6:8) = A;
                % Correct branching/topology information
                CylData(cyls,13) = (1:1:m)';
                CylData(cyls,10) = [cyls(2:end); 0];
                CylData(cyls(2:end),9) = cyls(1:end-1);
                par = CylData(cyls(1),9);
                if par > 0 && ~Keep(par)
                    par0 = CylData(par,9);
                    if Keep(par0) && CylData(par0,10) == par
                        CylData(cyls(1),9) = par0;
                    end
                end
                CylData(cyls,13) = (1:1:m)';
                
                % Correct child branches
                chi = vertcat(CChi{Cyls});
                if ~isempty(chi)
                    par = CylData(chi,9);
                    J = Keep(par);
                    par = par(~J)-1;
                    CylData(chi(~J),9) = par;
                    
                    par = CylData(chi,9);
                    rp = CylData(par,1);
                    sp = CylData(par,3:5);
                    ap = CylData(par,6:8);
                    lc = CylData(chi,2);
                    sc = CylData(chi,3:5);
                    ac = CylData(chi,6:8);
                    ec = sc+[lc.*ac(:,1) lc.*ac(:,2) lc.*ac(:,3)];
                    m = length(chi);
                    for k = 1:m
                        [d,V,h,B] = distances_to_line(sc(k,:),ap(k,:),sp(k,:));
                        V = V/d;
                        sc(k,:) = sp(k,:)+rp(k)*V+B;
                    end
                    ac = ec-sc;
                    [ac,lc] = normalize(ac);
                    CylData(chi,2) = lc;
                    CylData(chi,3:5) = sc;
                    CylData(chi,6:8) = ac;
                end
                
            else % odd number of cylinders
                I = [1 2:2:t]'; % select 1., 2., 4., 6., ...
                % Correct radii, axes and lengths
                E = S(end,:)+L(end)*A(end,:);
                S = S(I,:);
                l = L(1);
                a = S;
                a(1,:) = A(1,:);
                if length(I) > 2
                    a(2:end,:) = [S(3:end,:); E]-S(2:end,:);
                else
                    a(2,:) = E-S(2,:);
                end
                A = a;
                L = sqrt(sum(A.*A,2));
                L(1) = l;
                A = [A(:,1)./L A(:,2)./L A(:,3)./L];
                cyls = Cyls(I);
                m = length(I);
                Keep(cyls) = true;
                V = pi*CylData(Cyls,1).^2.*CylData(Cyls,2);
                J = (3:2:t)';
                V = V(I(2:end))+V(J);
                R = sqrt(V./L(2:end)/pi);
                CylData(cyls(2:end),1) = R;
                CylData(cyls,2) = L;
                CylData(cyls,6:8) = A;
                % Correct branching/topology information
                CylData(cyls,13) = (1:1:m)';
                CylData(cyls,10) = [cyls(2:end); 0];
                CylData(cyls(2:end),9) = cyls(1:end-1);
                par = CylData(cyls(1),9);
                if par > 0 && ~Keep(par)
                    par0 = CylData(par,9);
                    if Keep(par0) && CylData(par0,10) == par
                        CylData(cyls(1),9) = par0;
                    end
                end
                CylData(cyls,13) = (1:1:m)';
                
                % Correct child branches
                chi = vertcat(CChi{Cyls});
                if ~isempty(chi)
                    par = CylData(chi,9);
                    J = Keep(par);
                    par = par(~J)-1;
                    CylData(chi(~J),9) = par;
                    
                    par = CylData(chi,9);
                    rp = CylData(par,1);
                    sp = CylData(par,3:5);
                    ap = CylData(par,6:8);
                    lc = CylData(chi,2);
                    sc = CylData(chi,3:5);
                    ac = CylData(chi,6:8);
                    ec = sc+[lc.*ac(:,1) lc.*ac(:,2) lc.*ac(:,3)];
                    m = length(chi);
                    for k = 1:m
                        [d,V,h,B] = distances_to_line(sc(k,:),ap(k,:),sp(k,:));
                        V = V/d;
                        sc(k,:) = sp(k,:)+rp(k)*V+B;
                    end
                    ac = ec-sc;
                    [ac,lc] = normalize(ac);
                    CylData(chi,2) = lc;
                    CylData(chi,3:5) = sc;
                    CylData(chi,6:8) = ac;
                end
            end
            i = i+t;
        end
        % Change topology (parent, extension) indexes
        m = nnz(Keep);
        Ind(Keep) = (1:1:m)';
        I = CylData(:,9) > 0;
        CylData(I,9) = Ind(CylData(I,9));
        I = CylData(:,10) > 0;
        CylData(I,10) = Ind(CylData(I,10));
        
        % Update/reduce CylData
        CylData = CylData(Keep,:);
        
        if j < ReplaceIterations
            % Redetermine child cylinders
            nc = size(CylData,1);
            CChi = cell(nc,1);
            for i = 1:nc
                P = CylData(i,9);
                if P > 0
                    PE = CylData(P,10);
                    if PE ~= i
                        CChi{P} = [CChi{P}; i];
                    end
                end
            end
        end
    end
end
ModelData{1} = CylData;

nc = size(CylData,1);
disp([num2str(nc),' cylinders after final simplication'])

%fig_to_save=figure
%plot_cylinder_model(CylData,2,20,1)
