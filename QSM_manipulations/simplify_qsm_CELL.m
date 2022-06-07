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
