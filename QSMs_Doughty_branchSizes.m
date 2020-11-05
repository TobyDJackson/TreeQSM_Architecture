%% This script loops over QSMs and calculates architectural measures
% It then does a PCA and plots it out nicely


%% Loop over QSMs and calculate architectures
path='C:\Users\tobyd\Documents\Data\TLS_Data\QSMs\Global_QSMs\All\';
files=dir(strcat(path,'*.mat'));
num_files=length(files);

% Set up empty arrays for vertical-branch-radii (vbr) etc
vbr_min=nan(num_files,10); vbr_max=nan(num_files,10); vbr_median=nan(num_files,10); vbr_mean=nan(num_files,10);
hbr_min=nan(num_files,10); hbr_max=nan(num_files,10); hbr_median=nan(num_files,10); hbr_mean=nan(num_files,10);
height=nan(num_files,1); crown_radius=nan(num_files,1);

% Load the first tree
load(strcat(path,files(1).name) );
%architectures=Calculate_architectures(CylData,0);
plot_cylinder_model(CylData,1,20,1) %Just an example of what the QSM looks like.
SaveNames={files(1).name(1:end-4)};
SaveNames2={files(1).name(1:3)};

for i=2:num_files % Loop over the trees
    i
    load(strcat(path,files(i).name) );
    if length(CylData)<50; continue; end % skip QSMs with less than 50 cylinders because they produce errors.
    %temp_arch=Calculate_architectures(CylData,0); %These are the tree level architectural indices
    %architectures=cat(1,architectures,temp_arch);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bd = QSM_branch_dimensions(CylData); % This is the key function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Just save the outputs into arrays
    vbr_min(i,:)=bd.vbr_min; vbr_max(i,:)=bd.vbr_max;
    vbr_median(i,:)=bd.vbr_median; vbr_mean(i,:)=bd.vbr_mean;
    hbr_min(i,:)=bd.hbr_min; hbr_max(i,:)=bd.hbr_max;
    hbr_median(i,:)=bd.hbr_median; hbr_mean(i,:)=bd.hbr_mean;
    height(i)=bd.height; crown_radius(i)=bd.crown_radius;
    
    SaveNames=cat(1,SaveNames,{files(i).name(1:end-4)});
    SaveNames2=cat(1,SaveNames2,{files(i).name(1:3)});
end
%architectures=cat(2,cell2table(SaveNames),architectures);

%% PLOT the results
dark2=brewermap(11,'dark2');
sn=categorical(SaveNames2);
snn=grp2idx(sn);
%%
for i=1:629%1496
    
    num_in_group=length(find(snn==snn(i)));
    subplot(2,2,1)
    plot((vbr_max(i,:)),height(i).*(0.1:0.1:1),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Max branch radius (m)'); ylabel('Height (m)')
    subplot(2,2,2)
    plot(vbr_max(i,:),10:10:100,'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Max branch radius'); ylabel('Height %')
    subplot(2,2,3)
    plot(vbr_max(i,:)./vbr_max(i,1),10:10:100,'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Max branch radius / stem radius'); ylabel('Height %')
    subplot(2,2,4)
    plot(vbr_max(i,:)./(height(i)),10:10:100,'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Max branch radius / tree height'); ylabel('Height %')
end

 %% Slightly different plots
for i=1:629%1496
    num_in_group=length(find(snn==snn(i)));
    subplot(2,2,1)
    plot((vbr_max(i,:)),(0.1:0.1:1),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    %plot(0:
    xlabel('Max branch radius (m)'); ylabel('Height (m)')
    subplot(2,2,2)
    plot((vbr_median(i,:)),(0.1:0.1:1),'.','color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Median branch radius (m)'); ylabel('Height (m)')
    subplot(2,2,3)
    plot((vbr_max(i,:)),(0.1:0.1:1),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Max branch radius (m)'); ylabel('Height (m)')
    xlim([0 0.4])
    subplot(2,2,4)
    plot((vbr_median(i,:)),(0.1:0.1:1),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    xlabel('Median branch radius (m)'); ylabel('Height (m)')
    xlim([0 0.1])
end

%%
for i=1:629%1496
    num_in_group=length(find(snn==snn(i)));
    subplot(2,2,1)
    plot(crown_radius(i)*(0.1:0.1:1),(hbr_max(i,:)),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    ylabel('Max branch radius (m)'); xlabel('Crown radius (m)')
    subplot(2,2,2)
    plot(crown_radius(i)*(0.1:0.1:1),(hbr_median(i,:)),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    ylabel('Median branch radius (m)'); xlabel('Crown radius (m)')
    subplot(2,2,3)
    plot(crown_radius(i)*(0.1:0.1:1),(hbr_max(i,:)),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    ylabel('Max branch radius (m)'); xlabel('Crown radius (m)')
    ylim([0 0.4])
    subplot(2,2,4)
    plot(crown_radius(i)*(0.1:0.1:1),(hbr_median(i,:)),'color',cat(2,dark2(snn(i),:),15/num_in_group)); hold on
    ylabel('Median branch radius (m)'); xlabel('Crown radius (m)')
    ylim([0 0.1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Force needed to break branches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

MOR=50e6; %Modulus of rupture in Pascal (50-150e6)
% MOR = 0.12*WD - 7.14, where WD is wood density in kg/m3
branch_radius=0.1:-0.01:0.01; % Branch radius in m
num_files=4; % Branch length in m
branch_angle=0; % Branch angle from the horizontal
% Breaking condition
mass_break=(pi.*MOR.*branch_radius.^3)/(4*9.81.*num_files.*cos(branch_angle)); 
% Mass (kg) of an animal required to break the branch
mass_break=@(branch_radius,position_along_branch,MOR)(pi.*MOR.*branch_radius.^3)/(4*9.81.*position_along_branch.*cos(branch_angle)); 
% Mass (kg) of an animal required to break the branch

%%
branching_angle=0;
mass_break=@(branch_radius,branch_length,MOR)(pi.*MOR.*branch_radius.^3)./(4*9.81.*branch_length.*cos(branching_angle)); 
% Mass (kg) of an animal required to break the branch
branching_angle=0;
subplot(1,3,1)
branch_radius_in=0.01:0.001:0.05;
plot(branch_radius_in,mass_break(branch_radius_in,4,50e6),'color','black','linewidth',1.2); 
ylabel('Mass to break branch (kg)'); xlabel('Branch radius (m)'); title('length=4m, MOR=50MPa');
subplot(1,3,2)
branch_length_in=0.5:0.1:10;
plot(branch_length_in,mass_break(0.02,branch_length_in,50e6),'color','black','linewidth',1.2); 
xlabel('Branch length (m)'); title('radius=0.02m, MOR=50MPa');
subplot(1,3,3)
MOR_in=50e6:1e6:150e6;
plot(MOR_in,mass_break(0.02,4,MOR_in),'color','black','linewidth',1.2); 
xlabel('Modulus of rupture (Pa)'); title('radius=0.02m, position=4m');
