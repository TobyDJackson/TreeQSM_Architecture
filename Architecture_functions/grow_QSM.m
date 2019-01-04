function [new_qsm] = grow_QSM(QSM,Hgrow,fit)
% the current parameterization for this is in FinalTableTests

%%
load('C:\Users\Toby\Dropbox\Tobys_Stuff\MATLAB\QSM_Simulations\Danum_PCs\QSMs\Unsimplified_BO\A1n.mat')
%load('C:\Users\Toby\Dropbox\Tobys_Stuff\MATLAB\QSM_Simulations\Danum_PCs\QSMs\Simplified_buttress_BO\A1n_0.035_1_buttress.mat')
%QSM=QSM_buttress;
%%
a=1.2;
c=QSM.cylinder;
new_start=nan(size(c.start));
new_start(1,:)=c.start(1,:);
for i=2:length(c.radius)
    new_start(i,:)=new_start(c.parent(i),:)+Hgrow.*c.axis(c.parent(i),:).*c.length(c.parent(i));
end
new_len=Hgrow.*c.length;
new_rad=Rgrow.*c.radius;
new_qsm=QSM;
new_qsm.cylinder.length=new_len;
new_qsm.cylinder.radius=new_rad;
new_qsm.cylinder.start=new_start;
%%
if 1==2
    subplot(1,2,1)
    plot(c.start(:,3)); hold on;
    plot(new_start(:,3))
    subplot(1,2,2)
    plot(c.length); hold on;
    plot(new_len)
end
%%
if 1==2
    subplot(1,2,1)
    plot_cone_model(QSM.cylinder,1,20,1)
    subplot(1,2,2)
    plot_cone_model(new_qsm.cylinder,1,20,1)
end
%%
end

