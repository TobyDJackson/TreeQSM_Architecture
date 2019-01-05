function [bo] = recalculate_branching_orders(c,limit_ratio, PLOT)

%9.       index (row number in this file) of the parent cylinder
%10.    index (row number in this file) of the extension cylinder
%11.    branch (row number in the branch data-file) of the cylinder
%12.    branch order of the branch the cylinder belongs
%13.    running number of the cylinder in the branch it belongs
twigs=find(c(:,10)==0);
num_twigs=length(twigs);
bases=find(c(:,13)==1); %bottoms of the branches
bp=unique(c(bases(2:end),9)); % these are the branching points - the cylinders which hold the branchings
bo=zeros(length(c),1); %define all bo as zero

if PLOT==1
    for i=1:bp(1)
        label=num2str((bo(i))');        
        hold on
        text(double(c(i,3)),double(c(i,4)),double(c(i,5)),label)
        axis([-10 10 -10 10 0 30])
    end
end


for i=1:length(bp) %for all branching points
    ec=find(c(:,9)==bp(i)); %These are the starts of new branches - so THERE SHOULD BE THIS MANY BRANCHES
    radius_ratios=c(ec,1)./c(bp(i),1);
    continuations=zeros(size(radius_ratios)); %assume none are continuations
    continuations(find(radius_ratios>limit_ratio))=1; %If the child radius has more than xx% radius as the parent it is a CONTINUATION!
    
    %pause
    for j=1:length(ec)
        if continuations(j)==1
            bo(ec(j))=bo(c(ec(j),9)); %branching order stays the same if it is a continuation
        else
            bo(ec(j))=bo(c(ec(j),9))+1; %add 1 to the parent bo for the extension cylinder
        end
        % This is where I could insert a radius ratio continuation test
        for cyl=ec(j):length(c) %follow from extension cylinder to the top
              ext=c(cyl,10); %re-define extension cylinder
              if ext==0; break; end % stop when we reach the final twig, and go onto the next ec
              bo(ext)=bo(ec(j)); %add 1 to the bo for the extension cylinder
              if ext==any(bp); break; end %stop when we reach a branching point?
              
              if PLOT==1
                  [cyl ext bo(ext)]
                  label=num2str((bo(ext))');        
                  hold on
                  text(double(c(ext,3)),double(c(ext,4)),double(c(ext,5)),label)
                  axis([-10 10 -10 10 0 30])
                  pause
              end
              
        end
        
    end
end
bo(1)=0;

end

