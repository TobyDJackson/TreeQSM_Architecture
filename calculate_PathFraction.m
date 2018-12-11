function [pf] = calculate_PathFraction(ext_index,cyl_length,parent_index,PLOT)

    twigs=find(ext_index==0);
    num_twigs=length(twigs);

    for twig=1:length(twigs)
        path_len=0;
        parent=parent_index(twig);
        while parent>=2  %follow the path from this twig downwards
              path_len=path_len+cyl_length(parent);  %add on this length to the total length
              parent=parent_index(parent); %reassign the parent to the next one down
              %pause
        end    
        total_length(twig)=path_len;    %This reads the total path length for twig (count) into that row of total_length  
    end
    if PLOT==1
        hist(total_length)
        title('Path lengths')
        pause
    end
    
    pf=mean(total_length)/max(total_length); 
    clearvars parent path_len 
end

