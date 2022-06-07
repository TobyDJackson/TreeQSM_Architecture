function [CylData_BO_updated] = add_branching_order_above(CylData,row_index)
    
    % find all extension cylinders from a given cylinder.
    % Use a path fraction - twigs first approach. 
    % For each twig, follow the path down, saving the row numbers of each
    % cylinder in that path. 
    % if cyl number ==row_index is reached, stop and save the row numbers
    % cat together all these row numbers and find the unique ones
    % Plot a color plot of this to check (USE THIS IN ARCHITECTURE WORK TOO)
    % Add 1 to the branching order of all these cyls.
    % This could be used for re-iterations as well.
    
    % However, for the purposes of the Borneo paper I am only trying to
    % distinguish the crown from the stem, so I can just add 1 to
    % everything above it. Also, it would be good if the more mechanistic
    % approach worked.
    
    CylData_temp=CylData;

end

