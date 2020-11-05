function n_out = find_vertical_segments(centres,height)
    % find the maximum number of equally spaced vertical segments containing 
    % at least 1 cylinder
    scaled_height=centres(:,3)/height;
    
    n=2;
    temp=ones(n,1);
    while min(temp)>0;
        temp=ones(n,1);
        for i=1:n % loop over segments
            % find cylinders within this height segment
            temp(i)=length(find(scaled_height>(i-1)/n & scaled_height<=i/n));
        end
        n=n+1;
        %plot(temp); title(num2str(n)); pause; close all;
    end
    n_out=n-1;
end