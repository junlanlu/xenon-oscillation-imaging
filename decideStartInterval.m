function [ind_start, ind_inter, ind_end] = decideStartInterval(mask)

    % determine the starting slice and the interval for the montage
    num_slice = 16;

    sum_line = squeeze(sum(mask, [1 2]));
    binary_arr = sum_line > 300;

    [ind_start, ind_end] = getIndexMaxOnes(binary_arr);
    %ind_start = 30;
    %ind_end = 40;
    flt_inter = (ind_end-ind_start)/num_slice;

    % use 0.4 as a threshold to decide interval number
    if mod(flt_inter,1) > 0.4
        ind_inter = ceil(flt_inter);
    else
        ind_inter = floor(flt_inter);
    end

 end


function [index_start, index_end] = getIndexMaxOnes(arr)
    %returns the starting index and ending index of the max consecutive ones
    %intitialize count
    cur_count = 0;
    cur_start = 0;

    max_count = 0;
    pre_state = 0;

    index_start = 0;
    index_end = 128;
    n = size(arr, 1);
    for i = 0:n-1
        if (arr(i+1) == 0)
            cur_count = 0;
            if((pre_state == 1) && (cur_start == index_start))
                index_end = i-1;
            end
            pre_state = 0;
        else
            if(pre_state == 0)
                cur_start = i;
                pre_state = 1;
            end
            cur_count = cur_count + 1;
            if(cur_count > max_count)
                max_count = cur_count;
                index_start = cur_start;
            end
        end
    end
   
end
