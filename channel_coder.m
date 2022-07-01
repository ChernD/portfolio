function out_seq = channel_coder(in_seq)
    
    s = fix(sqrt(length(in_seq)));
    in_mtx = reshape(in_seq, [s, s]);
    if s^2 ~= length(in_seq)
        disp(['missed symbols', length(in_seq)-s^2]);
    end
    [col_ch, row_ch, last_checked] = check_parity(in_mtx);
    out_seq = [in_seq(1:s^2) col_ch row_ch' last_checked];

end

function [col_ch, row_ch, last_checked] = check_parity(mtx)
    col_ch = mod(sum(mtx, 1), 2); % parity checking along columns. return row vector
    row_ch = mod(sum(mtx, 2), 2); % parity checking along rows. return column vector
    last_checked = mod(sum(col_ch), 2) + mod(sum(row_ch), 2);
    if last_checked == 2
        last_checked = 1;
    end
    
end



