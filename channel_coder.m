

function out_seqs = channel_coder(in_seq, codeOrder, dbg)
    %{
    params :
            in_seq :[vector]: 1-d sequence of zeros and ones. 
                Should have length m*codeOrder to be correctly reshaped to square matrix.
                incorrect length raise warning and missing some number of last bits

            codeOrder :[scalar]: number of data bits. 

            dbg :[optional: bool]: debug mode to checking results by hand. 
                    default false

        returns :
            out_seq :[vector]: 1-d sequence of zeros and ones.
                size: m * (square root of codeOrder + 1)^2
    %}


    if nargin < 3
        dbg = false;
    end

    check_input(in_seq);
    codeOrder = check_codeOrder(codeOrder);
    m = fix(length(in_seq)/codeOrder); % how many sequences will be encoded
    missed = length(in_seq) - m * codeOrder;

    if missed~=0
        warning('Input size should be multiple of codeOrder. \n %d missed symbols', missed);
    end

    out_seqs = [];

    for i = 0:m-1
        out = base_channel_coder(in_seq((1 + i*codeOrder):((i+1)*codeOrder)));
        out_seqs = [out_seqs out];

        if dbg
            disp(["seq", out(1:codeOrder), 'check bits', out(codeOrder+1:length(out))]);
        end
    end

end


function out_seq = base_channel_coder(in_seq) 
    %{
        params :
            in_seq :[vector]: 1-d sequence of zeros and ones. must have length to be correctly reshaped to square matrix
        returns :
            out_seq :[vector]: 1-d sequence of zeros and ones.
                size: (square root of length in_seq + 1) in square

            

        example:
            size in_seq is 16
            mtx shape = 4x4
            col_ch size = 4
            transposed row_ch size = 4
            last_checked size = 1

            out_seq size = 16 + 4 + 4 + 1 = 25

            location:  data  col_ch row_ch last_checked

    %}   

    s = fix(sqrt(length(in_seq)));
    in_mtx = reshape(in_seq, [s, s]);
    missed = length(in_seq)-s^2;
    if missed~=0
        error('missed symbols %s', missed);
    end
    [col_ch, row_ch, last_checked] = check_parity(in_mtx);
    out_seq = [in_seq col_ch row_ch' last_checked];

end

function [col_ch, row_ch, last_checked] = check_parity(mtx)
    %{
        function creates the check bits for matrix

        params :
            mtx :[matrix]: square matrix of zeros and ones
        returns : 
            col_ch :[vector]: XOR sum of elements each column. 
            row_ch :[column vector]: XOR sum of elements each row. 
            last_checked :[scalar]: XOR(col_ch) OR XOR(row_ch).
            
        example:
             mtx(1,1)  mtx(2,1)  mtx(3,1)  row_ch(1)
             mtx(2,1)  mtx(2,2)  mtx(3,2)  row_ch(2)
             mtx(3,1)  mtx(2,3)  mtx(3,3)  row_ch(3)
             col_ch(1) col_ch(2) col_ch(3) last_checked

    %}
    col_ch = mod(sum(mtx, 1), 2); 
    row_ch = mod(sum(mtx, 2), 2); 
    last_checked = mod(sum(col_ch), 2) + mod(sum(row_ch), 2);
    if last_checked == 2
        last_checked = 1;
    end
    
end


function check_input(seq)
    %{
    check for binary input

    params:
        seq :[vector]:
    %}

    for i = 1:length(seq)
        if ~ismember(seq(i), [0 1])
            error('data should be consist from zeros and ones');
        end
    end

end


function codeOrder_out = check_codeOrder(codeOrder)
    %{
    code order validation
    2 or 3 will be square
    codeOrder = int^2 is valid
    other raise error 
    %}


    if ismember(codeOrder, [2 3])
        codeOrder_out = codeOrder^2;
        warning('codeOrder has been modified from %d to %d \n', codeOrder, codeOrder_out);
    elseif fix(sqrt(codeOrder)) == sqrt(codeOrder)
        codeOrder_out = codeOrder;
    else  
        error('invalid CodeOrder. CodeOrder should be 2, 3 or any int^2');
    end

end

