function f_out = SL1DlocalWENO_vec(fn, xshift, xshift_pos, epsilon, caseinx)
% fn : 6 x N

    if caseinx == 1
        Q1 = compute_weno5_Qiu_vec(fn(2:end,:), xshift_pos, epsilon);
        Q0 = compute_weno5_Qiu_vec(fn(1:5,:),  xshift_pos, epsilon);
        f_out = fn(4,:) - xshift * (Q1 - Q0);

    else
        fnr = fn(end:-1:1,:); % 6 --> 1
        Q1 = compute_weno5_Qiu_vec(fnr(2:end,:), xshift_pos, epsilon); % 5 -->1
        Q0 = compute_weno5_Qiu_vec(fnr(1:5,:),  xshift_pos, epsilon); % 6 --> 2
        f_out = fn(3,:) - xshift * (Q0 - Q1);
    end
end
