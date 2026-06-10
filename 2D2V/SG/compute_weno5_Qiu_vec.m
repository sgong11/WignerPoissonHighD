function f_m = compute_weno5_Qiu_vec(f_values, xshift, epsilon)
% f_values : 5 x N
% f_m      : 1 x N

    if nargin < 3
        epsilon = 1e-6;
    end

    % constants
    C0 = [1/30 0 -1/24 0 1/120;
         -13/60 -1/24 1/4 1/24 -1/30;
          47/60 5/8 -1/3 -1/8 1/20;
          9/20 -5/8 1/12 1/8 -1/30;
         -1/20 1/24 1/24 -1/24 1/120];

    ci = [ 1/3 -7/6 11/6 0 0;
           0  -1/6  5/6 1/3 0;
           0   0    1/3 5/6 -1/6];

    di = [0.1; 0.6; 0.3];

    D1 = [1 -4 3 0 0;
          0 1 0 -1 0;
          0 0 3 -4 1];

    D2 = [1 -2 1 0 0;
          0 1 -2 1 0;
          0 0 1 -2 1];

    % smoothness indicators
    beta = (D1*f_values).^2 + 13/3*(D2*f_values).^2;   % 3 x N

    wi = di ./ (beta + epsilon).^2;                    % 3 x N
    wi = wi ./ sum(wi,1);                              % normalized
    
    C = wi.' * ci;  % 5 x N

    % polynomial basis
    p = [1; xshift; xshift^2; xshift^3; xshift^4];

    % reconstruction
    
    f_m = f_values.' * C0;
    f_m(:,1) = sum(f_values.' .* C, 2); % N x 1
    f_m = (f_m * p).';
end
