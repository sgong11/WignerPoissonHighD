  function [f1] = Etof2d(f, E1)
  
  
  % Wei Guo's function to apply E1 to f;
  % E1 is also in the htensor format
  
    if(ndims(f)~=4 || ndims(E1)~=2)
        error('dimensions of f and E does not match');
    end
    
    
    

    f1 = f;
    Ax_leaf = [2 3]; % leaf node for E1

    x_leaf = [4 5]; % leaf node for f

    

    for ii=1:2
        
        k_A = size(E1.U{Ax_leaf(ii)}, 2);

        AxU = cell(1, k_A);
        for jj=1:k_A


            %Ajj = reshape(A.U{ii}(:, jj), m(mu), n(mu));
            AxU{jj} = bsxfun(@times,  E1.U{Ax_leaf(ii)}(:,jj), f.U{x_leaf(ii)});


        end
%         field1 = 'type'; value1 = {'.', '{}'};
%         field2 = 'subs'; value2 = {'U', {x_leaf(ii)}};
%         s = struct(field1, value1, field2, value2);
%       
%         f1 = builtin( 'subsasgn', f1, s, cell2mat(AxU)); 
        tt = x_leaf(ii);
        f1.U{tt} = cell2mat(AxU);
    end

    ii = 1;
    jj = 2;

    sz_A = size(E1.B{ii});
    sz_A(end+1:3) = 1; % force the sz_A has dimension three
    sz_x = size(f.B{jj});
    sz_x(end+1:3) = 1;

    % "3-D Kronecker product"
    f1.B{jj} = zeros(sz_A.*sz_x);
    for ss=1:sz_A(3)
        for rr=1:sz_x(3)
            f1.B{jj}(:, :, rr+(ss-1)*sz_x(3)) = ...
                kron(E1.B{ii}(:, :, ss), f.B{jj}(:, :, rr));
        end
    end

    f1.is_orthog = false;
end