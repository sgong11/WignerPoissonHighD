function z = add(x, y)
%ADD Add two HTD tensors with identical tree structure.
%
%   z = ht.add(x,y)
%
%   Input:
%       x,y - HTD structure with fields .U, .B, .tree, .rank
%   Output:
%       z - HTD structure with fields .U, .B, .tree, .rank

    z.tree = x.tree;
    z.U = cell(1, x.tree.N_node);
    z.B = cell(1, x.tree.N_node);
    z.rank = zeros(1, x.tree.N_node);

    for d = 1: x.tree.orders(1) % leaf node: concatenate bases
        node = x.tree.dim2ind(d);
        Ux = x.U{node};  % [n x r1]
        Uy = y.U{node};  % [n x r2]
        z.U{node} = [Ux, Uy];  % [n x (r1 + r2)]
        z.rank(node) = size(z.U{node}, 2);
    end

    for i = x.tree.postorder_nonleaf
        c1 = x.tree.children(i, 1);
        c2 = x.tree.children(i, 2);
        r1x = x.rank(c1); r2x = x.rank(c2); rx = x.rank(i);
        r1y = y.rank(c1); r2y = y.rank(c2); ry = y.rank(i);

        if i == 1  % root node: sum is along mode-3
            Bx = x.B{i};  % [r1x x r2x x 1]
            By = y.B{i};  % [r1y x r2y x 1]
            Bz = zeros(r1x + r1y, r2x + r2y, 1);
            Bz(1:r1x, 1:r2x, 1) = Bx;
            Bz(r1x+1:end, r2x+1:end, 1) = By;
        else
            Bx = x.B{i};  % [r1x x r2x x rx]
            By = y.B{i};  % [r1y x r2y x ry]
            Bz = zeros(r1x + r1y, r2x + r2y, rx + ry);
            Bz(1:r1x, 1:r2x, 1:rx) = Bx;
            Bz(r1x+1:end, r2x+1:end, rx+1:end) = By;
        end

        z.B{i} = Bz;
        z.rank(i) = size(Bz, 3);
    end
end
