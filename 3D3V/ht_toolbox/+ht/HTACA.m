function htd = HTACA(data, tree, base_rel_tol, tol_ratio, min_rank, max_rank)
%Hierarchical Tensor Adaptive Cross Approximation (HTACA) for tensor decomposition

    [subtree_pool, par2lsub_pool, lsub2par_pool, par2rsub_pool, rsub2par_pool] = ht.get_subtree_pool(tree);
    [tol_ACA,tol_rounding,first_call,max_norm_app] = initialize_tolerance(subtree_pool,base_rel_tol,tol_ratio);
    [min_rank_HTACA, max_rank_HTACA] = get_rank_bounds_for_pivoting(min_rank,max_rank,tree);
    [min_rank_pool, max_rank_pool] = get_rank_limit_pool(tree);
    tree_orig = tree;

    index_lprefix = {};
    index_rprefix = {};

    htd = recursive_htaca(data, index_lprefix, index_rprefix, {}, 1);

    function htd = recursive_htaca(data, index_lprefix, index_rprefix, data_app_list, node)
        tree_now  = subtree_pool{node}; 
        depth     = tree_now.depth; 
        modesizes = tree_now.modesizes; 

        if depth == 0
            data_local = @(cell_data) data(index_lprefix{:},cell_data{:},index_rprefix{:});
            htd.U{1} = data_local({(1:modesizes(1))'}) - ht.evaluate_fiber_sum(data_app_list,1:modesizes(1));
            htd.B{1} = []; 
            return
        elseif depth == 1
            data_local = @(cell_data) data(index_lprefix{:},cell_data{:},index_rprefix{:});
            [htd.U{2}, htd.U{3}, htd.B{1}] = aca_matrix(data_local,data_app_list,tree_now.modesizes, 9, node);
            htd.tree = tree_now; 
            htd.rank = [1 size(htd.U{2},2) size(htd.U{3},2)]; 
        else
            node_l = tree_orig.children(node,1);
            node_r = tree_orig.children(node,2);

            pivot = inf;  
            I = []; J = []; k = 1; 

            while abs(pivot) > tol_ACA(node) || k < min_rank_HTACA(node)
                if k == 1
                    ht_app_list = data_app_list;
                else
                    ht_app_list = [data_app_list ht_preR];
                end

                [pivot, index_l, index_r, index , I, J] = recursive_pivot_search(0, data, index_lprefix, index_rprefix, {}, ht_app_list, I, J, node, 0);
                
                if node == 1 && k == 1
                    max_norm_app = abs(pivot); 
                end

                ht_app_l_list = get_left_subtensors_list(ht_app_list,index_r,node);
                lprefix_new = index_lprefix; 
                rprefix_new = [index_r index_rprefix]; 
                ht_l = recursive_htaca(data, lprefix_new, rprefix_new, ht_app_l_list, node_l);
                
                ht_app_r_list = get_right_subtensors_list(ht_app_list,index_l,node);
                lprefix_new = [index_lprefix index_l]; 
                rprefix_new = index_rprefix; 
                ht_r = recursive_htaca(data, lprefix_new, rprefix_new, ht_app_r_list, node_r);
    
                B_root = 1./(pivot + 1e-15); 

                if k == 1
                    ht_preR = assemble_parent_tensor_lr(B_root,ht_l, ht_r, tree_now, subtree_pool{node_l}, subtree_pool{node_r}, lsub2par_pool{node_l}, rsub2par_pool{node_r});
                else
                    ht_1upd = assemble_parent_tensor_lr(B_root,ht_l, ht_r, tree_now, subtree_pool{node_l}, subtree_pool{node_r}, lsub2par_pool{node_l}, rsub2par_pool{node_r});
                end

                if first_call(node) && k == 1
                    tol_ACA(node) = tol_ACA(node)*max_norm_app; 
                    fro_norm_app  = ht.norm(ht_preR); 
                    tol_rounding(node) = tol_rounding(node)*fro_norm_app; 
                    first_call(node) = false; 
                end

                if k > 1
                    ht_preR = ht.add(ht_preR,ht_1upd); 
                    ht_preR = ht.truncate(ht_preR,tol_rounding(node),min_rank_pool{node},max_rank_pool{node}); 
                end

                if k >= max_rank_HTACA(node)
                    break
                end
                k = k + 1; 
            end

            if node == 1
                htd = ht.truncate(ht_preR,base_rel_tol*fro_norm_app,min_rank_pool{1},max_rank_pool{1});
            else
                htd = ht_preR;
            end
        end
    end

    function [U1, U2, B] = aca_matrix(data_local, data_app_list, Adim, sample_size, node)
        CUR_info.Ej    = []; CUR_info.delta = []; CUR_info.Ei    = [];
        I = []; J = []; k = 1; pivot = inf;
        
        while abs(pivot) > tol_ACA(node) || k < min_rank_HTACA(node)
            I_sample = random_sample(1,Adim(1),I,sample_size);
            J_sample = random_sample(1,Adim(2),J,sample_size);
        
            indices = [{I_sample'} {J_sample'}];
            sample_errs = abs( data_local(indices) - ht.evaluate_index_sum(data_app_list,indices{:}) - get_cur_value(CUR_info,indices{:}) );
            
            [~,ind] = max(sample_errs);
            i_star = I_sample(ind); j_star = J_sample(ind);
        
            indices = [{i_star} {(1:Adim(2))'}];
            err_J   = data_local(indices) - ht.evaluate_fiber_sum(data_app_list,indices{:}) - get_cur_value(CUR_info,indices{:});
            [~,j_star] = max(abs(err_J));

            indices = [{(1:Adim(1))'} {j_star}];
            Ej_ = data_local(indices) - ht.evaluate_fiber_sum(data_app_list,indices{:}) -  get_cur_value(CUR_info,indices{:});
            [val,i_star] = max(abs(Ej_));

            indices = [{i_star} {(1:Adim(2))'}];
            Ei_ = ( data_local(indices) - ht.evaluate_fiber_sum(data_app_list,indices{:}) - get_cur_value(CUR_info,indices{:}) ).';
        
            indices = [{i_star} {j_star}];
            pivot = data_local(indices) - ht.evaluate_index_sum(data_app_list,indices{:}) - get_cur_value(CUR_info,indices{:}) + 1e-15;
            
            CUR_info.delta = [CUR_info.delta,1/pivot];
            CUR_info.Ej    = [CUR_info.Ej, Ej_];
            CUR_info.Ei    = [CUR_info.Ei; Ei_];
            I = [I, i_star]; J = [J, j_star];

            if first_call(node) && k == 1
                tol_ACA(node) = tol_ACA(node)*max_norm_app;
                first_call(node) = false;
            end
        
            if k >= max_rank_HTACA(node)
                break
            end
            k = k + 1;
        end
        
        U1 = CUR_info.Ej; U2 = CUR_info.Ei.'; B = diag(CUR_info.delta); 
    end

    function value = get_cur_value(CUR_info,i,j)
        if ~isempty(CUR_info.delta)
            value = sum(CUR_info.Ej(i,:).*(CUR_info.delta.*CUR_info.Ei(:,j).'),2);
        else
            value = 0; 
        end
    end

    function [pivot, index_l, index_r , index, I, J] = recursive_pivot_search(pivot, data, index_lprefix, index_rprefix, index, ht_app_list, I, J, node, level)
        data_local = @(cell_data) data(index_lprefix{:},cell_data{:},index_rprefix{:});
        tree_now   = subtree_pool{node}; 
        depth      = tree_now.depth; 
        modesizes  = tree_now.modesizes; 

        if depth == 0
            vals              = data_local({(1:modesizes(1))'}) - ht.evaluate_fiber_sum(ht_app_list,1:modesizes(1));
            [abs_p_candi,ind] = max(abs(vals)); 
    
            if abs(pivot) < abs_p_candi
                index = {ind}; pivot = vals(ind); 
            end

            index_l = {}; index_r = {}; I = []; J = [];
            return
        else
            ch_now    = tree_now.children(1,:); 
            modesizes = tree_now.modesizes; 
            D         = tree_now.orders(1); 
            ord_l     = tree_now.orders(ch_now(1)); 
            ord_r     = tree_now.orders(ch_now(2)); 
            j_l       = cell(1,ord_l); j_r       = cell(1,ord_r); 
        
            max_iter    = D - 1; 
            sample_size = ht.determine_sample_size(D); 

            m_size  = [tree_now.freedom(ch_now(1)) tree_now.freedom(ch_now(2))];
            samp_i  = random_sample(1,m_size(1),I,sample_size); 
            samp_j  = random_sample(1,m_size(2),J,sample_size); 
            [j_l{:}]  = ind2sub(modesizes(1:ord_l),samp_i'); 
            [j_r{:}]  = ind2sub(modesizes(ord_l+1:end),samp_j'); 
            
            vals    = data_local([j_l j_r]) - ht.evaluate_index_sum(ht_app_list,j_l{:},j_r{:});
            [abs_p_candi,ind]   = max(abs(vals)); 
    
            if abs(pivot) < abs_p_candi
                index_l = cellfun(@(x) x(ind), j_l, 'UniformOutput', false); 
                index_r = cellfun(@(x) x(ind), j_r, 'UniformOutput', false); 
                pivot   = vals(ind); 
            else
                if isempty(index)
                    disp('Index is empty — pausing.');
                    keyboard    
                end
                index_l = index(1:ord_l); index_r = index(ord_l+1:end); 
            end

            for k = 1: max_iter
                ht_app_l_list   = get_left_subtensors_list(ht_app_list,index_r,node); 
                lprefix_new = index_lprefix; rprefix_new = [index_r index_rprefix]; 
                node_l = tree_orig.children(node,1); 
                [pivot, index_l_l, index_l_r, index_l] = recursive_pivot_search(pivot, data, lprefix_new, rprefix_new, index_l, ht_app_l_list, [], [], node_l, level + 1);
                
                ht_app_r_list = get_right_subtensors_list(ht_app_list,index_l,node);  
                lprefix_new = [index_lprefix index_l]; rprefix_new = index_rprefix; 
                node_r = tree_orig.children(node,2); 
                [pivot, index_r_l, index_r_r, index_r] = recursive_pivot_search(pivot, data, lprefix_new, rprefix_new, index_r, ht_app_r_list, [], [], node_r, level + 1);

                index = [index_l index_r];
            end

            if level == 0
                if ord_l > 1, i_ = sub2ind(modesizes(1:ord_l),index_l{:}); else, i_ = index_l{1}; end
                if ord_r > 1, j_ = sub2ind(modesizes(ord_l+1:end),index_r{:}); else, j_ = index_r{1}; end
                I = [I i_]; J = [J j_]; 
            else
                I = []; J = []; 
            end
        end
    end

    % --- PARFOR OPTIMIZED LEFT SUBTENSOR ---
    function ht_l_list = get_left_subtensors_list(ht_list,index_r,node)
        if isempty(ht_list)
            ht_l_list = {};
            return;
        end
        tree_now    = subtree_pool{node};
        node_l      = tree_orig.children(node,1);
        tree_l      = subtree_pool{node_l};
        par2lsub    = par2lsub_pool{node};
        N_list      = length(ht_list);
        D           = tree_now.orders(1);
        D_l         = tree_l.orders(1);
        N_l         = tree_l.N_node;
        N_nonleaf_r = D - D_l - 1;
    
        ht_l_list = cell(1,N_list);

        parfor k = 1:N_list
            curr_ht = ht_list{k};
            U_temp = curr_ht.U; B_temp = curr_ht.B; rank_  = curr_ht.rank;
            loc_U = cell(1,N_l); loc_B = cell(1,N_l); loc_rank = zeros(1, length(rank_)); 

            for d = 1:D_l
                n_l = tree_l.dim2ind(d); n_l_par = tree_now.dim2ind(d); 
                loc_U{n_l} = U_temp{n_l_par}; loc_rank(n_l) = rank_(n_l_par); 
            end
            for d = D_l+1 : D
                U_temp{tree_now.dim2ind(d)} = U_temp{tree_now.dim2ind(d)}(index_r{d-D_l},:);
            end
            for i = tree_now.postorder_nonleaf(N_nonleaf_r+1:end-1)
                loc_B{par2lsub(i)} = curr_ht.B{i}; loc_rank(par2lsub(i)) = rank_(i);
            end
            for i = tree_now.postorder_nonleaf(1:N_nonleaf_r)
                l = tree_now.children(i,1); r = tree_now.children(i,2); 
                B = B_temp{i}; Ul = U_temp{l}; Ur = U_temp{r}; 
                r_l = rank_(l); r_r = rank_(r); r_i = rank_(i); 
                AB = pagemtimes(Ul,reshape(B,r_l,[])); CC = repmat(Ur, 1, r_i); 
                V_ = reshape(AB.*CC,1,r_r,r_i); U_temp{i} = reshape(sum(V_,2),1,r_i); 
            end
            if tree_l.depth == 0
                r = tree_now.children(1,2); B = B_temp{1}*U_temp{r}.'; l = tree_now.children(1,1); 
                loc_U{1} = U_temp{l}*B; loc_B{1} = []; 
                ht_l_list{k} = struct('U', {loc_U}, 'B', {loc_B}, 'rank', loc_rank);
                continue;
            end
            r = tree_now.children(1,2); B = curr_ht.B{1}*U_temp{r}.'; l = tree_now.children(1,1); 
            if size(B,1) > 1, loc_B{1} = ht.ttm(curr_ht.B{l},B.',3); else, loc_B{1} = curr_ht.B{l}*B; end
            loc_rank(1) = 1;
            ht_l_list{k} = struct('U', {loc_U}, 'B', {loc_B}, 'rank', loc_rank, 'tree', tree_l);
        end
    end

    % --- PARFOR OPTIMIZED RIGHT SUBTENSOR ---
    function ht_r_list = get_right_subtensors_list(ht_list,index_l,node)
        if isempty(ht_list)
            ht_r_list = {};
            return;
        end
        tree_now    = subtree_pool{node}; node_r      = tree_orig.children(node,2); 
        tree_r      = subtree_pool{node_r}; par2rsub    = par2rsub_pool{node}; 
        N_list      = length(ht_list); D = tree_now.orders(1); D_r = tree_r.orders(1); 
        D_l         = D - D_r; N_r = tree_r.N_node; N_nonleaf_r = D_r - 1; 
        ht_r_list = cell(1,N_list);

        parfor k = 1: N_list
            curr_ht = ht_list{k};
            U_temp = curr_ht.U; B_temp = curr_ht.B; rank_  = curr_ht.rank; 
            loc_U = cell(1,N_r); loc_B = cell(1,N_r); loc_rank = zeros(1, length(rank_));

            for d = 1: D_r
                n_r = tree_r.dim2ind(d); n_r_par = tree_now.dim2ind(D_l+d); 
                loc_U{n_r} = U_temp{n_r_par}; loc_rank(n_r) = rank_(n_r_par); 
            end
            for d = 1 : D_l
                U_temp{tree_now.dim2ind(d)} = U_temp{tree_now.dim2ind(d)}(index_l{d},:);
            end
            for i = tree_now.postorder_nonleaf(1: N_nonleaf_r)
                loc_B{par2rsub(i)} = curr_ht.B{i}; loc_rank(par2rsub(i)) = curr_ht.rank(i);
            end
            for i = tree_now.postorder_nonleaf(N_nonleaf_r+1:end-1)
                l = tree_now.children(i,1); r = tree_now.children(i,2); 
                B = B_temp{i}; Ul = U_temp{l}; Ur = U_temp{r}; 
                r_l = rank_(l); r_r = rank_(r); r_i = rank_(i); 
                AB = pagemtimes(Ul,reshape(B,r_l,[])); CC = repmat(Ur, 1, r_i); 
                V_ = reshape(AB.*CC,1,r_r,r_i); U_temp{i} = reshape(sum(V_,2),1,r_i); 
            end
            if tree_r.depth == 0
                l = tree_now.children(1,1); B = U_temp{l}*B_temp{1}; r = tree_now.children(1,2); 
                loc_U{1} = U_temp{r}*B.'; loc_B{1} = []; 
                ht_r_list{k} = struct('U', {loc_U}, 'B', {loc_B}, 'rank', loc_rank);
                continue;
            end
            l = tree_now.children(1,1); B = U_temp{l}*curr_ht.B{1}; r = tree_now.children(1,2); 
            if size(B,2) > 1, loc_B{1} = ht.ttm(curr_ht.B{r},B,3); else, loc_B{1} = curr_ht.B{r}*B; end
            loc_rank(1) = 1; 
            ht_r_list{k} = struct('U', {loc_U}, 'B', {loc_B}, 'rank', loc_rank, 'tree', tree_r);
        end
    end

    function htd = assemble_parent_tensor_lr(B_root,ht_l, ht_r, tree_par, tree_l, tree_r, lsub2par, rsub2par)
        htd = ht.zeros(tree_par); htd.B{1} = B_root; htd.rank = ones(1,tree_par.N_node); 
        for d = 1: tree_l.orders(1)
            node = tree_l.dim2ind(d); htd.U{lsub2par(node)} = ht_l.U{node}; htd.rank(lsub2par(node)) = size(ht_l.U{node},2); 
        end
        for i = tree_l.postorder_nonleaf
            htd.B{lsub2par(i)} = ht_l.B{i}; htd.rank(lsub2par(i)) = size(ht_l.B{i},3); 
        end
        for d = 1: tree_r.orders(1)
            node = tree_r.dim2ind(d); htd.U{rsub2par(node)} = ht_r.U{node}; htd.rank(rsub2par(node)) = size(ht_r.U{node},2); 
        end
        for i = tree_r.postorder_nonleaf
            htd.B{rsub2par(i)} = ht_r.B{i}; htd.rank(rsub2par(i)) = size(ht_r.B{i},3); 
        end
        htd.tree = tree_par; 
    end

    function [min_rank_pool, max_rank_pool] = get_rank_limit_pool(tree)
        min_rank_pool = cell(1,tree.N_node); max_rank_pool = cell(1,tree.N_node);
        min_rank_pool{1} = min_rank; max_rank_pool{1} = max_rank;
        l = tree.children(1,1); r = tree.children(1,2);
        dfs(l,lsub2par_pool{l},1); dfs(r,rsub2par_pool{r},1);

        function dfs(node,sub2par,par_node)
            node_l = tree.children(node,1); node_r = tree.children(node,2); tree_now = subtree_pool{node};
            min_rank_pool{node} = zeros(1,subtree_pool{node}.N_node); max_rank_pool{node} = zeros(1,subtree_pool{node}.N_node);
            for d = 1: tree_now.orders(1)
                node_d = tree_now.dim2ind(d);  
                min_rank_pool{node}(node_d) = min_rank_pool{par_node}(sub2par(node_d)); 
                max_rank_pool{node}(node_d) = max_rank_pool{par_node}(sub2par(node_d)); 
            end
            for i = tree_now.postorder_nonleaf(1:end-1)
                min_rank_pool{node}(i) = min_rank_pool{par_node}(sub2par(i)); 
                max_rank_pool{node}(i) = max_rank_pool{par_node}(sub2par(i)); 
            end
            min_rank_pool{node}(1) = 1; max_rank_pool{node}(1) = 1; 
            if node_l > 0, dfs(node_l,lsub2par_pool{node_l},node); end
            if node_l > 0, dfs(node_r,rsub2par_pool{node_r},node); end
        end
    end
end

function [tol_ACA,tol_rounding,first_call,max_norm_app] = initialize_tolerance(subtree_pool,base_rel_tol,tol_ratio)
    tree_orig  = subtree_pool{1}; depth_orig = tree_orig.depth; children   = tree_orig.children;       
    nonleaf    = tree_orig.nonleaf; tol_ACA = zeros(1, tree_orig.N_node); tol_rounding = zeros(1, tree_orig.N_node); 
    first_call   = true(1, tree_orig.N_node); tol_ACA(1) = base_rel_tol; tol_rounding(:) = 1e-14; max_norm_app = 1; 
    dfs(1)
    function dfs(node)
        ch = children(node,:); 
        if nonleaf(ch(1)) 
            depth  = subtree_pool{ch(1)}.depth; ratio_ = (tol_ratio)^(depth_orig - depth); 
            tol_ACA(ch(1)) = ratio_*tol_ACA(1); dfs(ch(1)); 
        end
        if nonleaf(ch(2))
            depth  = subtree_pool{ch(2)}.depth; ratio_ = (tol_ratio)^(depth_orig - depth);
            tol_ACA(ch(2)) = ratio_*tol_ACA(1); dfs(ch(2)) 
        end
    end
end

function [min_rank_HTACA, max_rank_HTACA] = get_rank_bounds_for_pivoting(min_rank,max_rank,tree)
    min_rank_HTACA = zeros(1,tree.N_node); max_rank_HTACA = zeros(1,tree.N_node); 
    for i = tree.postorder_nonleaf
        ch = tree.children(i,:); 
        min_rank_HTACA(i) = 2*max(min_rank(ch(1)),min_rank(ch(2)));
        max_rank_HTACA(i) = 2*max(max_rank(ch(1)),max_rank(ch(2)));
        min_rank_HTACA(i) = min([min_rank_HTACA(i),tree.freedom(ch(1)),tree.freedom(ch(2))]);
        max_rank_HTACA(i) = min([max_rank_HTACA(i),tree.freedom(ch(1)),tree.freedom(ch(2))]);
    end
end

function sample = random_sample(low, high, unavailable, sample_size)
    sample = [];  
    if isempty(unavailable)
        val = randi(high - low + 1) + low - 1; sample = [sample, val]; unavailable = [unavailable, val];     
        num_samples = min([sample_size, (high - low) + 1]);
    else
        num_samples = min([sample_size, (high - low) - length(unavailable) + 1]);
    end
    while length(sample) < num_samples
        val = randi(high - low + 1) + low - 1;  
        if ~ismember(unavailable, val)         
            sample = [sample, val]; unavailable = [unavailable, val];  
        end
    end
end