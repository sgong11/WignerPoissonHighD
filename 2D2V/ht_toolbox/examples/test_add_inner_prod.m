
% test_add.m
% Unit test for ht.add function
toolbox_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(toolbox_root);

addpath('htucker_1.2');

fprintf("===== Testing ht.add =====\n");

% Define a 4D balanced tree
children4 = [2, 3; 4, 5; 6, 7; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind4  = [4, 5, 6, 7];  % map mode i → node index
modesizes4 = [16 16 16 16];
tree = ht.make_tree(children4, dim2ind4, modesizes4);

% Set ranks for two tensors
rank1 = [1 4 4 6 8 6 10];
rank2 = [1 3 5 7 6 5 9];

% Generate two random HTD tensors
x = ht.rand(tree, rank1);
y = ht.rand(tree, rank2);

% Add them
z = ht.add(x, y);

% Generate random test indices
Ntest = 200;
id = randi([1, 16], Ntest, 4);

% Evaluate all three
vx = ht.evaluate_index(x, id);
vy = ht.evaluate_index(y, id);
vz = ht.evaluate_index(z, id);

% Check relative error
vxy = vx + vy;
relerr = norm(vz - vxy) / norm(vxy);
fprintf("Relative error: %.2e\n", relerr);

% Plot if needed
figure;
plot(real(vz), 'b.-'); hold on;
plot(real(vxy), 'r--');
legend('z = x + y', 'vx + vy');
title('HTD Addition Evaluation');
xlabel('Test index');
ylabel('Value');
grid on;


x_ref = htensor(tree.children,tree.dim2ind,x.U,x.B,false);
y_ref = htensor(tree.children,tree.dim2ind,y.U,y.B,false);

tic
for k = 1: 1000
    val1 = innerprod(x_ref,y_ref);
end
toc

tic
for k = 1: 1000
val2 = ht.inner_product(x, y);
end
toc
fprintf("Relative error: %.2e\n", abs(val1-val2)/abs(val1));

