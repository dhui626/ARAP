
%% ARAP参数化
clear;
visualize = true;       %显示参数化结果
iteration = 50;         %迭代次数

%% 计算Tutte参数化结果作为初值
[x, t] = readObj('obj\Balls');
B = findBoundary(x, t);
v_count = size(x,1);         %总共点的个数
f_count = size(t,1);          %总共边的个数

uv = tutte(x, t, B);

if (visualize)
    subplot(231); drawmesh(t, x, B);title('input');
    subplot(232); drawmesh(t, uv, B);title('initial(Tutte)');
end

%% 计算并存储面积At
At = zeros(1, f_count);
for i=1:f_count
    index = t(i,:);
    a = norm(x( index(1),: )-x( index(2),: ),2);
    b = norm(x( index(2),: )-x( index(3),: ),2);
    c = norm(x( index(3),: )-x( index(1),: ),2);
    
    temp = (a*a+c*c-b*b)/2/a/c;
    At(i) = a*c*sin( acos(temp) )/2;    %计算三角面片的面积
end
At2 = [At; At];
At2 = At2(:);

%% 将原始3D triangle等距参数化到平面

flattened = zeros(f_count, 6);	%存储等距参数化结果，每行分别是三个点的二维坐标

for i=1:f_count
    index = t(i,:);
    a = norm(x( index(1),: )-x( index(2),: ),2);
    b = norm(x( index(2),: )-x( index(3),: ),2);
    c = norm(x( index(3),: )-x( index(1),: ),2);
    
    temp = (a*a+c*c-b*b)/2/a/c;
    flattened(i,[3 4]) = [a 0];
    flattened(i,[5 6]) = [c*temp c*sqrt(1-temp*temp)];
end

%% 计算M
M = sparse(2*f_count, v_count);
for i=1:f_count
	index = t(i,:);
    e12 = [flattened(i,[1 2]) - flattened(i,[3 4]); flattened(i,[1 2]) - flattened(i,[5 6])];
    
    temp = zeros(2, v_count);
    temp(1, index(1) )=1; temp(1, index(2) )=-1;
    temp(2, index(1) )=1; temp(2, index(3) )=-1;
    
    M([2*i-1 2*i], :) = e12 \ temp;
end
M = At2.*M;

%% 生成稀疏矩阵并预分解
MT = M';
lhs = MT * M;
dLhs = decomposition(lhs);

for it = 1:iteration %迭代次数
    %% Local

    J = M*uv;
    L = zeros(2*f_count, 2);
    
    for i=1:f_count
        Jt = J([2*i-1 2*i], :);
        [U,~,V] = svd(Jt);
        L([2*i-1 2*i], :) = U*V';
    end
    L = At2.*L;

    %% Global
    rhs = MT*L;
    uv = dLhs \ rhs;
    
    if (it==1 && visualize)
        %hold off
        subplot(234); drawmesh(t, uv, B);title('ARAP(iteration=1)');
    end
    
    if (it==10 && visualize)
        subplot(235); drawmesh(t, uv, B);title('ARAP(iteration=10)');
    end
end

if (visualize)
    subplot(236); drawmesh(t, uv, B);title('ARAP(final result)');
end

