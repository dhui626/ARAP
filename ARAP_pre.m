
%% ARAP 2D 网格变形 初始化
visualize = false;
%enableservice('AutomationServer', true);

%% 计算Tutte参数化结果作为初值
%[x, t] = readObj('obj\Balls');
B = findBoundary(x, t);
v_count = size(x,1);         %总共点的个数
f_count = size(t,1);          %总共边的个数

uv = tutte(x, t, B);

if (visualize)
    subplot(231); drawmesh(t, x, B);title('input');
    subplot(232); drawmesh(t, uv, B);title('initial(Tutte)');
end

%% 计算并存储cot
cot_data = zeros( v_count, v_count );
for i=1:f_count
    index = t(i,:);
    a = norm(x( index(1),: )-x( index(2),: ),2);
    b = norm(x( index(2),: )-x( index(3),: ),2);
    c = norm(x( index(3),: )-x( index(1),: ),2);
    
    temp = (a*a+c*c-b*b)/2/a/c;
    cot_data(index(2), index(3)) = cot( acos(temp) );
    temp = (a*a+b*b-c*c)/2/a/b;
    cot_data(index(3), index(1)) = cot( acos(temp) );
    temp = (b*b+c*c-a*a)/2/b/c;
    cot_data(index(1), index(2)) = cot( acos(temp) );
end

%% 记录(i,j)对应的点的指标（半边数据结构）
he_vdata = zeros( v_count, v_count );                %(i,j)元素的值是对应第三个点的指标
for i=1:f_count
    index = t(i,:);
    he_vdata(index(1), index(2)) = index(3);
    he_vdata(index(2), index(3)) = index(1);
    he_vdata(index(3), index(1)) = index(2);
end

%% 记录(i,j)对应的面的指标（半边数据结构）
he_fdata = zeros( v_count, v_count );                %(i,j)元素的值是对应面的指标
for i=1:f_count
    index = t(i,:);
    he_fdata(index(1), index(2)) = i;
    he_fdata(index(2), index(3)) = i;
    he_fdata(index(3), index(1)) = i;
end

%% 记录每个顶点的邻接点指标，防止每次都做重复计算
vNeighber = zeros( v_count, v_count );                %(i,j)元素的值是对应第三个点的指标
for i=1:f_count
    index = t(i,:);
    vNeighber(index(1), index(2)) = 1;
    vNeighber(index(1), index(3)) = 1;
    
    vNeighber(index(2), index(3)) = 1;
    vNeighber(index(2), index(1)) = 1;
    
    vNeighber(index(3), index(1)) = 1;
    vNeighber(index(3), index(2)) = 1;
end

%% 生成稀疏矩阵
lhs = zeros(v_count, v_count);

for i=1:v_count
    lhs(i,:) = -( cot_data(i,:) + cot_data(:,i)' );
    lhs(i,i) = -sum(lhs(i,:),2);
end

lhs(1,1) = lhs(1,1) + 1;

dLhs = decomposition(lhs);


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

for iteration = 1:10 %迭代次数
    %% Local

    %计算每个三角形对应变换的Jacobi矩阵 进行svd分解 计算并记录Lt
    Lt = zeros(2,2,f_count);

    for i=1:f_count
        index = t(i,:);
        Jt = cot_data(index(1),index(2)) * (uv(index(1),:) - uv(index(2),:))' * (flattened(i,[1 2]) - flattened(i,[3 4]))...
            + cot_data(index(2),index(3)) * (uv(index(2),:) - uv(index(3),:))' * (flattened(i,[3 4]) - flattened(i,[5 6]))...
            + cot_data(index(3),index(1)) * (uv(index(3),:) - uv(index(1),:))' * (flattened(i,[5 6]) - flattened(i,[1 2]));
        
%         Jt = [uv(index(1),:) - uv(index(2),:); uv(index(2),:) - uv(index(3),:)]' / [flattened(i,[1 2]) - flattened(i,[3 4]); flattened(i,[3 4]) - flattened(i,[5 6])]'; 
        
        [U,S,V] = svd(Jt);
        if ( det(Jt)>0 )
            Lt(:,:,i) = U*V';
        else
            Lt(:,:,i) = U*diag([1 -1])*V';
        end
    end

    %% Global

    rhs = zeros (v_count, 2);
    for i=1:v_count
        vj = find(vNeighber(i,:) == 1);
        vj_size = size(vj,2);
        sum2d = [0; 0];
        
        for m=1:vj_size
            f_index = he_fdata( i, vj(m) );
            if (f_index ~=0)
                if (t( f_index , 1) == i)
                    sum2d = sum2d + cot_data( i,vj(m) ) * Lt(:,:,f_index) * ( flattened(f_index,[1 2]) - flattened(f_index,[3 4]) )';
                elseif (t( f_index, 2) == i)
                    sum2d = sum2d + cot_data( i,vj(m) ) * Lt(:,:,f_index) * ( flattened(f_index,[3 4]) - flattened(f_index,[5 6]) )';
                elseif (t( f_index, 3) == i)
                    sum2d = sum2d + cot_data( i,vj(m) ) * Lt(:,:,f_index) * ( flattened(f_index,[5 6]) - flattened(f_index,[1 2]) )';
                end
            end
            
            f_index = he_fdata( vj(m), i );
            if (f_index ~=0)
                if (t( f_index, 1) == i)
                    sum2d = sum2d + cot_data( vj(m), i ) * Lt(:,:,f_index) * ( flattened(f_index,[1 2]) - flattened(f_index,[5 6]) )';
                elseif (t( f_index, 2) == i)
                    sum2d = sum2d + cot_data( vj(m), i ) * Lt(:,:,f_index) * ( flattened(f_index,[3 4]) - flattened(f_index,[1 2]) )';
                elseif (t( f_index, 3) == i)
                    sum2d = sum2d + cot_data( vj(m), i ) * Lt(:,:,f_index) * ( flattened(f_index,[5 6]) - flattened(f_index,[3 4]) )';
                end
            end
            
            rhs(i,:) = sum2d';
        end
    end

    uv = dLhs\ rhs;
    
end

lhs(1,1) = lhs(1,1) - 1; %前面为了使矩阵满秩固定了第一个点，现在把这个固定点的约束去除变为一般的cot Laplace阵。
LTL = lhs' * lhs;

%drawmesh(t, uv, B);


