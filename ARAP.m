
%% ARAP参数化
clear;
visualize = true;       %显示参数化结果
iteration = 20;         %迭代次数

%% 计算Tutte参数化结果作为初值
[x, t] = readObj('obj\Beetle_ABF');
B = findBoundary(x, t);
v_count = size(x,1);         %总共点的个数
f_count = size(t,1);          %总共边的个数

uv = tutte(x, t, B);

if (visualize)
    subplot(231); drawmesh(t, x, B);title('input');
    subplot(232); drawmesh(t, uv, B);title('initial(Tutte)');
end

%% 计算并存储cot
cot_data = sparse( v_count, v_count );
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
he_vdata = sparse( v_count, v_count );                %(i,j)元素的值是对应第三个点的指标
for i=1:f_count
    index = t(i,:);
    he_vdata(index(1), index(2)) = index(3);
    he_vdata(index(2), index(3)) = index(1);
    he_vdata(index(3), index(1)) = index(2);
end

%% 记录(i,j)对应的面的指标（半边数据结构）
he_fdata = sparse( v_count, v_count );                %(i,j)元素的值是对应面的指标
for i=1:f_count
    index = t(i,:);
    he_fdata(index(1), index(2)) = i;
    he_fdata(index(2), index(3)) = i;
    he_fdata(index(3), index(1)) = i;
end

%% 生成稀疏矩阵并预分解
lhs = sparse(v_count, v_count);

lhs(1,1) = 1; %固定第一个点的uv坐标到(0,0)
for i=2:v_count
    lhs(i,:) = -( cot_data(i,:) + cot_data(:,i)' );
    lhs(i,i) = -sum(lhs(i,:),2);
end

dLhs = decomposition(lhs);
%lhs_inv = pinv(lhs);

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

for it = 1:iteration %迭代次数
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
    for i=2:v_count         %第一行是由固定点构成的方程
        vj = FindNeighbor(i, t);
        vj_size = size(vj,1);
        sum = [0; 0];
        
        for m=1:vj_size
            f_index = he_fdata( i, vj(m) );
            if (f_index ~=0)
                if (t( f_index , 1) == i)
                    sum = sum + cot_data( i,vj(m) ) * Lt(:,:,f_index) * ( flattened(f_index,[1 2]) - flattened(f_index,[3 4]) )';
                elseif (t( f_index, 2) == i)
                    sum = sum + cot_data( i,vj(m) ) * Lt(:,:,f_index) * ( flattened(f_index,[3 4]) - flattened(f_index,[5 6]) )';
                elseif (t( f_index, 3) == i)
                    sum = sum + cot_data( i,vj(m) ) * Lt(:,:,f_index) * ( flattened(f_index,[5 6]) - flattened(f_index,[1 2]) )';
                end
            end
            
            f_index = he_fdata( vj(m), i );
            if (f_index ~=0)
                if (t( f_index, 1) == i)
                    sum = sum + cot_data( vj(m), i ) * Lt(:,:,f_index) * ( flattened(f_index,[1 2]) - flattened(f_index,[5 6]) )';
                elseif (t( f_index, 2) == i)
                    sum = sum + cot_data( vj(m), i ) * Lt(:,:,f_index) * ( flattened(f_index,[3 4]) - flattened(f_index,[1 2]) )';
                elseif (t( f_index, 3) == i)
                    sum = sum + cot_data( vj(m), i ) * Lt(:,:,f_index) * ( flattened(f_index,[5 6]) - flattened(f_index,[3 4]) )';
                end
            end
            
            rhs(i,:) = sum';
        end
    end

    uv = dLhs\ rhs;
    
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


%% 找到邻接点序号的函数
function vj = FindNeighbor(vi, t)
    temp = t( mod( find(t == vi)-1,size(t,1) )+1 ,:);  %找到有vi邻接点的所有行
    vj = unique(temp);  %删除重复元素
    vj( vj==vi ) = [];  %删除vi自身
end

