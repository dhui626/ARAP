
%% ARAP 2D 网格变形 每改变一次拖动点位置

local_time = 0;
global_time = 0;
solve_time = 0;

start_time = clock;

P2Pcounts = size(P2PDsts,1);
lhs_cal = LTL;
rhs = zeros (v_count, 2);

for i=1:P2Pcounts
    id = P2PVtxIds(i);
    lhs_cal(id, id) = lhs_cal(id, id) + 1;
end

dLhs = decomposition(lhs_cal);

end_time = clock;
precompute_time = etime(end_time,start_time);


for it = 1:iteration %迭代次数
    %% Local
    
    start_time = clock;
    
    %计算每个三角形对应变换的Jacobi矩阵 进行svd分解 计算并记录Lt
    Lt = zeros(2,2,f_count);

    for i=1:f_count
        index = t(i,:);
%         Jt = cot_data(index(1),index(2)) * (uv(index(1),:) - uv(index(2),:))' * (flattened(i,[1 2]) - flattened(i,[3 4]))...
%             + cot_data(index(2),index(3)) * (uv(index(2),:) - uv(index(3),:))' * (flattened(i,[3 4]) - flattened(i,[5 6]))...
%             + cot_data(index(3),index(1)) * (uv(index(3),:) - uv(index(1),:))' * (flattened(i,[5 6]) - flattened(i,[1 2]));
        
        Jt = [uv(index(1),:) - uv(index(2),:); uv(index(2),:) - uv(index(3),:)]' / [flattened(i,[1 2]) - flattened(i,[3 4]); flattened(i,[3 4]) - flattened(i,[5 6])]'; 
        
        [U,S,V] = svd(Jt);
        if ( det(Jt)>0 )
            Lt(:,:,i) = U*V';
        else
            Lt(:,:,i) = U*diag([1 -1])*V';
        end
    end
    
    end_time = clock;
    local_time = local_time + etime(end_time,start_time);
    
    %% Global
    start_time = clock;
    
    for i=1:v_count
        vj = find(vNeighber(i,:) == 1);  %找到邻接点序号
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
    
    rhs = lhs' * rhs;
    
    for n=1:P2Pcounts
        id = P2PVtxIds(n);
        rhs(id, :) = rhs(id, :) + P2PDsts(n, [1 2]);
    end
    
    end_time = clock;
    global_time = global_time + etime(end_time,start_time);
    
    start_time = clock;
    uv = dLhs\ rhs;
    end_time = clock;

    solve_time = solve_time + etime(end_time,start_time);
    
    
end

fprintf('precompute step costs %f s\n', precompute_time);
fprintf('local step costs %f s\n', local_time);
fprintf('global step costs %f s\n', global_time);
fprintf('solve step costs %f s\n', solve_time);

