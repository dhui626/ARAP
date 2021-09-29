
%% 初始化 使用Tutte方法
function xnew1 = tutte(x, t, B)

    N = size(x,1);          %总共点的个数
    n = N - size(B,2);      %内点的个数
    X = (1:N);
    X(B(:)) = [];   % B存储边界点的指标，X存储内点的指标

    %subplot(131); drawmesh(t, x, B);title('input');

    % 将边界点映到平面（圆上）
    b = x(B(:),:);  %边界点的集合
    dist = zeros(N-n,1);
    for i=1:N-n-1
        dist(i) = norm(b(i,:)-b(i+1,:),2);  
    end
    dist(N-n) = norm(b(N-n,:)-b(1,:),2);
    dist = dist/sum(dist);
    for i=1:N-n
        %b(i,:) = [cos( 2*pi*sum(dist(1:i)) ),sin( 2*pi*sum(dist(1:i)) ),0];
        b(i,:) = [cos( 2*pi*i/(N-n) ),sin( 2*pi*i/(N-n) ),0];
    end

    % 构建矩阵Wleft,Wright         %均匀权重

    Wleft = sparse(n,n);
    Wright = sparse(n,N-n);

    for i=1:n
        vj = FindNeighbor(X(i),t);
        di = size(vj,1);    %有多少邻接点

        vj_left = [];
        vj_right = [];
        for j=1:di
            vj_left = [vj_left,find(vj(j)==X)];
            vj_right = [vj_right,find(vj(j)==B)];
        end

        Wleft(i,vj_left) = -1/di;
        Wleft(i,i) = 1;
        Wright(i,vj_right) = 1/di;
    end

    % 求解稀疏方程组 Wleft * x = Wright * b

    xnew = decomposition(Wleft)\Wright*b;

    % 显示参数化结果

    xnew1 = zeros(N,2);
    xnew1(X,[1 2]) = xnew((1:n),[1 2]);
    xnew1(B,[1 2]) = b((1:N-n),[1 2]);
    %subplot(132); drawmesh(t, xnew1, B);title('initial(Tutte)');

end


%% 找到邻接点序号的函数
function vj = FindNeighbor(vi, t)
    temp = t( mod( find(t == vi)-1,size(t,1) )+1 ,:);  %找到有vi邻接点的所有行
    vj = unique(temp);  %删除重复元素
    vj( vj==vi ) = [];  %删除vi自身
end
    