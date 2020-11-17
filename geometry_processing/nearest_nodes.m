function [dist qi,pi] = nearest_nodes(q1,q2,p1,p2)
    dists_2p = ([q1 -ones(length(q1),1)]*[ones(1,length(p1)); p1']).^2 + ...
        ([q2 -ones(length(q2),1)]*[ones(1,length(p2)); p2']).^2;
    [dist,idxl] = min(dists_2p(:));
    [qi,pi] = ind2sub(size(dists_2p),idxl);
end
