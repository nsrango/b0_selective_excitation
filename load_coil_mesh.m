function [coil] = load_coil(coil_mesh,Xv,Yv,Zv)
    gamma = 42.57747892e6; %%Hz/T

    [L,areas_v,areas] = cotLaplacian(coil_mesh.Points,coil_mesh.ConnectivityList);  
    areas_v = 3*areas_v;
    rk = zeros(size(coil_mesh.Points,1),6);
    rk(:,1:3) = coil_mesh.Points(:,:);
    rk(:,4:6) = coil_mesh.vertexNormal;
    
    
    loopBasisX = sparse(coil_mesh.size(1),size(coil_mesh.Points,1));
    loopBasisY = sparse(coil_mesh.size(1),size(coil_mesh.Points,1));
    loopBasisZ = sparse(coil_mesh.size(1),size(coil_mesh.Points,1));
    
    E = coil_mesh.edges;
    boundary_vertex = [];
    for i = 1:length(E)
        vu = E(i,:);
        tu = cell2mat(coil_mesh.edgeAttachments(vu(1),vu(2)));
        if numel(tu) == 1
            boundary_vertex = [boundary_vertex,vu];
        elseif numel(tu) ~= 1
            ns = coil_mesh.faceNormal(tu');
            ev = diff(coil_mesh.Points(vu,:),1,1);
            ev = ev./norm(ev,2);
            as = areas(tu);
            fs = cross(ns,[ev;ev])./as;
            loopBasisX(tu,vu(1)) = loopBasisX(tu,vu(1))-fs(:,1);
            loopBasisY(tu,vu(1)) = loopBasisY(tu,vu(1))-fs(:,2);
            loopBasisZ(tu,vu(1)) = loopBasisZ(tu,vu(1))-fs(:,3);
            loopBasisX(tu,vu(2)) = loopBasisX(tu,vu(2))+fs(:,1);
            loopBasisY(tu,vu(2)) = loopBasisY(tu,vu(2))+fs(:,2);
            loopBasisZ(tu,vu(2)) = loopBasisZ(tu,vu(2))+fs(:,3);
        end
    end
    boundary_vertex = unique(boundary_vertex);
    loopBasis = reshape(ndSparse([loopBasisX, loopBasisY, loopBasisZ]),coil_mesh.size(1),[],3);
    [L,areas_v,areas] = cotLaplacian(coil_mesh.Points,coil_mesh.ConnectivityList);  
    
    rkp = rk;
    rkp(:,1:3) = rkp(:,1:3)*1e-3;
    
    
    
    %%
    T = construct_unitTransforms(coil_mesh);
    [x0,w0,~] = strang9();
    np = size(x0,1);
    x0 = [x0,zeros(np,1),ones(np,1)];
    %%
    pts = zeros(np,6,size(coil_mesh.Points,3));
    for i = 1:size(T,3)
        tmpy = x0*T(:,:,i)';
        pts(:,1:3,i) = 1e-3*tmpy(:,1:3);
        pts(:,4:6,i) = repmat(coil_mesh.faceNormal(i),np,1);
    end
    %%
    bfull2 = zeros(numel(Xv),coil_mesh.size(1),3);
    tic;
    for i = 1:size(bfull2,2)
        b = surface_basis2(pts(:,:,i),Xv*1e-3,Yv*1e-3,Zv*1e-3);
        w = 1e-6*areas(i)*w0;
        v = [x0(:,1),x0(:,2),(1-x0(:,1)).*(1-x0(:,2))];
        bfull2(:,i,:)= permute(b*(w.*v),[1,3,2]);
    end
    toc;
    %%
    b = zeros(numel(Xv),size(coil_mesh.Points,1));
    for i = 1:coil_mesh.size(1)
        for j = 1:3
            b(:,coil_mesh.ConnectivityList(i,j)) = b(:,coil_mesh.ConnectivityList(i,j))+bfull2(:,i,j);
        end
    end
   
    coil.mesh = coil_mesh;
    coil.b = gamma/(4*pi)*1.256637e-6*b;
    coil.areas = areas;
    coil.laplacian = L;
%     coil.boundary_vertex = setdiff(boundary_vertex,unique(reshape(coil_mesh.ConnectivityList(msh.TRIANGLES(:,4)==2,:),1,[])));
    coil.loopBasis = loopBasis;
%     coil.nLabel = unique(reshape(coil_mesh.ConnectivityList(msh.TRIANGLES(:,4)==2,:),1,[]));
end

