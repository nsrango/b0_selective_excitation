 %% Depends on gptoolbox, rf_tools, and matGeom (could be removed....)

clear
close all

gamma = 42.57747892e6; %%Hz/T

ex_band = 2000; % 1000 - 2000 Hz

down_rate = 4;

coil_mask = niftiread('data/040716_0012-label.nii.gz');

coil_mask = coil_mask(1:down_rate:end,1:down_rate:end,1:down_rate:end);

%coil_mask = padarray(coil_mask,[6 6 6],0,'both');

coil_m = zeros(size(coil_mask));
% 
coil_m(2:end -2,2:end-2,2:end-2) = coil_mask(2:end -2,2:end-2,2:end-2);


[m,n,p] = size(coil_mask);
[Xv,Yv,Zv] = meshgrid(1:n,1:m,1:p); 


D = smooth3(coil_m,'box',5);
mask_union_boundary_t = isosurface(Xv,Yv,Zv,D,0.5);
top_trimesh = triangulation(mask_union_boundary_t.faces,mask_union_boundary_t.vertices);

top_trimesh = triangulation(mask_union_boundary_t.faces,mask_union_boundary_t.vertices-2*top_trimesh.vertexNormal);


U = laplacian_smooth(mask_union_boundary_t.vertices,mask_union_boundary_t.faces,'uniform',[],.001);
top_trimesh = triangulation(mask_union_boundary_t.faces,U);


%%
V = top_trimesh.Points;
F = top_trimesh.ConnectivityList;
% Extract offset at minus 3% of bounind box diagonal length
iso = -0.05;
% Resolution grid → resolution of extracted offset surface
side = 60;
% Amount of smoothing to apply to distance field
sigma = 5;
bbd = norm(max(V)-min(V));
% Pad generously for positive iso values 
[BC,side,r] = voxel_grid([V;max(V)+10;min(V)-10],side);
D = signed_distance(BC,V,F);
D = reshape(D,side([2 1 3]));
% Smooth signed distance field
D = imfilter(D,fspecial('gaussian',9,sigma),'replicate');
BC3 = reshape(BC,[side([2 1 3]) 3]);
% Use matlab's built-in marching cubes iso-surface mesher (works most of the time)
surf = isosurface(BC3(:,:,:,1),BC3(:,:,:,2),BC3(:,:,:,3),D,iso*bbd);
SV = surf.vertices;
SF = surf.faces;
SVp = laplacian_smooth(SV,SF,'uniform',[],.001);


[W,G] = decimate_libigl(SVp,SF,0.2);
cut = W(:,1)<2 | W(:,1)>28;
cut_idx = find(cut);

cut_T = any(any(G==permute(cut_idx,[3,2,1]),3),2);
G(cut_T,:) = [];

[RV,IM,J,IMF] = remove_unreferenced(W,G);

top_trimesh = triangulation(IMF,RV);
edge = edges(IMF);
boundary_edges_mask = is_boundary_facet(edge,IMF);
%%
boundary_edges=edge(boundary_edges_mask,:);
boundary_edges_mutable = boundary_edges;

boundary_edges_grouped = {};
while(any(any(boundary_edges_mutable~=0,1),2))
start = [find(boundary_edges_mutable(:,1)~=0,1),1];
idx = [start(1),2];
first = true;
edge_connection_list = [];
while(idx(1)~=start(1) || first)
    first = false;
    [r,c] = find(boundary_edges_mutable==boundary_edges_mutable(idx(1),idx(2)));
    potential_list = [r,c];
    [~,row_match] = ismember(idx,potential_list,'rows');
    potential_list(row_match,:) = [];
    idx=potential_list;
    idx(2) = 3-idx(2);
    edge_connection_list = [edge_connection_list;idx(1)];
end
boundary_edges_mutable(edge_connection_list,:) = 0;
boundary_edges_grouped{end+1} = edge_connection_list;
end

boudary_nodes_grouped = {};
for boundary_edges_group = boundary_edges_grouped
    boudary_nodes_grouped{end+1} = unique(boundary_edges(boundary_edges_group{1},:));
end

n_boundary_nodes = sum(cellfun(@(x)(size(x,1)),boudary_nodes_grouped));
n_boundaries = numel(boudary_nodes_grouped);

%%
boundary_matrix = zeros(size(RV,1),size(RV,1)+n_boundaries);
interior_V = setdiff(1:size(RV,1),vertcat(boudary_nodes_grouped{:}));


boundary_matrix(sub2ind(size(boundary_matrix),interior_V,interior_V)) = 1;
for idx = numel(boundary_edges_grouped)
    boundary_matrix(boudary_nodes_grouped{idx},size(RV,1)+idx) = 0;
end
boundary_matrix( :, ~any(boundary_matrix,1) ) = [];

%%
[L,areas_v,areas] = cotLaplacian(top_trimesh.Points,top_trimesh.ConnectivityList);  
areas_v = 3*areas_v;

T = construct_unitTransforms(top_trimesh);
[x0,w0,~] = strang9();
np = size(x0,1);
x0 = [x0,zeros(np,1),ones(np,1)];

pts = zeros(np,6,size(top_trimesh.Points,3));
for i = 1:size(T,3)
    tmpy = x0*T(:,:,i)';
    pts(:,1:3,i) = 1e-3*tmpy(:,1:3);
    pts(:,4:6,i) = repmat(top_trimesh.faceNormal(i),np,1);
end

bfull2 = zeros(numel(Xv),top_trimesh.size(1),3);
tic;
for i = 1:size(bfull2,2)
    b = surface_basis2(pts(:,:,i),Xv*1e-3,Yv*1e-3,Zv*1e-3);
    w = 1e-6*areas(i)*w0;
    v = [x0(:,1),x0(:,2),(1-x0(:,1)).*(1-x0(:,2))];
    bfull2(:,i,:)= permute(b*(w.*v),[1,3,2]);
end
toc;
%%
b = zeros(numel(Xv),size(top_trimesh.Points,1));
for i = 1:top_trimesh.size(1)
    for j = 1:3
        b(:,top_trimesh.ConnectivityList(i,j)) = b(:,top_trimesh.ConnectivityList(i,j))+bfull2(:,i,j);
    end
end
%%


%%
% top_trimesh = triangulation(t_helmet_d.faces,t_helmet_d.vertices);
%%


figure(1)
p1 = trimesh(top_trimesh);
axis vis3d

p1.FaceColor = 'red';
p1.EdgeColor = 'none';
p1.FaceLighting = 'gouraud';
p1.FaceAlpha = '0.1';
camlight headlight
hold on

qq = quiver3(m/2,n/2,p/2,-15* 4 /down_rate,0,0,0);
qq.LineWidth = 5;
text(0,n/2+1,p/2+1,'Superior/Inferior','Color','black','FontSize',16);
qq = quiver3(m/2,n/2,p/2,0,-15* 4 /down_rate,0,0);
qq.LineWidth = 5;
text(m/2+1,0,p/2+1,'Front/Back','Color','red','FontSize',16);
qq = quiver3(m/2,n/2,p/2,0,0,10* 4 /down_rate,0);
qq.LineWidth = 5;
text(m/2+1,n/2+1,20 * 4 /down_rate ,'Left/Right','Color','blue','FontSize',16);


%%
% mask_brain = niftiread('MAP-C404_0012.nii.gz');
% 
% mask_brain = mask_brain(1:down_rate:end,1:down_rate:end,1:down_rate:end);
% 
% 
% move_lib = [-5,3,0; -5,10,0; -5,17,0; -5,3,2; -5,10,2; -5,17,2; -5,3,-2; -5,10,-2; -5,17,-2] * 4 / down_rate;
% shift =move_lib(2,:) ;%[-5,7,0];%[0,0,5]; % move to right
% %shift = [-5,7,0];
% s = size(mask_brain);
% mask_shift = padarray(mask_brain,[20 20 20]* 4 /down_rate,0,'both');
% 
% mask_brain = mask_shift(shift(1)+ 1 + 20 * 4 /down_rate:shift(1)+20* 4 /down_rate+s(1),shift(2)+1 + 20* 4 /down_rate :shift(2)+20* 4 /down_rate ...
% +s(2),shift(3)+1+ 20* 4 /down_rate :shift(3)+20* 4 /down_rate+s(3));

mask_brain = zeros(size(D));

mask_brain = double(sqrt(((Xv-15).^2+(Yv-15).^2+(Zv/10*15-15).^2)) <3.3);



D1 = smooth3(mask_brain,'box',1);

ccc = isosurface(Xv,Yv,Zv,D1,0.5);
head = triangulation(ccc.faces,ccc.vertices);

p2 = trimesh(head);
p2.FaceColor = 'green';
p2.EdgeColor = 'none';
p2.FaceLighting = 'gouraud';
p2.FaceAlpha = '0.3';
camlight headlight
hold on

mask_brain(mask_brain == 0) = 3;

%%%%
[x_mask, y_mask, z_mask] = ind2sub(size(mask_brain), find(mask_brain ==1));
x_min_m = min(x_mask);
x_max_m = max(x_mask);
y_min_m = min(y_mask);
y_max_m = max(y_mask);
z_min_m = min(z_mask);
z_max_m = max(z_mask);
%%%%

D = smooth3(coil_m,'box',1);
D(D<0.99) = 0  ; D(D>0.89) = 1;
D = double(D);

mask_brain = mask_brain.*D;


%%
% for i = 1:size(mask_brain,3)
% imshow(squeeze(mb_revise3(:,:,10)),[])
% pause(2)
% end

%%

ptstop = [
          6,9.18,15
          9,7.62,15
          11.4,6.58,15
          14,4.978,15
          17,3.42,15
          19,2.78,15
          22,2.22,15
          25,2.02,15
          
          6,8.83,11
          9,7.18,11
          11.23,6,11
          14,4.82,11
          17,3.86,11
          19,3.42,11
          22,3.1,11
          25,3.457,11
          
          6,9.717,8
          9,8.196,8
          11,7.27,8
          14,6.312,8
          17,5.625,8
          19,5.275,8
          22,5.375,8
          25,6.1,8
          
          6,11.98,5
          9,10.56,5
          11,9.94,5
          14,9.477,5
          17,9.225,5
          19,9.075,5
          22,9.025,5
          25,9.325,5
          %%
          %%
                ];
            
            
ptstop2 = ptstop + [0,30,0];

ptstop = [ptstop ; ptstop2];

side = [
        23,11,18.5
        16,11,18.5
        8,11,18.5
        23,16,18.5
        16,16,18.5
        8,16,18.5
        23,21,18.5
        16,21,18.5
        8,21,18.5       
        23,26,18.5
        16,26,18.5
        8,26,18.5   
        %%%%
        22,14,2.18
        22,19,2.14
        22,24,2.5
        16,14,1.5
        16,19,1.5
        16,24,2.02
        10,14,1.62
        10,19,1.5
        10,23.94,2
        ];
%ptstop = side;
ptstop = [ptstop; side] * 4 / down_rate;    

 side1 = side + [5,0,0];
 ptstop = [ptstop; side1];
 side1 = side + [0,3,0];
 ptstop = [ptstop; side1];

% ptstop1 = ptstop + [5,0,0];
% ptstop = [ptstop; ptstop1];

% ptstop = [];
% 
% for i = 3:3:27%[2,5,10,19,24,29]
%     for j = 1:3:30%[5,15,22,28]
%         for z = 1:4:20
%             ptstop = cat(2,ptstop,[i,j,z]');
%         end
%     end
% end
% ptstop = ptstop';

coil_loop.segments = {};


for i = 1:size(ptstop,1)
    [~,I] = min(vecnorm((top_trimesh.Points-ptstop(i,:)),2,2)); % caluculate the closet distance of the triangle point and chose point.
    n = top_trimesh.vertexNormal(I);
    p = top_trimesh.Points(I,:);
    phi = acos(n(3)/norm(n,2))*180/pi;
    theta = atan2(n(2),n(1))*180/pi;

    %circle = [p(1),p(2),p(3),47.5,phi,theta,0];
    circle = [p(1),p(2),p(3),2 * 4 / down_rate,phi,theta,0];
    c = circle3dPoint(circle,linspace(0,360,50));
    coil_loop.segments{i} = c';
    %plot3(c(:,1),c(:,2),c(:,3),'LineWidth',5); hold on;grid on;
end


segmentsI = cellfun(@(x)x*3 * 4*1e-3 * 4 / down_rate,coil_loop.segments,'UniformOutput',false);
coil_loop.b = zeros(numel(Xv),length(segmentsI));
for i = 1:length(segmentsI)
        fields = integrate_wirepath2(...
            {segmentsI{i}},...
            Xv*3 * 4 *1e-3 * 4 / down_rate,Yv*3* 4*1e-3 * 4 / down_rate,Zv*3* 4*1e-3 * 4 / down_rate,1); % shows the field points volume
        coil_loop.b(:,i) = gamma*fields(3,:); %hz
end

%%
figure(3);
trisurf(top_trimesh,'FaceColor','none'); axis equal; hold on;
trisurf(head,'FaceColor','none'); axis equal; hold on;
%h = slice(Xv,Yv,Zv,reshape(gamma*fields(3,:),size(Xv)),15 * 4 /down_rate,15 * 4 /down_rate,10 * 4 /down_rate);
%set(h,'edgecolor','none');
%set(gca,'CLim',[-50,50]);
for i = 1:length(coil_loop.segments)
ct = coil_loop.segments{i}';
plot3(ct(:,1),ct(:,2),ct(:,3),'LineWidth',5); hold on; axis equal;
end
[m,n,p] = size(coil_mask);
qq = quiver3(m/2,n/2,p/2,-20* 4 /down_rate,0,0,0);
qq.LineWidth = 5;
text(-15,n/2+1,p/2+1,'Superior/Inferior','Color','red','FontSize',16);
qq = quiver3(m/2,n/2 ,p/2,0,-20* 4 /down_rate,0,0);
qq.LineWidth = 5;
text(m/2+1,0,p/2+1,'Front/Back','Color','red','FontSize',16);
qq = quiver3(m/2,n/2,p/2,0,0,15* 4 /down_rate,0);
qq.LineWidth = 5;
text(m/2+1,n/2+1,20 * 4 /down_rate ,'Left/Right','Color','white','FontSize',16);
hold off;

%%
idxs = find(mask_brain(:,:,:)==1);
[I,J,K] = ind2sub(size(mask_brain),idxs);
      
brainbounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];

brainbounds_e = brainbounds+[-1,1
                           -1,1
                           -1,1];
      
%%      
mb_extend = zeros(size(mask_brain));
mb_extend(D==1) = 2;
mb_extend(brainbounds(1,1):brainbounds(1,2),1:brainbounds(2,1),:) = 3;
mb_extend(brainbounds(1,1):brainbounds(1,2),brainbounds(2,2):end,:) = 3;

mb_extend(brainbounds(1,1):brainbounds(1,2),brainbounds(2,1):brainbounds(2,2),brainbounds(3,1):brainbounds(3,2)) = 2;
mb_extend((mask_brain==1)) = 1 ;
mb_extend = (mb_extend).*D;

mb_extend_con = zeros(size(mask_brain));
mb_extend_con(D==1) = 2;
mb_extend_con(brainbounds(1,1):brainbounds(1,2),:,:) = 3;
mb_extend_con((mask_brain==1)) = 1 ;
mb_extend_con = (mb_extend_con).*D;


unshimmed = zeros(size(Xv));

%%
mb_patch = zeros(size(mask_brain));
mb_patch(D==1) = 2;
ds = diff(brainbounds_e,1,2)+1;
val = 3;
for i = 0
    for j = -2:2
        for k = -2:2
            if k + j <= 0
                val = 3;
            else 
                val = 3;
            end
            if j>=0
                tmp = padarray(mask_brain(:,1:(end-(j*ds(2))),:)==1,[0,min(j*ds(2),30),0],'pre');
            else
                tmp = padarray(mask_brain(:,1-(j*ds(2)):end,:)==1,[0,max(-j*ds(2),0),0],'post');
            end
            
            if k>=0
                tmp = padarray(tmp(:,:,1:(end-(k*ds(3))))==1,[0,0,min(k*ds(3),20)],'pre');
            else
                tmp = padarray(tmp(:,:,1-(k*ds(3)):end)==1,[0,0,max(-k*ds(3),0)],'post');
            end
             
            mb_patch(tmp==1) = val;
            imagesc(squeeze(mb_patch(15,:,:)))

                      
%             mb_patch(min(max((brainbounds(1,1):brainbounds(1,2))+i*ds(1),1),30),...
%                 min(max((brainbounds(2,1):brainbounds(2,2))+j*ds(2),1),30),...
%                 min(max((brainbounds(3,1):brainbounds(3,2))+k*ds(3),1),20)) = val;
%             val = val+1;
            
        end
    end
end

mb_patch(brainbounds(1,1):brainbounds(1,2),brainbounds(2,1):brainbounds(2,2),brainbounds(3,1):brainbounds(3,2)) = 2;


mb_patch(mask_brain==1) = 1;
mb_patch(D==0) = 0;


for s = 1:15
    imagesc([mb_patch(:,:,s) ; mb_extend(:,:,s)]); colorcet('L1'); axis equal;
    pause(0.1);
end

%%
% tic;
% [shimmed,amps,db0,mask_used,filt] ...
%            = perform_shim_ncc_fetal(unshimmed,mb_patch,reshape(b*boundary_matrix,[size(Xv),size(b*boundary_matrix,2)]),...
%                                            50,200000,{100,1e7,90/.66,0,0,0.999, 1000/.6, inf});
% toc;

tic;
[shimmed,amps,db0,mask_used,filt] ...
            = perform_shim_ncc_2s_fetal_yalmip(unshimmed,mb_patch,reshape(b*boundary_matrix,[size(Xv),size(boundary_matrix,2)]),...
                                            50,20000,90/.66,1000/.66);
toc;

% tic;
% [shimmed,amps,db0,mask_used,filt] ...
%             = perform_shim_ncc_2s_fetal_yalmip_binary(unshimmed,mb_patch,reshape(b*boundary_matrix,[size(Xv),size(b*boundary_matrix,2)]),...
%                                             200,2000000000,90/.66,1000/.66);
% toc;

rf2 = dzrf(500, 12, 'ex', 'pm',0.0006, 0.012);
x = -20 *3.5:0.01:20 *3.5;
ab2  = abr(rf2, x);

CK_param = interp1(toHz(x,8)+760+mean(filt),ab2ex(ab2),shimmed(:),'linear',0);



%%
ex = reshape(CK_param,size(shimmed));
f1 = figure(1); axis off;
a1 = gca;
f2 = figure(2); axis off;
a2 = gca;
f3 = figure(3); axis off;
a3 = gca;
ex2 = ex;
masku = mask_brain.*D;
ex2(masku(:,:,:)==0) = nan;


idxs = find(((mb_patch(:,:,:)==1 | mb_patch(:,:,:)==3).*ex2)>0.02);
[I,J,K] = ind2sub(size(ex2),idxs);

bounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];
      
idxs = find(mask_brain(:,:,:)==1);
[I,J,K] = ind2sub(size(ex2),idxs);
      
brainbounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];
      
idxs = find(D(:,:,:)~=0);
[I,J,K] = ind2sub(size(ex2),idxs);
      
bodybounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];      

boxmask = zeros(size(mask_brain));
boxmask(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) = 1;
bodybox = zeros(size(mask_brain));
bodybox(bodybounds(1,1):bodybounds(1,2),bodybounds(2,1):bodybounds(2,2),bodybounds(3,1):bodybounds(3,2)) = 1;


for s = 10
    
    [B,L] = bwboundaries(mask_brain(:,:,s)==1,'noholes');
    [B2,L2] = bwboundaries(boxmask(:,:,s)==1,'noholes');
    [B3,L3] = bwboundaries(bodybox(:,:,s)==1,'noholes');

    imagesc((mask_brain(:,:,s)~=0).*abs(ex2(:,:,s)),'parent',a1,[0,1]);
        axis off;
    set(a1,'Colormap',colorcet('L2'))
    hold(a1,'on');
    for k = 1:length(B)
       boundary = B{k};
       plot(a1,boundary(:,2), boundary(:,1), 'red', 'LineWidth', 2); 
    end
    for k = 1:length(B2)
       boundary = B2{k};
       plot(a1,boundary(:,2), boundary(:,1), 'green', 'LineWidth', 2);
    end
    for k = 1:length(B3)
       boundary = B3{k};
       plot(a1,boundary(:,2), boundary(:,1), 'blue', 'LineWidth', 2);
    end
    hold(a1, 'off');
   
    imagesc((D(:,:,s)~=0).*shimmed(:,:,s),'parent',a2,3000*[-1,1]+mean(filt)+500);
    set(a2,'Colormap',colorcet('D1A'))
    axis off;
    
    imagesc(mb_patch(:,:,s),'parent',a3);
    set(a3,'Colormap',colorcet('L1'))
    axis off;
    
    
    pause(1);
end

%%
figure(4);
patch('vertices',RV,'Faces',IMF,'FaceColor','interp','CData',10*amps'*boundary_matrix','edgecolor','black'); axis vis3d;
colorcet('D1A');
set(gca,'Clim', (300*[-1,1]))
hold on;
ax = gca;
% 
for i = 15
h = slice(Xv,Yv,Zv,double(mb_patch)*100-300,[i],[15],[10]);
set(gca,'Clim', (300*[-1,1]))

pause(0.5)
% delete(h);
end
hold off;

figure(5);
patch('vertices',RV,'Faces',IMF,'FaceColor','interp','CData',4*amps'*boundary_matrix','edgecolor','black'); axis vis3d;
colorcet('D1A');
set(gca,'Clim', (300*[-1,1]))
hold on;
ax = gca;
% 
for i = 15
h = slice(Xv,Yv,Zv,(mask_brain~=0).*abs(ex2)*300,[i],[15],[10]);
set(gca,'Clim', (300*[-1,1]))

pause(0.5)
% delete(h);
end
hold off;

figure(6);
patch('vertices',RV,'Faces',IMF,'FaceColor','interp','CData',10*amps'*boundary_matrix','edgecolor','black'); axis vis3d;
colorcet('D1A');
set(gca,'Clim', (300*[-1,1]))
hold on;
ax = gca;
% 
for i = 15
h = slice(Xv,Yv,Zv,(mask_brain~=0).*(shimmed-mean(filt))/5,[i],[15],[10]);
set(gca,'Clim', (300*[-1,1]))

pause(0.5)
% delete(h);
end
hold off;



%%
figure(5);
ax = gca;
color = colorcet('D1','N',256);
for i = 1:length(coil_loop.segments)
ct = coil_loop.segments{i}';
fill3(ct(:,1),ct(:,2),ct(:,3),color(round(max(min((amps(i)/50+0.5),1),0)*255+1),:)); hold on; axis equal;
end

for i = 1:30
h = slice(Xv,Yv,Zv,shimmed,[i],[],[10]);
ax.CLim = 2000*[-1,1]+mean(filt);
colorcet('D1A');
pause(0.5)
delete(h);
end
hold off;

%%



function map = mapfunc(shim, center, ref)
    shim = shim - center ;
    xxx = floor((size(ref,1)-1)/2 /(2*pi*1000/3)  * shim);
    if xxx < 1
        map = 0;
    elseif xxx > size(ref)
        map = 0;
    else
        map = ref(xxx);

    end
end

function flip = flipb0(b0,filt)

    flip = (b0-filt(2))*1/(filt(1)-filt(2));
    flip(flip>1) = 1;
    flip(flip<0) = 0;
    
end

function H = toHz(x, duration)
    H = 1000 * x/duration; 
end
