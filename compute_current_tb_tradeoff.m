 %%
clear
close all

gamma = 42.57747892e6; %%Hz/T

ex_band = 2000; % 1000 - 2000 Hz

down_rate = 4;

coil_mask = niftiread('040716_0012-label.nii.gz');

coil_mask = coil_mask(1:down_rate:end,1:down_rate:end,1:down_rate:end);

%coil_mask = padarray(coil_mask,[6 6 6],0,'both');

coil_m = zeros(size(coil_mask));
% 
coil_m(2:end -2,2:end-2,2:end-2) = coil_mask(2:end -2,2:end-2,2:end-2);


[m,n,p] = size(coil_mask);
[Xv,Yv,Zv] = meshgrid(1:n,1:m,1:p); 


D = smooth3(coil_m,'box',1);
mask_union_boundary_t = isosurface(Xv,Yv,Zv,D,0.5);
top_trimesh = triangulation(mask_union_boundary_t.faces,mask_union_boundary_t.vertices);

top_trimesh = triangulation(mask_union_boundary_t.faces,mask_union_boundary_t.vertices-2*top_trimesh.vertexNormal);


%%
mask_brain = niftiread('MAP-C404_0012.nii.gz');

mask_brain = mask_brain(1:down_rate:end,1:down_rate:end,1:down_rate:end);

move_lib = [-5,3,0; -5,10,0; -5,17,0; -5,3,2; -5,10,2; -5,17,2; -5,3,-2; -5,10,-2; -5,17,-2] * 4 / down_rate;
shift =move_lib(1,:) ;%[-5,7,0];%[0,0,5]; % move to right
s = size(mask_brain);
mask_shift = padarray(mask_brain,[20 20 20]* 4 /down_rate,0,'both');

mask_brain = mask_shift(shift(1)+ 1 + 20 * 4 /down_rate:shift(1)+20* 4 /down_rate+s(1),shift(2)+1 + 20* 4 /down_rate :shift(2)+20* 4 /down_rate ...
+s(2),shift(3)+1+ 20* 4 /down_rate :shift(3)+20* 4 /down_rate+s(3));
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
D = int16(D);


mask_brain_projection1 = sum(mask_brain==1,1)>0;
mask_brain_projection2 = sum(mask_brain==1,2)>0;
mask_brain_projection3 = sum(mask_brain==1,3)>0;

% 
mask_brain_p_extrude1 = repmat(mask_brain_projection1,[size(mask_brain_projection1,1),1,1]);
mask_brain_p_extrude2 = repmat(mask_brain_projection2,[1,size(mask_brain_projection2,2),1]);
mask_brain_p_extrude3 = repmat(mask_brain_projection3,[1,1,size(mask_brain_projection3,3)]);

mb_extend1 = mask_brain;
mb_extend1((mask_brain_p_extrude1==1)&(mask_brain==3)) = 2 ;
mb_extend1 = mb_extend1.*D;

mb_extend2 = mask_brain;
mb_extend2((mask_brain_p_extrude2==1)&(mask_brain==3)) = 2 ;
mb_extend2 = mb_extend2.*D;

mb_extend3 = mask_brain;
mb_extend3((mask_brain_p_extrude3==1)&(mask_brain==3)) = 2 ;
mb_extend3 = mb_extend3.*D;


mb_revise1 = zeros(size(mask_brain));
mb_revise1(x_min_m:x_max_m,:,:) = 3;
mb_revise1((mask_brain==1)) = 1 ;
mb_revise1(mb_revise1 == 0) = 2 ;
mb_revise1((mb_revise1 == 3) & (mb_extend1 == 2)) = 2;
mb_revise1 = int16(mb_revise1).*D;
% 
mb_revise2 = zeros(size(mask_brain));
mb_revise2(:,y_min_m:y_max_m,:) = 3;
mb_revise2((mask_brain==1)) = 1 ;
mb_revise2(mb_revise2 == 0) = 2 ;
mb_revise2((mb_revise2 == 3) & (mb_extend2 == 2)) = 2;
mb_revise2 = int16(mb_revise2).*D;

mb_revise3 = zeros(size(mask_brain));
mb_revise3(:,:,z_min_m:z_max_m) = 3;
mb_revise3((mask_brain==1)) = 1 ;
mb_revise3(mb_revise3 == 0) = 2 ;
mb_revise3((mb_revise3 == 3) & (mb_extend3 == 2)) = 2;
mb_revise3 = int16(mb_revise3).*D;

mask_brain = mask_brain.*D;



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
ptstop = [ptstop; side] * 4 / down_rate;    

side1 = side + [5,0,0];
ptstop = [ptstop; side1];
side1 = side + [0,3,0];
ptstop = [ptstop; side1];

coil_loop.segments = {};


for i = 1:size(ptstop,1)
    [~,I] = min(vecnorm((top_trimesh.Points-ptstop(i,:)),2,2)); % caluculate the closet distance of the triangle point and chose point.
    n = top_trimesh.vertexNormal(I);
    p = top_trimesh.Points(I,:);
    phi = acos(n(3)/norm(n,2))*180/pi;
    theta = atan2(n(2),n(1))*180/pi;

    circle = [p(1),p(2),p(3),2 * 4 / down_rate,phi,theta,0];
    c = circle3dPoint(circle,linspace(0,360,50));
    coil_loop.segments{i} = c';
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
mask_brain_slice = mask_brain;
% mask_brain_slice(:,:,[1:8,12:20]) = 3*D(:,:,[1:8,12:20]);
% mask_brain_slice([1:17,20:end],:,:) = 3*D([1:17,20:end],:,:);

 
idxs = find(mask_brain_slice(:,:,:)==1);
[I,J,K] = ind2sub(size(mask_brain_slice),idxs);
      
brainbounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];
brainbounds = brainbounds+[-1,1
                           -0,0
                           -2,2]; 
                       
mb_extend = zeros(size(mask_brain_slice));
mb_extend(D==1) = 2;
mb_extend(:,brainbounds(2,1):brainbounds(2,2),:) = 3;
mb_extend(brainbounds(1,1):brainbounds(1,2),brainbounds(2,1):brainbounds(2,2),brainbounds(3,1):brainbounds(3,2)) = 2;
mb_extend((mask_brain_slice==1)) = 1 ;
mb_extend = int16(mb_extend).*D;

mb_extend_con = zeros(size(mask_brain_slice));
mb_extend_con(D==1) = 2;
mb_extend_con(:,brainbounds(2,1):brainbounds(2,2),:) = 3;
mb_extend_con((mask_brain_slice==1)) = 1 ;
mb_extend_con = int16(mb_extend_con).*D;


unshimmed = zeros(size(Xv));

tic;
p = optim_shim_ncc_fetal_pulse_param(unshimmed,mb_extend,reshape(coil_loop.b,[size(Xv),size(ptstop,1)]));
toc;
%%

tbs = logspace(log10(50),log10(400),10);
bws = logspace(log10(1000),log10(3000),10);

[X,Y] = ndgrid(tbs,bws);

res = [];
for tb = tbs
   p2 = p(tb,[]);
   res = [res;p2(bws)];
end
disp('done');
%%
figure();
plot(tbs,res/3.5,'LineWidth',3);
xlabel('tb');
ylabel('number of turns requried')
legend(string(bws))
set(gca,'fontsize',22)
ylim([0,500]);
xlim([50,400]);

figure();
imagesc(mb_extend(:,:,10))

%%
figure(2);
surf(X,Y,res/3.5)
axis vis3d
xlabel('transition bandwidth');
ylabel('excitation bandwidth');
zlabel('current');

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
