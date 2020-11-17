 %%
clear
close all

gamma = 42.57747892e6; %%Hz/T

ex_band = 2000; % 1000 - 2000 Hz

down_rate = 4;

coil_mask = niftiread('D:\Molin\fetal\shim_coil_zoom\040716_0012-label.nii.gz');

coil_mask = coil_mask(1:down_rate:end,1:down_rate:end,1:down_rate:end);


coil_m = zeros(size(coil_mask));
% 
coil_m(2:end -2,2:end-2,2:end-2) = coil_mask(2:end -2,2:end-2,2:end-2);


[m,n,p] = size(coil_mask);
[Xv,Yv,Zv] = meshgrid(1:n,1:m,1:p); 
%[Xv,Yv,Zv] = meshgrid(1:n+4,1:m+4,1:p+4); 
%coil_m = padarray(coil_m,[2,2,2],0,'both');


D = smooth3(coil_m,'box',7);
%D = imdilate(D,strel('sphere',4));
mask_union_boundary_t = isosurface(Xv,Yv,Zv,D,0.5);
top_trimesh = triangulation(mask_union_boundary_t.faces,mask_union_boundary_t.vertices);

%%
mask_brain = niftiread('D:\Molin\fetal\result\segmentation\New folder\pred_intensity\MAP-C404\MAP-C404_0012.nii.gz');

mask_brain = mask_brain(1:down_rate:end,1:down_rate:end,1:down_rate:end);

move_lib = [-5,3,0; -5,10,0; -5,17,0; -5,3,2; -5,10,2; -5,17,2; -5,3,-2; -5,10,-2; -5,17,-2] * 4 / down_rate;
shift =move_lib(2,:) ;%[-5,7,0];%[0,0,5]; % move to right
%shift = [-5,7,0];
s = size(mask_brain);
mask_shift = padarray(mask_brain,[20 20 20]* 4 /down_rate,0,'both');

mask_brain = mask_shift(shift(1)+ 1 + 20 * 4 /down_rate:shift(1)+20* 4 /down_rate+s(1),shift(2)+1 + 20* 4 /down_rate :shift(2)+20* 4 /down_rate ...
+s(2),shift(3)+1+ 20* 4 /down_rate :shift(3)+20* 4 /down_rate+s(3));
D1 = smooth3(mask_brain,'box',5);

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

D = smooth3(coil_m,'box',7);
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
%ptstop = side;
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
%for dur = 3%[3,4,5,6]
%    for factor = 1.5%[1.5, 2]
dur = 6;
factor = 2;
DDD = 2.583; %1.46; %2.583

TB = dur*factor;

PP = load(['D:\Molin\fetal\meeting\7.23\results\TB_' num2str(TB) '_t_' num2str(dur) '.mat']);
profile = PP.profile;
bandw = 1000 * TB/dur;
transb = 1000 * DDD/ dur;

unshimmed = zeros(size(Xv));
tic;
[shimmed,amps,db0,mask_used,filt] ...
           = perform_shim_ncc_skull_linprog(unshimmed,mb_extend1,reshape(coil_loop.b,[size(Xv),size(ptstop,1)]),...
                                           50,50 * size(ptstop,1),{100,1e7,transb,0,0,0.99, bandw - transb, inf});
toc;

%%
figure(4);
set(gcf,'position',[0,0,1000,900])
tmpy = mask2vec(shimmed,mask_brain==1);
tmpy(tmpy>filt(1)+5000) = [];
h1 = histogram(tmpy,20); hold on;
h1.FaceAlpha = 0.5;
tmpy = mask2vec(shimmed,mb_extend1==3);
tmpy(tmpy<filt(2)-5000) = []; 
h2 = histogram(tmpy,100); hold on; 
h2.FaceAlpha = 0.5;
line([filt(1),filt(1)],[0,1000],'LineWidth',4,'color','black');
line([filt(2),filt(2)],[0,1000],'LineWidth',4,'color','black');
title(['high f: ' num2str(filt(1)) '      low f: ' num2str(filt(2)) '       transition bandwidth: ' num2str(transb) '      pass-bandwidth: ' num2str(bandw - transb)])
hold off;
saveas(gcf,['D:\Molin\fetal\meeting\7.23\results\rfshim_result\new' num2str(TB) '_' num2str(dur) '.png'])
%set(gca,'XLim',[-5000,5000]+mean(filt));

%% Plot on the same
for thred = [0.1, 0.5, 0.9]
close all
v = VideoWriter(['D:\Molin\fetal\meeting\7.28\result\idealrf_TB_' num2str(TB) 't_' num2str(dur) 'thred_' num2str(thred) '.avi']);
v.FrameRate = 2;
open(v)
shimmed_flip = flipb0(shimmed,filt);
b_m = zeros(size(D));
g_m = shimmed_flip.*(D>0);%zeros(size(v)); %kp_hm(:,:,:,5);
r_m = shimmed_flip.*(D>0) ;
flips = cat(4, r_m, g_m, b_m);
flips(flips >thred) = 1;
flips(flips <=thred) = 0;
figure(6);
set(gcf,'position',[0,0,1000,900])
set(gcf,'color','w');  
for i = 1 : size(D,3)
    cla
    h = imshow(squeeze(flips(:,:,i,:)),[0,1],'InitialMagnification','fit');
    %set(h, 'AlphaData', 0.1)
    hold on
    h2 = imshow(mask_brain(:,:,i),[0,2],'InitialMagnification','fit');
    set(h2, 'AlphaData', 0.3)
title(['127 coils, 20 turns, ideal rf' ' trans bw: ' num2str(transb) '  passbw: ' num2str(bandw - transb) '  thre: ' num2str(thred)],'Fontsize', 16)
frame = getframe(gcf);
writeVideo(v,frame);
drawnow
%pause(0.5)
end
close(v)
end
%%%
%%
for thred = [0.1, 0.5, 0.9]
ref = profile; 

[qq, ww, ee] = size(shimmed);
shimm_rf = shimmed;
center = (bandw - transb)/2 + filt(1) - 7*1000/dur;
for i = 1:qq
    for j = 1:ww
        for k = 1:ee
            shimm_rf(i,j,k) = mapfunc(shimmed(i,j,k), center, ref,dur);
        end
    end
end

%%%%plot on the same

v = VideoWriter(['D:\Molin\fetal\meeting\7.28\result\TB_' num2str(TB) 't_' num2str(dur) 'thred_' num2str(thred) 'center_' num2str(center) '.avi']);
v.FrameRate = 2;
open(v)
shimmed_flip = shimm_rf;
b_m = zeros(size(D));
g_m = shimmed_flip.*(D>0);%zeros(size(v)); %kp_hm(:,:,:,5);
r_m = shimmed_flip.*(D>0) ;
flips = cat(4, r_m, g_m, b_m);
flips (flips >thred) = 1;
flips (flips <=thred) = 0;
figure(5);
set(gcf,'position',[0,0,1000,900])
set(gcf,'color','w');  
for i = 1 : size(D,3)
    cla
    h = imshow(squeeze(flips(:,:,i,:)),[0,1],'InitialMagnification','fit');
    %set(h, 'AlphaData', 0.1)
    hold on
    h2 = imshow(mask_brain(:,:,i),[0,2],'InitialMagnification','fit');
    set(h2, 'AlphaData', 0.3)
title(['127 coils, 20 turns, ' ' trans bw: ' num2str(transb) '  passbw: ' num2str(bandw - transb) '  thre: ' num2str(thred)], 'Fontsize',16)
frame = getframe(gcf);
writeVideo(v,frame);
drawnow
%pause(0.5)
end
close(v)
end


%    end
%end

%% dif center location
for thred = [0.1, 0.5, 0.9]
    
    v = VideoWriter(['D:\Molin\fetal\meeting\7.28\result\newTB_' num2str(TB) 't_' num2str(dur) 'thred_' num2str(thred) '.avi']);
    v.FrameRate = 2;
    open(v)
    figure(5);
    set(gcf,'position',[0,0,1400,900])
    set(gcf,'color','w');  
    ref = profile;  
    shim_inter = [];
    cen_shift = [-30,0,30];
    for cen = cen_shift
    [qq, ww, ee] = size(shimmed);
    shimm_rf = shimmed;
    center = (bandw - transb)/2 + filt(1) - 7*1000/dur + cen;
    for i = 1:qq
        for j = 1:ww
            for k = 1:ee
                shimm_rf(i,j,k) = mapfunc(shimmed(i,j,k), center, ref,dur);
            end
        end
    end
    shim_inter = [shim_inter ; shimm_rf];
    end
%%%%plot on the same
    shim_inter = reshape(shim_inter, 3, size(shimmed,1),size(shimmed,2),size(shimmed,3));
    D1 = reshape([D; D; D], 3, size(shimmed,1),size(shimmed,2),size(shimmed,3));
shimmed_flip = shim_inter;
b_m = zeros(size(shim_inter));
g_m = shimmed_flip.*(D1>0);%zeros(size(v)); %kp_hm(:,:,:,5);
r_m = shimmed_flip.*(D1>0) ;
flips = cat(5, r_m, g_m, b_m);
flips (flips >thred) = 1;
flips (flips <=thred) = 0;

for i = 1 : size(D1,4)
    cla
    subplot(1,3,1)
    h = imshow(squeeze(flips(1,:,:,i,:)),[0,1],'InitialMagnification','fit');
    %set(h, 'AlphaData', 0.1)
    hold on
    h2 = imshow(mask_brain(:,:,i),[0,2],'InitialMagnification','fit');
    set(h2, 'AlphaData', 0.3)
    title(['censhift: ' num2str(cen_shift(1))])
    subplot(1,3,2)
    h = imshow(squeeze(flips(2,:,:,i,:)),[0,1],'InitialMagnification','fit');
    %set(h, 'AlphaData', 0.1)
    hold on
    h2 = imshow(mask_brain(:,:,i),[0,2],'InitialMagnification','fit');
    set(h2, 'AlphaData', 0.3)
    title([ 'censhift: ' num2str(cen_shift(2))])
    subplot(1,3,3)
    h = imshow(squeeze(flips(3,:,:,i,:)),[0,1],'InitialMagnification','fit');
    %set(h, 'AlphaData', 0.1)
    hold on
    h2 = imshow(mask_brain(:,:,i),[0,2],'InitialMagnification','fit');
    set(h2, 'AlphaData', 0.3)
    title(['censhift: ' num2str(cen_shift(3))])
frame = getframe(gcf);
writeVideo(v,frame);
drawnow
%pause(0.5)
end
close(v)
    
end

function map = mapfunc(shim, center, ref, dur)
    shim = shim - center ;
    xxx = floor((size(ref,1)-1)/2 /(7*1000/dur)  * shim);
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