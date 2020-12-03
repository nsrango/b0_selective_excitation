handles.fm = {};
handles.fm.mag = zeros(1,1);
handles.fm.mask = int16(zeros(1,1));
handles.fm.phase = zeros(1,1);
handles.fm.size = [1,1,1];
handles.fm.unshimmed = zeros(1,1);
handles.fm.predicted = zeros(1,1);
handles.fm.shimmed = zeros(1,1);
handles.fm.mask_used = zeros(1,1);
handles.fm.coil_fm = zeros(1,1,1);
handles.fm.coil_fm_raw = zeros(1,1,1);
handles.delta_TE = 0.001*2.46;
handles.bet_thresh = 0.3;
handles.slices = 1;
handles.opt_used = 1;
handles.parcel_used = 1;
%%
handles.fm.mag_path = 'data/GRE_FIELD_75SLI_DSHIM_0009';
handles = mag_import(handles);

handles.fm.unshimmed_path = 'data/GRE_FIELD_75SLI_DSHIM_0010';
handles = us_import(handles);


handles.fm.mask = mask_brain(handles.fm.mag,handles.bet_thresh,0,0,0,0,0);

handles.fm.coil_path = "field_maps_for_lincoln.mat";


load(handles.fm.coil_path);
handles.fm.coil_fm = data.coil.maps;
coilmask = 1|(handles.fm.coil_fm(:,:,:,1)~=handles.fm.coil_fm(1,1,1,1));
handles.fm.mask = handles.fm.mask;%%.*int16(coilmask);
handles.fm.unshimmed = unwrap_phase(handles.fm, handles.delta_TE);
%%
[nx,ny,nz] = size(handles.fm.unshimmed);
[X,Y,Z] = meshgrid(1:nx,1:ny,1:2:nz*2);

mask_excite = sqrt(((X-nx/2+20).^2+(Y-ny/2+10).^2+(Z-nz*2/2+25).^2))<=8;

idxs = find(mask_excite(:,:,:)==1 & handles.fm.mask(:,:,:)~=0);
[I,J,K] = ind2sub(size(mask_excite),idxs);
      
excite_bounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];
excite_bounds = excite_bounds+[-0,0
                           -5,5
                           -5,5]; 
                       
excite_bounds = max(min(excite_bounds, size(mask_excite)'),[1;1;1]);                 

D= handles.fm.mask~=0;

mb_extend = zeros(size(mask_excite));
mb_extend(D==1) = 2;
mb_extend(:,excite_bounds(2,1):excite_bounds(2,2),:) = 3;
mb_extend(excite_bounds(1,1):excite_bounds(1,2),excite_bounds(2,1):excite_bounds(2,2),excite_bounds(3,1):excite_bounds(3,2)) = 2;
mb_extend((mask_excite==1)&handles.fm.mask(:,:,:)~=0) = 1 ;
mb_extend = int16(mb_extend).*int16(D);

mb_extend_con = zeros(size(mask_excite));
mb_extend_con(D==1) = 2;
mb_extend_con(:,excite_bounds(2,1):excite_bounds(2,2),:) = 3;
mb_extend_con((mask_excite==1)) = 1 ;
mb_extend_con = int16(mb_extend_con).*int16(D);


imagesc(mb_extend(:,:,25))

%%

mb_patch = zeros(size(mask_excite));
mb_patch(D==1) = 2;
ds = diff(excite_bounds,1,2)+1;
val = 3;
for i = 0
    for j = -5:5
        for k = -5:5
            if k + j <= 0
                val = 3;
            else 
                val = 3;
            end
            if j>=0
                tmp = padarray(mask_excite(:,1:(end-(j*ds(2))),:)==1,[0,min(j*ds(2),size(mask_excite,2)),0],'pre');
            else
                tmp = padarray(mask_excite(:,1-(j*ds(2)):end,:)==1,[0,max(-j*ds(2),0),0],'post');
            end
            
            if k>=0
                tmp = padarray(tmp(:,:,1:(end-(k*ds(3))))==1,[0,0,min(k*ds(3),size(mask_excite,3))],'pre');
            else
                tmp = padarray(tmp(:,:,1-(k*ds(3)):end)==1,[0,0,max(-k*ds(3),0)],'post');
            end
             
            mb_patch(tmp==1) = val;
            imagesc(squeeze(mb_patch(:,:,37)))
        end
    end
end

mb_patch(excite_bounds(1,1):excite_bounds(1,2),excite_bounds(2,1):excite_bounds(2,2),excite_bounds(3,1):excite_bounds(3,2)) = 2;


mb_patch(mask_excite==1) = 1;
mb_patch(D==0) = 0;


for s = 25
    imagesc([mb_patch(:,:,s) ; mb_extend(:,:,s)]); colorcet('L1'); axis equal;
    pause(0.1);
end

%%

[shimmed,amps,db0,mask_used,filt] ...
           = perform_shim_ncc_2s_fetal_yalmip(handles.fm.unshimmed,mb_patch,handles.fm.coil_fm,...
                                           50,5000,200,1000);


N=1024;

t0=1e-3;  % width of a sidelobe  or half the width of the main lobe.  bandwidth of frequency response is roughly 1/t0

tau = 12e-3;
time = linspace(-tau/2,tau/2,N);
TBP = tau/t0;
alp=.46;  % 0.46 for Hanning window
b = sin(pi*time./t0)./(pi*time).*((1-alp)+alp*cos(pi*time/((TBP/2)*t0)));
b=b./max(abs(b));
b_all=[zeros(1,numel(b)) b zeros(1,numel(b))];

time_all = linspace(-1.5*tau,1.5*tau,3*N); 
BW=1/(time_all(2)-time_all(1));
freq = linspace(-BW/2,BW/2,numel(time_all));
inds=(numel(freq)/2-35):(numel(freq)/2+35);
bf = fftshift(fft(fftshift(b_all)));

CK_param = interp1(freq+1/t0/2+mean(filt),abs(bf)/max(abs(bf)),shimmed(:),'linear',0);


% plot(freq,abs(bf));xlim([-900,900])


% rf2 = dzrf(500, 12, 'ex', 'pm',0.0006, 0.012);
% x = -20 *3.5:0.01:20 *3.5;
% ab2  = abr(rf2, x);
% 
% CK_param = interp1(toHz(x,8)+760+mean(filt),ab2ex(ab2),shimmed(:),'linear',0);


disp('{')
for i = 1:size(amps,2) 
    allOneString = strcat('{',sprintf('%.2f,' , amps(:,i)),'},');
    disp(allOneString)
    allOneString = strcat('{',sprintf('%.2f,' , -amps(:,i)),'},');
    disp(allOneString)
    allOneString = strcat('{',sprintf('%.2f,' , 0*amps(:,i)),'},');
    disp(allOneString)
end
disp('};')
    
%%

ex = reshape(CK_param,size(shimmed));
f1 = figure(1); axis off;
a1 = gca;
f2 = figure(2); axis off;
a2 = gca;
f3 = figure(3); axis off;
a3 = gca;
ex2 = ex;
masku = mask_excite.*D;
ex2(masku(:,:,:)==0) = nan;


idxs = find(((mb_extend_con(:,:,:)==1 | mb_extend_con(:,:,:)==3).*ex)>0.1500);
[I,J,K] = ind2sub(size(ex2),idxs);

bounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];
      
idxs = find(mask_excite(:,:,:)==1);
[I,J,K] = ind2sub(size(ex2),idxs);
      
brainbounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];
      
idxs = find(D(:,:,:)~=0);
[I,J,K] = ind2sub(size(ex2),idxs);
      
bodybounds = [min(I) max(I);
          min(J) max(J);
          min(K) max(K);];      

boxmask = zeros(size(mask_excite));
boxmask(excite_bounds(1,1):excite_bounds(1,2),excite_bounds(2,1):excite_bounds(2,2),excite_bounds(3,1):excite_bounds(3,2)) = 1;
bodybox = zeros(size(mask_excite));
bodybox(bodybounds(1,1):bodybounds(1,2),bodybounds(2,1):bodybounds(2,2),bodybounds(3,1):bodybounds(3,2)) = 1;


for s = 1:75
    
   
    [B,L] = bwboundaries(mask_excite(:,:,s)==1,'noholes');
    [B2,L2] = bwboundaries(boxmask(:,:,s)==1,'noholes');
    [B3,L3] = bwboundaries(mb_patch(:,:,s)>=3,'noholes');

    imagesc((D(:,:,s)~=0).*abs(ex(:,:,s)),'parent',a1,[0,1]);
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
       plot(a1,boundary(:,2), boundary(:,1), 'red', 'LineWidth', 2);
    end
    hold(a1, 'off');
   
    imagesc((D(:,:,s)~=0).*shimmed(:,:,s),'parent',a2,500*[-1,1]+mean(filt));
    set(a2,'Colormap',colorcet('D1A'))
    axis off;
    
    imagesc(mb_extend(:,:,s),'parent',a3);
    set(a3,'Colormap',colorcet('L1'))
    axis off;
    
    
    pause(0.1);
end


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



