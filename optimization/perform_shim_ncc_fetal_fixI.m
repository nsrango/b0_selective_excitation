% function shimmed = perform_shim(unshimmed,mask,max_current)

function [shimmed,amps,db0,mask_used,filt] ...
            = perform_shim_ncc_fetal_fixI(unshimmed,mask_all,coil,...
                                            max_current,total_current,opts)
mask = mask_all == 1;
mask_sk = mask_all == 3;


testy = reshape(mask,[],1);

if nargin < 6
    tb = 70;
    b0 = 3;
else
    tb = opts{3};
    b0 = opts{4};
    if length(opts)<5
        excite_bw = 10000;
    else
        excite_bw = opts{5};
    end
end

f2naa = 0.7e-6;%ppm
gamma = 42.577478e6;


if sum(testy)<=50
    shimmed = unshimmed;
    amps = zeros(32,1);
    std_unshimmed = nan;
    std_shimmed = nan;
    return
end

unshimmed(isnan(unshimmed)) = 0;

nx = numel(coil(:,1,1,1));
ny = numel(coil(1,:,1,1));
nz = numel(coil(1,1,:,1));
nq = numel(coil(1,1,1,:));
unshimmed_vec = mask2vec(unshimmed,mask);
unshimmed_all = reshape(unshimmed,nx*ny*nz,1);
unshimmed2 = zeros(nx,ny,nz);
for i = 1:nz
    unshimmed2(:,:,i) = medfilt2(unshimmed(:,:,i),[2,2]);
end
unshimmed_vec_sk = mask2vec(unshimmed,mask_sk);


A = mask2vec(coil,mask);
A_sk = mask2vec(coil,mask_sk);
A_all = reshape(coil,[nx*ny*nz nq]);
n_b = size(A,1);
n_sk = size(A_sk,1);

maskt = mask;
maskt_sk = mask_sk;


A = mask2vec(coil,maskt);
A_sk = mask2vec(coil,maskt_sk);
n_b = size(A,1);
n_sk = size(A_sk,1);

unshimmed_vec = mask2vec(unshimmed,maskt);
unshimmed_vec_sk = mask2vec(unshimmed,maskt_sk);

n_us_b = length(unshimmed_vec);
n_us_sk = length(unshimmed_vec_sk);
A_ineq_b = [A, -1*ones(n_us_b,1), zeros(n_us_b,1)];
A_ineq_sk = [-1*A_sk, zeros(n_us_sk,1), ones(n_us_sk,1)];

%     sample_b = datasample(1:n_us_b,round(n_us_b*1/8),...
%                             'Replace',false,'Weights',weights_b);
%     sample_sk = datasample(1:n_us_sk,round(n_us_sk*1/8),...
%                             'Replace',false,'Weights',weights_sk);



amps_v = sdpvar(nq,1);
b_min = sdpvar(1,1);
b_max = sdpvar(1,1);

sk_max = sdpvar(1,1);
sk_min = sdpvar(1,1);
tb = sdpvar(1,1);
max_current = sdpvar(1,1);

amps_l1 = norm(amps_v,1);
constraint_b = b_min <= A*amps_v+unshimmed_vec;
constraint_b_max = b_max >= A*amps_v+unshimmed_vec;

constraint_sk = sk_max >= A_sk*amps_v+unshimmed_vec_sk;
constrain_min_sk = sk_min <= A_sk*amps_v+unshimmed_vec_sk;

conMaxI = [-max_current*(ones(nq,1)) <= amps_v <= max_current*(ones(nq,1))];
conL1I = amps_l1<=total_current;

con_noverlap = b_min-sk_max>=tb-gamma*f2naa*b0;

opts = sdpsettings('showprogress',1,'solver','mosek');
p = optimizer([conMaxI,conL1I,...
            constraint_b,constraint_sk,constrain_min_sk...
            constraint_b_max, sk_max-sk_min <= 2000,...
            con_noverlap],-tb,opts,max_current,tb);
        
        
amps = value(amps_v);
filt = [value((b_min)-246*(b0>0)),value((sk_max))-gamma*f2naa*b0-246*(b0>0)]; 
disp(value((b_min)-246*(b0>0)));
disp(value((sk_max))-gamma*f2naa*b0-246*(b0>0));
disp(value((b_min)-246*(b0>0))/2+(value((sk_max))-gamma*f2naa*b0-246*(b0>0))/2);

%%
% [Xm, Ym, Zm] = ndgrid(1:nx,1:ny,1:nz);
% dvar = dual(sconstraint_b);
% xm_ma_b = mask2vec(Xm,mask);
% ym_ma_b = mask2vec(Ym,mask);
% zm_ma_b = mask2vec(Zm,mask);
% 
% xm_ma_b(dvar<0.2) = [];
% ym_ma_b(dvar<0.2) = [];
% zm_ma_b(dvar<0.2) = [];
% 
% idxma_b = [xm_ma_b,ym_ma_b,zm_ma_b];
% 
% temp = zeros(nx,ny,nz);
% for i = 1:size(idxma_b,1) 
%     temp(idxma_b(i,1),idxma_b(i,2),idxma_b(i,3))=1; 
% end
% 
% [Xm, Ym, Zm] = ndgrid(1:nx,1:ny,1:nz);
% dvar = dual(sconstraint_sk);
% xm_ma_sk = mask2vec(Xm,mask_sk);
% ym_ma_sk = mask2vec(Ym,mask_sk);
% zm_ma_sk = mask2vec(Zm,mask_sk);
% 
% xm_ma_sk(dvar<0.2) = [];
% ym_ma_sk(dvar<0.2) = [];
% zm_ma_sk(dvar<0.2) = [];
% 
% idxma_sk = [xm_ma_sk,ym_ma_sk,zm_ma_sk];
% 
% temp_sk = zeros(nx,ny,nz);
% for i = 1:size(idxma_sk,1) 
%     temp_sk(idxma_sk(i,1),idxma_sk(i,2),idxma_sk(i,3))=1; 
% end

% figure(1); imagesc(temp(:,:,7)+2*temp_sk(:,:,7));

%%
s_all = reshape(A_all*amps(1:nq)+unshimmed_all,[nx ny nz]);
s_sk_max = max(mask2vec(s_all,mask_sk))-gamma*f2naa*b0+tb;

s_b = s_all;
s_b(~mask) = nan;
s_sk = s_all;
s_sk(~mask_sk) = nan;

maskt = zeros(size(mask));
maskt_sk = zeros(size(mask));

maskt(s_b>value((b_min))) = true;
maskt_sk(s_sk< value((sk_max))) = true;

amps = amps(1:nq);
shimmed = A_all*amps;

shimmed = reshape(shimmed,[nx ny nz]);

shimmed = (unshimmed + shimmed);

% flip polarity of current semaskt_sktting for shim_fun (double check on shim_fun_std!!!)
amps = double(amps);

db0 = 0;

mask_used = int16(maskt)+3*int16(maskt_sk);


end



