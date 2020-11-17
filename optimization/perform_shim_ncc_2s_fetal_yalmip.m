% function shimmed = perform_shim(unshimmed,mask,max_current)

function [shimmed,amps,db0,mask_used,filt] ...
            = perform_shim_ncc_2s_fetal_yalmip(unshimmed,mask_all,coil,...
                                            max_current,total_current,tb,excite_bw)
mask = mask_all == 1;
mask_d = mask_all == 3;
mask_u = mask_all ==4;

unshimmed(isnan(unshimmed)) = 0;

nx = numel(coil(:,1,1,1));
ny = numel(coil(1,:,1,1));
nz = numel(coil(1,1,:,1));
nq = numel(coil(1,1,1,:));
unshimmed_all = reshape(unshimmed,nx*ny*nz,1);
unshimmed2 = zeros(nx,ny,nz);
for i = 1:nz
    unshimmed2(:,:,i) = medfilt2(unshimmed(:,:,i),[2,2]);
end
unshimmed_vec_sk = mask2vec(unshimmed,mask_d);


A_all = reshape(coil,[nx*ny*nz nq]);

A = mask2vec(coil,mask);
A_d = mask2vec(coil,mask_d);
A_u = mask2vec(coil,mask_u);


unshimmed_vec = mask2vec(unshimmed,mask);
unshimmed_vec_d = mask2vec(unshimmed,mask_d);
unshimmed_vec_u = mask2vec(unshimmed,mask_u);


amps_v = sdpvar(nq,1);
e_min = sdpvar(1,1);
e_max = sdpvar(1,1);

d_max = sdpvar(1,1);
d_min = sdpvar(1,1);

u_min = sdpvar(1,1);
u_max = sdpvar(1,1);

slacks_e = sdpvar(length(unshimmed_vec),1);
slacks_d = sdpvar(length(unshimmed_vec_d),1);
slacks_u = sdpvar(length(unshimmed_vec_u),1);


amps_l1 = norm(amps_v,1);


sconstraint_min_e = e_min <= A*amps_v+unshimmed_vec;
sconstraint_max_e = e_max >= A*amps_v+unshimmed_vec;

sconstraint_max_d = d_max >= A_d*amps_v+unshimmed_vec_d;
sconstrain_min_d = d_min <= A_d*amps_v+unshimmed_vec_d;

sconstraint_max_u = u_max >= A_u*amps_v+unshimmed_vec_u;
sconstrain_min_u = u_min <= A_u*amps_v+unshimmed_vec_u;


conMaxI = [-max_current*(ones(nq,1)) <= amps_v <= max_current*(ones(nq,1))];
conL1I = amps_l1<=total_current;
slacPosb = slacks_e >= 0;
slacPossk = slacks_d >= 0;
con_noverlap = [e_min-d_max>=tb, u_min-e_max>=tb,e_max-e_min==excite_bw];

opts = sdpsettings('showprogress',1,'solver','mosek');
p = optimize([conMaxI,conL1I,...
             sconstraint_min_e,sconstraint_max_e...
            ,sconstraint_max_d,sconstrain_min_d...
            ,sconstraint_max_u,sconstrain_min_u
            con_noverlap],norm(amps_v,2),opts);
amps = value(amps_v);
filt = [value((e_min)),value((d_max))]; 



%%
s_all = reshape(A_all*amps(1:nq)+unshimmed_all,[nx ny nz]);

s_b = s_all;
s_b(~mask) = nan;
s_sk = s_all;
s_sk(~mask_d) = nan;

mask = zeros(size(mask));
mask_d = zeros(size(mask));

mask(s_b>value((e_min))) = true;
mask_d(s_sk< value((d_max))) = true;

amps = amps(1:nq);
shimmed = A_all*amps;

shimmed = reshape(shimmed,[nx ny nz]);

shimmed = (unshimmed + shimmed);

% flip polarity of current semaskt_sktting for shim_fun (double check on shim_fun_std!!!)
amps = double(amps);

db0 = 0;

mask_used = int16(mask)+3*int16(mask_d);


end


