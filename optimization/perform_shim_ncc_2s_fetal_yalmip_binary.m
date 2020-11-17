% function shimmed = perform_shim(unshimmed,mask,max_current)

function [shimmed,amps,db0,mask_used,filt] ...
            = perform_shim_ncc_2s_fetal_yalmip_binary(unshimmed,mask_all,coil,...
                                            max_current,total_current,tb,excite_bw)
mask = mask_all == 1;

regions = setdiff(unique(mask_all),[0,1,2]);
numr = numel(regions);

suppression_masks = bsxfun(@(x,y)(x==y),mask_all,permute(regions,[4,3,2,1]));


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

unshimmed_vec_1 = mask2vec(unshimmed,suppression_masks(:,:,:,1));


A_all = reshape(coil,[nx*ny*nz nq]);

A = mask2vec(coil,mask);
% A_d = mask2vec(coil,mask_d);
% A_u = mask2vec(coil,mask_u);
Ans = {};
unshimmed_vec_n = {};
for i = 1:numr
    Ans{i} = mask2vec(coil,suppression_masks(:,:,:,i));
    unshimmed_vec_n{i} = mask2vec(unshimmed,suppression_masks(:,:,:,i));
end

unshimmed_vec = mask2vec(unshimmed,mask);
% unshimmed_vec_d = mask2vec(unshimmed,mask_d);
% unshimmed_vec_u = mask2vec(unshimmed,mask_u);


amps_v = sdpvar(nq,1);
e_min = sdpvar(1,1);
e_max = sdpvar(1,1);

min_n = sdpvar(numr,1);
max_n = sdpvar(numr,1);

min_1 = sdpvar(1,1);
max_1 = sdpvar(1,1);

d_max = sdpvar(1,1);
d_min = sdpvar(1,1);

u_min = sdpvar(1,1);
u_max = sdpvar(1,1);

slacks_e = sdpvar(length(unshimmed_vec),1);
% slacks_d = sdpvar(length(unshimmed_vec_d),1);
% slacks_u = sdpvar(length(unshimmed_vec_u),1);


amps_l1 = norm(amps_v,1);


sconstraint_min_e = e_min <= A*amps_v+unshimmed_vec;
sconstraint_max_e = e_max >= A*amps_v+unshimmed_vec;

sconstraint_max_n=[];
sconstraint_min_n=[];
for i = 1:numr
    sconstraint_max_n = [sconstraint_max_n,(max_n(i) >= Ans{i}*amps_v+unshimmed_vec_n{i})];
    sconstraint_min_n = [sconstraint_min_n,(min_n(i) <= Ans{i}*amps_v+unshimmed_vec_n{i})];
end

% sconstraint_max_d = d_max >= A_d*amps_v+unshimmed_vec_d;
% sconstrain_min_d = d_min <= A_d*amps_v+unshimmed_vec_d;
% 
% sconstraint_max_u = u_max >= A_u*amps_v+unshimmed_vec_u;
% sconstrain_min_u = u_min <= A_u*amps_v+unshimmed_vec_u;


conMaxI = [-max_current*(ones(nq,1)) <= amps_v <= max_current*(ones(nq,1))];
conL1I = amps_l1<=total_current;
% con_noverlap = [e_min-d_max>=tb, u_min-e_max>=tb,e_max-e_min==excite_bw];

F = {};
t = {};
for i = 1:3
    F{i} = e_min-max_n(i)>=tb | min_n(i)-e_max>=tb;
    t{i} = [max_n(i)<=100000,min_n(i)>=-100000];
end



Fa = [F{:}];
ta = [t{:}];


opts = sdpsettings('showprogress',1,'solver','mosek');
p = optimize([conMaxI,conL1I,...
             sconstraint_min_e,sconstraint_max_e,...
%             ,sconstraint_max_d,sconstrain_min_d...
             sconstraint_max_n,sconstraint_min_n,...
%              sconstraint_max_u,sconstrain_min_u
            Fa,ta,e_max-e_min==excite_bw],[],opts);
amps = value(amps_v);
filt = [value((e_min)),value((e_min)-tb)]; 



%%
s_all = reshape(A_all*amps(1:nq)+unshimmed_all,[nx ny nz]);

amps = amps(1:nq);
shimmed = A_all*amps;

shimmed = reshape(shimmed,[nx ny nz]);

shimmed = (unshimmed + shimmed);

% flip polarity of current semaskt_sktting for shim_fun (double check on shim_fun_std!!!)
amps = double(amps);

db0 = 0;

mask_used = int16(mask);

end


