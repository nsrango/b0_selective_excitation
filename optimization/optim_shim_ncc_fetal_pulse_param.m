% function shimmed = perform_shim(unshimmed,mask,max_current)

function [p] ...
            = optim_shim_ncc_fetal_pulse_param(unshimmed,mask_all,coil)
mask = mask_all == 1;
mask_excite = mask_all == 3;

unshimmed(isnan(unshimmed)) = 0;

nx = numel(coil(:,1,1,1));
ny = numel(coil(1,:,1,1));
nz = numel(coil(1,1,:,1));
nq = numel(coil(1,1,1,:));

maskt_rej = mask;
maskt_ex = mask_excite;


A_rej = mask2vec(coil,maskt_rej);
A_ex = mask2vec(coil,maskt_ex);

unshimmed_vec_rej = mask2vec(unshimmed,maskt_rej);
unshimmed_vec_ex = mask2vec(unshimmed,maskt_ex);

amps_v = sdpvar(nq,1);
b_min = sdpvar(1,1);
b_max = sdpvar(1,1);

ex_max = sdpvar(1,1);
ex_min = sdpvar(1,1);
tb = sdpvar(1,1);
excite_bw = sdpvar(1,1);

max_current = sdpvar(1,1);

amps_l1 = norm(amps_v,1);
constraint_rej = b_min <= A_rej*amps_v+unshimmed_vec_rej;
constraint_rej_max = b_max >= A_rej*amps_v+unshimmed_vec_rej;

constraint_ex = ex_max >= A_ex*amps_v+unshimmed_vec_ex;
constrain_min_ex = ex_min <= A_ex*amps_v+unshimmed_vec_ex;

conMaxI = [-max_current*(ones(nq,1)) <= amps_v <= max_current*(ones(nq,1))];
% conL1I = amps_l1<=total_current;

con_noverlap = b_min-ex_max>=tb;

opts = sdpsettings('showprogress',0,'solver','mosek');
p = optimizer([conMaxI,...
            constraint_rej,constraint_ex,constrain_min_ex...
            constraint_rej_max, ex_max-ex_min <= excite_bw,...
            con_noverlap],max_current,opts,{tb,excite_bw},max_current);
        
        
end



