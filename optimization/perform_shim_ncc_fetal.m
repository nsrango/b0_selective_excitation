% function shimmed = perform_shim(unshimmed,mask,max_current)

function [shimmed,amps,db0,mask_used,filt] ...
            = perform_shim_ncc_fetal(unshimmed,mask_all,coil,...
                                            max_current,total_current,opts)
mask = mask_all == 1;
mask_sk = mask_all == 3;




testy = reshape(mask,[],1);

if nargin < 4
    tb = 70;
    b0 = 3;
else
    tb = opts{3};
    b0 = opts{4};
    
    if length(opts)<5
        f2naa = 0.7e-6;%ppm;
    else
        f2naa = opts{5};
    end
    if length(opts)<6
        alpha = 0.5;
    else
        alpha = opts{6};
    end
    if length(opts)<8
        highbw = inf;
        lowbw = inf;
    else
        highbw = opts{7};
        lowbw = opts{8};
    end
end

% f2naa = 0.7e-6;%ppm
gamma = 42.577478e6;


% if sum(testy)<=50
%     shimmed = unshimmed;
%     amps = zeros(32,1);
%     std_unshimmed = nan;
%     std_shimmed = nan;
%     return
% end

unshimmed(isnan(unshimmed)) = 0;

nx = numel(coil(:,1,1,1));
ny = numel(coil(1,:,1,1));
nz = numel(coil(1,1,:,1));
nq = numel(coil(1,1,1,:));
unshimmed_all = reshape(unshimmed,nx*ny*nz,1);

A = mask2vec(coil,mask);
A_sk = mask2vec(coil,mask_sk);
A_all = reshape(coil,[nx*ny*nz nq]);

maskt = mask;
maskt_sk = mask_sk;


A = mask2vec(coil,maskt);
A_sk = mask2vec(coil,maskt_sk);


unshimmed_vec = mask2vec(unshimmed,maskt);
unshimmed_vec_sk = mask2vec(unshimmed,maskt_sk);

nvar = 2*size(A,2)+4+size(A,1)+size(A_sk,1)+2;

sconstraint_b = sparse(size(A,1),nvar);
sconstraint_b(:,1:size(A,2)) = A;
sconstraint_b(:,1+2*size(A,2)+4:2*size(A,2)+4+size(A,1)) = 0;
sconstraint_b(:,1+2*size(A,2)+1) = -1;
sconstraint_b_RHS = -unshimmed_vec;

bmax_constraint = sparse(size(A,1),nvar);
bmax_constraint(:,1:size(A,2)) = A;
bmax_constraint(:,end-1) = -1;
bmax_constraint_RHS = -unshimmed_vec;

sconstraint_sk = sparse(size(A_sk,1),nvar);
sconstraint_sk(:,1:size(A_sk,2)) = A_sk;
sconstraint_sk(:,1+2*size(A,2)+4+size(A,1):end-2) = 0;
sconstraint_sk(:,1+2*size(A,2)+3) = -1;
sconstraint_sk_RHS = -unshimmed_vec_sk;

skmin_constraint = sparse(size(A_sk,1),nvar);
skmin_constraint(:,1:size(A_sk,2)) = A_sk;
skmin_constraint(:,end) = -1;
skmin_constraint_RHS = -unshimmed_vec_sk;

bwconstraint_high = sparse(1,nvar);
if(highbw < inf)
    bwconstraint_high(1,[1+2*size(A,2)+1,end-1]) = [-1,1];
    bwconstraint_high_RHS = highbw;
else
    bwconstraint_high_RHS = 0;
end

bwconstraint_low = sparse(1,nvar);
if(lowbw < inf)
    bwconstraint_low(1,[1+2*size(A,2)+3,end]) = [1,-1];
    bwconstraint_low_RHS = lowbw;
else
    bwconstraint_low_RHS = 0;
end

con_nooverlap = sparse(1,nvar);
con_nooverlap(1,1+2*size(A,2)+1) = 1;
con_nooverlap(1,1+2*size(A,2)+3) = -1;
con_nooverlap_RHS = tb-gamma*f2naa*b0;

abscon  = sparse(size(A,2),nvar);
abscon(1:size(A,2),1:size(A,2)) = speye(size(A,2));
abscon(1+size(A,2):2*size(A,2),1:size(A,2)) = -speye(size(A,2));
abscon(1:size(A,2),1+size(A,2):2*size(A,2)) = -speye(size(A,2));
abscon(1+size(A,2):2*size(A,2),1+size(A,2):2*size(A,2)) = -speye(size(A,2));
abscon_RHS = zeros(2*size(A,2),1);

abscon2 = sparse(1,nvar);
abscon2(:,1+size(A,2):2*size(A,2)) = 1;
abscon2_RHS = total_current;

Ap = [-sconstraint_b; 
       sconstraint_sk;
      -con_nooverlap; 
       abscon;
       abscon2;
       bmax_constraint;
       -skmin_constraint;
       bwconstraint_high;
       bwconstraint_low];
b = [-sconstraint_b_RHS;
      sconstraint_sk_RHS; 
     -con_nooverlap_RHS; 
      abscon_RHS;
      abscon2_RHS; 
      bmax_constraint_RHS;
      -skmin_constraint_RHS;
      bwconstraint_high_RHS;
      bwconstraint_low_RHS];
l = [-max_current* ones(2*size(A,2),1); -inf*ones(4,1); zeros(size(A,1),1); zeros(size(A_sk,1),1);-inf*ones(2,1)];
u = [max_current* ones(2*size(A,2),1); inf*ones(4,1); inf*ones(size(A,1),1); inf*ones(size(A_sk,1),1); inf*ones(2,1)];
f = [zeros(2*size(A,2),1); zeros(4,1); 2*alpha*zeros(size(A,1),1); 2*(1-alpha)*zeros(size(A_sk,1),1); zeros(2,1)];

            options = mskoptimset('');

options = mskoptimset(options,'Diagnostics','on','Display','iter');
[x, fval, exitflag, output, lambda] = linprog(f,Ap,b,[],[],l,u,options);
fprintf("%d\n",exitflag);
% amps_v = sdpvar(nq,1);
% b_min_slack = sdpvar(1,1);
% b_max = sdpvar(1,1);
% 
% sk_max_slack = sdpvar(1,1);
% sk_min = sdpvar(1,1);
% slacks_b = sdpvar(length(unshimmed_vec),1);
% slacks_sk = sdpvar(length(unshimmed_vec_sk),1);
% 
% 
% amps_l1 = norm(amps_v,1);
% sconstraint_b = b_min_slack <= A*amps_v+unshimmed_vec+slacks_b;
% sconstraint_b_max = b_max >= A*amps_v+unshimmed_vec+slacks_b;
% 
% sconstraint_sk = sk_max_slack >= A_sk*amps_v+unshimmed_vec_sk-slacks_sk;
% sconstrain_min_sk = sk_min <= A_sk*amps_v+unshimmed_vec_sk-slacks_sk;
% 
% conMaxI = [-max_current*(ones(nq,1)) <= amps_v <= max_current*(ones(nq,1))];
% conL1I = amps_l1<=total_current;
% slacPosb = slacks_b >= 0;
% slacPossk = slacks_sk >= 0;
% con_noverlap = b_min_slack-sk_max_slack>=tb-gamma*f2naa*b0;
% alpha = sdpvar(1,1);
% opts = sdpsettings('showprogress',1,'debug',1,'solver','mosek');
% % p = optimize([conMaxI,conL1I,...
% %             sconstraint_b,sconstraint_sk,sconstrain_min_sk...
% %             sconstraint_b_max, sk_max_slack-sk_min <= 10000,...
% %             slacPosb, slacPossk,...
% %             con_noverlap],+norm(slacks_b,1)+norm(slacks_sk,1),opts);

amps = x(1:size(A,2));
b_min_slack_v = x(1+2*size(A,2)+1);
sk_max_slack_v = x(1+2*size(A,2)+3);
filt = [((b_min_slack_v)-246*(b0>0)),((sk_max_slack_v))-gamma*f2naa*b0-246*(b0>0)]; 
% disp(((b_min_slack_v)-246*(b0>0)));
% disp(((sk_max_slack_v))-gamma*f2naa*b0-246*(b0>0));
% disp(((b_min_slack_v)-246*(b0>0))/2+(((sk_max_slack_v))-gamma*f2naa*b0-246*(b0>0))/2);


sum(x(size(A,2)+1:size(A,2)+1+size(A,1))>1e-5);
sum(x(size(A,2)+1+size(A,1)+1:size(A,2)+1+size(A,1)+1+size(A_sk,1))>1e-5);


%%
s_all = reshape(A_all*amps(1:nq)+unshimmed_all,[nx ny nz]);

s_b = s_all;
s_b(~mask) = nan;
s_sk = s_all;
s_sk(~mask_sk) = nan;

maskt = zeros(size(mask));
maskt_sk = zeros(size(mask));

maskt(s_b>((b_min_slack_v))) = true;
maskt_sk(s_sk< ((sk_max_slack_v))) = true;

amps = amps(1:nq);
shimmed = A_all*amps;

shimmed = reshape(shimmed,[nx ny nz]);

shimmed = (unshimmed + shimmed);

% flip polarity of current semaskt_sktting for shim_fun (double check on shim_fun_std!!!)
amps = double(amps);

db0 = 0;

mask_used = int16(maskt)+3*int16(maskt_sk);


end

