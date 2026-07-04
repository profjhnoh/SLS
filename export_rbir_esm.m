% export_rbir_esm.m
% Export the MATLAB RBIR effective-SINR-mapping tables (rbir.mat) for the C++
% SLS (rbir_esm.dat) plus an oracle CSV computed with the REAL
% wireless.internal.L2SM.calculateEffectiveSINR (alpha=beta=1).
%
% Usage:
%   matlab -batch "cd('D:/Dropbox/SLS/2026/Tbps_DLSLS_20251010'); export_rbir_esm"
%
% LICENSE NOTE: rbir.mat is MathWorks-proprietary data. The exported .dat and
% oracle CSV stay in this research repo ONLY — never in EMIMO_SLS_dist.
%
% Semantics (calculateEffectiveSINR.m, verified):
%   inputs/outputs in dB; sinr = sinr/beta; clamp to [minSNR, maxSNR];
%   scrbir = interp1(snr, rbir, sinr); avrbir = mean(scrbir);
%   snreff = alpha * interp1(rbir, snr, avrbir).

R = load(fullfile(matlabroot,'toolbox','shared','wireless','+wireless','+internal','+L2SM','rbir.mat'));
out_dat = fullfile(pwd, 'rbir_esm.dat');
out_dir = fullfile(pwd, 'bler_compare');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
oracle_csv = fullfile(out_dir, 'rbir_reference_tuples.csv');

MOD_KEEP = [4 16 64 256];   % QPSK/16/64/256QAM (Qm 2/4/6/8); we do not use 1024/4096

%% ---- 1) table export -----------------------------------------------------
fid = fopen(out_dat, 'w');
assert(fid > 0);
fprintf(fid, '# MATLAB rbir.mat export (802.11-14/1450r0 ext.) - MathWorks proprietary, DO NOT REDISTRIBUTE\n');
fprintf(fid, '# table key: modscheme (=2^Qm) ; rows: snr_dB rbir ; clamp to [min_snr, max_snr]\n');
fprintf(fid, 'MATLAB_RBIR_ESM v1\n');
fprintf(fid, 'tables %d\n', numel(MOD_KEEP));
for m = MOD_KEEP
    i = find(R.tableModScheme == m, 1);
    assert(~isempty(i), 'modscheme %d not in rbir.mat', m);
    T = double(R.rbirTable{i});
    fprintf(fid, 'table mod %d min_snr %.6f max_snr %.6f n %d\n', ...
        m, R.minSNR(i), R.maxSNR(i), size(T,1));
    fprintf(fid, '%.6f %.10e\n', T');
end
fprintf(fid, 'end\n');
fclose(fid);
d = dir(out_dat);
fprintf('Wrote %s (%.2f MB)\n', out_dat, d.bytes/1e6);

%% ---- 2) oracle: effective SINR for deterministic SINR-dB vectors ---------
% Vectors chosen to exercise: flat (ESM identity), ramps (dispersion),
% bimodal (worst-case spread), clamping at both table ends.
vecs = {
    'flat5',    5.0 * ones(1,10);
    'single3',  3.0;
    'ramp',     linspace(-5, 15, 21);
    'bimodal',  [ -8*ones(1,5), 12*ones(1,5) ];
    'wideramp', linspace(-15, 25, 41);
    'lowclamp', [ -25, -22, -18, -10 ];      % below minSNR -> clamp
    'highclamp',[ 20, 28, 35, 40 ];          % above maxSNR (mod-dep) -> clamp
    'mixed',    [ -3.2, 0.7, 4.4, 9.9, 2.1, -1.5, 7.3, 11.8 ];
};

rows = cell(0, 4);
for m = MOD_KEEP
    for v = 1:size(vecs,1)
        name = vecs{v,1};
        s    = vecs{v,2};
        eff  = wireless.internal.L2SM.calculateEffectiveSINR(s, m, 1.0, 1.0);
        rows(end+1,:) = {m, name, eff, strjoin(compose('%.4f', s), ' ')}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'modscheme','vec_id','eff_sinr_dB','sinr_dB_list'});
writetable(T, oracle_csv);
fprintf('Wrote %d oracle tuples to %s\n', height(T), oracle_csv);

% quick print
disp(T(1:8, {'modscheme','vec_id','eff_sinr_dB'}));
fprintf('Done.\n');
