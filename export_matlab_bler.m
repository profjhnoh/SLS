% export_matlab_bler.m
% Export MATLAB 5G Toolbox L2SM AWGN BLER tables to a consolidated text file
% for the C++ SLS (matlab_l2sm_bler.dat), plus a reference-tuple oracle CSV
% (bler_compare/matlab_reference_tuples.csv) for the C++ selftest.
%
% Usage:
%   matlab -batch "cd('D:/Dropbox/SLS/2026/Tbps_DLSLS_20251010'); export_matlab_bler"
%
% LICENSE NOTE: awgnTables is MathWorks-proprietary data. The exported .dat
% and oracle CSV stay in this research repo ONLY — never in the public
% EMIMO_SLS_dist distribution (ship this script instead).
%
% Table semantics replicated from (R2026a sources, verified):
%   nr5g.internal.L2SM getAWGNTable  : bucket by outer R ranges [Rlo Rhi],
%       inner: first tabulated R >= requested, exact Qm, exact Zc
%   nr5g.internal.L2SM sinrToCodeBLER: ECR clip [1/1024,1023/1024], R=round(ECR*1024),
%       nrDLSCHInfo(TBS, R/1024) -> BGN,Zc,C
%   wireless.internal.L2SM.interpolatePER: truncate at first BLER==0 row (inclusive),
%       10^interp1(snr, log10(per), 'linear','extrap'), per<1e-6->0, per>1->1
%   transport BLER = 1 - (1-cbBLER)^C

L2SM_path = fullfile(matlabroot,'toolbox','5g','netmod','+nr5g','+internal','L2SM.mat');
out_dat   = fullfile(pwd, 'matlab_l2sm_bler.dat');
out_dir   = fullfile(pwd, 'bler_compare');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
oracle_csv = fullfile(out_dir, 'matlab_reference_tuples.csv');

QM_KEEP = [2 4 6 8];   % drop 1024QAM (Qm=10): our MCS table tops out at 256QAM

fprintf('Loading %s\n', L2SM_path);
S = load(L2SM_path);
T = S.awgnTables;
fprintf('awgnTables fields: %s\n', strjoin(fieldnames(T), ', '));
fprintf('BGN values: [%s]\n', num2str(T.BGN'));

%% ---------------- Pass 1: count curves (after Qm filter / NaN trim) ------
nCurves = 0;
for b = 1:numel(T.BGN)
    G = T.data(b);
    fprintf('BGN=%d: outer R ranges = %d buckets\n', T.BGN(b), size(G.R,1));
    for g = 1:numel(G.data)
        sub = G.data(g);
        sz = size(sub.data);
        assert(numel(sz) >= 2, 'unexpected leaf dims');
        nR  = numel(sub.R); nQ = numel(sub.Qm); nZ = numel(sub.Zc);
        if g == 1
            fprintf('  bucket 1: inner R x Qm x Zc = %d x %d x %d, data dims [%s]\n', ...
                nR, nQ, nZ, num2str(sz));
            fprintf('  inner Qm values: [%s]\n', num2str(sub.Qm(:)'));
        end
        for j = 1:nQ
            if ~ismember(sub.Qm(j), QM_KEEP), continue; end
            nCurves = nCurves + nR * nZ;   % upper bound; empties skipped in pass 2
        end
    end
end
fprintf('Curve upper bound after Qm filter: %d\n', nCurves);

%% ---------------- Pass 2: write consolidated .dat -----------------------
fid = fopen(out_dat, 'w');
assert(fid > 0, 'cannot open %s', out_dat);
fprintf(fid, '# MATLAB R2026a 5G Toolbox L2SM awgnTables export - MathWorks proprietary, DO NOT REDISTRIBUTE\n');
fprintf(fid, '# curve key: bgn qm r1024 zc ; rows: snr_dB cbBLER ; TB BLER = 1-(1-cbBLER)^C\n');
fprintf(fid, 'MATLAB_L2SM_BLER v1\n');
countPos = ftell(fid);
fprintf(fid, 'curves %10d\n', 0);   % placeholder, patched after pass

nWritten = 0;
for b = 1:numel(T.BGN)
    bgn = T.BGN(b);
    G = T.data(b);
    for g = 1:numel(G.data)
        sub = G.data(g);
        D = double(sub.data);
        for i = 1:numel(sub.R)
            for j = 1:numel(sub.Qm)
                if ~ismember(sub.Qm(j), QM_KEEP), continue; end
                for k = 1:numel(sub.Zc)
                    c = D(:,:,i,j,k);
                    ok = ~any(isnan(c), 2);
                    c = c(ok, :);
                    if size(c,1) < 2, continue; end
                    fprintf(fid, 'curve bgn %d qm %d r1024 %g zc %d n %d\n', ...
                        bgn, sub.Qm(j), sub.R(i), sub.Zc(k), size(c,1));
                    fprintf(fid, '%.4f %.10e\n', c');
                    nWritten = nWritten + 1;
                end
            end
        end
    end
end
fprintf(fid, 'end\n');
fseek(fid, countPos, 'bof');
fprintf(fid, 'curves %10d\n', nWritten);
fclose(fid);
d = dir(out_dat);
fprintf('Wrote %d curves to %s (%.1f MB)\n', nWritten, out_dat, d.bytes/1e6);

%% ---------------- Oracle: reference tuples for C++ selftest --------------
% Our MCS table (must mirror receive_downlink.cpp tbs_info_bits_from_mcs):
%   mcs<=4 QPSK(2), <=10 16QAM(4), <=19 64QAM(6), else 256QAM(8)
CR_X1024 = [120,193,308,449,602, 378,434,490,553,616,658, ...
            466,517,567,616,666,719,772,822,873, ...
            682.5,711,754,797,841,885,916.5,948];
qm_of_mcs = @(m) 2*(m<=4) + 4*(m>4 & m<=10) + 6*(m>10 & m<=19) + 8*(m>19);

mcs_list  = [0 4 5 10 11 19 20 27];
sinr_list = [-8 -2 0 4 10 16 22];
tbs_list  = [292 300 3000 3824 3825 8000 50000 200000];   % incl. segmentation edges

rows = cell(0, 10);
for m = mcs_list
    qm = qm_of_mcs(m);
    cr = CR_X1024(m+1);
    ecr = min(max(cr/1024, 1/1024), 1023/1024);
    R = round(ecr * 1024);
    for tbs = tbs_list
        info = nrDLSCHInfo(tbs, R/1024);
        % --- curve selection (verbatim getAWGNTable logic) ---
        Gb  = T.data(T.BGN == info.BGN);
        sel = R >= Gb.R(:,1) & R <= Gb.R(:,2);
        sub = Gb.data(sel);
        i = find(sub.R - R >= 0, 1);
        j = find(sub.Qm == qm);
        k = find(sub.Zc == info.Zc);
        assert(~isempty(i) && ~isempty(j) && ~isempty(k), ...
            'no curve for mcs=%d tbs=%d (R=%d qm=%d Zc=%d BGN=%d)', m, tbs, R, qm, info.Zc, info.BGN);
        lut = double(sub.data(:,:,i,j,k));
        lut = lut(~any(isnan(lut),2), :);
        for s = sinr_list
            cb = wireless.internal.L2SM.interpolatePER(s, lut);
            tb = 1 - (1 - cb)^info.C;
            rows(end+1, :) = {s, m, qm, sub.R(i), tbs, info.BGN, info.Zc, info.C, cb, tb}; %#ok<SAGROW>
        end
    end
end
Tbl = cell2table(rows, 'VariableNames', ...
    {'esinr_dB','mcs','qm','r1024_curve','tbs','bgn','zc','C','cbBLER','tbBLER'});
writetable(Tbl, oracle_csv);
fprintf('Wrote %d oracle tuples to %s\n', height(Tbl), oracle_csv);

%% ---------------- Spot check vs legacy table -----------------------------
fprintf('\nSpot check (MCS 11, 64QAM R=466/1024, TBS=8000):\n');
for s = [0 4 8 12]
    r = Tbl(Tbl.mcs==11 & Tbl.tbs==8000 & Tbl.esinr_dB==s, :);
    if ~isempty(r)
        fprintf('  SINR=%+3d dB -> cbBLER=%.4e  tbBLER=%.4e (C=%d, Zc=%d)\n', ...
            s, r.cbBLER, r.tbBLER, r.C, r.zc);
    end
end
fprintf('Done.\n');
