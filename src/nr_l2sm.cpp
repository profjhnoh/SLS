#include "common.h"
#include "nr_l2sm.h"

#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>

// ============================================================================
// MCS -> (Qm, R*1024)  — the single source for the 28-entry MCS table
// (TS 38.214 Table 5.1.3.1-2 subset; values must stay identical to the legacy
//  CR_X1024 copies this replaces, including the fractional 682.5 / 916.5).
// ============================================================================

static const Real NR_CR_X1024[28] = {
	120, 193, 308, 449, 602,                      // MCS 0-4   QPSK
	378, 434, 490, 553, 616, 658,                 // MCS 5-10  16QAM
	466, 517, 567, 616, 666, 719, 772, 822, 873,  // MCS 11-19 64QAM
	682.5, 711, 754, 797, 841, 885, 916.5, 948 }; // MCS 20-27 256QAM

const int NR_CQI2MCS[16] = { 0,  0, 0, 2, 4, 6, 8, 11, 13, 16, 18, 21, 23, 25, 27, 27 };

void nr_mcs_to_qm_r(int mcs_idx, int& qm_bits, Real& r_x1024)
{
	if (mcs_idx < 0)  mcs_idx = 0;
	if (mcs_idx > 27) mcs_idx = 27;

	if      (mcs_idx <= 4)  qm_bits = QPSK;    // 2
	else if (mcs_idx <= 10) qm_bits = QAM16;   // 4
	else if (mcs_idx <= 19) qm_bits = QAM64;   // 6
	else                    qm_bits = QAM256;  // 8

	r_x1024 = NR_CR_X1024[mcs_idx];
}

// ============================================================================
// TS 38.214 §5.1.3.2 TBS quantization — shared by the per-layer TBS path in
// receive_downlink.cpp and by Regenerate_SINR_thresholds_from_matlab().
// The math replicates the legacy tbs_info_bits_from_mcs verbatim.
// ============================================================================

int nr_tbs_info_bits_from_mcs(int mcs_idx, int N_RE, int layer_factor)
{
	int  qm_bits;
	Real r_x1024;
	nr_mcs_to_qm_r(mcs_idx, qm_bits, r_x1024);
	Real coding_rate = r_x1024 / 1024.0;

	int N_info = (int)(N_RE * coding_rate * qm_bits * layer_factor * (1.0 - overhead));
	int info_bits = 0;
	if (N_info <= 3834)
	{
		int n = MAX(3, floor(log2(N_info)) - 6.);
		int N_info_a = MAX(24, pow(2, n) * floor(N_info / pow(2, n)));
		info_bits = N_info_a;
	}
	else
	{
		int n = floor(log2(N_info - 24)) - 5;
		int N_info_a = MAX(3840, pow(2, n) * round(((N_info - 24) / pow(2, n))));
		if (coding_rate <= (1. / 4.))
		{
			int C = ceil((N_info_a + 24) / 3816);
			info_bits = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
		}
		else if (N_info_a >= 8424)
		{
			int C = ceil((N_info_a + 24) / 8424);
			info_bits = 8 * C * ceil((N_info_a + 24) / (8 * C)) - 24;
		}
		else
		{
			info_bits = 8 * ceil((N_info_a + 24) / 8) - 24;
		}
	}
	return info_bits;
}

// ============================================================================
// TS 38.212 §5.2.2 code-block segmentation (nrDLSCHInfo equivalent)
// ============================================================================

NrSegmentInfo nr_dlsch_segment_info(int tbs_bits, double code_rate)
{
	NrSegmentInfo s;
	int A = tbs_bits;
	if (A < 24) A = 24;
	double R = code_rate;

	// TS 38.212 §7.2.2 base-graph selection
	s.bgn = (A <= 292 || (A <= 3824 && R <= 0.67) || R <= 0.25) ? 2 : 1;

	int       L   = (A > 3824) ? 24 : 16;   // TB CRC length
	long long B   = (long long)A + L;
	int       Kcb = (s.bgn == 1) ? 8448 : 3840;

	long long Bp;
	if (B <= Kcb) { s.C = 1;                                        Bp = B; }
	else          { s.C = (int)std::ceil((double)B / (Kcb - 24));   Bp = B + 24LL * s.C; }

	double Kp = (double)Bp / s.C;

	int Kb;
	if (s.bgn == 1) Kb = 22;
	else if (B > 640) Kb = 10;
	else if (B > 560) Kb = 9;
	else if (B > 192) Kb = 8;
	else              Kb = 6;

	// TS 38.212 Table 5.3.2-1 lifting sizes (ascending)
	static const int ZC[51] = {
		  2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
		 16,  18,  20,  22,  24,  26,  28,  30,  32,  36,  40,  44,  48,  52,
		 56,  60,  64,  72,  80,  88,  96, 104, 112, 120, 128, 144, 160, 176,
		192, 208, 224, 240, 256, 288, 320, 352, 384 };

	s.Zc = 384;
	for (int i = 0; i < 51; i++)
	{
		if ((double)Kb * ZC[i] >= Kp) { s.Zc = ZC[i]; break; }
	}
	s.K = ((s.bgn == 1) ? 22 : 10) * s.Zc;
	return s;
}

// ============================================================================
// Table storage + index. Immutable after Load (lookups run inside the OpenMP
// receive region — no lazy init, no mutable caches).
// ============================================================================

namespace {

struct L2smCurve
{
	int bgn, qm, r1024, zc;
	std::vector<double> snr;         // ascending
	std::vector<double> lg;          // log10(cbBLER); rows at/after firstZero unused
	int firstZero;                   // index of first row with cbBLER==0, or -1
};

struct RGroup
{
	int r1024;
	std::vector<std::pair<int,int> > zc_curves;  // (zc, curve index) sorted by zc
};

std::vector<L2smCurve>           g_curves;
std::map<int, std::vector<RGroup> > g_index;    // key = bgn*100 + qm; RGroups sorted by r1024
bool g_loaded = false;

// interpolatePER replica: truncate at first cbBLER==0 row (inclusive), linear
// interpolation of log10(BLER) vs SNR with two-sided extrapolation, then
// clamp <1e-6 -> 0 and >1 -> 1. Everything past the last non-zero point is 0.
double interp_per(const L2smCurve& c, double snr_dB)
{
	int nEnd = (c.firstZero >= 0) ? c.firstZero : (int)c.snr.size();  // rows [0, nEnd) are non-zero
	if (nEnd <= 0) return 0.0;                                        // curve starts at BLER 0
	if (c.firstZero >= 0 && snr_dB > c.snr[nEnd - 1]) return 0.0;     // interp toward log10(0) = -inf
	if (nEnd == 1) return (c.lg[0] >= 0.0) ? 1.0 : std::pow(10.0, c.lg[0]);

	int hi;
	if      (snr_dB <= c.snr[0])        hi = 1;             // left extrapolation
	else if (snr_dB >= c.snr[nEnd - 1]) hi = nEnd - 1;      // right extrapolation (no zero row)
	else
	{
		hi = (int)(std::upper_bound(c.snr.begin(), c.snr.begin() + nEnd, snr_dB) - c.snr.begin());
		if (hi < 1) hi = 1;
	}
	int lo = hi - 1;
	double t  = (snr_dB - c.snr[lo]) / (c.snr[hi] - c.snr[lo]);
	double lg = c.lg[lo] + t * (c.lg[hi] - c.lg[lo]);
	double per = std::pow(10.0, lg);
	if (per < 1e-6) per = 0.0;
	if (per > 1.0)  per = 1.0;
	return per;
}

// getAWGNTable replica on the flattened export: within (bgn, qm) pick the
// smallest tabulated R >= requested (fallback: largest tabulated R), then the
// exact Zc (fallback: nearest, tie -> larger — MATLAB would error here, but a
// missing Zc can only come from an incomplete export, so degrade gracefully).
const L2smCurve* select_curve(int bgn, int qm, int r1024_req, int zc_req, int* out_r1024)
{
	std::map<int, std::vector<RGroup> >::const_iterator it = g_index.find(bgn * 100 + qm);
	if (it == g_index.end() || it->second.empty()) return 0;
	const std::vector<RGroup>& groups = it->second;

	int gi = (int)groups.size() - 1;
	for (int i = 0; i < (int)groups.size(); i++)
	{
		if (groups[i].r1024 >= r1024_req) { gi = i; break; }
	}
	const RGroup& grp = groups[gi];
	if (out_r1024) *out_r1024 = grp.r1024;

	const std::vector<std::pair<int,int> >& zs = grp.zc_curves;
	int best = 0;
	int bestd = 1 << 30;
	for (int i = 0; i < (int)zs.size(); i++)
	{
		int d = zs[i].first - zc_req;
		if (d == 0) { best = i; break; }
		int ad = (d < 0) ? -d : d;
		if (ad < bestd || (ad == bestd && d > 0)) { bestd = ad; best = i; }
	}
	return &g_curves[zs[best].second];
}

} // namespace

bool Matlab_L2SM_loaded(void) { return g_loaded; }

bool Load_Matlab_L2SM_BLER(const char* path)
{
	std::ifstream in(path);
	if (!in.is_open())
	{
		cout << "=== matlab_l2sm_bler: cannot open " << path << " ===" << endl;
		return false;
	}

	g_curves.clear();
	g_index.clear();

	std::string line;
	long declared = -1;
	while (std::getline(in, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::istringstream ls(line);
		std::string tok;
		ls >> tok;
		if (tok == "MATLAB_L2SM_BLER") continue;
		if (tok == "curves") { ls >> declared; continue; }
		if (tok == "end") break;
		if (tok != "curve") continue;

		std::string kb, kq, kr, kz, kn;
		int bgn = 0, qm = 0, zc = 0, n = 0;
		double r = 0;
		ls >> kb >> bgn >> kq >> qm >> kr >> r >> kz >> zc >> kn >> n;
		if (n < 2 || bgn < 1 || bgn > 2) { cout << "=== matlab_l2sm_bler: bad curve header ===" << endl; return false; }

		L2smCurve c;
		c.bgn = bgn; c.qm = qm; c.zc = zc;
		c.r1024 = (int)(r + 0.5);
		c.snr.reserve(n); c.lg.reserve(n);
		c.firstZero = -1;
		for (int i = 0; i < n; i++)
		{
			double s, b;
			if (!(in >> s >> b)) { cout << "=== matlab_l2sm_bler: truncated curve data ===" << endl; return false; }
			c.snr.push_back(s);
			if (b <= 0.0)
			{
				if (c.firstZero < 0) c.firstZero = i;
				c.lg.push_back(-300.0);
			}
			else
			{
				c.lg.push_back(std::log10(b));
			}
		}
		std::getline(in, line);  // consume rest of last data line
		g_curves.push_back(c);
	}

	if (declared > 0 && (long)g_curves.size() != declared)
	{
		cout << "=== matlab_l2sm_bler: curve count mismatch (" << g_curves.size()
		     << " vs declared " << declared << ") ===" << endl;
		return false;
	}

	// Build (bgn,qm) -> sorted R groups -> Zc-sorted curve lists
	for (int idx = 0; idx < (int)g_curves.size(); idx++)
	{
		const L2smCurve& c = g_curves[idx];
		std::vector<RGroup>& groups = g_index[c.bgn * 100 + c.qm];
		int gi = -1;
		for (int i = 0; i < (int)groups.size(); i++)
			if (groups[i].r1024 == c.r1024) { gi = i; break; }
		if (gi < 0) { RGroup g; g.r1024 = c.r1024; groups.push_back(g); gi = (int)groups.size() - 1; }
		groups[gi].zc_curves.push_back(std::make_pair(c.zc, idx));
	}
	for (std::map<int, std::vector<RGroup> >::iterator it = g_index.begin(); it != g_index.end(); ++it)
	{
		std::sort(it->second.begin(), it->second.end(),
		          [](const RGroup& a, const RGroup& b) { return a.r1024 < b.r1024; });
		for (size_t i = 0; i < it->second.size(); i++)
			std::sort(it->second[i].zc_curves.begin(), it->second[i].zc_curves.end());
	}

	g_loaded = true;
	cout << "=== matlab_l2sm_bler: loaded " << g_curves.size() << " curves from " << path << " ===" << endl;
	return true;
}

// ============================================================================
// Lookups
// ============================================================================

double bler_lookup_matlab_core(double esinr_dB, int qm_bits, double r_x1024, int tbs_bits,
                               double* out_cb_bler, NrSegmentInfo* out_seg, int* out_curve_r1024)
{
	// ECR clip + integer R, exactly as sinrToCodeBLER does
	double ecr = r_x1024 / 1024.0;
	if (ecr < 1.0 / 1024.0)    ecr = 1.0 / 1024.0;
	if (ecr > 1023.0 / 1024.0) ecr = 1023.0 / 1024.0;
	int R = (int)(ecr * 1024.0 + 0.5);

	NrSegmentInfo seg = nr_dlsch_segment_info(tbs_bits, (double)R / 1024.0);
	if (out_seg) *out_seg = seg;

	const L2smCurve* c = select_curve(seg.bgn, qm_bits, R, seg.Zc, out_curve_r1024);
	if (!c)
	{
		if (out_cb_bler) *out_cb_bler = 1.0;
		return 1.0;  // no curve — pessimistic fallback
	}

	double cb = interp_per(*c, esinr_dB);
	if (out_cb_bler) *out_cb_bler = cb;
	if (cb <= 0.0) return 0.0;
	if (cb >= 1.0) return 1.0;
	return 1.0 - std::pow(1.0 - cb, (double)seg.C);
}

Real bler_lookup_matlab(Real esinr_linear, int mcs_idx, int tbs_bits)
{
	if (!g_loaded) return (Real)1.0;
	int  qm;
	Real r;
	nr_mcs_to_qm_r(mcs_idx, qm, r);
	if (tbs_bits < 24) tbs_bits = 24;
	double esinr_dB = 10.0 * std::log10((double)std::max(esinr_linear, (Real)1e-30));
	return (Real)bler_lookup_matlab_core(esinr_dB, qm, (double)r, tbs_bits);
}

// ============================================================================
// CQI threshold regeneration (Phase E): for each CQI's mapped MCS, find the
// SNR where TB BLER = 0.1 at a reference full-band single-layer allocation.
// ============================================================================

void Regenerate_SINR_thresholds_from_matlab(void)
{
	if (!g_loaded)
	{
		cout << "=== matlab_cqi_thresholds: table not loaded, keeping hardcoded thresholds ===" << endl;
		return;
	}

	const double TARGET = 0.1;
	int N_RE_a = num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms;
	int N_RE   = MIN(156, N_RE_a) * num_rb;

	cout << "=== matlab_cqi_thresholds: regenerating SINR_threshold_dB (ref N_RE=" << N_RE << ") ===" << endl;
	cout << "    CQI  MCS    TBS(ref)   old(dB)    new(dB)" << endl;

	Real prev = -1e9;
	for (int cqi = 1; cqi <= 15; cqi++)
	{
		int  m = NR_CQI2MCS[cqi];
		int  qm;
		Real r;
		nr_mcs_to_qm_r(m, qm, r);
		int tbs = nr_tbs_info_bits_from_mcs(m, N_RE, 1);

		// TB BLER is monotone non-increasing in SNR: bisect for BLER = TARGET
		double lo = -15.0, hi = 40.0;
		for (int it = 0; it < 60; it++)
		{
			double mid = 0.5 * (lo + hi);
			double b = bler_lookup_matlab_core(mid, qm, (double)r, tbs);
			if (b > TARGET) lo = mid; else hi = mid;
		}
		Real th = (Real)(0.5 * (lo + hi));

		// enforce strict monotonicity across CQI
		if (th <= prev + 0.01)
		{
			cout << "    [warn] CQI " << cqi << " threshold non-monotone, clamping" << endl;
			th = prev + 0.01;
		}

		printf("    %3d  %3d  %9d   %+7.2f    %+7.2f\n", cqi, m, tbs,
		       (double)SINR_threshold_dB[cqi - 1], (double)th);
		SINR_threshold_dB[cqi - 1] = th;
		prev = th;
	}
}

// ============================================================================
// RBIR effective-SINR mapping (calculateEffectiveSINR replica, alpha=beta=1)
// ============================================================================

namespace {

struct RbirTable
{
	double minSNR, maxSNR;
	std::vector<double> snr;   // dB, ascending
	std::vector<double> rbir;  // ascending (information bits/symbol)
};

std::map<int, RbirTable> g_rbir;   // key = modscheme (2^Qm: 4/16/64/256)
bool g_rbir_loaded = false;

double lin_interp(const std::vector<double>& x, const std::vector<double>& y, double xv)
{
	// x ascending; xv assumed within [x.front(), x.back()] (callers clamp)
	if (xv <= x.front()) return y.front();
	if (xv >= x.back())  return y.back();
	int hi = (int)(std::upper_bound(x.begin(), x.end(), xv) - x.begin());
	if (hi < 1) hi = 1;
	int lo = hi - 1;
	double t = (xv - x[lo]) / (x[hi] - x[lo]);
	return y[lo] + t * (y[hi] - y[lo]);
}

} // namespace

bool RBIR_ESM_loaded(void) { return g_rbir_loaded; }

bool Load_RBIR_ESM(const char* path)
{
	std::ifstream in(path);
	if (!in.is_open())
	{
		cout << "=== rbir_esm: cannot open " << path << " ===" << endl;
		return false;
	}

	g_rbir.clear();
	std::string line;
	while (std::getline(in, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::istringstream ls(line);
		std::string tok;
		ls >> tok;
		if (tok == "MATLAB_RBIR_ESM" || tok == "tables") continue;
		if (tok == "end") break;
		if (tok != "table") continue;

		std::string km, kmin, kmax, kn;
		int mod = 0, n = 0;
		double mn = 0, mx = 0;
		ls >> km >> mod >> kmin >> mn >> kmax >> mx >> kn >> n;
		if (n < 2) { cout << "=== rbir_esm: bad table header ===" << endl; return false; }

		RbirTable t;
		t.minSNR = mn; t.maxSNR = mx;
		t.snr.reserve(n); t.rbir.reserve(n);
		for (int i = 0; i < n; i++)
		{
			double s, r;
			if (!(in >> s >> r)) { cout << "=== rbir_esm: truncated table ===" << endl; return false; }
			t.snr.push_back(s);
			t.rbir.push_back(r);
		}
		std::getline(in, line);
		g_rbir[mod] = t;
	}

	g_rbir_loaded = !g_rbir.empty();
	cout << "=== rbir_esm: loaded " << g_rbir.size() << " modulation tables from " << path << " ===" << endl;
	return g_rbir_loaded;
}

double rbir_effective_sinr_dB(const double* sinr_dB, int n, int qm_bits)
{
	int mod = 1 << qm_bits;   // 2^Qm: 4/16/64/256
	std::map<int, RbirTable>::const_iterator it = g_rbir.find(mod);
	if (it == g_rbir.end() || n <= 0)
	{
		// no table — degrade to plain dB average (should not happen when loaded)
		double s = 0; for (int i = 0; i < n; i++) s += sinr_dB[i];
		return (n > 0) ? s / n : -30.0;
	}
	const RbirTable& t = it->second;

	double av = 0.0;
	for (int i = 0; i < n; i++)
	{
		double s = sinr_dB[i];
		if (s < t.minSNR) s = t.minSNR;
		if (s > t.maxSNR) s = t.maxSNR;
		av += lin_interp(t.snr, t.rbir, s);
	}
	av /= n;

	return lin_interp(t.rbir, t.snr, av);   // inverse map (rbir ascending)
}

Real rbir_effective_sinr_linear(const std::vector<Real>& sinr_linear, int qm_bits)
{
	int n = (int)sinr_linear.size();
	if (n <= 0) return (Real)0.0;
	std::vector<double> s_dB(n);
	for (int i = 0; i < n; i++)
		s_dB[i] = 10.0 * std::log10((double)std::max(sinr_linear[i], (Real)1e-30));
	double eff_dB = rbir_effective_sinr_dB(&s_dB[0], n, qm_bits);
	return (Real)std::pow(10.0, eff_dB / 10.0);
}

bool Matlab_RBIR_selftest(const char* tuples_csv)
{
	std::ifstream in(tuples_csv);
	if (!in.is_open())
	{
		cout << "=== rbir_selftest: cannot open " << tuples_csv << " ===" << endl;
		return false;
	}

	std::string line;
	std::getline(in, line);   // header

	int nRows = 0, nFail = 0;
	while (std::getline(in, line))
	{
		if (line.empty()) continue;
		// modscheme,vec_id,eff_sinr_dB,s1 s2 s3 ...   (list is comma-free)
		size_t c1 = line.find(',');
		size_t c2 = line.find(',', c1 + 1);
		size_t c3 = line.find(',', c2 + 1);
		if (c1 == std::string::npos || c2 == std::string::npos || c3 == std::string::npos) continue;
		int    mod   = atoi(line.substr(0, c1).c_str());
		double eff_o = atof(line.substr(c2 + 1, c3 - c2 - 1).c_str());
		std::istringstream vs(line.substr(c3 + 1));
		std::vector<double> s;
		double v;
		while (vs >> v) s.push_back(v);
		if (s.empty()) continue;
		nRows++;

		int qm = 0;
		while ((1 << (qm + 1)) <= mod) qm++;   // mod = 2^qm
		double eff_c = rbir_effective_sinr_dB(&s[0], (int)s.size(), qm);

		if (std::fabs(eff_c - eff_o) > 1e-3)
		{
			nFail++;
			if (nFail <= 5)
				printf("    [FAIL] mod=%d n=%d : got %.5f dB, want %.5f dB\n",
				       mod, (int)s.size(), eff_c, eff_o);
		}
	}

	cout << "=== rbir_selftest: " << (nRows - nFail) << "/" << nRows << " tuples pass ===" << endl;
	return (nRows > 0) && (nFail == 0);
}

// ============================================================================
// Throughput-maximizing MCS selection (precomputed dB-grid decision table)
// ============================================================================

namespace {

const double TPUT_GRID_LO   = -20.0;   // dB
const double TPUT_GRID_HI   =  40.0;   // dB
const double TPUT_GRID_STEP =  0.05;   // dB
const int    TPUT_GRID_N    = (int)((TPUT_GRID_HI - TPUT_GRID_LO) / TPUT_GRID_STEP) + 1;

std::vector<int> g_tput_mcs_grid;      // best MCS per grid point; empty = not built
int g_tput_nre_ref = 0;

// direct argmax at one effective SINR (dB)
int tput_argmax_direct(double esinr_dB, const int* tbs_per_mcs)
{
	int    best_m = 0;
	double best_s = -1.0;
	for (int m = 0; m < 28; m++)
	{
		int  qm;
		Real r;
		nr_mcs_to_qm_r(m, qm, r);
		double bler  = bler_lookup_matlab_core(esinr_dB, qm, (double)r, tbs_per_mcs[m]);
		double score = (double)tbs_per_mcs[m] * (1.0 - bler);
		if (score > best_s + 1e-9) { best_s = score; best_m = m; }   // tie -> lower MCS
	}
	return best_m;
}

} // namespace

bool TputMCS_Grid_ready(void) { return !g_tput_mcs_grid.empty(); }

void Build_TputMCS_Grid(int n_re_ref)
{
	if (!g_loaded)
	{
		cout << "=== tput_mcs: BLER table not loaded, grid not built ===" << endl;
		return;
	}
	g_tput_nre_ref = n_re_ref;

	// TBS per MCS at the reference allocation, single layer (matches the
	// per-layer receiver TBS convention: layer_factor = 1).
	int tbs_per_mcs[28];
	for (int m = 0; m < 28; m++)
		tbs_per_mcs[m] = nr_tbs_info_bits_from_mcs(m, n_re_ref, 1);

	g_tput_mcs_grid.assign(TPUT_GRID_N, 0);
	for (int i = 0; i < TPUT_GRID_N; i++)
	{
		double s = TPUT_GRID_LO + i * TPUT_GRID_STEP;
		g_tput_mcs_grid[i] = tput_argmax_direct(s, tbs_per_mcs);
	}

	// Sanity print: decision boundaries (where the best MCS switches)
	cout << "=== tput_mcs: grid built (N_RE_ref=" << n_re_ref << "). Decision boundaries: ===" << endl;
	int prev = g_tput_mcs_grid[0];
	int shown = 0;
	printf("    %.2f dB -> MCS %d\n", TPUT_GRID_LO, prev);
	for (int i = 1; i < TPUT_GRID_N && shown < 40; i++)
	{
		if (g_tput_mcs_grid[i] != prev)
		{
			prev = g_tput_mcs_grid[i];
			printf("    %+.2f dB -> MCS %d\n", TPUT_GRID_LO + i * TPUT_GRID_STEP, prev);
			shown++;
		}
	}
}

int nr_mcs_maxtput_from_dB(Real esinr_dB)
{
	if (g_tput_mcs_grid.empty()) return 0;
	double x = ((double)esinr_dB - TPUT_GRID_LO) / TPUT_GRID_STEP;
	int i = (int)(x + 0.5);
	if (i < 0) i = 0;
	if (i >= TPUT_GRID_N) i = TPUT_GRID_N - 1;
	return g_tput_mcs_grid[i];
}

// ============================================================================
// Selftest against the MATLAB-generated oracle
// ============================================================================

bool Matlab_BLER_selftest(const char* tuples_csv)
{
	std::ifstream in(tuples_csv);
	if (!in.is_open())
	{
		cout << "=== matlab_bler_selftest: cannot open " << tuples_csv << " ===" << endl;
		return false;
	}

	std::string line;
	std::getline(in, line);  // header

	int nRows = 0, nFail = 0;
	while (std::getline(in, line))
	{
		if (line.empty()) continue;
		for (size_t i = 0; i < line.size(); i++) if (line[i] == ',') line[i] = ' ';
		std::istringstream ls(line);
		double esinr, r1024_curve, cb_o, tb_o;
		int mcs, qm_o, tbs, bgn_o, zc_o, C_o;
		if (!(ls >> esinr >> mcs >> qm_o >> r1024_curve >> tbs >> bgn_o >> zc_o >> C_o >> cb_o >> tb_o)) continue;
		nRows++;

		int  qm;
		Real r;
		nr_mcs_to_qm_r(mcs, qm, r);

		double cb = 0;
		NrSegmentInfo seg;
		int curve_r = 0;
		double tb = bler_lookup_matlab_core(esinr, qm, (double)r, tbs, &cb, &seg, &curve_r);

		bool ok = (qm == qm_o) && (seg.bgn == bgn_o) && (seg.Zc == zc_o) && (seg.C == C_o)
		          && (curve_r == (int)(r1024_curve + 0.5));
		if (ok)
		{
			// clamp regions must match exactly; transition region in log10 domain
			if      (cb_o <= 0.0) ok = (cb <= 0.0);
			else if (cb_o >= 1.0) ok = (cb >= 1.0);
			else                  ok = (cb > 0.0) && (std::fabs(std::log10(cb) - std::log10(cb_o)) < 1e-3);
			if (ok)
			{
				if      (tb_o <= 0.0) ok = (tb <= 0.0);
				else if (tb_o >= 1.0) ok = (tb >= 1.0);
				else                  ok = (tb > 0.0) && (std::fabs(std::log10(tb) - std::log10(tb_o)) < 1e-3);
			}
		}

		if (!ok)
		{
			nFail++;
			if (nFail <= 5)
			{
				printf("    [FAIL] esinr=%+.1f mcs=%d tbs=%d : got bgn=%d zc=%d C=%d r=%d cb=%.4e tb=%.4e"
				       " / want bgn=%d zc=%d C=%d r=%d cb=%.4e tb=%.4e\n",
				       esinr, mcs, tbs, seg.bgn, seg.Zc, seg.C, curve_r, cb, tb,
				       bgn_o, zc_o, C_o, (int)(r1024_curve + 0.5), cb_o, tb_o);
			}
		}
	}

	cout << "=== matlab_bler_selftest: " << (nRows - nFail) << "/" << nRows << " tuples pass ===" << endl;
	return (nRows > 0) && (nFail == 0);
}
