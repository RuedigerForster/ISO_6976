// iso6976.cpp
//
// Clean-room implementation of ISO 6976:2016
// "Natural Gas — Calculation of calorific values, density, relative density
//  and Wobbe indices from composition"
//
// Formula numbers in comments refer to ISO 6976:2016.
// Annex B formula numbers refer to the uncertainty propagation annex.
//
// Copyright (C) 2026 Paul Lysakowski
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include <stdexcept>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

// ============================================================================
// Physical constants
// ============================================================================

static const double P_REF   = 101.325;    // reference pressure [kPa]
static const double T_ZERO  = 273.15;     // 0 °C in Kelvin
static const double R_MOL   = 8.3144621;  // molar gas constant [J/(mol·K)]
static const double U_R_MOL = 0.0000075;  // std. uncertainty of R

static const double M_AIR   = 28.96546;   // molar mass of dry air [kg/kmol], Table A.1
static const double U_M_AIR = 0.00017;

// ============================================================================
// ISO 6976:2016 Table A.2 — molar masses and atomic composition
//
// Columns: M [kg/kmol], nC, nH, nN, nO, nS, nHe, nNe, nAr
// Row order matches ISO 6976:2016 Table A.2, components 1–60
// ============================================================================

static const double TAB_M[60][9] = {
  { 16.04246,  1,  4, 0, 0, 0, 0, 0, 0 },  //  1 methane
  { 30.06904,  2,  6, 0, 0, 0, 0, 0, 0 },  //  2 ethane
  { 44.09562,  3,  8, 0, 0, 0, 0, 0, 0 },  //  3 propane
  { 58.12220,  4, 10, 0, 0, 0, 0, 0, 0 },  //  4 n-butane
  { 58.12220,  4, 10, 0, 0, 0, 0, 0, 0 },  //  5 isobutane
  { 72.14878,  5, 12, 0, 0, 0, 0, 0, 0 },  //  6 n-pentane
  { 72.14878,  5, 12, 0, 0, 0, 0, 0, 0 },  //  7 isopentane
  { 72.14878,  5, 12, 0, 0, 0, 0, 0, 0 },  //  8 neopentane
  { 86.17536,  6, 14, 0, 0, 0, 0, 0, 0 },  //  9 n-hexane
  { 86.17536,  6, 14, 0, 0, 0, 0, 0, 0 },  // 10 2-methylpentane
  { 86.17536,  6, 14, 0, 0, 0, 0, 0, 0 },  // 11 3-methylpentane
  { 86.17536,  6, 14, 0, 0, 0, 0, 0, 0 },  // 12 2,2-dimethylbutane
  { 86.17536,  6, 14, 0, 0, 0, 0, 0, 0 },  // 13 2,3-dimethylbutane
  { 100.20194, 7, 16, 0, 0, 0, 0, 0, 0 },  // 14 n-heptane
  { 114.22852, 8, 18, 0, 0, 0, 0, 0, 0 },  // 15 n-octane
  { 128.25510, 9, 20, 0, 0, 0, 0, 0, 0 },  // 16 n-nonane
  { 142.28168,10, 22, 0, 0, 0, 0, 0, 0 },  // 17 n-decane
  { 28.05316,  2,  4, 0, 0, 0, 0, 0, 0 },  // 18 ethylene
  { 42.07974,  3,  6, 0, 0, 0, 0, 0, 0 },  // 19 propylene
  { 56.10632,  4,  8, 0, 0, 0, 0, 0, 0 },  // 20 1-butene
  { 56.10632,  4,  8, 0, 0, 0, 0, 0, 0 },  // 21 cis-2-butene
  { 56.10632,  4,  8, 0, 0, 0, 0, 0, 0 },  // 22 trans-2-butene
  { 56.10632,  4,  8, 0, 0, 0, 0, 0, 0 },  // 23 isobutylene
  { 70.13290,  5, 10, 0, 0, 0, 0, 0, 0 },  // 24 1-pentene
  { 40.06386,  3,  4, 0, 0, 0, 0, 0, 0 },  // 25 propadiene
  { 54.09044,  4,  6, 0, 0, 0, 0, 0, 0 },  // 26 1,2-butadiene
  { 54.09044,  4,  6, 0, 0, 0, 0, 0, 0 },  // 27 1,3-butadiene
  { 26.03728,  2,  2, 0, 0, 0, 0, 0, 0 },  // 28 acetylene
  { 70.13290,  5, 10, 0, 0, 0, 0, 0, 0 },  // 29 cyclopentane
  { 84.15948,  6, 12, 0, 0, 0, 0, 0, 0 },  // 30 methylcyclopentane
  { 98.18606,  7, 14, 0, 0, 0, 0, 0, 0 },  // 31 ethylcyclopentane
  { 84.15948,  6, 12, 0, 0, 0, 0, 0, 0 },  // 32 cyclohexane
  { 98.18606,  7, 14, 0, 0, 0, 0, 0, 0 },  // 33 methylcyclohexane
  { 112.21264, 8, 16, 0, 0, 0, 0, 0, 0 },  // 34 ethylcyclohexane
  { 78.11184,  6,  6, 0, 0, 0, 0, 0, 0 },  // 35 benzene
  { 92.13842,  7,  8, 0, 0, 0, 0, 0, 0 },  // 36 toluene
  { 106.16500, 8, 10, 0, 0, 0, 0, 0, 0 },  // 37 ethylbenzene
  { 106.16500, 8, 10, 0, 0, 0, 0, 0, 0 },  // 38 o-xylene
  { 32.04186,  1,  4, 0, 1, 0, 0, 0, 0 },  // 39 methanol
  { 48.10746,  1,  4, 0, 0, 1, 0, 0, 0 },  // 40 methanethiol
  { 2.01588,   0,  2, 0, 0, 0, 0, 0, 0 },  // 41 hydrogen
  { 18.01528,  0,  2, 0, 1, 0, 0, 0, 0 },  // 42 water
  { 34.08088,  0,  2, 0, 0, 1, 0, 0, 0 },  // 43 hydrogen sulphide
  { 17.03052,  0,  3, 1, 0, 0, 0, 0, 0 },  // 44 ammonia
  { 27.02534,  1,  1, 1, 0, 0, 0, 0, 0 },  // 45 hydrogen cyanide
  { 28.0101,   1,  0, 0, 1, 0, 0, 0, 0 },  // 46 carbon monoxide
  { 60.0751,   1,  0, 0, 1, 1, 0, 0, 0 },  // 47 carbonyl sulphide
  { 76.1407,   1,  0, 0, 0, 2, 0, 0, 0 },  // 48 carbon disulphide
  { 4.002602,  0,  0, 0, 0, 0, 1, 0, 0 },  // 49 helium
  { 20.1797,   0,  0, 0, 0, 0, 0, 1, 0 },  // 50 neon
  { 39.948,    0,  0, 0, 0, 0, 0, 0, 1 },  // 51 argon
  { 28.0134,   0,  0, 2, 0, 0, 0, 0, 0 },  // 52 nitrogen
  { 31.9988,   0,  0, 0, 2, 0, 0, 0, 0 },  // 53 oxygen
  { 44.0095,   1,  0, 0, 2, 0, 0, 0, 0 },  // 54 carbon dioxide
  { 64.0638,   0,  0, 0, 2, 1, 0, 0, 0 },  // 55 sulphur dioxide
  { 156.30826,11, 24, 0, 0, 0, 0, 0, 0 },  // 56 n-undecane
  { 170.33484,12, 26, 0, 0, 0, 0, 0, 0 },  // 57 n-dodecane
  { 184.36142,13, 28, 0, 0, 0, 0, 0, 0 },  // 58 n-tridecane
  { 198.38800,14, 30, 0, 0, 0, 0, 0, 0 },  // 59 n-tetradecane
  { 212.41458,15, 32, 0, 0, 0, 0, 0, 0 },  // 60 n-pentadecane
};

// ============================================================================
// ISO 6976:2016 Table A.3 — summing factors s_i
//
// Columns: s at 0 °C, 15 °C, 15.55 °C, 20 °C, uncertainty u(s)
// Used for compression factor Z, Eq. (1)
// ============================================================================

static const double TAB_S[60][5] = {
  { 0.04886, 0.04452, 0.04437, 0.04317, 0.0005 },
  { 0.0997,  0.0919,  0.0916,  0.0895,  0.0011 },
  { 0.1465,  0.1344,  0.1340,  0.1308,  0.0016 },
  { 0.2022,  0.1840,  0.1834,  0.1785,  0.0039 },
  { 0.1885,  0.1722,  0.1717,  0.1673,  0.0031 },
  { 0.2586,  0.2361,  0.2354,  0.2295,  0.0107 },
  { 0.2458,  0.2251,  0.2244,  0.2189,  0.0088 },
  { 0.2245,  0.2040,  0.2033,  0.1979,  0.0060 },
  { 0.3319,  0.3001,  0.2990,  0.2907,  0.0271 },
  { 0.3114,  0.2826,  0.2816,  0.2740,  0.0221 },
  { 0.2997,  0.2762,  0.2754,  0.2690,  0.0234 },
  { 0.2530,  0.2350,  0.2344,  0.2295,  0.0173 },
  { 0.2836,  0.2632,  0.2625,  0.2569,  0.0207 },
  { 0.4076,  0.3668,  0.3654,  0.3547,  0.1001 },
  { 0.4845,  0.4346,  0.4329,  0.4198,  0.1002 },
  { 0.5617,  0.5030,  0.5010,  0.4856,  0.1006 },
  { 0.6713,  0.5991,  0.5967,  0.5778,  0.1006 },
  { 0.0868,  0.0799,  0.0797,  0.0778,  0.0010 },
  { 0.1381,  0.1267,  0.1263,  0.1232,  0.0016 },
  { 0.1964,  0.1776,  0.1770,  0.1721,  0.0041 },
  { 0.2075,  0.1870,  0.1863,  0.1810,  0.0045 },
  { 0.2072,  0.1868,  0.1862,  0.1809,  0.0043 },
  { 0.1966,  0.1777,  0.1770,  0.1721,  0.0037 },
  { 0.2622,  0.2297,  0.2287,  0.2208,  0.0102 },
  { 0.1417,  0.1313,  0.1310,  0.1282,  0.0025 },
  { 0.2063,  0.1862,  0.1855,  0.1803,  0.0110 },
  { 0.1993,  0.1739,  0.1731,  0.1673,  0.0038 },
  { 0.0936,  0.0836,  0.0833,  0.0808,  0.0024 },
  { 0.2409,  0.2221,  0.2215,  0.2164,  0.0137 },
  { 0.2817,  0.2612,  0.2605,  0.2548,  0.0262 },
  { 0.4227,  0.3684,  0.3666,  0.3531,  0.1006 },
  { 0.2939,  0.2686,  0.2677,  0.2610,  0.0325 },
  { 0.3667,  0.3317,  0.3305,  0.3213,  0.0668 },
  { 0.5275,  0.4547,  0.4524,  0.4345,  0.1006 },
  { 0.2752,  0.2527,  0.2520,  0.2460,  0.0274 },
  { 0.3726,  0.3359,  0.3347,  0.3251,  0.1002 },
  { 0.4129,  0.3797,  0.3785,  0.3694,  0.1002 },
  { 0.4852,  0.4411,  0.4396,  0.4277,  0.1004 },
  { 0.5806,  0.4464,  0.4423,  0.4117,  0.0233 },
  { 0.1909,  0.1700,  0.1693,  0.1640,  0.0117 },
  { -0.01,   -0.01,   -0.01,   -0.01,   0.0250 },
  { 0.3093,  0.2562,  0.2546,  0.2419,  0.0150 },
  { 0.1006,  0.0923,  0.0920,  0.0898,  0.0023 },
  { 0.1230,  0.1100,  0.1096,  0.1062,  0.0021 },
  { 0.3175,  0.2765,  0.2751,  0.2644,  0.0076 },
  { 0.0258,  0.0217,  0.0215,  0.0203,  0.0010 },
  { 0.1211,  0.1114,  0.1110,  0.1084,  0.0054 },
  { 0.2182,  0.1958,  0.1951,  0.1894,  0.0098 },
  { -0.01,   -0.01,   -0.01,   -0.01,   0.0250 },
  { -0.01,   -0.01,   -0.01,   -0.01,   0.0250 },
  { 0.0307,  0.0273,  0.0272,  0.0262,  0.0010 },
  { 0.0214,  0.0170,  0.0169,  0.0156,  0.0010 },
  { 0.0311,  0.0276,  0.0275,  0.0265,  0.0010 },
  { 0.0821,  0.0752,  0.0749,  0.0730,  0.0020 },
  { 0.1579,  0.1406,  0.1400,  0.1356,  0.0035 },
  { 0.7228,  0.6402,  0.6374,  0.6159,  0.1006 },
  { 0.8567,  0.7615,  0.7583,  0.7335,  0.1006 },
  { 0.9129,  0.8061,  0.8026,  0.7748,  0.1006 },
  { 1.0135,  0.8940,  0.8900,  0.8589,  0.1006 },
  { 1.1176,  0.9849,  0.9804,  0.9459,  0.1006 },
};

// Temperature indices for summing-factor table (0 °C, 15 °C, 15.55 °C, 20 °C)
static const double T_S[4]  = { 0.0, 15.0, 15.55, 20.0 };
// Z_air at the same temperatures (Table A.1) and their uncertainties
static const double Z_AIR[4]   = { 0.999419, 0.999595, 0.999601, 0.999645 };
static const double U_Z_AIR[4] = { 0.000015, 0.000015, 0.000015, 0.000015 };

// ============================================================================
// ISO 6976:2016 Table A.4 — ideal-gas gross calorific values H_ch,i^o [kJ/mol]
//
// Columns: at 0 °C, 15 °C, 15.55 °C, 20 °C, 25 °C, uncertainty u(H_ch)
// ============================================================================

static const double TAB_HC[60][6] = {
  {  892.92,  891.51,  891.46,  891.05,  890.58, 0.19 },
  { 1564.35, 1562.14, 1562.06, 1561.42, 1560.69, 0.51 },
  { 2224.03, 2221.10, 2220.99, 2220.13, 2219.17, 0.51 },
  { 2883.35, 2879.76, 2879.63, 2878.58, 2877.40, 0.72 },
  { 2874.21, 2870.58, 2870.45, 2869.39, 2868.20, 0.72 },
  { 3542.91, 3538.60, 3538.45, 3537.19, 3535.77, 0.23 },
  { 3536.01, 3531.68, 3531.52, 3530.25, 3528.83, 0.23 },
  { 3521.75, 3517.44, 3517.28, 3516.02, 3514.61, 0.25 },
  { 4203.24, 4198.24, 4198.06, 4196.60, 4194.95, 0.32 },
  { 4195.64, 4190.62, 4190.44, 4188.97, 4187.32, 0.53 },
  { 4198.27, 4193.22, 4193.04, 4191.56, 4189.90, 0.53 },
  { 4185.86, 4180.83, 4180.65, 4179.17, 4177.52, 0.48 },
  { 4193.68, 4188.61, 4188.43, 4186.94, 4185.28, 0.46 },
  { 4862.88, 4857.18, 4856.98, 4855.31, 4853.43, 0.67 },
  { 5522.41, 5516.01, 5515.78, 5513.90, 5511.80, 0.76 },
  { 6182.92, 6175.82, 6175.56, 6173.48, 6171.15, 0.81 },
  { 6842.69, 6834.90, 6834.62, 6832.33, 6829.77, 0.87 },
  { 1413.55, 1412.12, 1412.07, 1411.65, 1411.18, 0.21 },
  { 2061.57, 2059.43, 2059.35, 2058.73, 2058.02, 0.34 },
  { 2721.57, 2718.71, 2718.60, 2717.76, 2716.82, 0.39 },
  { 2714.88, 2711.94, 2711.83, 2710.97, 2710.00, 0.50 },
  { 2711.09, 2708.26, 2708.16, 2707.33, 2706.40, 0.47 },
  { 2704.88, 2702.06, 2701.96, 2701.13, 2700.20, 0.42 },
  { 3381.32, 3377.76, 3377.63, 3376.59, 3375.42, 0.73 },
  { 1945.26, 1943.97, 1943.92, 1943.54, 1943.11, 0.60 },
  { 2597.15, 2595.12, 2595.05, 2594.46, 2593.79, 0.40 },
  { 2544.14, 2542.11, 2542.03, 2541.44, 2540.77, 0.41 },
  { 1301.86, 1301.37, 1301.35, 1301.21, 1301.05, 0.32 },
  { 3326.14, 3322.19, 3322.05, 3320.89, 3319.59, 0.36 },
  { 3977.05, 3972.46, 3972.29, 3970.95, 3969.44, 0.56 },
  { 4637.20, 4631.93, 4631.74, 4630.20, 4628.47, 0.71 },
  { 3960.68, 3956.02, 3955.85, 3954.49, 3952.96, 0.32 },
  { 4609.33, 4604.08, 4603.89, 4602.36, 4600.64, 0.71 },
  { 5272.76, 5266.90, 5266.69, 5264.97, 5263.05, 0.95 },
  { 3305.12, 3302.90, 3302.81, 3302.16, 3301.43, 0.27 },
  { 3952.77, 3949.83, 3949.72, 3948.86, 3947.89, 0.51 },
  { 4613.16, 4609.54, 4609.40, 4608.34, 4607.15, 0.66 },
  { 4602.18, 4598.64, 4598.52, 4597.48, 4596.31, 0.76 },
  {  766.60,  765.09,  765.03,  764.59,  764.09, 0.13 },
  { 1241.64, 1240.28, 1240.23, 1239.84, 1239.39, 0.32 },
  {  286.64,  286.15,  286.13,  285.99,  285.83, 0.02 },
  {   45.064,  44.431,  44.408,  44.222,  44.013, 0.004 },
  {  562.93,  562.38,  562.36,  562.19,  562.01, 0.23 },
  {  384.57,  383.51,  383.47,  383.16,  382.81, 0.18 },
  {  671.92,  671.67,  671.66,  671.58,  671.50, 1.26 },
  {  282.80,  282.91,  282.91,  282.95,  282.98, 0.06 },
  {  548.01,  548.14,  548.15,  548.19,  548.23, 0.24 },
  { 1104.05, 1104.32, 1104.33, 1104.40, 1104.49, 0.43 },
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 49 helium (no combustion)
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 50 neon
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 51 argon
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 52 nitrogen
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 53 oxygen
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 54 carbon dioxide
  {    0.0,    0.0,    0.0,    0.0,    0.0,  0.0 },  // 55 sulphur dioxide
  { 7502.22, 7493.73, 7493.42, 7490.93, 7488.14, 1.54 },
  { 8162.43, 8153.24, 8152.91, 8150.21, 8147.19, 1.13 },
  { 8821.88, 8811.99, 8811.63, 8808.73, 8805.48, 1.21 },
  { 9481.71, 9471.12, 9470.73, 9467.63, 9464.15, 1.32 },
  {10141.65,10130.23,10129.82,10126.52,10122.82, 1.44 },
};

// Temperature indices for calorific value table (0 °C, 15 °C, 15.55 °C, 20 °C, 25 °C)
static const double T_HC[5] = { 0.0, 15.0, 15.55, 20.0, 25.0 };

// Standard enthalpy of vaporisation of water at T_HC temperatures [kJ/mol]
// and uncertainty (ISO 6976:2016 Table A.4 footnote / used in Eq. (3))
static const double L_VAP[5]   = { 45.064, 44.431, 44.408, 44.222, 44.013 };
static const double U_L_VAP[5] = { 0.004,  0.004,  0.004,  0.004,  0.004  };

// ============================================================================
// Molar mass uncertainties u(M_i) and inter-component correlations r(M_i,M_j)
// ISO 6976:2016 Eqs. (24)–(25)
//
// Computed once via double-checked locking on first call.
// ============================================================================

static double u_M[60];
static double r_M[60][60];
static bool   tables_ready = false;

static void init_M_tables() {
  if (tables_ready) return;

  // Atomic mass standard uncertainties [kg/kmol] from ISO 6976:2016 Eq. (25)
  const double u_aC  = 0.0004,   u_aH  = 0.000035, u_aN  = 0.0001;
  const double u_aO  = 0.00015,  u_aS  = 0.0025;
  const double u_aHe = 0.000001, u_aNe = 0.0003,   u_aAr = 0.0005;

  // u(M_i) — Eq. (25)
  for (int i = 0; i < 60; i++) {
    double nc  = TAB_M[i][1], nh  = TAB_M[i][2], nn  = TAB_M[i][3];
    double no  = TAB_M[i][4], ns  = TAB_M[i][5];
    double nhe = TAB_M[i][6], nne = TAB_M[i][7], nar = TAB_M[i][8];
    u_M[i] = std::sqrt(
      nc*nc*u_aC*u_aC + nh*nh*u_aH*u_aH + nn*nn*u_aN*u_aN +
      no*no*u_aO*u_aO + ns*ns*u_aS*u_aS +
      nhe*nhe*u_aHe*u_aHe + nne*nne*u_aNe*u_aNe + nar*nar*u_aAr*u_aAr
    );
  }

  // r(M_i, M_j) — Eq. (24)
  for (int i = 0; i < 60; i++) {
    for (int j = 0; j < 60; j++) {
      if (u_M[i] == 0.0 || u_M[j] == 0.0) { r_M[i][j] = 0.0; continue; }
      double cov =
        TAB_M[i][1]*TAB_M[j][1]*u_aC*u_aC +
        TAB_M[i][2]*TAB_M[j][2]*u_aH*u_aH +
        TAB_M[i][3]*TAB_M[j][3]*u_aN*u_aN +
        TAB_M[i][4]*TAB_M[j][4]*u_aO*u_aO +
        TAB_M[i][5]*TAB_M[j][5]*u_aS*u_aS +
        TAB_M[i][6]*TAB_M[j][6]*u_aHe*u_aHe +
        TAB_M[i][7]*TAB_M[j][7]*u_aNe*u_aNe +
        TAB_M[i][8]*TAB_M[j][8]*u_aAr*u_aAr;
      r_M[i][j] = cov / (u_M[i] * u_M[j]);
    }
  }
  tables_ready = true;
}

// ============================================================================
// Table index helpers
// ============================================================================

// Returns column index for summing-factor and Z_air tables (t2 in {0,15,15.55,20})
static int idx_s(double t) {
  if (t > 15.54 && t < 15.56) t = 15.55;
  for (int i = 0; i < 4; i++) if (t == T_S[i]) return i;
  return -1;
}

// Returns column index for calorific value table (t1 in {0,15,15.55,20,25})
static int idx_hc(double t) {
  if (t > 15.54 && t < 15.56) t = 15.55;
  for (int i = 0; i < 5; i++) if (t == T_HC[i]) return i;
  return -1;
}

// ============================================================================
// Core property calculations
// ============================================================================

// Eq. (1): compression factor Z
static double Z_val(const NumericVector& x, double p2, int ki) {
  double sum = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum += x[i] * TAB_S[i][ki];
  }
  return 1.0 - p2 / P_REF * sum * sum;  // Eq. (1)
}

// Eq. (2): ideal-gas gross calorific value H_ch^o [kJ/mol]
static double Hc_o_G(const NumericVector& x, int kh) {
  double sum = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum += x[i] * TAB_HC[i][kh];
  }
  return sum;
}

// Eq. (3): ideal-gas net calorific value H_cn^o [kJ/mol]
static double Hc_o_N(const NumericVector& x, int kh) {
  double H_gross = Hc_o_G(x, kh);
  double L = L_VAP[kh];
  double sum_nH = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum_nH += x[i] * TAB_M[i][2];  // nH atoms
  }
  return H_gross - sum_nH / 2.0 * L;
}

// Eq. (5): molar mass M [kg/kmol]
static double M_val(const NumericVector& x) {
  double sum = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum += x[i] * TAB_M[i][0];
  }
  return sum;
}

// Eq. (8): ideal molar volume V_o [m^3/mol]  (used internally; = R*T/p)
static double V_o(double t2, double p2) {
  return R_MOL * (t2 + T_ZERO) / p2;
}

// ============================================================================
// Uncertainty calculations (ISO 6976:2016 Annex B)
// ============================================================================

// Intermediate: sqrt((1-Z)*p2/p_ref), Eq. (B.7)
static double sigma_Z(double Z, double p2) {
  return std::sqrt((1.0 - Z) * p2 / P_REF);
}

// Annex B, Eq. (B.4): u(H_ch^o)
static double u_Hc_o_G(const NumericVector& x, const NumericVector& u,
                        const NumericMatrix& r, int kh) {
  // Composition-driven term
  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      s1 += TAB_HC[i][kh] * u[i] * r(i,j) * TAB_HC[j][kh] * u[j];
    }
  }
  // Table-value uncertainty term
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  return std::sqrt(s1 + s2);
}

// Annex B, Eq. (B.5): u(H_cm^o) — mass basis gross CV
static double u_Hm_o_G(const NumericVector& x, const NumericVector& u,
                        const NumericMatrix& r, int kh) {
  init_M_tables();
  double Hg = Hc_o_G(x, kh);
  double M  = M_val(x);
  double Hm = Hg / M;

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = TAB_HC[i][kh] / Hg - TAB_M[i][0] / M;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = TAB_HC[j][kh] / Hg - TAB_M[j][0] / M;
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (x[j] == 0.0 || r_M[i][j] == 0.0) continue;
      s3 += x[i]*u_M[i]*r_M[i][j]*x[j]*u_M[j];
    }
  }
  return std::sqrt(s1 + s2 / (Hg*Hg) + s3 / (M*M)) * Hm;
}

// Annex B, Eq. (B.6): u(H_cv^o) — volumetric basis gross CV (ideal gas)
static double u_Hv_o_G(const NumericVector& x, const NumericVector& u,
                        const NumericMatrix& r, int kh, int ks,
                        double Z, double p2) {
  double Hg  = Hc_o_G(x, kh);
  double Vox = V_o(T_S[ks], p2);
  double Hv  = Hg / Vox;        // ideal-gas volumetric GCV = Hc_o_G / V_o
  double sg  = sigma_Z(Z, p2);

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = TAB_HC[i][kh] / Hg + 2.0 * TAB_S[i][ks] * sg / Z;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = TAB_HC[j][kh] / Hg + 2.0 * TAB_S[j][ks] * sg / Z;
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s3 += x[i]*x[i] * TAB_S[i][4]*TAB_S[i][4];
  }
  double tmp = s1 + s2/(Hg*Hg) + 4.0*sg*sg*s3/(Z*Z) +
               (U_R_MOL/R_MOL)*(U_R_MOL/R_MOL);
  return std::sqrt(tmp) * Hv;
}

// Annex B, Eq. (B.8): u(H_cn^o)
static double u_Hc_o_N(const NumericVector& x, const NumericVector& u,
                        const NumericMatrix& r, int kh) {
  double L    = L_VAP[kh];
  double uL   = U_L_VAP[kh];

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = TAB_HC[i][kh] - L / 2.0 * TAB_M[i][2];
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = TAB_HC[j][kh] - L / 2.0 * TAB_M[j][2];
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double sum_nH = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum_nH += x[i] * TAB_M[i][2];
  }
  return std::sqrt(s1 + s2 + (sum_nH/2.0)*(sum_nH/2.0)*uL*uL);
}

// Annex B, Eq. (B.9): u(H_cm^n) — mass basis net CV
static double u_Hm_o_N(const NumericVector& x, const NumericVector& u,
                        const NumericMatrix& r, int kh) {
  init_M_tables();
  double Hn  = Hc_o_N(x, kh);
  double M   = M_val(x);
  double Hmn = Hn / M;
  double L   = L_VAP[kh];
  double uL  = U_L_VAP[kh];

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = (TAB_HC[i][kh] - L/2.0*TAB_M[i][2]) / Hn - TAB_M[i][0] / M;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = (TAB_HC[j][kh] - L/2.0*TAB_M[j][2]) / Hn - TAB_M[j][0] / M;
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (x[j] == 0.0 || r_M[i][j] == 0.0) continue;
      s3 += x[i]*u_M[i]*r_M[i][j]*x[j]*u_M[j];
    }
  }
  double sum_nH = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum_nH += x[i] * TAB_M[i][2];
  }
  double tmp = s1 + s2/(Hn*Hn) + s3/(M*M) +
               (sum_nH/(2.0*Hn))*(sum_nH/(2.0*Hn))*uL*uL;
  return std::sqrt(tmp) * Hmn;
}

// Annex B, Eq. (B.10): u(H_cv^n) — volumetric net CV (ideal gas)
static double u_Hv_o_N(const NumericVector& x, const NumericVector& u,
                        const NumericMatrix& r, int kh, int ks,
                        double Z, double p2) {
  double Hn  = Hc_o_N(x, kh);
  double Vox = V_o(T_S[ks], p2);
  double Hvn = Hn / Vox;
  double L   = L_VAP[kh];
  double uL  = U_L_VAP[kh];
  double sg  = sigma_Z(Z, p2);

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = (TAB_HC[i][kh] - L/2.0*TAB_M[i][2]) / Hn + 2.0*TAB_S[i][ks]*sg/Z;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = (TAB_HC[j][kh] - L/2.0*TAB_M[j][2]) / Hn + 2.0*TAB_S[j][ks]*sg/Z;
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s3 += x[i]*x[i] * TAB_S[i][4]*TAB_S[i][4];
  }
  double sum_nH = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum_nH += x[i] * TAB_M[i][2];
  }
  double tmp = s1 + s2/(Hn*Hn) + 4.0*sg*sg*s3/(Z*Z) +
               (U_R_MOL/R_MOL)*(U_R_MOL/R_MOL) +
               (sum_nH/(2.0*Hn))*(sum_nH/(2.0*Hn))*uL*uL;
  return std::sqrt(tmp) * Hvn;
}

// Annex B, Eq. (B.11): u(rho) — density (real gas)
static double u_D(const NumericVector& x, const NumericVector& u,
                  const NumericMatrix& r, int ks,
                  double Z, double D, double p2) {
  init_M_tables();
  double M  = M_val(x);
  double sg = sigma_Z(Z, p2);

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = TAB_M[i][0] / M + 2.0 * TAB_S[i][ks] * sg / Z;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = TAB_M[j][0] / M + 2.0 * TAB_S[j][ks] * sg / Z;
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (x[j] == 0.0 || r_M[i][j] == 0.0) continue;
      s2 += x[i]*u_M[i]*r_M[i][j]*x[j]*u_M[j];
    }
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s3 += x[i]*x[i] * TAB_S[i][4]*TAB_S[i][4];
  }
  double tmp = s1 + s2/(M*M) + 4.0*sg*sg*s3/(Z*Z) +
               (U_R_MOL/R_MOL)*(U_R_MOL/R_MOL);
  return std::sqrt(tmp) * D;
}

// Annex B, Eq. (B.12): u(d) — relative density (real gas)
static double u_G(const NumericVector& x, const NumericVector& u,
                  const NumericMatrix& r, int ks,
                  double Z, double G, double p2) {
  init_M_tables();
  double M       = M_val(x);
  double Z_air_k = Z_AIR[ks];
  double uZa     = U_Z_AIR[ks];
  double sg      = sigma_Z(Z, p2);

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = TAB_M[i][0] / M + 2.0 * TAB_S[i][ks] * sg / Z;
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = TAB_M[j][0] / M + 2.0 * TAB_S[j][ks] * sg / Z;
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (x[j] == 0.0 || r_M[i][j] == 0.0) continue;
      s2 += x[i]*u_M[i]*r_M[i][j]*x[j]*u_M[j];
    }
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s3 += x[i]*x[i] * TAB_S[i][4]*TAB_S[i][4];
  }
  double tmp = s1 + s2/(M*M) + 4.0*sg*sg*s3/(Z*Z) +
               (uZa/Z_air_k)*(uZa/Z_air_k);
  return std::sqrt(tmp) * G;
}

// Annex B, Eq. (B.13): u(W_gross) — Wobbe index gross
static double u_W_G(const NumericVector& x, const NumericVector& u,
                    const NumericMatrix& r, int kh, int ks,
                    double Z, double W_G, double p2) {
  init_M_tables();
  double Hg      = Hc_o_G(x, kh);
  double M       = M_val(x);
  double Z_air_k = Z_AIR[ks];
  double uZa     = U_Z_AIR[ks];
  double sg      = sigma_Z(Z, p2);

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = TAB_HC[i][kh]/Hg + TAB_S[i][ks]*sg/Z - TAB_M[i][0]/(2.0*M);
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = TAB_HC[j][kh]/Hg + TAB_S[j][ks]*sg/Z - TAB_M[j][0]/(2.0*M);
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s3 += x[i]*x[i] * TAB_S[i][4]*TAB_S[i][4];
  }
  double s4 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (x[j] == 0.0 || r_M[i][j] == 0.0) continue;
      s4 += x[i]*u_M[i]*r_M[i][j]*x[j]*u_M[j];
    }
  }
  double tmp = s1
    + s2 / (Hg*Hg)
    + sg*sg*s3 / (Z*Z)
    + s4 / (4.0*M*M)
    + (U_R_MOL/R_MOL)*(U_R_MOL/R_MOL)
    + (U_M_AIR/(2.0*M_AIR))*(U_M_AIR/(2.0*M_AIR))
    + (uZa/(2.0*Z_air_k))*(uZa/(2.0*Z_air_k));
  return std::sqrt(tmp) * W_G;
}

// Annex B, Eq. (B.14): u(W_net)
static double u_W_N(const NumericVector& x, const NumericVector& u,
                    const NumericMatrix& r, int kh, int ks,
                    double Z, double W_N, double p2) {
  init_M_tables();
  double Hn      = Hc_o_N(x, kh);
  double M       = M_val(x);
  double L       = L_VAP[kh];
  double uL      = U_L_VAP[kh];
  double Z_air_k = Z_AIR[ks];
  double uZa     = U_Z_AIR[ks];
  double sg      = sigma_Z(Z, p2);

  double s1 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (u[i] == 0.0) continue;
    double ci = (TAB_HC[i][kh] - L/2.0*TAB_M[i][2])/Hn + TAB_S[i][ks]*sg/Z - TAB_M[i][0]/(2.0*M);
    for (int j = 0; j < 60; j++) {
      if (u[j] == 0.0 || r(i,j) == 0.0) continue;
      double cj = (TAB_HC[j][kh] - L/2.0*TAB_M[j][2])/Hn + TAB_S[j][ks]*sg/Z - TAB_M[j][0]/(2.0*M);
      s1 += ci * u[i] * r(i,j) * cj * u[j];
    }
  }
  double s2 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s2 += x[i]*x[i] * TAB_HC[i][5]*TAB_HC[i][5];
  }
  double s3 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    s3 += x[i]*x[i] * TAB_S[i][4]*TAB_S[i][4];
  }
  double s4 = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    for (int j = 0; j < 60; j++) {
      if (x[j] == 0.0 || r_M[i][j] == 0.0) continue;
      s4 += x[i]*u_M[i]*r_M[i][j]*x[j]*u_M[j];
    }
  }
  double sum_nH = 0.0;
  for (int i = 0; i < 60; i++) {
    if (x[i] == 0.0) continue;
    sum_nH += x[i] * TAB_M[i][2];
  }
  double tmp = s1
    + s2 / (Hn*Hn)
    + sg*sg*s3 / (Z*Z)
    + s4 / (4.0*M*M)
    + (U_R_MOL/R_MOL)*(U_R_MOL/R_MOL)
    + (U_M_AIR/(2.0*M_AIR))*(U_M_AIR/(2.0*M_AIR))
    + (uZa/(2.0*Z_air_k))*(uZa/(2.0*Z_air_k))
    + (sum_nH/(2.0*Hn))*(sum_nH/(2.0*Hn))*uL*uL;
  return std::sqrt(tmp) * W_N;
}

// ============================================================================
// Single exported Rcpp function
// ============================================================================

//' Calculate all ISO 6976:2016 properties and their uncertainties
//'
//' @param x       Numeric vector of length 60: mole fractions [mol/mol]
//' @param u_x     Numeric vector of length 60: standard uncertainties
//' @param r_x     Numeric matrix 60x60: correlation matrix
//' @param t1      Combustion temperature [°C]: 0, 15, 15.55, 20, or 25
//' @param t2      Volume reference temperature [°C]: 0, 15, 15.55, or 20
//' @param p2      Reference pressure [kPa]: 90–110
//' @param k       Coverage factor (default 1)
//' @return Named list of properties and uncertainties
// [[Rcpp::export]]
Rcpp::List iso6976_calc(
    Rcpp::NumericVector x,
    Rcpp::NumericVector u_x,
    Rcpp::NumericMatrix r_x,
    double t1, double t2, double p2,
    double k = 1.0
) {
  // Validate
  if (x.size() != 60 || u_x.size() != 60)
    Rcpp::stop("x and u_x must have length 60");
  if (r_x.nrow() != 60 || r_x.ncol() != 60)
    Rcpp::stop("r_x must be a 60x60 matrix");
  if (p2 < 90.0 || p2 > 110.0)
    Rcpp::stop("p2 out of application range (90–110 kPa)");

  int kh = idx_hc(t1);
  if (kh < 0) Rcpp::stop("t1 must be 0, 15, 15.55, 20, or 25 °C");

  int ks = idx_s(t2);
  if (ks < 0) Rcpp::stop("t2 must be 0, 15, 15.55, or 20 °C");

  init_M_tables();

  // ---- values ----
  double M   = M_val(x);
  double Z   = Z_val(x, p2, ks);
  if (Z <= 0.9)
    Rcpp::stop("Computed Z <= 0.9: outside application range");

  double Hcg = Hc_o_G(x, kh);
  double Hcn = Hc_o_N(x, kh);
  double Hmg = Hcg / M;
  double Hmn = Hcn / M;

  double Vo  = V_o(t2, p2);           // ideal molar volume [m^3/mol]
  double V   = Z * Vo;                 // real molar volume [m^3/mol]  Eq.(11)

  double Hvg_o = Hcg / Vo;            // ideal-gas vol. GCV  Eq.(7)
  double Hvn_o = Hcn / Vo;            // ideal-gas vol. NCV  Eq.(9)
  double Hvg   = Hcg / V;             // real-gas  vol. GCV  Eq.(10)
  double Hvn   = Hcn / V;             // real-gas  vol. NCV  Eq.(12)

  double G_o = M / M_AIR;             // ideal rel. density  Eq.(13)
  double D_o = M / Vo;                // ideal density       Eq.(14)

  // Eq.(17): real relative density (Eq.18 for Z_air)
  double Z_air = Z_AIR[ks];
  double Z_air_real = 1.0 - p2/P_REF * (1.0 - Z_air);  // Eq.(18)
  double G = G_o * Z_air_real / Z;                       // Eq.(17)
  double D = D_o / Z;                                    // Eq.(19)

  double Wg_o = Hvg_o / std::sqrt(G_o);  // Eq.(15)
  double Wn_o = Hvn_o / std::sqrt(G_o);  // Eq.(16)
  double Wg   = Hvg   / std::sqrt(G);    // Eq.(20)
  double Wn   = Hvn   / std::sqrt(G);    // Eq.(21)

  // ---- uncertainties ----
  double uHcg   = u_Hc_o_G(x, u_x, r_x, kh);
  double uHcn   = u_Hc_o_N(x, u_x, r_x, kh);
  double uHmg   = u_Hm_o_G(x, u_x, r_x, kh);
  double uHmn   = u_Hm_o_N(x, u_x, r_x, kh);
  double uHvg_o = u_Hv_o_G(x, u_x, r_x, kh, ks, Z, p2);
  double uHvn_o = u_Hv_o_N(x, u_x, r_x, kh, ks, Z, p2);
  // For real-gas volumetric, the standard only gives ideal-gas formulas in
  // Annex B; u(Hv_G) = u(Hv_o_G) to first order (Z cancels in relative unc.)
  double uHvg   = uHvg_o;
  double uHvn   = uHvn_o;

  double uD  = u_D(x, u_x, r_x, ks, Z, D, p2);
  double uG  = u_G(x, u_x, r_x, ks, Z, G, p2);
  double uWg = u_W_G(x, u_x, r_x, kh, ks, Z, Wg, p2);
  double uWn = u_W_N(x, u_x, r_x, kh, ks, Z, Wn, p2);

  // ---- assemble result ----
  return Rcpp::List::create(
    Rcpp::Named("M")       = M,
    Rcpp::Named("Z")       = Z,
    Rcpp::Named("G_o")     = G_o,
    Rcpp::Named("D_o")     = D_o,
    Rcpp::Named("G")       = G,      Rcpp::Named("u_G")   = k * uG,
    Rcpp::Named("D")       = D,      Rcpp::Named("u_D")   = k * uD,
    Rcpp::Named("Hcg")     = Hcg,    Rcpp::Named("u_Hcg") = k * uHcg,
    Rcpp::Named("Hcn")     = Hcn,    Rcpp::Named("u_Hcn") = k * uHcn,
    Rcpp::Named("Hmg")     = Hmg,    Rcpp::Named("u_Hmg") = k * uHmg,
    Rcpp::Named("Hmn")     = Hmn,    Rcpp::Named("u_Hmn") = k * uHmn,
    Rcpp::Named("Hvg_o")   = Hvg_o,  Rcpp::Named("u_Hvg_o") = k * uHvg_o,
    Rcpp::Named("Hvn_o")   = Hvn_o,  Rcpp::Named("u_Hvn_o") = k * uHvn_o,
    Rcpp::Named("Hvg")     = Hvg,    Rcpp::Named("u_Hvg") = k * uHvg,
    Rcpp::Named("Hvn")     = Hvn,    Rcpp::Named("u_Hvn") = k * uHvn,
    Rcpp::Named("Wg_o")    = Wg_o,
    Rcpp::Named("Wn_o")    = Wn_o,
    Rcpp::Named("Wg")      = Wg,     Rcpp::Named("u_Wg")  = k * uWg,
    Rcpp::Named("Wn")      = Wn,     Rcpp::Named("u_Wn")  = k * uWn
  );
}
