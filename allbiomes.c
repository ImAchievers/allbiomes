#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <pthread.h>
#include <inttypes.h>
#include <stdatomic.h>
#include <math.h>
#include <float.h>

#include "generator.h"
#include "biomenoise.h"

#define MCV             MC_NEWEST
#define START_SEED      1000000000000000000ULL
#define NUM_THREADS     18
#define MAX_AREA        4203264              /* max continent area in blocks^2 */
#define MAX_MISS        8                    /* Skip if more biomes than this missing at 1:16 */
#define RING_R          ((int)(0.6 * sqrt((double)MAX_AREA)))  /* escape-ring radius, auto-scaled to MAX_AREA */
#define RING_MIN        8                    /* min ocean hits out of 16 on the ring */
#define C2_STEP         48                   /* stage-C2 climate scan step in blocks */
#define COARSE          1                    /* coarse pre-flood size reject before stage B */
#define COARSE_MARGIN   1.5
#define VERIFY_ISOLATION 1                   /* 1:4 re-flood to drop land-bridged hits */
#define ISO_MARGIN      1.2
#define SEARCH_RADIUS   2000                 /* candidate-core search radius in blocks */
#define A0_STEP         1000                 /* A0 grid step in blocks */
#define LOW_OCT         4                    /* low-frequency octaves kept in the low-pass */
#define A_MARGIN        0.06                 /* stage-A low-pass threshold slack */
#define FILL_FACTOR     1.15                 /* flood budget factor over MAX_AREA */
#define MIN_REPORT      36                   /* min biomes present to log a near-miss */

#define LAND_C          (-0.19)
#define TEMP_HI         0.55
#define TEMP_LO         (-0.45)
#define EROS_HI         0.5501
#define EROS_LO         (-0.375)

#define CHUNK           1024
#define STATUS_INTERVAL 20
#define CKPT_INTERVAL   30
#define RESULTS_FILE    "allbiomes_results.txt"
#define CHECKPOINT_FILE "allbiomes_checkpoint.txt"
#define MEASURE_MAX     15000000             /* single <seed> <x> <z>: max continent area to measure */

/* exact climate min/max over region */
static double getParaDescent(const DoublePerlinNoise *para, double factor,
    int x, int z, int w, int h, int i0, int j0, int maxrad,
    int maxiter, double alpha, void *data, int (*func)(void*,int,int,double))
{
    int dirx = 0, dirz = 0, dira;
    int k, i, j;
    double v, vd, va;
    v = factor * sampleDoublePerlin(para, x+i0, 0, z+j0);
    if (func) {
        if (func(data, x+i0, z+j0, factor < 0 ? -v : v))
            return nan("");
    }
    i = i0; j = j0;
    for (k = 0; k < maxiter; k++) {
        if (dirx == 0) dirx = +1;
        if (i+dirx >= 0 && i+dirx < w)
            vd = factor * sampleDoublePerlin(para, x+i+dirx, 0, z+j);
        else vd = v;
        if (vd >= v) {
            dirx *= -1;
            if (i+dirx >= 0 && i+dirx < w)
                vd = factor * sampleDoublePerlin(para, x+i+dirx, 0, z+j);
            else vd = v;
            if (vd >= v)
                dirx = 0;
        }
        if (dirx) {
            dira = (int)(dirx * alpha * (v - vd));
            if (abs(dira) > 2 && i+dira >= 0 && i+dira < w) {
                va = factor * sampleDoublePerlin(para, x+i+dira, 0, z+j);
                if (va < vd) {
                    i += dira;
                    v = va;
                    goto L_x_end;
                }
            }
            v = vd;
            i += dirx;
        L_x_end:
            if (func) {
                if (func(data, x+i, z+j, factor < 0 ? -v : v))
                    return nan("");
            }
        }
        if (dirz == 0) dirz = +1;
        if (j+dirz >= 0 && j+dirz < h)
            vd = factor * sampleDoublePerlin(para, x+i, 0, z+j+dirz);
        else vd = v;
        if (vd >= v) {
            dirz *= -1;
            if (j+dirz >= 0 && j+dirz < h)
                vd = factor * sampleDoublePerlin(para, x+i, 0, z+j+dirz);
            else vd = v;
            if (vd >= v)
                dirz = 0;
        }
        if (dirz) {
            dira = (int)(dirz * alpha * (v - vd));
            if (abs(dira) > 2 && j+dira >= 0 && j+dira < h) {
                va = factor * sampleDoublePerlin(para, x+i, 0, z+j+dira);
                if (va < vd) {
                    j += dira;
                    v = va;
                    goto L_z_end;
                }
            }
            j += dirz;
            v = vd;
        L_z_end:
            if (func) {
                if (func(data, x+i, z+j, factor < 0 ? -v : v))
                    return nan("");
            }
        }
        if (dirx == 0 && dirz == 0) {
            int c;
            for (c = 0; c < 4; c++) {
                dirx = (c & 1) ? -1 : +1;
                dirz = (c & 2) ? -1 : +1;
                if (i+dirx < 0 || i+dirx >= w || j+dirz < 0 || j+dirz >= h)
                    continue;
                vd = factor * sampleDoublePerlin(para, x+i+dirx, 0, z+j+dirz);
                if (vd < v) {
                    v = vd;
                    i += dirx;
                    j += dirz;
                    break;
                }
            }
            if (c >= 4)
                break;
        }
        if (abs(i - i0) > maxrad || abs(j - j0) > maxrad)
            break;
    }
    return v;
}

static int getParaRange(const DoublePerlinNoise *para, double *pmin, double *pmax,
    int x, int z, int w, int h, void *data, int (*func)(void*,int,int,double))
{
    const double beta = 1.5;
    const double factor = 10000;
    const double perlin_grad = 2.0 * 1.875;
    double v, lmin, lmax, dr, vdif, small_regime;
    char *skip = NULL;
    int i, j, step, ii, jj, ww, hh, skipsiz;
    int maxrad, maxiter;
    int err = 1;

    if (pmin) *pmin = DBL_MAX;
    if (pmax) *pmax = -DBL_MAX;

    lmin = DBL_MAX, lmax = 0;
    for (i = 0; i < para->octA.octcnt; i++) {
        double lac = para->octA.octaves[i].lacunarity;
        if (lac < lmin) lmin = lac;
        if (lac > lmax) lmax = lac;
    }

    small_regime = 1e3 * sqrt(lmax);
    if (w*h < small_regime) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                v = factor * sampleDoublePerlin(para, x+i, 0, z+j);
                if (func) {
                    err = func(data, x+i, z+j, v);
                    if (err)
                        return err;
                }
                if (pmin && v < *pmin) *pmin = v;
                if (pmax && v > *pmax) *pmax = v;
            }
        }
        return 0;
    }

    step = (int) (0.5 / lmin - FLT_EPSILON) + 1;

    dr = lmax / lmin * beta;
    for (j = 0; j < h; j += step) {
        for (i = 0; i < w; i += step) {
            if (pmin) {
                v = getParaDescent(para, +factor, x, z, w, h, i, j,
                    step, step, dr, data, func);
                if (v != v) goto L_end;
                if (v < *pmin) *pmin = v;
            }
            if (pmax) {
                v = -getParaDescent(para, -factor, x, z, w, h, i, j,
                    step, step, dr, data, func);
                if (v != v) goto L_end;
                if (v > *pmax) *pmax = v;
            }
        }
    }

    step = (int) (1.0 / (perlin_grad * lmax + FLT_EPSILON)) + 1;

    vdif = 0;
    for (i = 0; i < para->octA.octcnt; i++) {
        const PerlinNoise *p = para->octA.octaves + i;
        double contrib = step * p->lacunarity * 1.0;
        if (contrib > 1.0) contrib = 1;
        vdif += contrib * p->amplitude;
    }
    for (i = 0; i < para->octB.octcnt; i++) {
        const double lac_factB = 337.0 / 331.0;
        const PerlinNoise *p = para->octB.octaves + i;
        double contrib = step * p->lacunarity * lac_factB;
        if (contrib > 1.0) contrib = 1;
        vdif += contrib * p->amplitude;
    }
    vdif = fabs(factor * vdif * para->amplitude);

    maxrad = step;
    maxiter = step*2;
    ww = (w+step-1) / step;
    hh = (h+step-1) / step;
    skipsiz = (ww+1) * (hh+1) * sizeof(*skip);
    skip = (char*) malloc(skipsiz);

    if (pmin) {
        memset(skip, 0, skipsiz);
        for (jj = 0; jj <= hh; jj++) {
            j = jj * step; if (j >= h) j = h-1;
            for (ii = 0; ii <= ww; ii++) {
                i = ii * step; if (i >= w) i = w-1;
                if (skip[jj*ww+ii]) continue;
                v = factor * sampleDoublePerlin(para, x+i, 0, z+j);
                if (func) {
                    int e = func(data, x+i, z+j, v);
                    if (e) { err = e; goto L_end; }
                }
                if (pmax && v > *pmax) *pmax = v;
                dr = beta * (v - *pmin) / vdif;
                if (dr > 1.0) {
                    int a, b, r = (int) dr;
                    for (b = 0; b < r; b++) {
                        if (b+jj < 0 || b+jj >= hh) continue;
                        for (a = -r+1; a < r; a++) {
                            if (a+ii < 0 || a+ii >= ww) continue;
                            skip[(b+jj)*ww + (a+ii)] = 1;
                        }
                    }
                    continue;
                }
                v = getParaDescent(para, +factor, x, z, w, h, i, j,
                    maxrad, maxiter, dr, data, func);
                if (v != v) goto L_end;
                if (v < *pmin) *pmin = v;
            }
        }
    }

    if (pmax) {
        memset(skip, 0, skipsiz);
        for (jj = 0; jj <= hh; jj++) {
            j = jj * step; if (j >= h) j = h-1;
            for (ii = 0; ii <= ww; ii++) {
                i = ii * step; if (i >= w) i = w-1;
                if (skip[jj*ww+ii]) continue;
                v = -factor * sampleDoublePerlin(para, x+i, 0, z+j);
                if (func) {
                    int e = func(data, x+i, z+j, -v);
                    if (e) { err = e; goto L_end; }
                }
                dr = beta * (v + *pmax) / vdif;
                if (dr > 1.0) {
                    int a, b, r = (int) dr;
                    for (b = 0; b < r; b++) {
                        if (b+jj < 0 || b+jj >= hh) continue;
                        for (a = -r+1; a < r; a++) {
                            if (a+ii < 0 || a+ii >= ww) continue;
                            skip[(b+jj)*ww + (a+ii)] = 1;
                        }
                    }
                    continue;
                }
                v = -getParaDescent(para, -factor, x, z, w, h, i, j,
                    maxrad, maxiter, dr, data, func);
                if (v != v) goto L_end;
                if (v > *pmax) *pmax = v;
            }
        }
    }

    err = 0;
L_end:
    if (skip)
        free(skip);
    return err;
}

/* Low-pass continentalness */
static void lowpassContInit(BiomeNoise * bn, uint64_t seed, int K) {
    static const double ampfull[9] = { 1,1,2,2,2,1,1,1,1 };
    double amp[9] = { 0 };
    for (int i = 0; i < K && i < 9; i++)
        amp[i] = ampfull[i];
    Xoroshiro pxr;
    xSetSeed(&pxr, seed);
    uint64_t xlo = xNextLong(&pxr), xhi = xNextLong(&pxr);
    pxr.lo = xlo ^ 0x83886c9d0ae3a662;
    pxr.hi = xhi ^ 0xafa638a61b42e8ad;
    xDoublePerlinInit(&bn->climate[NP_CONTINENTALNESS], &pxr, bn->oct, amp, -9, 9, -1);
    bn->climate[NP_CONTINENTALNESS].amplitude = 1.5;
}

/* required biomes */
static const int REQ[] = {
    plains, desert, windswept_hills, forest, taiga, swamp,
    river, frozen_river, snowy_plains, beach, jungle,
    sparse_jungle, stony_shore, snowy_beach,
    birch_forest, dark_forest, snowy_taiga, old_growth_pine_taiga,
    windswept_forest, savanna, savanna_plateau, badlands,
    wooded_badlands, sunflower_plains,
    windswept_gravelly_hills, flower_forest, ice_spikes,
    old_growth_birch_forest, old_growth_spruce_taiga,
    windswept_savanna, eroded_badlands, bamboo_jungle, meadow,
    grove, snowy_slopes, jagged_peaks, frozen_peaks, stony_peaks,
    mangrove_swamp, cherry_grove, pale_garden, dappled_forest,
};
#define NREQ ((int)(sizeof(REQ) / sizeof(REQ[0])))
static const char * REQNAME[] = {
    "plains","desert","windswept_hills","forest","taiga","swamp","river",
    "frozen_river","snowy_plains","beach","jungle","sparse_jungle","stony_shore",
    "snowy_beach","birch_forest","dark_forest","snowy_taiga","old_growth_pine_taiga",
    "windswept_forest","savanna","savanna_plateau","badlands","wooded_badlands",
    "sunflower_plains","windswept_gravelly_hills","flower_forest","ice_spikes",
    "old_growth_birch_forest","old_growth_spruce_taiga","windswept_savanna",
    "eroded_badlands","bamboo_jungle","meadow","grove","snowy_slopes",
    "jagged_peaks","frozen_peaks","stony_peaks","mangrove_swamp","cherry_grove",
    "pale_garden","dappled_forest"
};

/* climate boxes */
#define NF (-99.0f)
#define PF (99.0f)
typedef struct {
    int id;
    float t0, t1, h0, h1, c0, c1, e0, e1;
    float wa0, wa1, wb0, wb1;
} BBox;
static const BBox BOXES[] = {
    { 131, NF,-.1440f, NF,-.0929f, -.1903f,PF, .4415f,.5411f, .2603f,PF, NF,-.2603f },
    { 182, .1963f,.5433f, NF,PF, -.1903f,PF, NF,-.3728f, .4f,PF, NF,-.4f },
    { 186, -.1508f,.2031f, .3f,PF, .0301f,PF, NF,.0524f, .2603f,PF, NF,-.2603f },
    { 165, .5501f,PF, NF,.1071f, -.1903f,PF, NF,.0524f, .0444f,PF, NF,-.0444f },
    { 140, NF,-.4502f, NF,-.3501f, -.1977f,PF, NF,PF, .0444f,PF, NF,-.0444f },
    { 185, -.4434f,.2031f, NF,-.0929f, .0228f,PF, NF,.0524f, .2603f,PF, NF,-.2603f },
    { 38,  .5501f,PF, .1071f,PF, -.1903f,PF, NF,.0524f, NF,PF, 0,0 },
    { 180, NF,.2031f, NF,PF, -.1903f,PF, NF,-.3728f, .4f,PF, NF,-.4f },
    { 181, NF,.2031f, NF,PF, -.1903f,PF, NF,-.3728f, .4f,PF, NF,-.4f },
    { 184, .1963f,PF, NF,PF, -.1169f,PF, .5501f,PF, NF,PF, 0,0 },
    { 132, -.1576f,.2031f, NF,-.3501f, -.1903f,PF, NF,PF, NF,-.0444f, .0444f,PF },
    { 129, -.1576f,.2031f, NF,-.3501f, -.1903f,PF, NF,PF, .0444f,PF, NF,-.0444f },
};
#define NBOX ((int)(sizeof(BOXES) / sizeof(BOXES[0])))

static volatile sig_atomic_t g_stop = 0;
static void on_sigint(int sig) { (void)sig; g_stop = 1; }
static atomic_ullong g_next, g_seeds, g_hits;
static _Atomic int g_bestbio;
static _Atomic int64_t g_bestarea;
static _Atomic uint64_t g_bestseed;
static pthread_mutex_t g_outmtx = PTHREAD_MUTEX_INITIALIZER;
static FILE * g_out;
static int64_t g_areacap = MAX_AREA;

#define FCELL 96        /* stage-B fill cell size in blocks */
#define FHALF (FCELL/2)
#define FR 72           /* fill radius in FCELL cells */
#define FW (2 * FR + 1)
#define DMAX 700        /* max 1:16 grid width for stage D */
#define CG (2 * FCELL)  /* coarse pre-flood cell size in blocks */
#define CHALF (CG / 2)
#define FRc (FR / 2)    /* coarse pre-flood radius in CG cells */
#define FWc (2 * FRc + 1)
#define ISO_W 1536      /* isolation re-flood window, 1:4 cells */
#define ISO_R (ISO_W / 2)

typedef struct {
    BiomeNoise bnL; /* continentalness, low-pass */
    BiomeNoise bnC; /* continentalness, full */
    BiomeNoise bnT; /* temperature, full */
    BiomeNoise bnE; /* erosion, full */
    int haveC, haveT, haveE;
    Generator g;
    int haveG;
    uint64_t seed;
    int8_t * lc;
    uint8_t * vis;
    int32_t * q;
    int * ids;
    int32_t * q16;
    uint8_t * comp;
    size_t idsz;
    int gx0, gz0, gw, gh;
    uint8_t * iso_vis;
    int32_t * iso_q;
} Ctx;

static inline double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static inline int fdiv(int a, int b) {
    return a >= 0 ? a / b : -((-a + b - 1) / b);
}
static inline double sampC(const BiomeNoise * bn, int np, double bx, double bz) {
    return sampleDoublePerlin(bn->climate + np, bx * 0.25, 0, bz * 0.25);
}

static void ctx_seed(Ctx * cx, uint64_t seed) {
    cx->seed = seed;
    cx->haveC = cx->haveT = cx->haveE = cx->haveG = 0;
    lowpassContInit(&cx->bnL, seed, LOW_OCT);
}
static BiomeNoise * needC(Ctx * cx) {
    if (!cx->haveC) { setClimateParaSeed(&cx->bnC, cx->seed, 0, NP_CONTINENTALNESS, -1); cx->haveC = 1; }
    return &cx->bnC;
}
static BiomeNoise * needT(Ctx * cx) {
    if (!cx->haveT) { setClimateParaSeed(&cx->bnT, cx->seed, 0, NP_TEMPERATURE, -1); cx->haveT = 1; }
    return &cx->bnT;
}
static BiomeNoise * needE(Ctx * cx) {
    if (!cx->haveE) { setClimateParaSeed(&cx->bnE, cx->seed, 0, NP_EROSION, -1); cx->haveE = 1; }
    return &cx->bnE;
}
static Generator * needG(Ctx * cx) {
    if (!cx->haveG) { applySeed(&cx->g, DIM_OVERWORLD, cx->seed); cx->haveG = 1; }
    return &cx->g;
}

/* budgeted flood fill */
typedef struct {
    int cells;
    int bx1, bz1, bx2, bz2;
    int cxblk, czblk;
} Fill;

static inline int landcell(Ctx * cx, BiomeNoise * bn, int cxc, int czc) {
    int idx = (czc + FR) * FW + (cxc + FR);
    int8_t v = cx->lc[idx];
    if (v) return v == 1;
    double c = sampC(bn, NP_CONTINENTALNESS, cxc * FCELL + FHALF, czc * FCELL + FHALF);
    cx->lc[idx] = (c >= LAND_C) ? 1 : 2;
    return c >= LAND_C;
}

static int stageB(Ctx * cx, int sx, int sz, int budget, Fill * out) {
    BiomeNoise * bn = &cx->bnL;
    memset(cx->lc, 0, FW * FW);
    memset(cx->vis, 0, FW * FW);
    int c0x = fdiv(sx, FCELL), c0z = fdiv(sz, FCELL);
    if (abs(c0x) >= FR - 1 || abs(c0z) >= FR - 1) return 0;
    int fx = 0, fz = 0, found = 0;
    for (int r = 0; r <= 3 && !found; r++)
        for (int dz = -r; dz <= r && !found; dz++)
            for (int dx = -r; dx <= r && !found; dx++)
                if (landcell(cx, bn, c0x + dx, c0z + dz)) { fx = c0x + dx; fz = c0z + dz; found = 1; }
    if (!found) return 0;
    int qh = 0, qt = 0;
    int32_t * q = cx->q;
    q[qt++] = (fz + FR) * FW + (fx + FR);
    cx->vis[q[0]] = 1;
    int n = 0;
    long sx_sum = 0, sz_sum = 0;
    int x1 = fx, x2 = fx, z1 = fz, z2 = fz;
    while (qh < qt) {
        int p = q[qh++];
        int px = p % FW - FR, pz = p / FW - FR;
        if (++n > budget) return 0;
        if (abs(px) >= FR - 1 || abs(pz) >= FR - 1) return 0;
        if (px < x1) x1 = px;
        if (px > x2) x2 = px;
        if (pz < z1) z1 = pz;
        if (pz > z2) z2 = pz;
        sx_sum += px; sz_sum += pz;
        static const int DX[4] = { 1,-1,0,0 }, DZ[4] = { 0,0,1,-1 };
        for (int k = 0; k < 4; k++) {
            int nx = px + DX[k], nz = pz + DZ[k];
            int ni = (nz + FR) * FW + (nx + FR);
            if (cx->vis[ni]) continue;
            cx->vis[ni] = 1;
            if (landcell(cx, bn, nx, nz)) q[qt++] = ni;
        }
    }
    out->cells = n;
    out->bx1 = x1 * FCELL; out->bz1 = z1 * FCELL;
    out->bx2 = x2 * FCELL + (FCELL-1); out->bz2 = z2 * FCELL + (FCELL-1);
    out->cxblk = (int)(sx_sum * FCELL / n) + FHALF;
    out->czblk = (int)(sz_sum * FCELL / n) + FHALF;
    return 1;
}

/* coarse pre-flood */
static inline int landcellC(Ctx * cx, int cxc, int czc) {
    int idx = (czc + FRc) * FWc + (cxc + FRc);
    int8_t v = cx->lc[idx];
    if (v) return v == 1;
    double c = sampC(&cx->bnL, NP_CONTINENTALNESS, cxc * CG + CHALF, czc * CG + CHALF);
    cx->lc[idx] = (c >= LAND_C) ? 1 : 2;
    return c >= LAND_C;
}
static int coarseTooBig(Ctx * cx, int sx, int sz, int budget) {
    memset(cx->lc, 0, FWc * FWc);
    memset(cx->vis, 0, FWc * FWc);
    int c0x = fdiv(sx, CG), c0z = fdiv(sz, CG);
    if (abs(c0x) >= FRc - 1 || abs(c0z) >= FRc - 1) return 0;
    int fx = 0, fz = 0, found = 0;
    for (int r = 0; r <= 3 && !found; r++)
        for (int dz = -r; dz <= r && !found; dz++)
            for (int dx = -r; dx <= r && !found; dx++)
                if (landcellC(cx, c0x + dx, c0z + dz)) { fx = c0x + dx; fz = c0z + dz; found = 1; }
    if (!found) return 0;
    int32_t * q = cx->q;
    int qh = 0, qt = 0;
    q[qt++] = (fz + FRc) * FWc + (fx + FRc);
    cx->vis[q[0]] = 1;
    int n = 0;
    while (qh < qt) {
        int p = q[qh++];
        int px = p % FWc - FRc, pz = p / FWc - FRc;
        if (++n > budget) return 1;
        if (abs(px) >= FRc - 1 || abs(pz) >= FRc - 1) return 0;
        static const int DX[4] = { 1,-1,0,0 }, DZ[4] = { 0,0,1,-1 };
        for (int k = 0; k < 4; k++) {
            int nx = px + DX[k], nz = pz + DZ[k];
            int ni = (nz + FRc) * FWc + (nx + FRc);
            if (cx->vis[ni]) continue;
            cx->vis[ni] = 1;
            if (landcellC(cx, nx, nz)) q[qt++] = ni;
        }
    }
    return 0;
}

/* exact temperature/erosion extremes over the bbox */
static int stageC(Ctx * cx, const Fill * f) {
    int x4 = f->bx1 >> 2, z4 = f->bz1 >> 2;
    int w4 = ((f->bx2 - f->bx1) >> 2) + 1, h4 = ((f->bz2 - f->bz1) >> 2) + 1;
    double lo, hi;
    BiomeNoise * bt = needT(cx);
    if (getParaRange(&bt->climate[NP_TEMPERATURE], &lo, &hi, x4, z4, w4, h4, NULL, NULL)) return 0;
    if (hi < TEMP_HI || lo > TEMP_LO) return 0;
    BiomeNoise * be = needE(cx);
    if (getParaRange(&be->climate[NP_EROSION], &lo, &hi, x4, z4, w4, h4, NULL, NULL)) return 0;
    if (hi < EROS_HI || lo > EROS_LO) return 0;
    return 1;
}

static inline int inbox(const BBox * b, float t, float h, float c, float e, float w, float d) {
    if (t < b->t0 - d || t > b->t1 + d) return 0;
    if (h < b->h0 - d || h > b->h1 + d) return 0;
    if (c < b->c0 - d || c > b->c1 + d) return 0;
    if (e < b->e0 - d || e > b->e1 + d) return 0;
    if (w >= b->wa0 - d && w <= b->wa1 + d) return 1;
    if (b->wb0 != b->wb1 && w >= b->wb0 - d && w <= b->wb1 + d) return 1;
    return 0;
}

/* rare biome climate pre-gate */
static int stageC2(Ctx * cx, const Fill * f, int verbose) {
    Generator * g = needG(cx);
    const BiomeNoise * bn = &g->bn;
    uint8_t seen[NBOX];
    memset(seen, 0, sizeof seen);
    int need = NBOX;
    const int step = C2_STEP;
    for (int bz = f->bz1; bz <= f->bz2; bz += step) {
        for (int bx = f->bx1; bx <= f->bx2; bx += step) {
            if (sampC(&cx->bnL, NP_CONTINENTALNESS, bx, bz) < LAND_C - 0.04) continue;
            float t = sampC(bn, NP_TEMPERATURE, bx, bz);
            float h = sampC(bn, NP_HUMIDITY, bx, bz);
            float c = sampC(bn, NP_CONTINENTALNESS, bx, bz);
            float e = sampC(bn, NP_EROSION, bx, bz);
            float w = sampC(bn, NP_WEIRDNESS, bx, bz);
            for (int b = 0; b < NBOX; b++) {
                if (seen[b]) continue;
                if (inbox(&BOXES[b], t, h, c, e, w, 0.08f)) { seen[b] = 1; need--; }
            }
            if (!need) break;
        }
        if (!need) break;
    }
    if (need) {
        if (verbose) {
            printf("  C2 gate: reject, climate-impossible biomes:");
            for (int b = 0; b < NBOX; b++)
                if (!seen[b]) printf(" %d", BOXES[b].id);
            printf("\n");
        }
        return 0;
    }
    return 1;
}

typedef struct {
    int64_t area;
    int nfound;
    int cxblk, czblk;
    int reached_bio16;
    char missing[512];
    int miss[NREQ];
    int nmiss;
} DRes;

/* 1:16 biome scan + exact continent area */
static int stageD1(Ctx * cx, const Fill * f, DRes * dr, int verbose) {
    Generator * g = needG(cx);
    int gx0 = fdiv(f->bx1, 16) - 4, gz0 = fdiv(f->bz1, 16) - 4;
    int gw = (f->bx2 - f->bx1) / 16 + 9, gh = (f->bz2 - f->bz1) / 16 + 9;
    if (gw > DMAX || gh > DMAX) return 0;
    Range r = { 16, gx0, gz0, gw, gh, 64, 1 };
    if (genBiomes(g, cx->ids, r)) return 0;
    cx->gx0 = gx0; cx->gz0 = gz0; cx->gw = gw; cx->gh = gh;

    memset(cx->comp, 0, (size_t)gw * gh);
    int sxc = fdiv(f->cxblk, 16) - gx0, szc = fdiv(f->czblk, 16) - gz0;
    int found = 0, fx = 0, fz = 0;
    for (int rr = 0; rr <= 8 && !found; rr++)
        for (int dz = -rr; dz <= rr && !found; dz++)
            for (int dx = -rr; dx <= rr && !found; dx++) {
                int x = sxc + dx, z = szc + dz;
                if (x < 0 || z < 0 || x >= gw || z >= gh) continue;
                if (!isOceanic(cx->ids[z * gw + x])) { fx = x; fz = z; found = 1; }
            }
    if (!found) return 0;

    int32_t * q = cx->q16;
    int qh = 0, qt = 0;
    q[qt++] = fz * gw + fx;
    cx->comp[q[0]] = 1;
    long n = 0, sxs = 0, szs = 0;
    uint8_t pres[256];
    memset(pres, 0, sizeof pres);
    while (qh < qt) {
        int p = q[qh++];
        int px = p % gw, pz = p / gw;
        n++; sxs += px; szs += pz;
        int id = cx->ids[p];
        if (id >= 0 && id < 256) pres[id] = 1;
        static const int DX[4] = { 1,-1,0,0 }, DZ[4] = { 0,0,1,-1 };
        for (int k = 0; k < 4; k++) {
            int nx = px + DX[k], nz = pz + DZ[k];
            if (nx < 0 || nz < 0 || nx >= gw || nz >= gh) continue;
            int ni = nz * gw + nx;
            if (cx->comp[ni]) continue;
            cx->comp[ni] = 1;
            if (!isOceanic(cx->ids[ni])) q[qt++] = ni;
        }
        if (n * 256 > g_areacap) return 0;
    }
    dr->area = n * 256;
    dr->cxblk = (int)((gx0 + sxs / n) * 16) + 8;
    dr->czblk = (int)((gz0 + szs / n) * 16) + 8;

    dr->nmiss = 0;
    for (int i = 0; i < NREQ; i++)
        if (!pres[REQ[i]]) dr->miss[dr->nmiss++] = REQ[i];
    if (verbose) {
        printf("  D1: area=%lld (%ld cells) center=(%d,%d) missing@16=%d:",
            (long long)dr->area, n, dr->cxblk, dr->czblk, dr->nmiss);
        for (int i = 0; i < dr->nmiss; i++) printf(" %d", dr->miss[i]);
        printf("\n");
    }
    if (dr->nmiss > MAX_MISS) { dr->nfound = NREQ - dr->nmiss; return 1; }
    dr->reached_bio16 = 1;
    return 1;
}

/* 1:4 fallback for biomes missing at 1:16 */
static void stageD2(Ctx * cx, DRes * dr, int verbose) {
    Generator * g = needG(cx);
    int gx0 = cx->gx0, gz0 = cx->gz0, gw = cx->gw, gh = cx->gh;
    int *miss = dr->miss; int nmiss = dr->nmiss;
    uint8_t pres[256]; memset(pres, 0, sizeof pres);

    if (nmiss > 0) {
        const BBox * mb[NREQ];
        int anyboxless = 0;
        for (int i = 0; i < nmiss; i++) {
            mb[i] = NULL;
            for (int b = 0; b < NBOX; b++)
                if (BOXES[b].id == miss[i]) { mb[i] = &BOXES[b]; break; }
            if (!mb[i]) anyboxless = 1;
        }
        const BiomeNoise * bn = &g->bn;
        for (int p = 0; p < gw * gh && nmiss > 0; p++) {
            if (!cx->comp[p] || isOceanic(cx->ids[p])) continue;
            int px = p % gw + gx0, pz = p / gw + gz0;
            double x4c = px * 4 + 2, z4c = pz * 4 + 2;
            float t = sampleDoublePerlin(bn->climate + NP_TEMPERATURE, x4c, 0, z4c);
            float h = sampleDoublePerlin(bn->climate + NP_HUMIDITY, x4c, 0, z4c);
            float c = sampleDoublePerlin(bn->climate + NP_CONTINENTALNESS, x4c, 0, z4c);
            float e = sampleDoublePerlin(bn->climate + NP_EROSION, x4c, 0, z4c);
            float w = sampleDoublePerlin(bn->climate + NP_WEIRDNESS, x4c, 0, z4c);
            int gate = anyboxless;
            for (int i = 0; i < nmiss && !gate; i++)
                if (mb[i] && inbox(mb[i], t, h, c, e, w, 0.08f)) gate = 1;
            if (!gate) continue;
            for (int sz4 = 0; sz4 < 4 && nmiss > 0; sz4++)
                for (int sx4 = 0; sx4 < 4 && nmiss > 0; sx4++) {
                    int id = getBiomeAt(g, 4, px * 4 + sx4, 64, pz * 4 + sz4);
                    if (id < 0 || id >= 256 || pres[id]) continue;
                    pres[id] = 1;
                    for (int i = 0; i < nmiss; i++)
                        if (miss[i] == id) {
                            for (int j = i + 1; j < nmiss; j++) miss[j - 1] = miss[j];
                            nmiss--; anyboxless = 0;
                            for (int k = 0; k < nmiss; k++) {
                                mb[k] = NULL;
                                for (int b = 0; b < NBOX; b++)
                                    if (BOXES[b].id == miss[k]) mb[k] = &BOXES[b];
                                if (!mb[k]) anyboxless = 1;
                            }
                            break;
                        }
                }
        }
        if (verbose && nmiss) {
            printf("  D2: still missing after 1:4 fallback:");
            for (int i = 0; i < nmiss; i++) printf(" %d", miss[i]);
            printf("\n");
        }
    }
    dr->nmiss = nmiss;
    dr->nfound = NREQ - nmiss;
}

static const int RAY_LO[]   = { 1800, 1300, 2300, 2800, 3300 };
static const int RAY_FULL[] = { 1800, 2400, 1500, 2100, 2700, 1200, 3000, 3300 };
static const int DIRX[8] = { 1,-1,0,0,1,1,-1,-1 };
static const int DIRZ[8] = { 0,0,1,-1,1,-1,1,-1 };
static const double RINGX[16] = {
    1.00000, 0.92388, 0.70711, 0.38268, 0.00000, -0.38268, -0.70711, -0.92388,
    -1.00000, -0.92388, -0.70711, -0.38268, 0.00000, 0.38268, 0.70711, 0.92388
};
static const double RINGZ[16] = {
    0.00000, 0.38268, 0.70711, 0.92388, 1.00000, 0.92388, 0.70711, 0.38268,
    0.00000, -0.38268, -0.70711, -0.92388, -1.00000, -0.92388, -0.70711, -0.38268
};

static inline int is_full_result(const DRes * dr) {
    return dr->nfound == NREQ;
}

/* 1:4 re-flood to drop land-bridges */
static int isolationOK(Ctx * cx, const DRes * dr) {
    Generator * g = needG(cx);
    memset(cx->iso_vis, 0, (size_t)ISO_W * ISO_W);
    int c0a = dr->cxblk >> 2, c0b = dr->czblk >> 2;
    int fa = c0a, fb = c0b, found = 0;
    for (int r = 0; r <= 8 && !found; r++)
        for (int da = -r; da <= r && !found; da++)
            for (int db = -r; db <= r && !found; db++) {
                int bid = getBiomeAt(g, 4, c0a + da, 64, c0b + db);
                if (bid >= 0 && !isOceanic(bid)) { fa = c0a + da; fb = c0b + db; found = 1; }
            }
    if (!found) return 1;
    int budget = (int)((double)MAX_AREA / 16.0 * ISO_MARGIN) + 1;
    int32_t * q = cx->iso_q;
    int qh = 0, qt = 0;
    q[qt++] = (fb - c0b + ISO_R) * ISO_W + (fa - c0a + ISO_R);
    cx->iso_vis[q[0]] = 1;
    int n = 0;
    static const int DA[4] = { 1,-1,0,0 }, DB[4] = { 0,0,1,-1 };
    while (qh < qt) {
        int p = q[qh++];
        int pa = p % ISO_W - ISO_R + c0a, pb = p / ISO_W - ISO_R + c0b;
        if (++n > budget) return 0;
        for (int k = 0; k < 4; k++) {
            int na = pa + DA[k], nb = pb + DB[k];
            int wa = na - c0a + ISO_R, wb = nb - c0b + ISO_R;
            if (wa < 0 || wa >= ISO_W || wb < 0 || wb >= ISO_W) return 0;
            int ni = wb * ISO_W + wa;
            if (cx->iso_vis[ni]) continue;
            cx->iso_vis[ni] = 1;
            int bid = getBiomeAt(g, 4, na, 64, nb);
            if (bid >= 0 && !isOceanic(bid)) q[qt++] = ni;
        }
    }
    return 1;
}

static void build_missing(DRes * dr) {
    int * miss = dr->miss, nmiss = dr->nmiss;
    dr->missing[0] = 0;
    for (int i = 0; i < nmiss - 1; i++)
        for (int j = i + 1; j < nmiss; j++) {
            int ii = -1, jj = -1;
            for (int k = 0; k < NREQ; k++) {
                if (REQ[k] == miss[i]) ii = k;
                if (REQ[k] == miss[j]) jj = k;
            }
            if (ii >= 0 && jj >= 0 && strcmp(REQNAME[ii], REQNAME[jj]) > 0) {
                int t = miss[i]; miss[i] = miss[j]; miss[j] = t;
            }
        }
    for (int i = 0; i < nmiss; i++)
        for (int k = 0; k < NREQ; k++)
            if (REQ[k] == miss[i]) {
                if (i) strncat(dr->missing, "+", sizeof(dr->missing) - strlen(dr->missing) - 1);
                strncat(dr->missing, REQNAME[k], sizeof(dr->missing) - strlen(dr->missing) - 1);
                break;
            }
}

static void report(Ctx * cx, DRes * dr, int verbose) {
    if (dr->nfound < MIN_REPORT || dr->area > MAX_AREA) return;
    if (VERIFY_ISOLATION && !isolationOK(cx, dr)) {
        if (verbose) printf("  isolation: reject (land bridge to more land)\n");
        return;
    }
    int full = is_full_result(dr);
    if (!full) build_missing(dr);
    pthread_mutex_lock(&g_outmtx);
    if (dr->nfound > atomic_load(&g_bestbio) ||
        (dr->nfound == atomic_load(&g_bestbio) &&
         (atomic_load(&g_bestarea) == 0 || dr->area < atomic_load(&g_bestarea)))) {
        atomic_store(&g_bestbio, dr->nfound);
        atomic_store(&g_bestarea, dr->area);
        atomic_store(&g_bestseed, cx->seed);
    }
    atomic_fetch_add(&g_hits, 1);
    fprintf(g_out, "%llu | %d | %lld | %s\n",
        (unsigned long long)cx->seed, dr->nfound, (long long)dr->area, full ? "ALL" : dr->missing);
    fflush(g_out);
    pthread_mutex_unlock(&g_outmtx);
    printf("Seed: %llu | Biomes: %d | Size: %lld | [%d, %d]\n",
        (unsigned long long)cx->seed, dr->nfound, (long long)dr->area, dr->cxblk, dr->czblk);
    fflush(stdout);
}

static void test_seed(Ctx * cx, uint64_t seed, int verbose) {
    ctx_seed(cx, seed);
    atomic_fetch_add(&g_seeds, 1);

    double candT = 0.2 - A_MARGIN;
    double oceanT = LAND_C + A_MARGIN;
    struct { int x, z; double c; } cand[64];
    int nc = 0;
    double cmin = 99;
    for (int z = -SEARCH_RADIUS; z <= SEARCH_RADIUS; z += A0_STEP)
        for (int x = -SEARCH_RADIUS; x <= SEARCH_RADIUS; x += A0_STEP) {
            double c = sampC(&cx->bnL, NP_CONTINENTALNESS, x, z);
            if (c < cmin) cmin = c;
            if (c >= candT && nc < 64) { cand[nc].x = x; cand[nc].z = z; cand[nc].c = c; nc++; }
        }
    if (!nc || cmin > oceanT) return;

    for (int i = 1; i < nc; i++)
        for (int j = i; j > 0 && cand[j].c > cand[j - 1].c; j--) {
            __typeof__(cand[0]) t = cand[j]; cand[j] = cand[j - 1]; cand[j - 1] = t;
        }
    int keep[4], nk = 0;
    for (int i = 0; i < nc && nk < 4; i++) {
        int ok = 1;
        for (int j = 0; j < nk; j++) {
            int dx = cand[i].x - cand[keep[j]].x, dz = cand[i].z - cand[keep[j]].z;
            if (dx * dx + dz * dz < 1200 * 1200) { ok = 0; break; }
        }
        if (ok) keep[nk++] = i;
    }

    int budget = (int)(MAX_AREA * FILL_FACTOR / (FCELL * FCELL)) + 1;
    int donebb[4][4];
    int ndone = 0;

    for (int k = 0; k < nk; k++) {
        int px = cand[keep[k]].x, pz = cand[keep[k]].z;
        if (verbose) printf("candidate %d at (%d,%d) c_low=%.3f\n", k, px, pz, cand[keep[k]].c);
        int skip = 0;
        for (int j = 0; j < ndone; j++)
            if (px >= donebb[j][0] && px <= donebb[j][2] && pz >= donebb[j][1] && pz <= donebb[j][3]) skip = 1;
        if (skip) { if (verbose) printf("  inside previous fill, skip\n"); continue; }

        /* low-pass ocean rays */
        int dir;
        for (dir = 0; dir < 8; dir++) {
            int ok = 0;
            for (int ri = 0; ri < (int)(sizeof RAY_LO / sizeof *RAY_LO); ri++) {
                double c = sampC(&cx->bnL, NP_CONTINENTALNESS,
                    px + DIRX[dir] * RAY_LO[ri], pz + DIRZ[dir] * RAY_LO[ri]);
                if (c <= oceanT) { ok = 1; break; }
            }
            if (!ok) break;
        }
        if (dir < 8) { if (verbose) printf("  A1 fail dir %d\n", dir); continue; }

        /* ring pre-fill reject */
        if (RING_MIN > 0) {
            int oc = 0;
            for (int a = 0; a < 16; a++)
                if (sampC(&cx->bnL, NP_CONTINENTALNESS,
                        px + RING_R * RINGX[a], pz + RING_R * RINGZ[a]) < LAND_C) oc++;
            if (oc < RING_MIN) {
                if (verbose) printf("  ring reject (ocean %d/16 < %d at r=%d)\n", oc, RING_MIN, RING_R);
                continue;
            }
        }

        /* full-octave ocean rays */
        BiomeNoise * bnc = needC(cx);
        for (dir = 0; dir < 8; dir++) {
            int ok = 0;
            for (int ri = 0; ri < (int)(sizeof RAY_FULL / sizeof *RAY_FULL); ri++) {
                double c = sampC(bnc, NP_CONTINENTALNESS,
                    px + DIRX[dir] * RAY_FULL[ri], pz + DIRZ[dir] * RAY_FULL[ri]);
                if (c <= LAND_C + 0.02) { ok = 1; break; }
            }
            if (!ok) break;
        }
        if (dir < 8) { if (verbose) printf("  A2 fail dir %d\n", dir); continue; }

        /* coarse pre-flood size reject */
        if (COARSE) {
            int cbud = (int)(MAX_AREA * COARSE_MARGIN / (double)(CG * CG)) + 1;
            if (coarseTooBig(cx, px, pz, cbud)) {
                if (verbose) printf("  B0 coarse: confidently > max-area, skip fine flood\n");
                continue;
            }
        }

        /* flood fill of the low-pass field */
        Fill f;
        if (!stageB(cx, px, pz, budget, &f)) {
            if (verbose) printf("  B fail (unbounded or > %d cells)\n", budget);
            continue;
        }
        if (ndone < 4) {
            donebb[ndone][0] = f.bx1; donebb[ndone][1] = f.bz1;
            donebb[ndone][2] = f.bx2; donebb[ndone][3] = f.bz2; ndone++;
        }
        if (verbose)
            printf("  B: cells=%d (~%d blocks^2) bbox=(%d,%d)..(%d,%d) centroid=(%d,%d)\n",
                f.cells, f.cells * FCELL * FCELL, f.bx1, f.bz1, f.bx2, f.bz2, f.cxblk, f.czblk);

        /* exact climate extremes over bbox */
        if (!stageC(cx, &f)) { if (verbose) printf("  C fail (temp/erosion range)\n"); continue; }
        if (verbose) printf("  C: temp/erosion extremes present\n");

        /* rarest-biome climate pre-gate */
        if (!stageC2(cx, &f, verbose)) continue;

        /* biome presence at 1:16 + exact area */
        DRes dr;
        memset(&dr, 0, sizeof dr);
        if (!stageD1(cx, &f, &dr, verbose)) {
            if (verbose) printf("  D1 fail (area > max or degenerate)\n");
            continue;
        }
        if (!dr.reached_bio16) {
            if (verbose) printf("  D1: too many biomes missing at 1:16 (%d)\n", dr.nmiss);
            report(cx, &dr, verbose);
            continue;
        }

        /* 1:4 fallback for biomes missed at 1:16 */
        stageD2(cx, &dr, verbose);
        report(cx, &dr, verbose);
    }
}

static Ctx * ctx_new(void) {
    Ctx * cx = calloc(1, sizeof * cx);
    initBiomeNoise(&cx->bnL, MCV);
    initBiomeNoise(&cx->bnC, MCV);
    initBiomeNoise(&cx->bnT, MCV);
    initBiomeNoise(&cx->bnE, MCV);
    setupGenerator(&cx->g, MCV, 0);
    cx->lc = malloc(FW * FW);
    cx->vis = malloc(FW * FW);
    cx->q = malloc(FW * FW * sizeof(int32_t));
    cx->idsz = (size_t)DMAX * DMAX;
    cx->ids = malloc(cx->idsz * sizeof(int));
    cx->q16 = malloc(cx->idsz * sizeof(int32_t));
    cx->comp = malloc(cx->idsz);
    cx->iso_vis = malloc((size_t)ISO_W * ISO_W);
    cx->iso_q = malloc((size_t)ISO_W * ISO_W * sizeof(int32_t));
    return cx;
}

static void * worker(void * arg) {
    (void)arg;
    Ctx * cx = ctx_new();
    while (!g_stop) {
        uint64_t s0 = atomic_fetch_add(&g_next, CHUNK);
        uint64_t s1 = s0 + CHUNK;
        for (uint64_t s = s0; s < s1 && !g_stop; s++)
            test_seed(cx, s, 0);
    }
    return NULL;
}

static void write_checkpoint(void) {
    FILE * f = fopen(CHECKPOINT_FILE, "w");
    if (f) { fprintf(f, "%llu\n", (unsigned long long)atomic_load(&g_next)); fclose(f); }
}

static void * ckpt(void * a) {
    (void)a;
    while (!g_stop) {
        sleep(CKPT_INTERVAL);
        write_checkpoint();
    }
    return NULL;
}

static void * statmon(void * a) {
    (void)a;
    double t0 = now_s(), tl = t0;
    uint64_t last = 0;
    while (!g_stop) {
        sleep(STATUS_INTERVAL);
        double tn = now_s(), dt = tn - tl;
        if (dt <= 0) dt = 1;
        uint64_t seeds = atomic_load(&g_seeds);
        fprintf(stderr, "[%.0fs] %.2fM  %.0f/s | Hits=%llu | Best=%d/%d %lld %llu\n",
            tn - t0, seeds / 1e6, (seeds - last) / dt,
            (unsigned long long)atomic_load(&g_hits),
            (int)atomic_load(&g_bestbio), NREQ,
            (long long)atomic_load(&g_bestarea), (unsigned long long)atomic_load(&g_bestseed));
        last = seeds; tl = tn;
    }
    return NULL;
}

static void put_biome_name(char * dst, const char * name) {
    int cap = 1;
    while (*name) {
        char c = *name++;
        if (c == '_') { *dst++ = ' '; cap = 1; }
        else { *dst++ = (cap && c >= 'a' && c <= 'z') ? c - 32 : c; cap = 0; }
    }
    *dst = 0;
}

/* measure the continent at (bx,bz): area, biome count, missing biomes */
static int measure_seed(Ctx * cx, uint64_t seed, int bx, int bz) {
    ctx_seed(cx, seed);
    g_areacap = MEASURE_MAX;
    int budget = (int)(MEASURE_MAX * FILL_FACTOR / (FCELL * FCELL)) + 1;
    Fill f;
    DRes dr;
    memset(&dr, 0, sizeof dr);
    if (!stageB(cx, bx, bz, budget, &f) || !stageD1(cx, &f, &dr, 0)) {
        printf("No continent under %d Blocks^2 at (%d, %d)\n", MEASURE_MAX, bx, bz);
        return 0;
    }
    stageD2(cx, &dr, 0);

    char miss[1536];
    miss[0] = 0;
    for (int i = 0; i < dr.nmiss; i++)
        for (int k = 0; k < NREQ; k++)
            if (REQ[k] == dr.miss[i]) {
                if (i) strcat(miss, ", ");
                put_biome_name(miss + strlen(miss), REQNAME[k]);
                break;
            }
    if (!dr.nmiss) strcpy(miss, "NONE");
    printf("Size: %lld Blocks^2 | Biome Count: %d/%d | Missing: %s\n",
        (long long)dr.area, dr.nfound, NREQ, miss);
    return 0;
}

static int runMeasure(const char * ss, const char * sx, const char * sz) {
    char * e;
    unsigned long long seed = strtoull(ss, &e, 10);
    if (e == ss || *e) { fprintf(stderr, "bad seed\n"); return 1; }
    long x = strtol(sx, &e, 10);
    if (e == sx || *e) { fprintf(stderr, "bad x\n"); return 1; }
    long z = strtol(sz, &e, 10);
    if (e == sz || *e) { fprintf(stderr, "bad z\n"); return 1; }
    Ctx * cx = ctx_new();
    return measure_seed(cx, seed, (int)x, (int)z);
}

static int runSingle(const char * ss) {
    char * e;
    unsigned long long v = strtoull(ss, &e, 10);
    if (e == ss || *e) { fprintf(stderr, "bad seed\n"); return 1; }
    g_out = stdout;
    Ctx * cx = ctx_new();
    printf("testing seed %llu (max-area %lld)\n", v, (long long)MAX_AREA);
    test_seed(cx, v, 1);
    return 0;
}

int main(int argc, char ** argv) {
    if (argc >= 2 && !strcmp(argv[1], "single")) {
        if (argc == 3) return runSingle(argv[2]);
        if (argc == 5) return runMeasure(argv[2], argv[3], argv[4]);
        fprintf(stderr, "usage: %s single <seed> [<x> <z>]\n", argv[0]);
        return 1;
    }
    if (argc > 1) { fprintf(stderr, "usage: %s | %s single <seed> [<x> <z>]\n", argv[0], argv[0]); return 1; }

    uint64_t start = START_SEED;
    int resumed = 0;
    FILE * cf = fopen(CHECKPOINT_FILE, "r");
    if (cf) {
        unsigned long long v;
        if (fscanf(cf, "%llu", &v) == 1 && v > start) { start = v; resumed = 1; }
        fclose(cf);
    }
    atomic_store(&g_next, start);

    g_out = fopen(RESULTS_FILE, "a");
    if (!g_out) { perror(RESULTS_FILE); return 1; }

    if (resumed)
        fprintf(stderr, "All Biome Islands | Resume from %llu | %d threads | Size Threshold: %lld\n",
            (unsigned long long)start, NUM_THREADS, (long long)MAX_AREA);
    else
        fprintf(stderr, "All Biome Islands | start %llu | %d threads | Size Threshold: %lld\n",
            (unsigned long long)start, NUM_THREADS, (long long)MAX_AREA);

    signal(SIGINT, on_sigint);
    double t0 = now_s();
    pthread_t th[NUM_THREADS], cth, sth;
    for (int t = 0; t < NUM_THREADS; t++)
        pthread_create(&th[t], NULL, worker, NULL);
    pthread_create(&cth, NULL, ckpt, NULL);
    pthread_create(&sth, NULL, statmon, NULL);
    for (int t = 0; t < NUM_THREADS; t++) pthread_join(th[t], NULL);
    g_stop = 1;
    write_checkpoint();

    unsigned long long seeds = atomic_load(&g_seeds);
    double dt = now_s() - t0;
    fprintf(stderr, "Done: %llu seeds in %.0fs\n", seeds, dt);
    fclose(g_out);
    return 0;
}
