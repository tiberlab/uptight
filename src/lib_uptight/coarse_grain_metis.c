/* Optional METIS bridge.  Keeping METIS types in C avoids ABI differences in
 * idx_t between METIS builds and the Fortran library.
 *
 * Also provides weak default implementations of the host-log callbacks that
 * Fortran (coarse_grain / JD) expects.  Inside tiberCAD, UptWrapper.cpp supplies
 * strong definitions that override these and route into Messages::info. */
#include <stdio.h>
#include <stdlib.h>

__attribute__((weak)) void upt_log_message(const char *message, int length)
{
  if (message == NULL || length <= 0) {
    fputs("(icgn)\n", stdout);
    return;
  }
  fputs("(icgn) ", stdout);
  fwrite(message, 1, (size_t)length, stdout);
  fputc('\n', stdout);
}

__attribute__((weak)) void upt_jd_log_message(const char *message, int length)
{
  if (message == NULL || length <= 0) {
    fputs("(jd)\n", stdout);
    return;
  }
  fputs("(jd) ", stdout);
  fwrite(message, 1, (size_t)length, stdout);
  fputc('\n', stdout);
}

#ifdef UPT_HAVE_METIS
#include <metis.h>
#endif

int upt_icgn_metis_partition(int nvtxs, const int *xadj, const int *adjncy,
                           const int *vwgt, const int *adjwgt, int nparts,
                           int ufactor, int seed, int *part)
{
#ifdef UPT_HAVE_METIS
  idx_t *xa = NULL, *ad = NULL, *vw = NULL, *aw = NULL, *pt = NULL;
  idx_t nv = (idx_t)nvtxs, nc = (idx_t)nparts, edgecut = 0;
  idx_t options[METIS_NOPTIONS];
  int i, nedge = xadj[nvtxs];
  xa = malloc((size_t)(nvtxs + 1) * sizeof(*xa));
  ad = malloc((size_t)nedge * sizeof(*ad));
  vw = malloc((size_t)nvtxs * sizeof(*vw));
  aw = malloc((size_t)nedge * sizeof(*aw));
  pt = malloc((size_t)nvtxs * sizeof(*pt));
  if (!xa || !ad || !vw || !aw || !pt) goto fail;
  for (i=0;i<=nvtxs;i++) xa[i] = (idx_t)xadj[i];
  for (i=0;i<nedge;i++) { ad[i] = (idx_t)adjncy[i]; aw[i] = (idx_t)adjwgt[i]; }
  for (i=0;i<nvtxs;i++) vw[i] = (idx_t)vwgt[i];
  idx_t ncon = 1;
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_UFACTOR] = (idx_t)ufactor;
  options[METIS_OPTION_SEED] = (idx_t)seed;
  if (METIS_PartGraphKway(&nv, &ncon, xa, ad, vw, NULL, aw, &nc, NULL,
                          NULL, options, &edgecut, pt) != METIS_OK) goto fail;
  for (i=0;i<nvtxs;i++) part[i] = (int)pt[i];
  free(xa); free(ad); free(vw); free(aw); free(pt);
  return 0;
fail:
  free(xa); free(ad); free(vw); free(aw); free(pt);
  return -1;
#else
  (void)nvtxs; (void)xadj; (void)adjncy; (void)vwgt; (void)adjwgt;
  (void)nparts; (void)ufactor; (void)seed; (void)part;
  return -1;
#endif
}
