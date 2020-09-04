#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP fitglmm_brent(SEXP Y_in, SEXP X_in, SEXP Phi_in, SEXP sqrtW_in, SEXP method_in, SEXP dispersion_in, SEXP tau_in, SEXP fixtau_in, SEXP tol_in, SEXP taumin_in, SEXP taumax_in, SEXP tauregion_in);

SEXP fitglmm_nm(SEXP Y_in, SEXP X_in, SEXP q_in, SEXP Phi_in, SEXP W_in, SEXP method_in, SEXP dispersion_in, SEXP tau_in, SEXP fixtau_in, SEXP maxiter_in, SEXP tol_in);

SEXP fitglmm_ai(SEXP Y_in, SEXP X_in, SEXP q_in, SEXP Phi_in, SEXP ng_in, SEXP group_in, SEXP W_in, SEXP tau_in, SEXP fixtau_in);

SEXP glmm_score_text(SEXP res_in, SEXP P_in, SEXP infile_in, SEXP outfile_in, SEXP tol_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP p_in, SEXP infile_nrow_skip_in, SEXP infile_sep_in, SEXP infile_na_in, SEXP infile_ncol_skip_in, SEXP infile_ncol_print_in, SEXP infile_header_print_in, SEXP nperbatch_in, SEXP select_in);

SEXP glmm_score_text_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP infile_in, SEXP outfile_in, SEXP tol_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP p_in, SEXP infile_nrow_skip_in, SEXP infile_sep_in, SEXP infile_na_in, SEXP infile_ncol_skip_in, SEXP infile_ncol_print_in, SEXP infile_header_print_in, SEXP nperbatch_in, SEXP select_in);

SEXP glmm_score_bed(SEXP res_in, SEXP P_in, SEXP bimfile_in, SEXP bedfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in);

SEXP glmm_score_bed_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP bimfile_in, SEXP bedfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in);

SEXP glmm_wald_text(SEXP n_in, SEXP snp_in, SEXP infile_in, SEXP tol_in, SEXP center_in, SEXP miss_method_in, SEXP p_in, SEXP infile_nrow_skip_in, SEXP infile_sep_in, SEXP infile_na_in, SEXP infile_ncol_skip_in, SEXP infile_ncol_print_in, SEXP snp_col_in, SEXP select_in);

SEXP glmm_wald_bed(SEXP n_in, SEXP snp_in, SEXP bimfile_in, SEXP bedfile_in, SEXP center_in, SEXP miss_method_in, SEXP select_in);

SEXP glmm_score_bgen13(SEXP res_in, SEXP P_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

SEXP glmm_score_bgen11(SEXP res_in, SEXP P_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

SEXP glmm_score_bgen13_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

SEXP glmm_score_bgen11_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

SEXP bgenHeader(SEXP bgenfile_in);

SEXP getVariantPos(SEXP bgenfile_in, SEXP offset_in, SEXP mbgen_in, SEXP nbgen_in, SEXP compression_in, SEXP layout_in, SEXP cores_in);
  
static const R_CallMethodDef R_CallDef[]  = {
  {"fitglmm_brent", (DL_FUNC) &fitglmm_brent, 12},
  {"fitglmm_nm", (DL_FUNC) &fitglmm_nm, 11},
  {"fitglmm_ai", (DL_FUNC) &fitglmm_ai, 9},
  {"glmm_score_text", (DL_FUNC) &glmm_score_text, 19},
  {"glmm_score_text_sp", (DL_FUNC) &glmm_score_text_sp, 21},
  {"glmm_score_bed", (DL_FUNC) &glmm_score_bed, 12},
  {"glmm_score_bed_sp", (DL_FUNC) &glmm_score_bed_sp, 14},
  {"glmm_wald_text", (DL_FUNC) &glmm_wald_text, 14},
  {"glmm_wald_bed", (DL_FUNC) &glmm_wald_bed, 7},
  {"glmm_score_bgen13", (DL_FUNC) &glmm_score_bgen13, 17},
  {"glmm_score_bgen11", (DL_FUNC) &glmm_score_bgen11, 17},
  {"glmm_score_bgen13_sp", (DL_FUNC) &glmm_score_bgen13_sp, 19},
  {"glmm_score_bgen11_sp", (DL_FUNC) &glmm_score_bgen11_sp, 19},
  {"bgenHeader", (DL_FUNC) &bgenHeader, 1},
  {"getVariantPos", (DL_FUNC) &getVariantPos, 7},
  {NULL, NULL, 0}
};

void R_init_GMMAT(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
