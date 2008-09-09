/*
 *  ratlas_cblas.h
 *  Ruby Numerical Library - RNum
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#ifndef RATLAS_CBLAS_H
#define RATLAS_CBLAS_H
#include "ratlas.h"
VALUE ratlas_gesv_bang(VALUE self, VALUE astor, VALUE bstor);
VALUE ratlas_getrf_bang(VALUE self, VALUE astor);
VALUE ratlas_getrs_bang(VALUE self, VALUE trans, VALUE astor, VALUE rowperm,
        VALUE bstor);
VALUE ratlas_getri_bang(VALUE self, VALUE astor, VALUE rowperm);
VALUE ratlas_scal_bang(VALUE self, VALUE alpha, VALUE xstorage);
VALUE ratlas_axpy_bang(VALUE self, VALUE alpha, VALUE xstorage,  VALUE ystorage);
VALUE ratlas_dot(VALUE self, VALUE xstor, VALUE ystor);
VALUE ratlas_dotc(VALUE self, VALUE xstor, VALUE ystor);
VALUE ratlas_nrm2(VALUE self, VALUE xstor);
VALUE ratlas_asum(VALUE self, VALUE xstor);
VALUE ratlas_iamax(VALUE self, VALUE xstor);
VALUE ratlas_gemv_bang(VALUE self, VALUE transa, VALUE alpha, VALUE astor,
        VALUE xstor, VALUE beta, VALUE ystor);
VALUE ratlas_symv_bang(VALUE self, VALUE alpha, VALUE astor, VALUE xstor,
        VALUE beta, VALUE ystor);
VALUE ratlas_trmv_bang(VALUE self, VALUE astor, VALUE xstor, VALUE uplo,
        VALUE trans, VALUE diag);
VALUE ratlas_trsv_bang(VALUE self, VALUE astor, VALUE xstor, VALUE uplo,
        VALUE trans, VALUE diag);
VALUE ratlas_ger_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE ystor,
        VALUE astor);
VALUE ratlas_gerc_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE ystor,
        VALUE astor);
VALUE ratlas_syr_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE astor);
VALUE ratlas_syr2_bang(VALUE self, VALUE alpha, VALUE xstor, VALUE ystor,
        VALUE astor);
VALUE ratlas_gemm_bang(VALUE self, VALUE transa, VALUE transb, VALUE alpha,
        VALUE astor, VALUE bstor, VALUE beta, VALUE cstor);
VALUE ratlas_symm_bang(VALUE self, VALUE side, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor);
VALUE ratlas_hemm_bang(VALUE self, VALUE side, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor);
VALUE ratlas_syrk_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE beta, VALUE cstor);
VALUE ratlas_herk_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE beta, VALUE cstor);
VALUE ratlas_syr2k_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor);
VALUE ratlas_her2k_bang(VALUE self, VALUE trans, VALUE alpha, VALUE astor,
        VALUE bstor, VALUE beta, VALUE cstor);
#endif
