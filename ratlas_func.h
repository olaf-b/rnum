/*
 *  ratlas_func.h
 *  Ruby Numeric module 
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#ifndef RATLAS_FUNC_H
#define RATLAS_FUNC_H
#include "ratlas.h"
int ratlas_imax(int a, int b);
int ratlas_imin(int a, int b);
char ratlas_blasflag2lapack(int cblasflag);
VALUE ratlas_rb_complex_new(double re, double im);
VALUE ratlas_get_complex(Ratlas_Complex *c, int offset);
void ratlas_set_complex(Ratlas_Complex *c, int offset, VALUE newval);
void ratlas_store_complex_row(int argc, VALUE *argv, Ratlas_Matrix *rmat,
        int i);
void ratlas_store_complex_column(int argc, VALUE *argv, Ratlas_Matrix *rmat,
        int j);
void ratlas_store_float_row(int argc, VALUE *argv, Ratlas_Matrix *rmat, int i);
void ratlas_store_float_column(int argc, VALUE *argv, Ratlas_Matrix *rmat,
                     int j);
Ratlas_Matrix* ratlas_matrix_alloc (int type, int matrixtype, int nrow,
        int ncol, void *dataptr);
Ratlas_Matrix* ratlas_tripack_alloc (int type, int matrixtype, int nrow,
        int ncol, int upper, void *dataptr);
void ratlas_matrix_free(Ratlas_Matrix *rmat);
void ratlas_matrix_ptr_free(Ratlas_Matrix *rmat);
int ratlas_isnumeric(VALUE obj);
void ratlas_assertnum(VALUE obj);
int ratlas_get_type(VALUE obj);
int ratlas_check_type(VALUE obj, int ratlastype);
VALUE ratlas_matrix_new(VALUE self, VALUE arg);
VALUE ratlas_memadr(VALUE self, VALUE stor);
VALUE ratlas_memsize(VALUE self, VALUE stor);
VALUE ratlas_matrix_from_memadr(VALUE self, VALUE nrow, VALUE ncol,
               VALUE memadr);
VALUE ratlas_cmatrix_new(VALUE self, VALUE arg);
VALUE ratlas_vector_new(VALUE self, VALUE arg);
VALUE ratlas_cvector_new(VALUE self, VALUE arg);
VALUE ratlas_zeros_new(VALUE self, VALUE arg);
VALUE ratlas_czeros_new(VALUE self, VALUE arg);
VALUE ratlas_ones_new(VALUE self, VALUE arg);
VALUE ratlas_eye_new(VALUE self, VALUE arg);
VALUE ratlas_cones_new(VALUE self, VALUE arg);
VALUE ratlas_vec2diag(VALUE self, VALUE stor);
VALUE ratlas_diag2vec(VALUE self, VALUE stor);
VALUE ratlas_diag_new(VALUE self, VALUE stor);
VALUE ratlas_vector_to_a (VALUE self, VALUE storage);
VALUE ratlas_matrix_to_a (VALUE self, VALUE storage);
VALUE ratlas_get_column(VALUE self, VALUE storage, VALUE colidx);
VALUE ratlas_get_columns(VALUE self, VALUE storage, VALUE colidx);
VALUE ratlas_get_row(VALUE self, VALUE storage, VALUE rowidx);
VALUE ratlas_get_rows(VALUE self, VALUE storage, VALUE rowidx);
VALUE ratlas_set_column_bang(VALUE self, VALUE storage, VALUE colidx, VALUE newval);
VALUE ratlas_map_bang(VALUE self, VALUE stor, VALUE block);
VALUE ratlas_each(VALUE self, VALUE stor, VALUE block);
VALUE ratlas_each_with_index(VALUE self, VALUE stor, VALUE block);
VALUE ratlas_each_with_ijindex(VALUE self, VALUE stor, VALUE block);
VALUE ratlas_zip_bang(VALUE self, VALUE stor, VALUE argv, VALUE block);
VALUE ratlas_set_columns_bang(VALUE self, VALUE storage, VALUE colidx, VALUE newval);
VALUE ratlas_set_row_bang(VALUE self, VALUE storage, VALUE rowidx, VALUE newval);
VALUE ratlas_set_rows_bang(VALUE self, VALUE storage, VALUE rowidx, VALUE newval);
VALUE ratlas_size(VALUE self, VALUE storage);
VALUE ratlas_get_one(VALUE self, VALUE storage, VALUE idx);
VALUE ratlas_get_many(VALUE self, VALUE storage, VALUE idx);
VALUE ratlas_get_by_range(VALUE self, VALUE storage, VALUE range);
VALUE ratlas_get2d_one(VALUE self, VALUE storage, VALUE row, VALUE col);
VALUE ratlas_get2d_many(VALUE self, VALUE storage, VALUE rows, VALUE cols);
VALUE ratlas_get2d_by_range(VALUE self, VALUE storage, VALUE rowrange,
               VALUE colrange);
VALUE ratlas_set_one_bang(VALUE self, VALUE storage, VALUE idx, VALUE newval);
VALUE ratlas_set_many_bang(VALUE self, VALUE storage, VALUE idx, VALUE newval);
VALUE ratlas_set_by_range_bang(VALUE self, VALUE storage, VALUE range, 
        VALUE newvals);
VALUE ratlas_set2d_one_bang(VALUE self, VALUE storage, VALUE row, VALUE col,
        VALUE newval);
VALUE ratlas_set2d_many_bang(VALUE self, VALUE storage, VALUE row, VALUE col,
        VALUE newval);
VALUE ratlas_set2d_by_range_bang(VALUE self, VALUE storage, VALUE rowrange, 
        VALUE colrange, VALUE newvals);
VALUE ratlas_concat(VALUE self, VALUE stor1, VALUE stor2);
VALUE ratlas_hcat(VALUE self, VALUE stor1, VALUE stor2);
VALUE ratlas_vcat(VALUE self, VALUE stor1, VALUE stor2);
VALUE ratlas_indgen_bang(VALUE self, VALUE stor, VALUE start);
Ratlas_Matrix* ratlas_matrix_clone(Ratlas_Matrix *orig);
Ratlas_Matrix* ratlas_matrix_dup(Ratlas_Matrix *orig);
VALUE ratlas_storage_clone(VALUE self, VALUE oldstore);
VALUE ratlas_storage_dup(VALUE self, VALUE oldstore);
VALUE ratlas_storage_alloc(VALUE self, VALUE arg);
VALUE ratlas_complex_storage_alloc(VALUE self, VALUE arg);
VALUE ratlas_re(VALUE self, VALUE stor);
VALUE ratlas_im(VALUE self, VALUE stor);
VALUE ratlas_reshape_bang(VALUE self, VALUE stor, VALUE arg);
VALUE ratlas_mul(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_mul_bang(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_div(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_div_bang(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_add_bang(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_add(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_rem_bang(VALUE self, VALUE stor, VALUE denom);
VALUE ratlas_madd(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_madd_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_msub(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_msub_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_mmul_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_mmul(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_mdiv_bang(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_mdiv(VALUE self, VALUE transa, VALUE astor, VALUE transc,
        VALUE cstor);
VALUE ratlas_sub(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_sub_bang(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_pow_bang(VALUE self, VALUE stor, VALUE alpha);
VALUE ratlas_sqrt_bang(VALUE self, VALUE stor);
VALUE ratlas_log_bang(VALUE self, VALUE stor);
VALUE ratlas_exp_bang(VALUE self, VALUE stor);
VALUE ratlas_sum(VALUE self, VALUE stor);
VALUE ratlas_rowsum(VALUE self, VALUE stor);
VALUE ratlas_colsum(VALUE self, VALUE stor);
VALUE ratlas_colon(VALUE self, VALUE trans1, VALUE stor1, VALUE trans2,
        VALUE stor2);
VALUE ratlas_admul_bang(VALUE self, VALUE astor, VALUE dstor);
VALUE ratlas_damul_bang(VALUE self, VALUE dstor, VALUE astor);
VALUE ratlas_dadd_bang(VALUE self, VALUE astor, VALUE dstor);
VALUE ratlas_dsub_bang(VALUE self, VALUE astor, VALUE dstor);
VALUE ratlas_colpiv_bang(VALUE self, VALUE stor, VALUE piv);
VALUE ratlas_transpose(VALUE self, VALUE stor);
void ratlas_matrix_sym2dense(int uplo, Ratlas_Matrix *mat);
VALUE ratlas_storage_sym2dense_bang(VALUE self, VALUE uplo, VALUE stor);
void ratlas_matrix_uplo_nullify(int unitdiag, int uplo, Ratlas_Matrix *mat);
VALUE ratlas_storage_diag2dense(VALUE self, VALUE dstor);
VALUE ratlas_storage_real2complex(VALUE self, VALUE restor);
#endif
