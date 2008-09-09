/*
 *  ratlas.c
 *  Ruby Numerical Library - RNum
 *    (C) Copyright 2006- by Olaf Trygve BERGLIHN
 *
 *  This program is free software.  You can distribute/modify this program
 *  under the same terms as Ruby itself.  There is absolutely no warranty with
 *  this software.
*/

#include "ratlas.h"
#include "ratlas_lapack.h"
#include "ratlas_cblas.h"

#define PROPAGATESYM 1 /* undef when full support for symmetric matrices is
              implemented */

ID id_Complex, id_Range, id_atre, id_atim, id_imag, id_call;
ID id_first, id_last, id_exclude_end;
VALUE complex_class, range_class, ratlas_module, blas_module, lapack_module;
VALUE ratlas_storage_class;

/* Module init stuff */
void Init_ratlas() {
    rb_require("complex");
    id_Complex = rb_intern("Complex");
    id_Range = rb_intern("Range");
    complex_class = rb_const_get(rb_cObject, id_Complex);
    range_class = rb_const_get(rb_cObject, id_Range);
    id_atre = rb_intern("@real");
    id_atim = rb_intern("@image");
    id_imag = rb_intern("imag");
    id_call = rb_intern("call");
    id_first = rb_intern("first");
    id_last = rb_intern("last");
    id_exclude_end = rb_intern("exclude_end?");
    ratlas_module = rb_define_module("RAtlas");
    blas_module = rb_define_module("Blas");
    lapack_module = rb_define_module("Lapack");

    ratlas_storage_class = rb_define_class_under(ratlas_module, 
            "Storage", rb_cObject);

    /* BLAS/ATLAS constants */
    rb_define_const(ratlas_module, "ROWMAJOR", INT2FIX(CblasRowMajor));
    rb_define_const(ratlas_module, "COLMAJOR", INT2FIX(CblasColMajor));
    rb_define_const(ratlas_module, "NOTRANS", INT2FIX(CblasNoTrans));
    rb_define_const(ratlas_module, "TRANS", INT2FIX(CblasTrans));
    rb_define_const(ratlas_module, "CONJTRANS", INT2FIX(CblasConjTrans));
    rb_define_const(ratlas_module, "UPPER", INT2FIX(CblasUpper));
    rb_define_const(ratlas_module, "LOWER", INT2FIX(CblasLower));
    rb_define_const(ratlas_module, "NONUNITDIAG", INT2FIX(CblasNonUnit));
    rb_define_const(ratlas_module, "UNITDIAG", INT2FIX(CblasUnit));
    rb_define_const(ratlas_module, "LEFT", INT2FIX(CblasLeft));
    rb_define_const(ratlas_module, "RIGHT", INT2FIX(CblasRight));


    /* Matrix creation and manipulation */
    rb_define_singleton_method(ratlas_module, "matrix", 
            ratlas_matrix_new, 1);
    rb_define_singleton_method(ratlas_module, "memadr", ratlas_memadr, 1);
    rb_define_singleton_method(ratlas_module, "memsize", ratlas_memsize, 1);
    rb_define_singleton_method(ratlas_module, "matrix_from_memadr",
            ratlas_matrix_from_memadr, 3);
    rb_define_singleton_method(ratlas_module, "cmatrix",
            ratlas_cmatrix_new, 1);
    rb_define_singleton_method(ratlas_module, "vector",
            ratlas_vector_new, 1);
    rb_define_singleton_method(ratlas_module, "cvector",
            ratlas_cvector_new, 1);
    rb_define_singleton_method(ratlas_module, "zeros", ratlas_zeros_new, 1);
    rb_define_singleton_method(ratlas_module, "czeros",
            ratlas_czeros_new, 1);
    rb_define_singleton_method(ratlas_module, "ones", ratlas_ones_new, 1);
    rb_define_singleton_method(ratlas_module, "eye", ratlas_eye_new, 1);
    rb_define_singleton_method(ratlas_module, "indgen!", ratlas_indgen_bang,
            2);
    rb_define_singleton_method(ratlas_module, "cones", ratlas_cones_new, 1);
    rb_define_singleton_method(ratlas_module, "clone", 
            ratlas_storage_clone, 1);
    rb_define_singleton_method(ratlas_module, "dup", ratlas_storage_dup, 1);
    rb_define_singleton_method(ratlas_module, "alloc", 
            ratlas_storage_alloc, 1);
    rb_define_singleton_method(ratlas_module, "calloc", 
            ratlas_complex_storage_alloc, 1);
    
    rb_define_singleton_method(ratlas_module, "reshape!", 
            ratlas_reshape_bang, 2);
    
    rb_define_singleton_method(ratlas_module, "matrix_to_a", 
            ratlas_matrix_to_a, 1);
    rb_define_singleton_method(ratlas_module, "vector_to_a", 
            ratlas_vector_to_a, 1);
    rb_define_singleton_method(ratlas_module, "get_one", ratlas_get_one, 2);
    rb_define_singleton_method(ratlas_module, "get_many", ratlas_get_many,
                   2);
    rb_define_singleton_method(ratlas_module, "get_by_range", 
            ratlas_get_by_range, 2);
    rb_define_singleton_method(ratlas_module, "get2d_one", 
            ratlas_get2d_one, 3);
    rb_define_singleton_method(ratlas_module, "get2d_many", 
            ratlas_get2d_many, 3);
    rb_define_singleton_method(ratlas_module, "get2d_by_range", 
            ratlas_get2d_by_range, 3);

    rb_define_singleton_method(ratlas_module, "size", ratlas_size, 1);
    rb_define_singleton_method(ratlas_module, "set_one!", ratlas_set_one_bang, 3);
    rb_define_singleton_method(ratlas_module, "set_many!", 
            ratlas_set_many_bang, 3);
    rb_define_singleton_method(ratlas_module, "set_by_range!",
                   ratlas_set_by_range_bang, 3);
    rb_define_singleton_method(ratlas_module, "set2d_one!",
                   ratlas_set2d_one_bang, 4);
    rb_define_singleton_method(ratlas_module, "set2d_many!",
                   ratlas_set2d_many_bang, 4);
    rb_define_singleton_method(ratlas_module, "set2d_by_range!",
                   ratlas_set2d_by_range_bang, 4);
    rb_define_singleton_method(ratlas_module, "concat", ratlas_concat, 2);
    rb_define_singleton_method(ratlas_module, "hcat", ratlas_hcat, 2);
    rb_define_singleton_method(ratlas_module, "vcat", ratlas_vcat, 2);
    rb_define_singleton_method(ratlas_module, "map!", ratlas_map_bang, 2);
    rb_define_singleton_method(ratlas_module, "zip!", ratlas_zip_bang, 3);
    rb_define_singleton_method(ratlas_module, "column", 
            ratlas_get_column, 2);
    rb_define_singleton_method(ratlas_module, "columns", 
            ratlas_get_columns, 2);
    rb_define_singleton_method(ratlas_module, "row", ratlas_get_row, 2);
    rb_define_singleton_method(ratlas_module, "rows", ratlas_get_rows, 2);
    rb_define_singleton_method(ratlas_module, "set_column!", 
            ratlas_set_column_bang, 3);
    rb_define_singleton_method(ratlas_module, "set_columns!", 
            ratlas_set_columns_bang, 3);
    rb_define_singleton_method(ratlas_module, "set_row!", 
            ratlas_set_row_bang, 3);
    rb_define_singleton_method(ratlas_module, "set_rows!",
            ratlas_set_rows_bang, 3);
    rb_define_singleton_method(ratlas_module, "re", ratlas_re, 1);
    rb_define_singleton_method(ratlas_module, "im", ratlas_im, 1);
    rb_define_singleton_method(ratlas_module, "sum", ratlas_sum, 1);
    rb_define_singleton_method(ratlas_module, "rowsum", ratlas_rowsum, 1);
    rb_define_singleton_method(ratlas_module, "colsum", ratlas_colsum, 1);
    

    rb_define_singleton_method(ratlas_module, "vec2diag", ratlas_vec2diag, 
            1);
    rb_define_singleton_method(ratlas_module, "diag2vec", ratlas_diag2vec, 
            1);
    rb_define_singleton_method(ratlas_module, "diag", ratlas_diag_new, 1);
    rb_define_singleton_method(ratlas_module, "colpiv!", 
            ratlas_colpiv_bang, 2);
    rb_define_singleton_method(ratlas_module, "transpose", 
            ratlas_transpose, 1);
    rb_define_singleton_method(ratlas_module, "sym2dense!", 
            ratlas_storage_sym2dense_bang, 2);
    rb_define_singleton_method(ratlas_module, "diag2dense", 
            ratlas_storage_diag2dense, 1);
    rb_define_singleton_method(ratlas_module, "real2complex", 
            ratlas_storage_real2complex, 1);
        

    /* ratlas scalar operations */
    rb_define_singleton_method(ratlas_module, "add!", ratlas_add_bang, 2);
    rb_define_singleton_method(ratlas_module, "add", ratlas_add, 2);
    rb_define_singleton_method(ratlas_module, "madd!", ratlas_madd_bang, 4);
    rb_define_singleton_method(ratlas_module, "madd", ratlas_madd, 4);
    rb_define_singleton_method(ratlas_module, "sub!", ratlas_sub_bang, 2);
    rb_define_singleton_method(ratlas_module, "sub", ratlas_sub, 2);
    rb_define_singleton_method(ratlas_module, "msub!", ratlas_msub_bang, 4);
    rb_define_singleton_method(ratlas_module, "msub", ratlas_msub, 4);
    rb_define_singleton_method(ratlas_module, "mul!", ratlas_mul_bang, 2);
    rb_define_singleton_method(ratlas_module, "mul", ratlas_mul, 2);
    rb_define_singleton_method(ratlas_module, "mmul!", ratlas_mmul_bang, 4);
    rb_define_singleton_method(ratlas_module, "mmul", ratlas_mmul, 4);
    rb_define_singleton_method(ratlas_module, "mdiv!", ratlas_mdiv_bang, 4);
    rb_define_singleton_method(ratlas_module, "mdiv", ratlas_mdiv, 4);

    rb_define_singleton_method(ratlas_module, "div!", ratlas_div_bang, 2);
    rb_define_singleton_method(ratlas_module, "div", ratlas_div, 2);
    rb_define_singleton_method(ratlas_module, "pow!", ratlas_pow_bang, 2);
    rb_define_singleton_method(ratlas_module, "sqrt!", ratlas_sqrt_bang, 1);
    rb_define_singleton_method(ratlas_module, "log!", ratlas_log_bang, 1);
    rb_define_singleton_method(ratlas_module, "exp!", ratlas_exp_bang, 1);
    rb_define_singleton_method(ratlas_module, "colon", ratlas_colon, 4);
    rb_define_singleton_method(ratlas_module, "rem!", ratlas_rem_bang, 2);
    

    /* ratlas diagonal operations */
    rb_define_singleton_method(ratlas_module, "damul!", ratlas_damul_bang,
            2);
    rb_define_singleton_method(ratlas_module, "admul!", ratlas_admul_bang,
            2);
    rb_define_singleton_method(ratlas_module, "dadd!", ratlas_dadd_bang, 2);
    rb_define_singleton_method(ratlas_module, "dsub!", ratlas_dsub_bang, 2);


    /* Lapack functions */
    rb_define_singleton_method(lapack_module, "gesv!", ratlas_gesv_bang, 2);
    rb_define_singleton_method(lapack_module, "posv!", ratlas_posv_bang, 3);
    rb_define_singleton_method(lapack_module, "sysv!", ratlas_sysv_bang, 3);
    rb_define_singleton_method(lapack_module, "getrf!", ratlas_getrf_bang,
            1);
    rb_define_singleton_method(lapack_module, "potrf!", ratlas_potrf_bang,
            2);
    rb_define_singleton_method(lapack_module, "sytrf!", ratlas_sytrf_bang,
            2);
    rb_define_singleton_method(lapack_module, "getrs!", ratlas_getrs_bang,
            4);
    rb_define_singleton_method(lapack_module, "potrs!", ratlas_potrs_bang,
            3);
    rb_define_singleton_method(lapack_module, "sytrs!", ratlas_sytrs_bang,
            4);
    rb_define_singleton_method(lapack_module, "getri!", ratlas_getri_bang,
            2);
    rb_define_singleton_method(lapack_module, "potri!", ratlas_potri_bang,
            2);
    rb_define_singleton_method(lapack_module, "sytri!", ratlas_sytri_bang,
            3);
    rb_define_singleton_method(lapack_module, "gesvd!", ratlas_gesvd_bang, 1);

    /* Atlas Blas level 1 */
    rb_define_singleton_method(blas_module, "scal!", ratlas_scal_bang, 2);
    rb_define_singleton_method(blas_module, "axpy!", ratlas_axpy_bang, 3);
    rb_define_singleton_method(blas_module, "dot", ratlas_dot, 2);
    rb_define_singleton_method(blas_module, "dotc", ratlas_dotc, 2);
    rb_define_singleton_method(blas_module, "nrm2", ratlas_nrm2, 1);
    rb_define_singleton_method(blas_module, "asum", ratlas_asum, 1);
    rb_define_singleton_method(blas_module, "iamax", ratlas_iamax, 1);

    /* Atlas Blas level 2 */
    rb_define_singleton_method(blas_module, "gemv!", ratlas_gemv_bang, 6);
    rb_define_singleton_method(blas_module, "symv!", ratlas_symv_bang, 5);
    rb_define_singleton_method(blas_module, "trmv!", ratlas_trmv_bang, 5);
    rb_define_singleton_method(blas_module, "trsv!", ratlas_trsv_bang, 5);
    rb_define_singleton_method(blas_module, "ger!", ratlas_ger_bang, 4);
    rb_define_singleton_method(blas_module, "gerc!", ratlas_gerc_bang, 4);
    rb_define_singleton_method(blas_module, "syr!", ratlas_syr_bang, 3);
    rb_define_singleton_method(blas_module, "her!", ratlas_syr_bang, 3);
    rb_define_singleton_method(blas_module, "syr2!", ratlas_syr2_bang, 4);
    rb_define_singleton_method(blas_module, "her2!", ratlas_syr2_bang, 4);
    
    /* Atlas Blas level 3 */
    rb_define_singleton_method(blas_module, "gemm!", ratlas_gemm_bang, 7);
    rb_define_singleton_method(blas_module, "symm!", ratlas_symm_bang, 6);
    rb_define_singleton_method(blas_module, "hemm!", ratlas_hemm_bang, 6);
    rb_define_singleton_method(blas_module, "syrk!", ratlas_syrk_bang, 5);
    rb_define_singleton_method(blas_module, "herk!", ratlas_herk_bang, 5);
    rb_define_singleton_method(blas_module, "syr2k!", ratlas_syr2k_bang,
            6);
    rb_define_singleton_method(blas_module, "her2k!", ratlas_her2k_bang,
            6);
}
