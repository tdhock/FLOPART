/*-*- compile-command: "R CMD INSTALL .." -*- */
  
#include "funPieceListLog.h"
#include "FLOPART.h"
#include <R.h>
#include <R_ext/Rdynload.h>
  

  void FLOPART_interface
    (int *data_ptr, double *weight_ptr,
     int *data_count, double *penalty,
     double *cost_mat, int *end_vec,
     double *mean_vec, int *intervals_mat,
     int *label_starts, int *label_ends,
     int *label_types, int *label_count){
    int status = FLOPART
    (data_ptr, weight_ptr,
     *data_count, *penalty,
     cost_mat, end_vec, mean_vec, intervals_mat, label_starts,
     label_ends, label_types, *label_count);
    if(status == ERROR_MIN_MAX_SAME){
      error("data[i]=%d for all i", data_ptr[0]);
    }
  }
  
  R_CMethodDef cMethods[] = {
    {"FLOPART_interface",
     (DL_FUNC) &FLOPART_interface, 12
      //,{INTSXP, REALSXP, REALSXP, INTSXP}
    },
    {NULL, NULL, 0} 
  };
  
  extern "C" {
    void R_init_FLOPART(DllInfo *info) {
      R_registerRoutines(info, cMethods, NULL, NULL, NULL);
      //R_useDynamicSymbols call says the DLL is not to be searched for
      //entry points specified by character strings so .C etc calls will
      //only find registered symbols.
      R_useDynamicSymbols(info, FALSE);
    }
    
  }
  