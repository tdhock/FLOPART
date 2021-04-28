#include <Rcpp.h>
#include <R.h>
#include "FLOPART.h"

//' Interface to FLOPART C++ code
//'
//' @param data_vec Integer vector of non-negative count data
//' @param weight_vec Numeric vector of positive weights (same size as data_vec)
//' @param penalty non-negative real-valued penalty (larger for fewer peaks)
//' @param label_type_vec Integer vector of label types
//' @param label_start_vec Integer vector of label starts
//' @param label_end_vec Integer vector of label ends
//' @return List with named elements: cost_mat and intervals_mat (one row for each data point, first column up, second down), segments_df (one row for each segment in the optimal model)
// [[Rcpp::export]]
Rcpp::List FLOPART_interface
(const Rcpp::IntegerVector data_vec,
 const Rcpp::NumericVector weight_vec,
 const double penalty,
 const Rcpp::IntegerVector label_type_vec,
 const Rcpp::IntegerVector label_start_vec,
 const Rcpp::IntegerVector label_end_vec){
  int data_count = data_vec.size();
  if(data_count < 2){
    Rcpp::stop("need at least two data points");
  }
  if(data_count != weight_vec.size()){
    Rcpp::stop("data_vec and weight_vec should be same size");
  }
  int label_count = label_type_vec.size();
  if(label_count != label_start_vec.size()){
    Rcpp::stop("label_start_vec and label_type_vec should be same size");
  }
  if(label_count != label_end_vec.size()){
    Rcpp::stop("label_end_vec and label_type_vec should be same size");
  }
  Rcpp::NumericMatrix cost_mat(data_count, 2);
  Rcpp::IntegerMatrix intervals_mat(data_count, 2);
  Rcpp::IntegerVector rev_end_vec(data_count);
  Rcpp::NumericVector rev_mean_vec(data_count);
  int status = FLOPART
    (&data_vec[0], &weight_vec[0],
     data_count, penalty,
     &label_type_vec[0], &label_start_vec[0], &label_end_vec[0], label_count,
     &cost_mat[0], &rev_end_vec[0], &rev_mean_vec[0], &intervals_mat[0]
     );
  if(status == ERROR_MIN_MAX_SAME){
    Rcpp::stop("data[i]=%d for all i", data_vec[0]);
  }
  //convert to segments data table.
  int seg_count=1;
  while(seg_count < data_count && 0 <= rev_end_vec[seg_count-1]){
    seg_count++;
  }
  Rcpp::NumericVector seg_mean_vec(seg_count);
  Rcpp::IntegerVector seg_start_vec(seg_count);
  Rcpp::IntegerVector seg_end_vec(seg_count);
  for(int seg_i=0; seg_i < seg_count; seg_i++){
    int mean_index = seg_count-1-seg_i;
    seg_mean_vec[seg_i] = rev_mean_vec[mean_index];
    if(mean_index==0){
      seg_end_vec[seg_i] = data_count;
    }else{
      seg_end_vec[seg_i] = rev_end_vec[mean_index-1]+1;
    }
    if(seg_i==0){
      seg_start_vec[seg_i] = 1;
    }else{
      seg_start_vec[seg_i] = rev_end_vec[mean_index]+2;
    }
  }
  return Rcpp::List::create
    (Rcpp::Named("cost_mat", cost_mat),
     Rcpp::Named("intervals_mat", intervals_mat),
     Rcpp::Named("segments_df", Rcpp::DataFrame::create
		 (Rcpp::Named("mean", seg_mean_vec),
		  Rcpp::Named("start", seg_start_vec),
		  Rcpp::Named("end", seg_end_vec)))

     );
}
