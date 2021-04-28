#define ERROR_MIN_MAX_SAME 1
#define LABEL_NOPEAKS 0
#define LABEL_PEAKSTART 1
#define LABEL_PEAKEND -1
#define LABEL_UNLABELED -2

int FLOPART
(const int *data_vec,
 const double *weight_vec,
 const int data_count,
 const double penalty,
 const int *label_starts,
 const int *label_ends,
 const int *label_types,
 const int label_count,
 // the following matrices are for output.
 double *cost_mat,
 int *end_vec,
 double *mean_vec,
 int *intervals_mat
 );

