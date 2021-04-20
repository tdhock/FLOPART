#define LABEL_NOPEAKS 0
#define LABEL_PEAKSTART 1
#define LABEL_PEAKEND -1
#define LABEL_UNLABELED -2


int FLOPART
(int *data_vec, double *weight_vec,
 int data_count, double penalty,
 // the following matrices are for output.
 double *cost_mat,
 int *end_vec,
 double *mean_vec,
 int *intervals_mat,
 int *label_starts,
 int *label_ends,
 int *label_types,
 int label_count);

