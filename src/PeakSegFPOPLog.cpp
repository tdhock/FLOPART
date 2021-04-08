/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>
#include <R.h> // for Rprintf
#include "PeakSegFPOPLog.h"

int PeakSegFPOPLog
  (int *data_vec, double *weight_vec, int data_count,
   double penalty,
   // the following matrices are for output.
   // cost_mat and intervals_mat store the optimal cost and number of intervals
   // at each time point, for the up and down cost models.
   // end_vec and mean_vec store the best model up to and including the
   // last data point. 
   double *cost_mat, //data_count x 2.
   int *end_vec,//data_count
   double *mean_vec,//data_count
   int *intervals_mat,//data_count x 2
   int *label_starts,
   int *label_ends,
   int *label_types,
   int label_count){ //matrix where rows are (start, end, labelType)
  
  bool at_beginning = false;
  bool at_end = false;  
  int curr_label_index = 0;
  int curr_label_type = label_types[0]; //1 = peakStart, 0 = noPeak, -1 = peakEnd
  int label_to_check; //used to see if current data point is in a label
  bool label_found;
  
  
  double min_log_mean=INFINITY, max_log_mean=-INFINITY;
  for(int data_i=0; data_i<data_count; data_i++){
    
    double log_data = log( (double)(data_vec[data_i]) );
    if(log_data < min_log_mean){
      min_log_mean = log_data;
    }
    if(max_log_mean < log_data){
      max_log_mean = log_data;
    }
  }
  if(min_log_mean == max_log_mean){
    return ERROR_MIN_MAX_SAME;
  }
  std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count * 2);
  PiecewisePoissonLossLog *up_cost, *down_cost, *up_cost_prev, *down_cost_prev;
  PiecewisePoissonLossLog min_prev_cost, cost_of_change_up, cost_of_change_down;
  int verbose=0;
  double cum_weight_i = 0.0, cum_weight_prev_i;
  for(int data_i=0; data_i<data_count; data_i++){
    cum_weight_i += weight_vec[data_i];
    up_cost = &cost_model_mat[data_i];
    down_cost = &cost_model_mat[data_i + data_count];
    
    label_to_check = curr_label_index;
    label_found = false;
    at_beginning = false;
    at_end = false;
    while(data_i >= label_starts[label_to_check] && !label_found){
      if(data_i <= label_ends[label_to_check]){
        if(data_i==label_starts[label_to_check]){
          at_beginning = true;
Rprintf("data point %d is at beginning of label \n", data_i);
        } 
        if(data_i==label_ends[label_to_check]){
          at_end = true;
Rprintf("data point %d is at end of label \n", data_i);
        }
        
        label_found = true; 
        curr_label_type = label_types[label_to_check];
        curr_label_index = label_to_check;
        
      }
      label_to_check++;
    }
    if(!label_found){
      curr_label_type = LABEL_UNLABELED;
    }
    if(curr_label_type==LABEL_PEAKSTART){
      Rprintf("data point %d is in peak start \n", data_i);
    }
    if(curr_label_type==LABEL_PEAKEND){
      Rprintf("data point %d is in peak end \n", data_i);
    }
    if(curr_label_type==LABEL_NOPEAKS){
      Rprintf("data point %d is in no peaks \n", data_i);
    }
    if(curr_label_type==LABEL_UNLABELED){
      Rprintf("data point %d is unlabeled \n", data_i);
    }
    
    
    
    
    //---UP COST CALCULATIONS---
    
    //if in middle of peak end
    if(data_i == 0){
      up_cost->piece_list.emplace_back
      (1.0, -data_vec[0], 0.0,
       min_log_mean, max_log_mean, -1, false);
Rprintf("initializing up cost! \n");      
    }
    
    else if(curr_label_type == LABEL_PEAKEND && !at_beginning && !at_end){
Rprintf("UP: in middle of peak end \n");
      //up cost is just previous up cost
      up_cost = up_cost_prev;
    }
    
    //if in no peaks or beginning of peak start or end of peak end
    else if(curr_label_type == LABEL_NOPEAKS || (curr_label_type == LABEL_PEAKSTART 
                                                   && !at_beginning)
              || (curr_label_type == LABEL_PEAKEND && at_end)){
Rprintf("UP: in no peaks or not beginning of peak start or end of peak end \n");   
Rprintf("SETTING UP INF \n");
              //up cost is infinite
              up_cost->set_infinite();
    }
    
    //if unlabeled or not beginning of peak start or beginning of peak end  
    else if(curr_label_type == LABEL_UNLABELED 
              || (curr_label_type== LABEL_PEAKSTART && !at_beginning) 
              || (curr_label_type == LABEL_PEAKEND && at_beginning)){
Rprintf("UP: unlabeled or not beginning of peakStart or beginning of peakEnd \n");              
              //up cost is min of prev_up, min_less(prev+down + penalty)
              cost_of_change_up.set_to_min_less_of(down_cost_prev, verbose);
      //seg end  
      cost_of_change_up.set_prev_seg_end(data_i-1);
      cost_of_change_up.addPenalty(penalty, cum_weight_prev_i);
      up_cost->set_to_min_env_of(up_cost_prev, &cost_of_change_up, verbose);
      
    }
     
    //---DOWN COST CALCULATIONS--- 
    if(data_i==0){
      // initialization Cdown_1(m)=gamma_1(m)/w_1
      down_cost->piece_list.emplace_back
      (1.0, -data_vec[0], 0.0,
       min_log_mean, max_log_mean, -1, false);
Rprintf("initializing down cost! \n");      
      
    }
    
    
    //if at beginning of peak end or end of peak start, infinite down
    else if((curr_label_type == LABEL_PEAKEND && at_beginning) ||
       (curr_label_type == LABEL_PEAKSTART && at_end)){
Rprintf("DOWN: at beginning of peak end of end of peak start \n");
Rprintf("SETTING DOWN INF \n");
      down_cost->set_infinite();
    }
    
    //if not at beginning of noPeaks or middle of peak start
    else if((curr_label_type == LABEL_NOPEAKS && !at_beginning)
         || (curr_label_type == LABEL_PEAKSTART && !at_beginning && !at_end)){
Rprintf("DOWN: at not beginning of no peaks or middle of peak start \n");     
      down_cost = down_cost_prev;
      
    }
    
    //if unlabeled or first or no peaks/peakStart or not at beginning of peak end
    else if(curr_label_type == LABEL_UNLABELED 
              || (curr_label_type == LABEL_PEAKSTART && at_beginning)
              || (curr_label_type == LABEL_NOPEAKS && at_beginning)
              || (curr_label_type == LABEL_PEAKEND && !at_beginning)){
Rprintf("DOWN: in unlabeled, beginning peak start, beginning no peaks,  or not beginning of peak end! doing change down or beginning no peaks \n");              
              //down cost is min(prev_down, min_more(prev_up + penalty))
      cost_of_change_down.set_to_min_more_of(up_cost_prev, verbose);
      cost_of_change_down.set_prev_seg_end(data_i-1);
      cost_of_change_down.addPenalty(penalty, cum_weight_prev_i);
      down_cost->set_to_min_env_of(down_cost_prev, &cost_of_change_down, verbose);
    } 
    
    
    up_cost->adjustWeights(cum_weight_prev_i, cum_weight_i, weight_vec,
                           data_i, data_vec);
Rprintf("before adjusting weights, is down cost infinite? %d \n", down_cost->is_infinite());
    
    down_cost->adjustWeights(cum_weight_prev_i, cum_weight_i, weight_vec,
                             data_i, data_vec);
Rprintf("after adjusting weights, is down cost infinite? %d \n", down_cost->is_infinite());
    
    cum_weight_prev_i = cum_weight_i;
    up_cost_prev = up_cost;
    down_cost_prev = down_cost;
Rprintf("leaving for loop, is down cost infinite? %d \n", down_cost->is_infinite());
  }
  // Decoding the cost_model_vec, and writing to the output matrices.
  double best_cost, best_log_mean, prev_log_mean;
  int prev_seg_end=data_count;
  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  for(int i=0; i<2*data_count; i++){
    up_cost = &cost_model_mat[i];
    intervals_mat[i] = up_cost->piece_list.size();
    up_cost->Minimize
      (cost_mat+i, &best_log_mean,
       &prev_seg_end, &prev_log_mean);
  }
  // last segment is down (offset N) so the second to last segment is
  // up (offset 0).
  int prev_seg_offset = 0;
  down_cost = &cost_model_mat[data_count*2-1];
  down_cost->Minimize
    (&best_cost, &best_log_mean,
     &prev_seg_end, &prev_log_mean);
  mean_vec[0] = exp(best_log_mean);
  end_vec[0] = prev_seg_end;
  int out_i=1;
  while(0 <= prev_seg_end){
    // up_cost is actually either an up or down cost.
    up_cost = &cost_model_mat[prev_seg_offset + prev_seg_end];
    Rprintf("decoding out_i=%d prev_seg_end=%d prev_seg_offset=%d\n", out_i, prev_seg_end, prev_seg_offset);
    up_cost->print();
    if(prev_log_mean != INFINITY){
      //equality constraint inactive
      best_log_mean = prev_log_mean;
    }
    up_cost->findMean
      (best_log_mean, &prev_seg_end, &prev_log_mean);
    mean_vec[out_i] = exp(best_log_mean);
Rprintf("setting mean_vec[out_i] to %d \n", mean_vec[out_i]);
    end_vec[out_i] = prev_seg_end;
    // change prev_seg_offset and out_i for next iteration.
    if(prev_seg_offset==0){
      //up_cost is actually up
      prev_seg_offset = data_count;
    }else{
      //up_cost is actually down
      prev_seg_offset = 0;
    }
    out_i++;
  }
return 0;
}

