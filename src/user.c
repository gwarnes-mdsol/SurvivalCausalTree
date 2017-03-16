/*
 * split.Rule = user (set temporarily as CT)
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;

int
userinit(int n, double *y[], int maxcat, char **error,
		int *size, int who, double *wt, double *treatment,
		int bucketnum, int bucketMax, double *train_to_est_ratio)
{
	if (who == 1 && maxcat > 0) {
		graycode_init0(maxcat);
		countn = (int *) ALLOC(2 * maxcat, sizeof(int));
		tsplit = countn + maxcat;
		treatment_effect = (double *) ALLOC(8 * maxcat, sizeof(double));
		wts = treatment_effect + maxcat;
		trs = wts + maxcat;
		sums = trs + maxcat;
		wtsums = sums + maxcat;
		trsums = wtsums + maxcat;
		wtsqrsums = trsums + maxcat;
		trsqrsums = wtsqrsums + maxcat;
	}
	*size = 1;
	*train_to_est_ratio = n * 1.0 / ct.NumHonest;
	return 0;
}



void
userss(int n, double *y[], double *value,  double *con_mean, double *tr_mean,
		double *risk, double *wt, double *treatment, double max_y,
		double alpha, double train_to_est_ratio, double *propensity, double *censoringProb)
{
  //Rprintf("%f\n", censoringProb[2]);
	int i;
	double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */
	double ttreat = 0.,tcon = 0.;
	double effect;
	double tr_var, con_var;
	double con_sqr_sum = 0., tr_sqr_sum = 0.;
	double pp = 0.;
	for (i = 0; i < n; i++) {
		temp1 += *y[i] * wt[i] * treatment[i] / propensity[i];
		temp0 += *y[i] * wt[i] * (1 - treatment[i]) / (1-propensity[i]);
		twt += wt[i];
		ttreat += wt[i] * treatment[i] / propensity[i];
		tcon += wt[i] * (1-treatment[i]) / (1-propensity[i]);

		tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
	}

	effect = temp1 / ttreat - temp0 / tcon;
	tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
	con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));

	*tr_mean = temp1 / ttreat;
	*con_mean = temp0 / tcon;
	*value = effect;
	*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect;
		//(1 - alpha) * (1 + train_to_est_ratio) * twt * (tr_var /ttreat  + con_var / (twt - ttreat));
}
int compare_double(const void *a,const void *b) {
  double *x = (double *) a;
  double *y = (double *) b;
  // return *x - *y; // this is WRONG...
  if (*x < *y) return -1;
  else if (*x > *y) return 1; return 0;
}
int compare_int( const void* a, const void* b )
{
  if( *(int*)a == *(int*)b ) return 0;
  return *(int*)a < *(int*)b ? -1 : 1;
}

void user(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split,
		int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
		double train_to_est_ratio, double *propensity, double *censoringProb)
{
  //Rprintf("%f\t",propensity[2]);
  Rprintf("%f\t",censoringProb[3]);
	int i, j;
	double temp;
	double left_sum, right_sum;
	double left_tr_sum, right_tr_sum;
	double left_tr, right_tr;
	double left_wt, right_wt;
	int left_n, right_n;
	double best;
	int direction = LEFT;
	int where = 0;
	double node_effect, left_effect, right_effect;
	double left_temp, right_temp;
	int min_node_size = minsize;

	double tr_var, con_var;
	double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
	double left_tr_var, left_con_var, right_tr_var, right_con_var;

	//added definition
	double left_tr_sum2,right_tr_sum2;

	right_wt = 0.;
	right_tr = 0.;
	right_sum = 0.;
	right_tr_sum = 0.;
	right_sqr_sum = 0.;
	right_tr_sqr_sum = 0.;
	right_n = n;
	for (i = 0; i < n; i++) {
		right_wt += wt[i];
		right_tr += wt[i] * treatment[i];
		right_sum += *y[i] * wt[i];
		right_tr_sum += *y[i] * wt[i] * treatment[i];
		//added

		//right_tr_sum2 += *y[i]*wt[i]*treatment[i]*alpha[i];
		//
		//right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
		//right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
	}
  //TODO: For every sample, estimate K(Y) with kaplan-meier estimator
  //qsort( elements, num_elem, sizeof(int), compare_int );


	temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
	/*
	 tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
	con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
		- (right_sum - right_tr_sum) * (right_sum - right_tr_sum)
		/ ((right_wt - right_tr) * (right_wt - right_tr));
	 */
	node_effect =   temp * temp * right_wt ;

	if (nclass == 0) {
		/* continuous predictor */
		left_wt = 0;
		left_tr = 0;
		left_n = 0;
		left_sum = 0;
		left_tr_sum = 0;
		left_sqr_sum = 0;
		left_tr_sqr_sum = 0;
		best = 0;

		for (i = 0; right_n > edge; i++) {
			left_wt += wt[i];
			right_wt -= wt[i];
			left_tr += wt[i] * treatment[i];
			right_tr -= wt[i] * treatment[i];
			left_n++;
			right_n--;
			temp = *y[i] * wt[i] * treatment[i];
			left_tr_sum += temp;
			right_tr_sum -= temp;
			left_sum += *y[i] * wt[i];
			right_sum -= *y[i] * wt[i];
			//temp = (*y[i]) *  (*y[i]) * wt[i];
			//left_sqr_sum += temp;
			//right_sqr_sum -= temp;
			//temp = (*y[i]) * (*y[i]) * wt[i] * treatment[i];
			//left_tr_sqr_sum += temp;
			//right_tr_sqr_sum -= temp;
      //printf("%f",propensity);
			if (x[i + 1] != x[i] && left_n >= edge &&
					(int) left_tr >= min_node_size &&
					(int) left_wt - (int) left_tr >= min_node_size &&
					(int) right_tr >= min_node_size &&
					(int) right_wt - (int) right_tr >= min_node_size) {

				left_temp = left_tr_sum / left_tr -
					(left_sum - left_tr_sum) / (left_wt - left_tr);

				right_temp = right_tr_sum / right_tr -
					(right_sum - right_tr_sum) / (right_wt - right_tr);

				left_effect = left_temp * left_temp * left_wt;
				/*
				 right_effect = alpha * right_temp * right_temp * right_wt;
				*/
				//printf("%f",alpha[1]);
				right_effect = right_temp * right_temp * right_wt;

				temp = left_effect + right_effect - node_effect;
				//temp=fabs(right_temp-left_temp);
				if (temp > best) {
					best = temp;
					where = i;
					if (left_temp < right_temp)
						direction = LEFT;
					else
						direction = RIGHT;
				}
			}
		}

		*improve = best;
		if (best > 0) {         /* found something */
			csplit[0] = direction;
			*split = (x[where] + x[where + 1]) / 2;
		}
	}
}


double
userpred(double *y, double wt, double treatment, double *yhat, double propensity)
{
	double ystar;
	double temp;

	ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
	temp = ystar - *yhat;
	return temp * temp * wt;
}
