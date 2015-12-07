#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "datamodel.h"

float find_correlation(float *vals, int start, int end, int offset);

int main(int argc, char**argv) {

	int len = 1000;
	float *vals1 = (float *)malloc(1000*sizeof(float));
	float *vals2 = (float *)malloc(1000*sizeof(float));
	int i = 0;
	
	while (i < len) {
		// I want 100 to be 2*pi
		float val = sin( i * 2 * 3.14159 / 100);
		vals1[i] = val;
		val = cos(i * 2 * 3.14159 / 100);
		vals2[i] = val;
		//printf("Val %d %f\n",i,val);
		i++;
	}
	
	i = 0;
	
	while (i < len-200) {
		int l = 1;
		while (l < 100) {
			float corr = find_correlation(vals2,i,i+100,l);
			printf("Corr %d\t%d\t%d\t%f\n",i,l,(i/100),corr);
			l++;
		}
		i++;
	}
  return 0;
}

	float find_correlation(float *vals, int start, int end, int offset) {
	
	float sum_sq_x = 0.0;
	float sum_sq_y = 0.0;
	float sum_coproduct = 0.0;
	float mean_x = vals[start];
	float mean_y = vals[start+offset];
		
		
	int len = (end-start+1);

	int i = 1;
		
	while (i < len) {
			
		float sweep = (i - 1.0) / i;
		float delta_x = vals[i+start] - mean_x;
		float delta_y = vals[i+start+offset] - mean_y;
		sum_sq_x += delta_x * delta_x * sweep;
		sum_sq_y += delta_y * delta_y * sweep;
		sum_coproduct += delta_x * delta_y * sweep;
		mean_x += delta_x / i;
		mean_y += delta_y / i;
		i++;
	}
		
	float pop_sd_x = sqrt( sum_sq_x / len);
	float pop_sd_y = sqrt( sum_sq_y / len );
	float cov_x_y  = sum_coproduct / len;
		
	if (pop_sd_x > 0 && pop_sd_y > 0) {
		return cov_x_y / (pop_sd_x * pop_sd_y);
	} else {
		return 0.0;
	}
}

