/*

	This builds a kdtree as an array, optimized for coalesced access on GPU
	[x1, x2, ..., xn, y1, y2, ... yn, z1, z2, ... zn] for n points

*/

#include <stdlib.h>
#include <stdio.h>
#include "kdtree.h"

int kdtree_partition(double* points, int l, int r, int dim, int num_points) {
	
	// the index into the points array where each dimension starts
	int x_idx = num_points*0;
	int y_idx = num_points*1;
	int z_idx = num_points*2;
	int this_idx = num_points*dim;
	int swap = l+1; // the element to the right of the pivot

	// iterate through all the elements except the pivot value (which is always the left-most element)
	for(int i=l+1; i<=r; i++) { 	
		if(points[this_idx + i] < points[this_idx + l]) { // compare the current element with the pivot ; if element is in the wrong spot, switch places 
			// save old values from each dimension of this point
			double x = points[x_idx + i];
			double y = points[y_idx + i];
			double z = points[z_idx + i];	

			points[x_idx + i] = points[x_idx + swap];
			points[y_idx + i] = points[y_idx + swap];
			points[z_idx + i] = points[z_idx + swap];

			points[x_idx + swap] = x;
			points[y_idx + swap] = y;
			points[z_idx + swap] = z;
	
			swap++;
		}		
	}

	// swap the pivot
//	double* p = points[l]; // pivot value
	double x = points[x_idx + l];
	double y = points[y_idx + l];
	double z = points[z_idx + l];

	points[x_idx + l] = points[x_idx + (swap - 1)];
	points[y_idx + l] = points[y_idx + (swap - 1)];
	points[z_idx + l] = points[z_idx + (swap - 1)];

	points[x_idx + (swap - 1)] = x;
	points[y_idx + (swap - 1)] = y;
	points[z_idx + (swap - 1)] = z;

	return swap-1; // the partition point
}

void kdtree_quicksort(double* points, int num_points, int l, int r, int dim) { 

	if(l > r) return;

	int p;
	p = kdtree_partition(points, l, r, dim, num_points);
	
	kdtree_quicksort(points, l, p-1, dim, num_points);
	kdtree_quicksort(points, p+1, r, dim, num_points);
}

void kdtree_build_r(double* points, int num_points, double* node, int l, int r, int dim) {

	//quicksort(points, num_points, l, r, dim);		

}

// assumes data contains points groups of 3 doubles: x, y, z
double* kdtree_build(double* points, int num_points) {

	double* root = malloc(num_points * 3);
	kdtree_build_r(points, num_points, root, 0, num_points-1, 0);
	
	return root;
}

