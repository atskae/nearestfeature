/*

	this builds a kdtree as a Structure of Array (SoA), optimized for coalesced access on gpu
	[x1, x2, ..., xn, y1, y2, ... yn, z1, z2, ... zn] for n points

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // log2()

#include "kdtree.h"

// lat range: -90 to 90 ; longt range: -180 to 180
// invalid point; (lat=-150, longt=-250) = (x=1073.585280926619, y=4323.631369885733, z=4554.477733167408) 
#define INVALID_X 1073 

int partition(double** points, int l, int r, int dim) {
	
	int swap = l+1;
	for(int i=l+1; i<=r; i++) { // iterate through all the elements except the pivot value (which is always the left-most element)
	
		if(points[i][dim] < points[l][dim]) { // compare the current element with the pivot ; if element is in the wrong spot, switch places 
			double* p = points[swap];
			points[swap] = points[i];
			points[i] = p;
			swap++;
		}		
	}

	// swap the pivot
	double* p = points[l]; // pivot value
	points[l] = points[swap-1];
	points[swap-1] = p;

	return swap-1; // the partition point
}

void quicksort(double** points, int l, int r, int dim) { 

	if(l > r) return;

	int sep = partition(points, l, r, dim);
	
	quicksort(points, l, sep-1, dim);
	quicksort(points, sep+1, r, dim);
}

void kdtree_print(kdtree* kdt) {

	// [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	// new_line: 1, 3, 7
	printf("--kdtree: %i points, %i nodes\n", kdt->num_points, kdt->num_nodes);
	int new_line = 1;
	for(int i=0; i<kdt->num_nodes; i++) {
		if(i == new_line) {
			printf("\n");
			new_line = (new_line * 2) + 1;		
		}
		if( !(fabs(kdt->x[i]- INVALID_X) < 1e-32) ) printf("%i (x=%.3f, y=%.3f, z=%.3f), ", i, kdt->x[i], kdt->y[i], kdt->z[i]);
		else printf("%i (null) ", i);
	}
	printf("\n--\n");
}

void kdtree_build_r(double** points, int axis, kdtree* kdt, int l, int r, int index) {

	// leaf node ; point
	if(r == l || l > r || r-l == 0) {

		// if need to allocate more memory
		if(kdt->num_nodes ==  kdt->num_nodes_limit) {
			kdt->num_nodes_limit = kdt->num_nodes_limit * 2; // double the size
			int num_bytes = kdt->num_nodes_limit * sizeof(double); 

			kdt->x = realloc(kdt->x, num_bytes);
			for(int i=kdt->num_nodes; i<kdt->num_nodes_limit; i++) {
				kdt->x[i] = INVALID_X;
			}
			kdt->y = realloc(kdt->y, num_bytes);
			kdt->z = realloc(kdt->z, num_bytes);	
			if(!kdt->x || !kdt->y || !kdt->z) {
				printf("kdtree_build_r() realloc failed for new buffer size %i\n", num_bytes);
				exit(1);
			}
		}	

		kdt->x[index] = points[r][0]; // x
		kdt->y[index] = points[r][1]; // y
		kdt->z[index] = points[r][2]; // z 
		kdt->num_nodes++;
		
		// set children to null
		kdt->x[index*2 + 1] = INVALID_X;
		kdt->x[index*2 + 2] = INVALID_X;
		kdt->num_nodes+=2;
		
		return;
	}
	
	int split_index = l + (r-l)/2; // median index
	double split;
	if( ((r-l)+1)%2 == 0) split = (points[split_index][axis] + points[split_index+1][axis])/2; // even number of elements 
	else split = points[split_index][axis]; // odd number of elments

	// if need to allocate more memory
	if(kdt->num_nodes == kdt->num_nodes_limit) {
		kdt->num_nodes_limit = kdt->num_nodes_limit * 2; // double the size
		int num_bytes = kdt->num_nodes_limit * sizeof(double); 
		
		kdt->x = realloc(kdt->x, num_bytes);
		for(int i=kdt->num_nodes; i<kdt->num_nodes_limit; i++) {
			kdt->x[i] = INVALID_X;
		}
		kdt->y = realloc(kdt->y, num_bytes);
		kdt->z = realloc(kdt->z, num_bytes);	
		if(!kdt->x || !kdt->y || !kdt->z) {
			printf("kdtree_build_r() realloc failed for new buffer size %i\n", num_bytes);
			exit(1);
		}
	}

	switch(axis) {
		case 0:
			kdt->x[index] = split;
			kdt->y[index] = 0;
			kdt->z[index] = 0;
			break;
		case 1:
			kdt->x[index] = 0;
			kdt->y[index] = split;
			kdt->z[index] = 0;
			break;
		case 2:
			kdt->x[index] = 0;
			kdt->y[index] = 0;
			kdt->z[index] = split;
			break;
	
	}		
	kdt->num_nodes++;

	// left child	
	quicksort(points, l, split_index, (axis+1)%3); // sort in the next axis 
	kdtree_build_r(points, (axis+1)%3, kdt, l, split_index, 2*index+1);	

	// right child
	quicksort(points, split_index+1, r, (axis+1)%3); // sort in the next axis 
	kdtree_build_r(points, (axis+1)%3, kdt, split_index+1, r, 2*index+2);	

}

kdtree* kdtree_build(double** points, int num_points) {
	
	kdtree* kdt = (kdtree*) malloc(sizeof(kdtree));
	kdt->num_points = num_points;
	kdt->num_nodes = 0;
	kdt->num_nodes_limit = (num_points * 3) * 2;
	kdt->x = malloc(kdt->num_nodes_limit * sizeof(double));
	for(int i=0; i<kdt->num_nodes_limit; i++) {
		kdt->x[i] = INVALID_X;
	}
	kdt->y = malloc(kdt->num_nodes_limit * sizeof(double));
	kdt->z = malloc(kdt->num_nodes_limit * sizeof(double));

	kdtree_build_r(points, 0, kdt, 0, num_points-1, 0);	
//	kdtree_print(kdt);

	return kdt;
}

