/*

	this builds a kdtree as a Structure of Array (SoA), optimized for coalesced access on gpu
	[x1, x2, ..., xn, y1, y2, ... yn, z1, z2, ... zn] for n points

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // log2()
#include <string.h>

#include "kdtree.h"

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

void print_point(kdtree* kdt, int idx) { 

	if(idx > kdt->array_lim) {
		printf("Print point %i out of bounds %i.\n", idx, kdt->num_points);
		return;
	}
	printf("idx=%i, ", idx);
	if(kdt->emptys[idx]) printf("null\n");
	else printf("(%.12f, %.12f, %.12f)\n", kdt->x[idx], kdt->y[idx], kdt->z[idx]);
}

void kdtree_print(kdtree* kdt) {

	// [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	// new_line: 1, 3, 7
	printf("--kdtree: %i points, %i nodes\n", kdt->num_points, kdt->num_nodes);
	int new_line = 1;
	int level = 0;
	printf("level 0, axis 0\n");
	for(int i=0; i<kdt->array_lim; i++) {
		if(i == new_line) {
			level++;
			printf("level %i, axis %i\n", level, level%3);
			new_line = (new_line * 2) + 1;		
		}
		if(!kdt->emptys[i]) {
			if(kdt->emptys[2*i + 1] && kdt->emptys[2*i + 2]) printf("**leaf** ");
			printf("%i (x=%.12f, y=%.12f, z=%.12f)\n", i, kdt->x[i], kdt->y[i], kdt->z[i]);
		}
		else printf("%i (null)\n", i);
	}
	printf("--\n");
}

void kdtree_build_r(double** points, int axis, kdtree* kdt, int l, int r, int index) {

//	printf("valid %i, index=%i, max_index=%i\n", index < kdt->array_lim, index, kdt->array_lim);

	// leaf node ; point
	if(r == l || l > r || r-l == 0) {
	
		// if need to allocate more memory
		if(index*2 + 2 > kdt->array_lim ) {			
			int prev_lim = kdt->array_lim;
			kdt->array_lim = index*2 + 2; // the index to the right child
			kdt->emptys = realloc(kdt->emptys, (kdt->array_lim + 1) * sizeof(char));
			kdt->axes = realloc(kdt->axes, (kdt->array_lim + 1) * sizeof(int));

			if(!kdt->emptys || !kdt->axes) {
				printf("Realloc %lu bytes kdt->emptys or kdt->axes failed\n", (kdt->array_lim) * sizeof(char) + 1);
				return;
			}	
			for(int i=prev_lim; i<=kdt->array_lim; i++) {
				kdt->emptys[i] = 1;
				kdt->axes[i] = -1;
			}	
			int num_bytes = (kdt->array_lim + 1) * sizeof(double); 			
			kdt->x = realloc(kdt->x, num_bytes);
			kdt->y = realloc(kdt->y, num_bytes);
			kdt->z = realloc(kdt->z, num_bytes);	
			if(!kdt->x || !kdt->y || !kdt->z) {
				printf("kdtree_build_r() realloc failed for new buffer size %i\n", num_bytes);
				return;
			}
		}
		kdt->x[index] = points[r][0]; // x
		kdt->y[index] = points[r][1]; // y
		kdt->z[index] = points[r][2]; // z 
		kdt->axes[index] = axis;
		kdt->num_nodes++;
		// mark this node as not empty
		kdt->emptys[index] = 0;
	
		// mark the children as empty
		kdt->emptys[index*2 + 1] = 1;
		kdt->emptys[index*2 + 2] = 1;

		return;
	}
	
	int split_index = l + (r-l)/2; // median index of the sorted points
	double split;
	if( ((r-l)+1)%2 == 0) split = (points[split_index][axis] + points[split_index+1][axis])/2; // even number of elements 
	else split = points[split_index][axis]; // odd number of elments

	// if need to allocate more memory
	if(index*2 + 2 > kdt->array_lim) {			
		int prev_lim = kdt->array_lim;
		kdt->array_lim = index*2 + 2; // the index to the right child
		kdt->emptys = realloc(kdt->emptys, (kdt->array_lim + 1) * sizeof(char));
		kdt->axes = realloc(kdt->axes, (kdt->array_lim + 1) * sizeof(int));
	
		if(!kdt->emptys || !kdt->axes) {
			printf("Realloc %lu bytes kdt->emptys or axes failed\n", (kdt->array_lim) * sizeof(char) + 1);
			return;
		}	
		for(int i=prev_lim; i<=kdt->array_lim; i++) {
			kdt->emptys[i] = 1;
			kdt->axes[i] = -1;
		}	
		int num_bytes = (kdt->array_lim + 1) * sizeof(double); 			
		kdt->x = realloc(kdt->x, num_bytes);
		kdt->y = realloc(kdt->y, num_bytes);
		kdt->z = realloc(kdt->z, num_bytes);	
		if(!kdt->x || !kdt->y || !kdt->z) {
			printf("kdtree_build_r() realloc failed for new buffer size %i\n", num_bytes);
			return;
		}
	}
	
	kdt->x[index] = 0;
	kdt->y[index] = 0;
	kdt->z[index] = 0;
	
	switch(axis) {
		case 0:
			kdt->x[index] = split;
			break;
		case 1:
			kdt->y[index] = split;
			break;
		case 2:
			kdt->z[index] = split;
			break;
		default:
			printf("This should never happen.\n");
			break;	
	}
	kdt->axes[index] = axis;		
	kdt->num_nodes++;
	// mark this node as not empty
	kdt->emptys[index] = 0;

	// left child	
	quicksort(points, l, split_index, (axis+1)%3); // sort in the next axis 
	kdtree_build_r(points, (axis+1)%3, kdt, l, split_index, 2*index + 1);	

	// right child
	quicksort(points, split_index+1, r, (axis+1)%3); // sort in the next axis 
	kdtree_build_r(points, (axis+1)%3, kdt, split_index+1, r, 2*index + 2);	

}

kdtree* kdtree_build(double** points, int num_points) {
	
	kdtree* kdt = (kdtree*) malloc(sizeof(kdtree));
	kdt->num_points = num_points;
	kdt->num_nodes = 0;
	kdt->array_lim = num_points * 2; // the highest index in the array ; aka size of the array ; this is NOT the same as the number of nodes in the tree
	kdt->emptys = (char*) malloc( (kdt->array_lim + 1) * sizeof(char)); // must use this to mark the gaps in the array as empty nodes
	kdt->axes = (int*) malloc( (kdt->array_lim + 1) * sizeof(int));
	if(!kdt->emptys || !kdt->axes) {
		printf("malloc() kdt->emptys or kdt->axes failed\n");
		return NULL;
	}
	memset(kdt->emptys, 1, (kdt->array_lim + 1) * sizeof(char)); // all the nodes are empty initially
	memset(kdt->axes, -1, (kdt->array_lim + 1) * sizeof(int));
	
	kdt->x = malloc( (kdt->array_lim + 1) * sizeof(double));
	kdt->y = malloc( (kdt->array_lim + 1) * sizeof(double));
	kdt->z = malloc( (kdt->array_lim + 1)* sizeof(double));
	if(!kdt->x || !kdt->y || !kdt->z) {
		printf("malloc() kdt points failed\n");
		return NULL;
	}

	quicksort(points, 0, num_points-1, 0); // sort in the first axis 
	printf("kdtree_build_r\n");
	kdtree_build_r(points, 0, kdt, 0, num_points-1, 0);	
	//	kdtree_print(kdt);

	return kdt;
}
