/*

	this builds a kdtree as a Structure of Array (SoA), optimized for coalesced access on gpu
	[x1, x2, ..., xn, y1, y2, ... yn, z1, z2, ... zn] for n points

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // log2()
#include <string.h>

#include "kdtree.h"

/*
	
	GPU

*/

int max_leaf_idx = 0;
double INVALID_X = 0.0; // to indicate inner nodes that are empty ;

double getRadians(double degrees) { return (degrees * PI)/180; }

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

void kdtree_build_gpu_r(double** points, int axis, kdtree* kdt, int l, int r, int index, double** pts_x, double** pts_y, double** pts_z) {

	// base case: leaf node ; point
	if(r == l || l > r || r-l == 0) {
	
		// if need to allocate more memory
		if(index*2 + 2 > kdt->len_nodes-1) {			
			int prev_lim = kdt->len_nodes - 1;
			kdt->len_nodes = (index*2 + 2) + 1; // index of the right child
			kdt->leaves = realloc(kdt->leaves, kdt->len_nodes/8 + 1);			
			if(!kdt->leaves) {
				printf("Realloc kdt->emptys failed\n");
				return;
			}
			// set newly allocated nodes to not leaves	
			for(int i=prev_lim; i<kdt->len_nodes/8 + 1; i++) {
				kdt->leaves[i] = 0;
			}
	
			int num_bytes = kdt->len_nodes * sizeof(double); 			
			*pts_x = realloc(*pts_x, num_bytes);
			*pts_y = realloc(*pts_y, num_bytes);
			*pts_z = realloc(*pts_z, num_bytes);	
			if(!(*pts_x) || !(*pts_y) || !(*pts_z) ) {
				printf("kdtree_build_r() realloc failed for new buffer size %i\n", num_bytes);
				return;
			}
		}
		(*pts_x)[index] = points[r][0]; // x
		(*pts_y)[index] = points[r][1]; // y
		(*pts_z)[index] = points[r][2]; // z 
		kdt->num_nodes++;

	//	printf("leaf is built!: idx=%i ", index);
	//	printf("(%.2f, %.2f, %.2f)\n", (*pts_x)[index], (*pts_y)[index], (*pts_z)[index]);
		// set this bit to indicate that this node is a leaf
		kdt->leaves[index/8] |= 0x1<<(index%8);
		max_leaf_idx = (index > max_leaf_idx) ? index : max_leaf_idx;
		
		// mark that the leaf nodes are empty
		(*pts_x)[index*2 + 1] = INVALID_X;
		(*pts_x)[index*2 + 2] = INVALID_X;	

		return;
	} // if leaf ; end
	
	int split_index = l + (r-l)/2; // median index of the sorted points
	double split;
	if( ((r-l)+1)%2 == 0) split = (points[split_index][axis] + points[split_index+1][axis])/2; // even number of elements 
	else split = points[split_index][axis]; // odd number of elments

	// if need to allocate more memory
	if(index*2 + 2 > kdt->len_nodes-1) {			
		int prev_lim = kdt->len_nodes - 1;
		kdt->len_nodes = (index*2 + 2) + 1; // index of the right child
		kdt->leaves = realloc(kdt->leaves, kdt->len_nodes/8 + 1);			
		if(!kdt->leaves) {
			printf("Realloc kdt->emptys failed\n");
			return;
		}
		// set newly allocated nodes to not leaves	
		for(int i=prev_lim; i<kdt->len_nodes/8 + 1; i++) {
			kdt->leaves[i] = 0;
		}

		int num_bytes = kdt->len_nodes * sizeof(double); 			
		*pts_x = realloc(*pts_x, num_bytes);
		*pts_y = realloc(*pts_y, num_bytes);
		*pts_z = realloc(*pts_z, num_bytes);	
		if(!(*pts_x) || !(*pts_y) || !(*pts_z) ) {
			printf("kdtree_build_r() realloc failed for new buffer size %i\n", num_bytes);
			return;
		}
	}
	
	(*pts_x)[index] = 0;
	(*pts_y)[index] = 0;
	(*pts_z)[index] = 0;
	
	switch(axis) {
		case 0:
			(*pts_x)[index] = split;
			break;
		case 1:
			(*pts_y)[index] = split;
			break;
		case 2:
			(*pts_z)[index] = split;
			break;
		default:
			printf("This should never happen.\n");
			break;	
	}
	kdt->num_nodes++;
	
	// left child	
	quicksort(points, l, split_index, (axis+1)%3); // sort in the next axis 
	kdtree_build_gpu_r(points, (axis+1)%3, kdt, l, split_index, 2*index + 1, pts_x, pts_y, pts_z);	

	// right child
	quicksort(points, split_index+1, r, (axis+1)%3); // sort in the next axis 
	kdtree_build_gpu_r(points, (axis+1)%3, kdt, split_index+1, r, 2*index + 2, pts_x, pts_y, pts_z);	

}

kdtree* kdtree_build_gpu(double** points, int num_points) {

	double invalid = getRadians(-1); // longitude value can never be -1
	INVALID_X = EARTH_RADIUS * cos(invalid) * cos(invalid); 
	
	kdtree* kdt = (kdtree*) malloc(sizeof(kdtree));
	kdt->len_nodes = num_points * 2; // the length of the nodes array
	kdt->leaves = (char*) malloc( (kdt->len_nodes/8 + 1) * sizeof(char)); // a bitmap of the leaf nodes ; 8 bits/char 
	if(!kdt->len_nodes || !kdt->leaves) {
		printf("malloc() kdt->len_nodes or kdt->leaves failed\n");
		return NULL;
	}
	memset(kdt->leaves, 0, kdt->len_nodes/8 + 1); // bitmap that indicates which nodes are leaves
	
	int num_bytes = kdt->len_nodes * sizeof(double);	
	double* pts_x = malloc(num_bytes);
	double* pts_y = malloc(num_bytes);
	double* pts_z = malloc(num_bytes);
	if(!pts_x || !pts_y || !pts_z) {
		printf("malloc() pts failed\n");
		return NULL;
	}

	quicksort(points, 0, num_points-1, 0); // sort in the first axis 
	printf("kdtree_build_gpu_r\n");
	kdtree_build_gpu_r(points, 0, kdt, 0, num_points-1, 0, &pts_x, &pts_y, &pts_z);	

	kdt->len_nodes = max_leaf_idx + 1;
	num_bytes = kdt->len_nodes * sizeof(double);	
	kdt->nodes = malloc(num_bytes * 3); // for each dimension x,y,z
	memcpy(kdt->nodes, pts_x, num_bytes);
	memcpy(kdt->nodes + kdt->len_nodes, pts_y, num_bytes);
	memcpy(kdt->nodes + kdt->len_nodes*2, pts_z, num_bytes);

	free(pts_x);
	free(pts_y);
	free(pts_z);

	kdt->leaves = realloc(kdt->leaves, kdt->len_nodes/8 + 1);

//	printf("kdtree in 1D array format\n");
//	for(int i=0; i<kdt->len_nodes; i++) {
//		if( fabs(kdt->nodes[i] - INVALID_X) < 1e-32) {
//			printf("idx=%i) null\n", i);
//			continue;
//		}
//		if(kdt->leaves[i/8] & (0x1 << (i%8))) printf("**leaf** ");
//		printf("idx=%i) (%.2f, %.2f, %.2f)\n", i, kdt->nodes[kdt->len_nodes*0 + i], kdt->nodes[kdt->len_nodes*1 + i], kdt->nodes[kdt->len_nodes*2 + i]);
//	}

	return kdt;
}

/*

	CPU

*/

void kdtree_build_cpu_r(double** points, node* root, int axis, int l, int r) { 

	// leaf node	
	if(r == l || l > r || r-l == 0) {
		root->p = r;
		root->leaf = 1; 
		return;
	}

	// calculate the split value	
	int split_index = l + (r-l)/2; // median index
	if(((r-l)+1)%2 == 0) root->split = (points[split_index][axis] + points[split_index+1][axis])/2; // even number of elements
	else root->split = points[split_index][axis]; // odd number of elments; clear median
	
	root->axis = axis;	
	root->left = NULL; 
	root->right = NULL;
	root->leaf = 0;	

	node* new_node;
	int next_axis = (axis+1) % 3; // 3 dimensions x,y,z
	
	// left child
	new_node = (node*) malloc(sizeof(node));	
	root->left = new_node;
	// sort all the values to the left of split index
	quicksort(points, l, split_index, next_axis); 
	kdtree_build_cpu_r(points, root->left, next_axis, l, split_index);	

	// right child	
	new_node = (node*) malloc(sizeof(node));	
	root->right = new_node;	

	// sort all the values the right of the split index
	quicksort(points, split_index+1, r, next_axis); // sort in the next axis 
	kdtree_build_cpu_r(points, root->right, next_axis, split_index+1, r);		
}

node* kdtree_build_cpu(double**points, int num_points) {

	printf("kdtree_build_cpu\n");	
	node* root = (node*) malloc(sizeof(node));	

	quicksort(points, 0, num_points-1, 0);
	kdtree_build_cpu_r(points, root, 0, 0, num_points-1);
	
	return root;	
}
