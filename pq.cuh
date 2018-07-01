#ifndef PQ_H
#define PQ_H

typedef struct pqueue {
	int max_size;
	int num_elems;
	int* elems; // kdtree indicies
	double* dists; // distances from kdtree point to query point	
} pqueue;

// these were removed in order for kernels to return void
//void pq_print(pqueue* q);
//pqueue* pq_build(int size);

__device__ void pq_insert(pqueue* q, int new_idx, double new_dist); // returns -1 if the heap is full
__device__ void pq_extract(pqueue* q, int* value, double* dist); // sets the index in the kdtree to value
 
#endif
