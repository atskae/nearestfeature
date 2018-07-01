#include "pq.cuh"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//typedef struct pqueue {
//	int max_size;
//	int num_elems;
//	int* elems; // kdtree indicies
//	double* dists; // distances from kdtree point to query point	
//} pqueue;

// has_left_child(pqueue* q, int index) { return index*2 + 1 < q->num_elems }
//

//__device__ void swap(pqueue* q, int a, int b) {
//	
//	// save the current values at index a
//	int idx = q->elems[a];
//	int dist = q->dists[a];
//	// replace
//	q->elems[a] = q->elems[b];
//	q->dists[a] = q->dists[b];
//
//	q->elems[b] = idx;
//	q->dists[b] = dist;
//}

//
//__device__ void pq_print(pqueue* q) {
//
//	// [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
//	// new_line: 1, 3, 7
//	printf("--priority queue: num_elems=%i, max_size=%i\n", q->num_elems, q->max_size);
//	int new_line = 1;
//	int level = 0;
//	printf("level 0\n");
//	for(int i=0; i<q->num_elems; i++) {
//		if(i == new_line) {
//			level++;
//			printf("level %i\n", level);
//			new_line = (new_line * 2) + 1;		
//		}
//		printf("%i kdt_idx=%i, dist=%.12f\n", i, q->elems[i], q->dists[i]);
//	}
//	printf("--\n");
//}
//
//__device__ pqueue* pq_build(int size) { // should this be the number of registers/thread?
////	printf("pq_build()\n");
//	pqueue* q = (pqueue*) malloc(sizeof(pqueue));
//	q->max_size = size;
//	q->num_elems = 0;
//	q->elems = (int*) malloc(sizeof(int) * q->max_size);
//	q->dists = (double*) malloc(sizeof(double) * q->max_size);
//	
//	memset(q->elems, 0, size * sizeof(int));
//	memset(q->dists, 0, size * sizeof(double));
//
//	return q;
//}

//
//__device__ void heapify_up(pqueue* q) {	
//	
////	printf("heapify_up()\n");
//
//	// start at the last element in the heap
//	int idx = q->num_elems - 1;
//
////	printf("idx=%i, current dist %.12f, parent's distance %.12f (at idx=%i)\n", idx, q->dists[idx], q->dists[(idx-1)/2], (idx-1)/2);
//	while(idx > 0 && q->dists[idx] < q->dists[(idx - 1)/2]) { // check if the current node is smaller than its parent; if so, the current node is out of place!
////		printf("swapping %i and %i\n",idx, (idx-1)/2 );
//		swap(q, idx, (idx-1)/2); // swap the parent and child
//		idx = (idx-1)/2; // go up to the next element (the parent of the current node)
//	}	
//}

//
//__device__ void heapify_down(pqueue* q) {	
//
//	// start at the first element of the heap
//	int idx = 0;
//	
//	while((idx*2+1) < q->num_elems - 1) { // while this node has a left child ; no need to check for the right child since the heap is always filled from the left
//		// get the index of the smaller of the two children
//		int smaller_child = idx*2 + 1; // default to the left child
//		if(q->dists[idx*2+2] < q->dists[idx*2 + 1]) smaller_child = idx*2 + 2; // if right child is actually smaller than left child, set smaller to right
//
//		if(q->dists[idx] < q->dists[smaller_child]) break; // if this node is already smaller than its children, then we are done
//
//		// swap the current node with the smaller child
//		swap(q, idx, smaller_child);
//		
//		// move down the heap to the smaller_child
//		idx = smaller_child;	
//	}
//}

// the thread must check whether the pq is full BEFORE the function call
__device__ void pq_insert(pqueue* q, int new_idx, double new_dist) {
//	printf("pq_insert()\n");

	// add the newest element to the last spot
	q->elems[q->num_elems] = new_idx;
	q->dists[q->num_elems] = new_dist;
	// increase queue size
	q->num_elems++;

//	heapify_up(q);
	// start at the last element in the heap
	int idx = q->num_elems - 1;
	int parent_idx = (idx-1)/2;

//	printf("idx=%i, current dist %.12f, parent's distance %.12f (at idx=%i)\n", idx, q->dists[idx], q->dists[(idx-1)/2], (idx-1)/2);
	while(idx > 0 && q->dists[idx] < q->dists[parent_idx]) { // check if the current node is smaller than its parent; if so, the current node is out of place!
//		printf("swapping %i and %i\n",idx, (idx-1)/2 );
		
		//swap(q, idx, (idx-1)/2); // swap the parent and child
		// save the current values at index a
		int val = q->elems[idx];
		int dist = q->dists[idx];
		// replace
		q->elems[idx] = q->elems[parent_idx];
		q->dists[idx] = q->dists[parent_idx];
	
		q->elems[parent_idx] = val;
		q->dists[parent_idx] = dist;	
	
		idx = parent_idx; // go up to the next element (the parent of the current node)
	}		
}

// returns the index of the kdtree which contains the shortest distance from the query point so far
__device__ void pq_extract(pqueue* q, int* result, double* dist) {
//	printf("pq_extract()\n");
	if(q->num_elems == 0) {
		*result = -1;
		return;
	}

	// extract the minimum element
	int min = q->elems[0];
	double d = q->dists[0];
	q->num_elems--;

	// move the last element to the front
	q->elems[0] = q->elems[q->num_elems];
	q->dists[0] = q->dists[q->num_elems];	
	
	//	heapify_down(q);	
	// start at the first element of the heap
	int idx = 0;
	
	while((idx*2+1) < q->num_elems - 1) { // while this node has a left child ; no need to check for the right child since the heap is always filled from the left
		// get the index of the smaller of the two children
		int smaller_child = idx*2 + 1; // default to the left child
		if(q->dists[idx*2+2] < q->dists[idx*2 + 1]) smaller_child = idx*2 + 2; // if right child is actually smaller than left child, set smaller to right

		if(q->dists[idx] < q->dists[smaller_child]) break; // if this node is already smaller than its children, then we are done

		// swap the current node with the smaller child
//		swap(q, idx, smaller_child);
		// save the current values at index a
		int val = q->elems[idx];
		int dist = q->dists[idx];
		// replace
		q->elems[idx] = q->elems[smaller_child];
		q->dists[idx] = q->dists[smaller_child];
	
		q->elems[smaller_child] = val;
		q->dists[smaller_child] = dist;
	
		// move down the heap to the smaller_child
		idx = smaller_child;	
	}

	*result = min;
	*dist = d; 

	return;
}

