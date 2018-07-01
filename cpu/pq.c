#include "pq.h"
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

void swap(pqueue* q, int a, int b) {
	
	// save the current values at index a
	int idx = q->elems[a];
	int dist = q->dists[a];
	// replace
	q->elems[a] = q->elems[b];
	q->dists[a] = q->dists[b];

	q->elems[b] = idx;
	q->dists[b] = dist;
}

void pq_print(pqueue* q) {

	// [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	// new_line: 1, 3, 7
	printf("--priority queue: num_elems=%i, max_size=%i\n", q->num_elems, q->max_size);
	int new_line = 1;
	int level = 0;
	printf("level 0\n");
	for(int i=0; i<q->num_elems; i++) {
		if(i == new_line) {
			level++;
			printf("level %i\n", level);
			new_line = (new_line * 2) + 1;		
		}
		printf("%i kdt_idx=%i, dist=%.12f\n", i, q->elems[i], q->dists[i]);
	}
	printf("--\n");
}

pqueue* pq_build(int size) { // should this be the number of registers/thread?
//	printf("pq_build()\n");
	pqueue* q = malloc(sizeof(pqueue));
	q->max_size = size;
	q->num_elems = 0;
	q->elems = malloc(sizeof(int) * q->max_size);
	q->dists = malloc(sizeof(double) * q->max_size);
	
	memset(q->elems, 0, size * sizeof(int));
	memset(q->dists, 0, size * sizeof(double));

	return q;
}

void heapify_up(pqueue* q) {	
	
//	printf("heapify_up()\n");

	// start at the last element in the heap
	int idx = q->num_elems - 1;

//	printf("idx=%i, current dist %.12f, parent's distance %.12f (at idx=%i)\n", idx, q->dists[idx], q->dists[(idx-1)/2], (idx-1)/2);
	while(idx > 0 && q->dists[idx] < q->dists[(idx - 1)/2]) { // check if the current node is smaller than its parent; if so, the current node is out of place!
//		printf("swapping %i and %i\n",idx, (idx-1)/2 );
		swap(q, idx, (idx-1)/2); // swap the parent and child
		idx = (idx-1)/2; // go up to the next element (the parent of the current node)
	}	
}


void heapify_down(pqueue* q) {	

	// start at the first element of the heap
	int idx = 0;
	
	while((idx*2+1) < q->num_elems - 1) { // while this node has a left child ; no need to check for the right child since the heap is always filled from the left
		// get the index of the smaller of the two children
		int smaller_child = idx*2 + 1; // default to the left child
		if(q->dists[idx*2+2] < q->dists[idx*2 + 1]) smaller_child = idx*2 + 2; // if right child is actually smaller than left child, set smaller to right

		if(q->dists[idx] < q->dists[smaller_child]) break; // if this node is already smaller than its children, then we are done

		// swap the current node with the smaller child
		swap(q, idx, smaller_child);
		
		// move down the heap to the smaller_child
		idx = smaller_child;	
	}
}


int pq_insert(pqueue* q, int idx, double dist) {
//	printf("pq_insert()\n");

	if(q->num_elems == q->max_size) return -1; // check if queue is full

	// add the newest element to the last spot
	q->elems[q->num_elems] = idx;
	q->dists[q->num_elems] = dist;
	// increase queue size
	q->num_elems++;

	heapify_up(q);

	return 0;	
}

// returns the index of the kdtree
int pq_extract(pqueue* q) {
//	printf("pq_extract()\n");
	if(q->num_elems == 0) return -1;

	// extract the minimum element
	int min = q->elems[0];
	q->num_elems--;

	// move the last element to the front
	q->elems[0] = q->elems[q->num_elems];
	q->dists[0] = q->dists[q->num_elems];	
	
	// heapify down
	heapify_down(q);	

	return min;
}

