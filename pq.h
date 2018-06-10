#ifndef PQ_H
#define PQ_H

typedef struct pqueue {
	int max_size;
	int num_elems;
	int* elems; // kdtree indicies
	double* dists; // distances from kdtree point to query point	
} pqueue;

//void swap(pqueue* q, int a, int b);
void pq_print(pqueue* q);
pqueue* pq_build(int size);
//void heapify_up(pqueue* q);
//void heapify_down(pqueue* q);
int pq_insert(pqueue* q, int idx, double dist); // returns -1 if the heap is full
int pq_extract(); // returns the index in the kd-tree

#endif
