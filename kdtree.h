#ifndef KDTREE_H
#define KDTREE_H

typedef struct kdtree {
	double* x;
	double* y;
	double* z;
	int num_points;
	int num_nodes;
	int num_nodes_limit; // used only for building tree
} kdtree;

int partition(double** points, int l, int r, int dim); 
void quicksort(double** points, int l, int r, int dim); 
void kdtree_print(kdtree* kdt);
void kdtree_build_r(double** points, int axis, kdtree* kdt, int l, int r, int index);
kdtree* kdtree_build(double** points, int num_points);

#endif
