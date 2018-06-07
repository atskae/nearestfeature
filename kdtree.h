#ifndef KDTREE_H
#define KDTREE_H

typedef struct kdtree {
	double* x;
	double* y;
	double* z;
	int num_points;
	int num_nodes;
	int array_lim;
	char* emptys; // keeps track of the empty nodes ; this is used to check if the node is a leaf
} kdtree;

int partition(double** points, int l, int r, int dim); 
void quicksort(double** points, int l, int r, int dim); 
void print_point(kdtree* kdt, int idx);
void kdtree_print(kdtree* kdt);
void kdtree_build_r(double** points, int axis, kdtree* kdt, int l, int r, int index);
kdtree* kdtree_build(double** points, int num_points);

#endif
