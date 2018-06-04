#ifdef KDTREE_H
#define KDTREE_H

int kdtree_partition(double* points, int num_points, int l, int r, int dim);
void kdtree_quicksort(double* points, int num_points, int l, int r, int dim); 
void kdtree_build_r(double* points, int num_points);
double* kdtree_build(double* points, int num_points);

#endif
