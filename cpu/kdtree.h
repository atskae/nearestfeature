#ifndef KDTREE_H
#define KDTREE_H

#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371 // in km

//double INVALID_X; // to indicate inner nodes that are empty ;

typedef struct kdtree {
	double* x;
	double* y;
	double* z;
	int num_points;
	int num_nodes;
	int array_lim;
	char* emptys; // keeps track of the empty nodes ; this is used to check if the node is a leaf
	int* axes;

	// GPU stuff	
	int len_nodes; // the length of nodes array`
	double* nodes; // 1D array of all split nodes and x,y,z coordinates at the leaves
	char* leaves;
	double invalid_x;
} kdtree;

double getRadians(double degrees);

int partition(double** points, int l, int r, int dim); 
void quicksort(double** points, int l, int r, int dim); 
void print_point(kdtree* kdt, int idx);
void kdtree_print(kdtree* kdt);
void kdtree_build_r(double** points, int axis, kdtree* kdt, int l, int r, int index);
kdtree* kdtree_build(double** points, int num_points);

// GPU stuff
kdtree* kdtree_build_gpu(double** points, int num_points);
void kdtree_build_gpu_r(double** points, int axis, kdtree* kdt, int l, int r, int index, double** pts_x, double** pts_y, double** pts_z);


/*

	CPU stuff

*/

typedef struct node {
	int p; // index into the points array
	int axis;
	double split; // the value that determines which nodes go on the left/right
	char leaf;
	
	struct node* left;
	struct node* right;
} node;

void kdtree_build_cpu_r(double** points, node* root, int axis, int l, int r);
node* kdtree_build_cpu(double** points, int num_points);

#endif
