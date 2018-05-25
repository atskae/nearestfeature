//#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h> 	
#include <string.h>
#include "jsmn.h" // JSON file parser
#include <math.h>
/*
	Nearest Feature
	Given a dataset of tokenStr and a given point p, find the closest feature from p
*/

#define BUFFER_SIZE 160000 

// Error handling macro
#define CHECK(function) {															\
	const cudaError_t error = function;												\
	if(error != cudaSuccess) {														\
		printf("ERROR in %s:%d\n", __FILE__, __LINE__); 							\
		printf("error code:%d, reason %s\n", error, cudaGetErrorString(error)); 	\
		exit(1);																	\
	}																				\
}

int p[2] = {0, 0}; // find the closest feature from this point
char* file = "rivers-small.geojson"; // file containing coordinates of the tokenStr
char* jsonString; 
jsmntok_t* tokens; 
int numTokens = 0;
int numPoints = 0;
int numNodes = 0;
int duplicates = 0; // lines overlap points so there will be duplicate points in the dataset ; duplicates are not placed into the tree
int uniquePoints = 0;

typedef struct point {
	long double x;
	long double y;
} point;

// stucture of arrays
//typedef struct points {
//	float* x;
//	float* y;
//} points;

typedef enum axis_t {x_axis, y_axis} axis_t;

typedef struct node {
	axis_t axis; // either x_axis or y_axis
	long double split; // the value that determines which nodes go on the left/right
	point* p; // only contains point if this node has no children
	char visited;
	int level;
	struct node* left;
	struct node* right;
} node;

// returns the time elapsed in seconds 
void elapsed(time_t start, time_t end) {
	double t = ( ((double) (end - start)) / CLOCKS_PER_SEC );	
	printf("%.5f ms elapsed.\n", t*1000);
}

void getTokenString(int index, char* ptr) {
	jsmntok_t t = tokens[index];
	int strPtr = 0;

	for(int i=t.start; i<t.end; i++) {
		ptr[strPtr] = jsonString[i];
		strPtr++;
	}
	ptr[strPtr] = '\0';
}

void printToken(int index) {
	jsmntok_t t = tokens[index];
	// add a nicer switch statement to print the type of token	
	printf("-----\nTOKEN %i: type %i, start %i, end %i, size %i, ", index, t.type, t.start, t.end, t.size);
	for(int i=t.start; i<t.end; i++) {
		printf("%c", jsonString[i]);
	}
	printf("\n-----\n");
}

void parseJSON(char* file) {

	printf("parseJSON() will parse %s\n", file);
	jsonString = (char*) malloc(BUFFER_SIZE);
	// obtain the JSON string
    FILE *f;
    char c;
    int index = 0;
    f = fopen(file, "rt");
	if(!f) {
		printf("Failed to open file %s\n", file);
		exit(1);
	}
    while((c = fgetc(f)) != EOF){
        jsonString[index] = c;
        index++;
		if(index == BUFFER_SIZE) {
			jsonString = realloc(jsonString, index+BUFFER_SIZE); // allocates more memory
		}    
	}
	fclose(f);

    jsonString[index] = '\0';
	index++; // index saves the size of file (in bytes) for now
	jsonString = realloc(jsonString, index); // adjusts buffer size
	printf("%i bytes\n", index);
	
	jsmn_parser parser;
	numTokens = jsmn_parse(&parser, jsonString, strlen(jsonString), NULL, 0); // 1 token contains all the tokenStr
	printf("%i tokens expected\n", numTokens);
	tokens = (jsmntok_t*) malloc(numTokens * sizeof(jsmntok_t));
	jsmn_init(&parser);

 	numTokens = jsmn_parse(&parser, jsonString, strlen(jsonString), tokens, numTokens); // 1 token contains all the tokenStr
	if(numTokens < 0) {
		printf("Failed to parse JSON file; %i\n", numTokens);
		exit(1);
	} else printf("%i tokens parsed.\n", numTokens);
	
	free(jsonString);
}

int partition_x(point* points, int l, int r) {
	
	int swap = l+1;
	for(int i=l+1; i<=r; i++) { // iterate through all the elements except the pivot value (which is always the left-most element)
	
		if(points[i].x < points[l].x) { // compare the current element with the pivot ; if element is in the wrong spot, switch places 
				point p = points[swap];
				points[swap] = points[i];
				points[i] = p;
				swap++;
		}		

	}

	// swap the pivot
	point p = points[l]; // pivot value
	points[l] = points[swap-1];
	points[swap-1] = p;

	return swap-1; // the partition point
}

int partition_y(point* points, int l, int r) {
	
	int swap = l+1;
	for(int i=l+1; i<=r; i++) { // iterate through all the elements except the pivot value (which is always the left-most element)
	
		if(points[i].y < points[l].y) { // compare the current element with the pivot ; if element is in the wrong spot, switch places 
			point p = points[swap];
			points[swap] = points[i];
			points[i] = p;
			swap++;
		}		
	}

	// swap the pivot
	point p = points[l]; // pivot value
	points[l] = points[swap-1];
	points[swap-1] = p;

	return swap-1; // the partition point
}

void quicksort(point* points, int l, int r, char dim) { 

	if(l > r) return;

	int sep = 0;

	if(dim == 'x') sep = partition_x(points, l, r);
	else sep = partition_y(points, l, r);

	quicksort(points, l, sep-1, dim);
	quicksort(points, sep+1, r, dim);
}

double max(double a, double b) {
	if(a > b) return a;
	else return b;
}

//typedef struct node {
//	axis_t axis; // either x_axis or y_axis
//	double split; // the value that determines which nodes go on the left/right
//	point* p; // only contains point if this node has no children
//	char visited;
//	struct node* left;
//	struct node* right;
//} node;

void printTree(node* queue[], int* head, int* tail, int count) {
	
	if(count == numNodes) return;

	// take a node out of the queue
	node* n = queue[ (*tail) ];
	(*tail) = ((*tail) + 1)%numNodes;
	
	printf("axis %i, level %i ", n->axis, n->level);
	if(n->p) printf("point (%.3Lf,%.3Lf)\n", n->p->x, n->p->y);
	else printf("split %.3Lf\n", n->split);

	// put children in queue
	if(n->left) {
		queue[*head] = n->left;
		(*head)=( *head + 1)%numNodes;
	}
	if(n->right) {
		queue[*head] = n->right;
		(*head)=( *head + 1)%numNodes;
	}

	printTree(queue, head, tail, ++count);	

}

void buildTree_r(node* root, point* p_x, axis_t axis, int l, int r) { // p_x = sorted points in x_axis, p_y = sorted points in y_axis

	char dup = p_x[l].x == p_x[r].x && p_x[l].y == p_x[r].y && l != r; // checks if the points are duplicates
	if(r == l || l > r || dup) { // a single point OR there are duplicate points
		if(dup){
			duplicates++;
		//	printf("duplicate ");
		}
	//	printf("leaf node (%.3Lf, %.3Lf)\n", p_x[r].x, p_x[r].y);
		root->p = &p_x[r];
		root->axis = (++axis)%2; // obtain the previous axis
		root->visited = 1;
		uniquePoints++;

		return;
	}
	int split_index = l;

	if(root->visited == 0) {	
		//float split value;
		long double split;
		if(axis == x_axis) {
			split = p_x[l].x + (p_x[r].x - p_x[l].x)/2;
	//		printf("x-axis calculating split: %.3Lf + (%.3Lf - %.3Lf)/2\n", p_x[l].x, p_x[r].x, p_x[l].x);		
		} else {
			split = p_x[l].y + (p_x[r].y - p_x[l].y)/2;
	//		printf("y-axis calculating split: %.3Lf + (%.3Lf - %.3Lf)/2\n", p_x[l].y, p_x[r].y, p_x[l].y);	
		}
		// either the range contains all the same nodes....? Can be merged into one node I think..
	
		root->split = split;
		root->axis = axis;
		root->p = NULL;
		root->left = NULL; // point
		root->right = NULL;
		root->visited = 1;
	}
	
	// find the split index
	if(axis == x_axis) {
		for(int i=l; i<=r; i++) {
			if(p_x[i].x > root->split) {
				split_index = max(l, i-1);
				break;
			}
		}	
	} else {
		for(int i=l; i<=r; i++) {
			if(p_x[i].y > root->split) {
				split_index = max(l, i-1);
				break;
			}		
		}
	}
//	printf("axis %i, l=%i, r=%i, split %.3Lf, split_index %i, level %i\n", axis, l, r, root->split, split_index, root->level);
//	for(int i=0; i<numPoints; i++) {	
//		printf("%i: (%f, %f)\n", i, p_x[i].x, p_x[i].y);
//	}
//
	node* n_l;
	if(!root->left) {
		n_l = (node*) malloc(sizeof(node));
		n_l->visited = 0;
		root->left = n_l;
		root->left->level = root->level + 1;
		numNodes++;
	} else {
		printf("This should never happen\n");
	}	
	// sort all the values to the left of split index
	if(axis == x_axis) quicksort(p_x, l, split_index, 'y'); // sort in the opposite direction
	else quicksort(p_x, l, split_index, 'x');
//	printf("axis %i l=%i, r=%i split %.3Lf, visiting empty left node next\n", root->axis, l, r, root->split);
	buildTree_r(root->left, p_x, (axis+1)%2, l, split_index);	

	node* n_r;
	if(!root->right) {
		n_r = (node*) malloc(sizeof(node));
		n_r->visited = 0;
		root->right = n_r;
		root->right->level = root->level + 1;
		numNodes++;
	}
	if(axis == x_axis) quicksort(p_x, split_index+1, r, 'y'); // sort in the opposite direction
	else quicksort(p_x, split_index+1, r, 'x');
//	printf("axis %i l=%i, r=%i split %.3Lf, visiting empty right node next\n", root->axis, l, r, root->split);
	buildTree_r(root->right, p_x, (axis+1)%2, split_index+1, r);			
}

node* buildTree() {
	
	printf("buildTree()\n");
	char tokenStr[BUFFER_SIZE];
	int features_index = 0; // the index into the token array
	
	for(int i=1; i<numTokens; i++) { 
		if(tokens[i].type == JSMN_ARRAY && (tokens[i-1].end - tokens[i-1].start) == 8 && jsonString[ tokens[i-1].start ] == 'f') {	
			getTokenString(i-1, tokenStr);
			if(strcmp(tokenStr, "features") == 0) {
				features_index = i;
				break;
			}
		}
	}

	int index = features_index;
	numPoints = tokens[features_index].size * 2; // each feature has 2 points
	point* points_x = (point*) malloc(numPoints * sizeof(point));
	printf("%d features\n", tokens[features_index].size);
	
	// structure of arrays
	// sorted by x coordinate
//	points ps_x;
//	ps_x.x = (double*) malloc(numPoints * sizeof(double));
//	ps_x.y = (double*) malloc(numPoints * sizeof(double));
//	// sorted by y coordinate	
//	points ps_y;
//	ps_y.x = (double*) malloc(numPoints * sizeof(double));
//	ps_y.y = (double*) malloc(numPoints * sizeof(double));
//
	int points_index = 0; // index into ps
	memset(tokenStr, '\0', BUFFER_SIZE);	
	
	// to normalize the points
	long double total_x = 0;
	long double total_y = 0;
	for(int i=0; i<tokens[features_index].size; i++) { // for each feature
		while(strcmp(tokenStr, "coordinates") != 0) {
			if(tokens[index].type == JSMN_STRING) getTokenString(index, tokenStr);
			index++;
		}
		index+=3; // goes to the token where the x-coordinate starts
		point p;

		for(int j=0; j<2; j++) { // for each point; 2 points per line
			getTokenString(index, tokenStr);
			char* end = &tokenStr[tokens[index].size - 1]; // get a pointer to the last digit?
			p.x = strtold(tokenStr, &end);
			total_x += p.x;
//			printf("str %s, p.x obtained:%.3Lf\n", tokenStr, p.x);
//			printf("p.x string %s\n", tokenStr);	
			getTokenString(index+1, tokenStr);
			end = &tokenStr[tokens[index+1].size - 1]; // get a pointer to the last digit?
			p.y = strtold(tokenStr, &end);
			total_y += p.y;
//			printf("p.y string %s\n", tokenStr);							
			points_x[i*2 + j] = p; 	
			index+=3; // move three tokens ahead to get to next point
			
		//	points_index++;
		}		
	}

	printf("%i points in dataset\n", numPoints);	
	
	point* points_y = (point*) malloc(numPoints * sizeof(point));
	points_y = memcpy(points_y, points_x, sizeof(point)*numPoints);	
	if(!points_y) {
		printf("points_y malloc failed\n");
		exit(1);
	}

	quicksort(points_x, 0, numPoints-1, 'x');
	quicksort(points_y, 0, numPoints-1, 'y');	

//	memcpy(ps_y.x, ps_x.x, numPoints*sizeof(float));
//	memcpy(ps_y.y, ps_x.y, numPoints*sizeof(float));
//
//	quicksort_p(ps_x, 0, numPoints, 'x');
//	quicksort_p(ps_y, 0, numPoints, 'y');
	
//	typedef struct node {
//		axis_t axis; // either x_axis or y_axis
//		float split; // the value that determines which nodes go on the left/right
//		point p*; // only contains point if this node has no children
//		struct node* left;
//		struct node* right;
//	} node;
	
	node* root = (node*) malloc(sizeof(node));	
	root->visited = 0;
	root->level = 0;
	numNodes++;

	double x_mean = total_x/numPoints;
	double y_mean = total_y/numPoints;
//	printf("normalized with mean (%f,%f)\n ", x_mean, y_mean);
//	for(int i=0; i<numPoints; i++) {
//		points_x[i].x -= x_mean;
//		points_x[i].y -= y_mean;
//		printf("%i (%.3Lf, %.3Lf) (%lu bytes)\n", i, points_x[i].x, points_x[i].y, sizeof(points_x[i].x));
//	}

	printf("Build tree! los gehts!\n");
	buildTree_r(root, points_x, x_axis, 0, numPoints-1);
	printf("Build done. %i nodes, %i unique, %i duplicates\n", numNodes, uniquePoints, duplicates);
	node* queue[numNodes];
	queue[0] = root;	
	int head = 1;
	int tail = 0;
//	printTree(queue, &head, &tail, 0);
	
	return root;	
}
//
//point* findBetter(node* n, point* best, point* p) {
//
//	// base case
//	if(n->p) {
//	
//	}
//	
//	// else, traverse the tree
//	if(axis == x_axis) {
//		if(p.x > p_node->x) findNearestPoint(n->right, p);
//		else best = findNearestPoint(n->left, p);
//		
//		// check if there is any point better than best
//		best = findBetter(n->right, best, p);		
//
//	} else {
//		if(p.y > p_node->y) best = findNearestPoint(n->right, p);
//		else best = findNearestPoint(n->left, p);
//	}
//
//}
//
//point* findBetter(node* n, point* best, point* p) {
//	//base case: leaf node
//	
//	// check if there is any point better than best
//	if( fabsl(p->x - p->x) > fabsl(best->x - p->x)) {
//		findBetter();
//	}	
//}
//
//point* findNearestPoint(node* n, point p) {
//		
//	printf("axis %i, level %i ", n->axis, n->level);
//	if(n->p) printf("point (%.3Lf,%.3Lf)\n", n->p->x, n->p->y);
//	else printf("split %.3Lf\n", n->split);
//
//	// base case: leaf node
//	if(n->p) return n->p;
//	
//	axis_t axis = n->axis;
//	point* p_node = n->p;	
//
//	point* best;
//	if(axis == x_axis) {
//		if(p.x > p_node->x) best = findNearestPoint(n->right, p);
//		else best = findNearestPoint(n->left, p);
//	
//		// checks if the opposite side has a closer point	
//		best = findBetter(n->right, best, p);	
//	
//	} else {
//		if(p.y > p_node->y) best = findNearestPoint(n->right, p);
//		else best = findNearestPoint(n->left, p);
//	}		
//
//	return best;
//}

int main(int argc, char* argv[]) {

	time_t start, end;
		
	if(argv[1]) file = argv[1];	
	if(argc == 4) {
		p[0] = atoi(argv[2]); // x-coordinate
		p[1] = atoi(argv[3]); // y coordinate
	}

	printf("Finding the nearest neighbor from point (%i, %i) in dataset %s\n", p[0], p[1], file);
	parseJSON(file);
	start = clock();
	node* kdtree = buildTree();
	end = clock();
	printf("tree build in ");
	elapsed(start, end);

	point a;
	a.x = 31.192;
	a.y = -22.371;
	printf("findNearestPoint() from point (%.3Lf, %.3Lf)\n", a.x, a.y);
//	start = clock();	
//	point*b = findNearestPoint(kdtree, a);
//	end = clock();
//	printf("query in ");
//	elapsed(start, end);
//	printf("the closest point from point (%.3Lf, %.3Lf) is (%.3Lf, %.3Lf)\n", a.x, a.y, b->x, b->y);
	return 0;		
}
