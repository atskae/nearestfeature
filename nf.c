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

#define BUFFER_SIZE 1600000

int DIM = 0; // dimension of the coordinates
double* ref_pt; // the point where we want to find the closest feature

char* file; // file containing coordinates of the tokenStr
char* jsonString; 
jsmntok_t* tokens; 
int numTokens = 0;
int numPoints = 0;
int numNodes = 0;
int duplicates = 0; // lines overlap points so there will be duplicate points in the dataset ; duplicates are not placed into the tree
int uniquePoints = 0;

typedef struct node {
	int axis;
	double split; // the value that determines which nodes go on the left/right
	double* p; // only contains point if this node has no children
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
	printf("%i numTokens expected.\n", numTokens);
	tokens = (jsmntok_t*) malloc(numTokens * sizeof(jsmntok_t));
	jsmn_init(&parser);

 	numTokens = jsmn_parse(&parser, jsonString, strlen(jsonString), tokens, numTokens); // 1 token contains all the tokenStr
	if(numTokens < 0) {
		printf("Failed to parse JSON file; %i\n", numTokens);
		exit(1);
	} else printf("%i tokens parsed.\n", numTokens);
	
	free(jsonString);
}

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

	int sep = 0;
	sep = partition(points, l, r, dim);
	
	quicksort(points, l, sep-1, dim);
	quicksort(points, sep+1, r, dim);
}

double max(double a, double b) {
	if(a > b) return a;
	else return b;
}

void printPoint(double* p) {
	printf("(");
	for(int i=0; i<DIM; i++) {
		if(i == DIM - 1) printf("%.5f)\n", p[i]);
		else printf("%.5f, ", p[i]);
	}
}

void printTree(node* queue[], int* head, int* tail, int count) {
	
	if(count == numNodes) return;

	// take a node out of the queue
	node* n = queue[ (*tail) ];
	(*tail) = ((*tail) + 1) % numNodes;
	
	printf("axis %i, level %i ", n->axis, n->level);
	if(n->p) {
		printf("point (");
		for(int i=0; i<DIM; i++) {
			if(i == DIM - 1) printf("%.7f)\n", n->p[i]);
			else printf("%.7f, ", n->p[i]);
		}
	}
	else printf("split %.7f\n", n->split);

	// put children in queue
	if(n->left) {
		queue[*head] = n->left;
		(*head)=( *head + 1) % numNodes;
	}
	if(n->right) {
		queue[*head] = n->right;
		(*head)=( *head + 1) % numNodes;
	}

	printTree(queue, head, tail, ++count);	

}

void buildTree_r(node* root, double** points, int axis, int l, int r) { // points = sorted points in x_axis, p_y = sorted points in y_axis
	
//	printf("buildTree_r l=%i, r=%i\n", l, r);

	char dup = 1;
	for(int i=0; i<DIM; i++) {
		if(points[r][i] != points[l][i]) { // don't compare floating point like this!!!!! fix
			dup = 0;
			break;	
		}
	}
	if(r == l || l > r || dup) { // a single point OR there are duplicate points
		if(dup) duplicates++;		
		root->p = points[r];
//		printf("leaf node (");
//		for(int i=0; i<DIM; i++) {
//			if(i == DIM - 1) printf("%.7f)\n", root->p[i]);
//			else printf("%.7f, ", root->p[i]);
//		}

		root->axis = (axis + (DIM-1)) % DIM; // obtain the previous axis
		root->visited = 1;
		uniquePoints++;
		return;
	}
	int split_index = l;

	if(root->visited == 0) {	
		root->split = points[l][axis] + (points[r][axis] - points[l][axis])/2;	
		root->axis = axis;
		root->left = NULL; 
		root->right = NULL;
		root->visited = 1;
	}
	
	// find the split index
	for(int i=l; i<=r; i++) {
		if(points[i][axis] > root->split) {
			split_index = max(l, i-1);
			break;
		}
	}	

//	printf("axis %i, l=%i, r=%i, split %.3f, split_index %i, level %i\n", axis, l, r, root->split, split_index, root->level);
//	for(int i=0; i<numPoints; i++) {	
//		printf("%i: (%f, %f)\n", i, points[i].x, points[i].y);
//	}
//
	node* n_l;
	if(!root->left) {
		n_l = (node*) malloc(sizeof(node));
		n_l->visited = 0;
		n_l->p = NULL;
		n_l->level = root->level + 1;
		root->left = n_l;
		numNodes++;
	} else {
		printf("This should never happen\n");
	}
	
	// sort all the values to the left of split index
	quicksort(points, l, split_index, (axis+1)%DIM); // sort in the next axis 
	buildTree_r(root->left, points, (axis+1)%DIM, l, split_index);	

	node* n_r;
	if(!root->right) {
		n_r = (node*) malloc(sizeof(node));
		n_r->visited = 0;
		n_r->p = NULL;
		n_r->level = root->level + 1;
		root->right = n_r;
		numNodes++;
	} else {
		printf("This should also never happen...\n");
	}
	// sort all the values the right of the split index
	quicksort(points, split_index+1, r, (axis+1)%DIM); // sort in the next axis 
	buildTree_r(root->right, points, (axis+1)%DIM, split_index+1, r);	
	
}

node* buildTree() {
	
	printf("buildTree()\n");
	char tokenStr[BUFFER_SIZE];
	int features_index = 0; // the index into the token array

	// obtain the index into the tokens array where the features start	
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
	double** points = (double**) malloc(numPoints*DIM * sizeof(double*));
	int p_index = 0; // index into points
	printf("%d features\n", tokens[features_index].size);
	memset(tokenStr, '\0', BUFFER_SIZE);	

	int pp_f = 0; // points per feature ; 2 for lines, 1 for points
	for(int i=0; i<tokens[features_index].size; i++) { // for each feature
		while(strcmp(tokenStr, "coordinates") != 0) {
			if(tokens[index].type == JSMN_STRING) getTokenString(index, tokenStr);
			index++;
		}
		getTokenString(index-2, tokenStr);
//		printf("type %s\n", tokenStr);
		if(strstr(tokenStr, "Line") != NULL) pp_f = 2; // points per feature
		else pp_f = 1;
	
		index+=3; // goes to the token where the x-coordinate starts
		char* end;
	
		for(int j=0; j<pp_f; j++) { // for each point in the feature
			points[p_index] = malloc(DIM * sizeof(double));
	
			for(int k=0; k<DIM; k++) { 
				getTokenString(index+k, tokenStr);
				end = &tokenStr[tokens[index+k].size - 1]; // get a pointer to the last digit
				points[p_index][k] = strtod(tokenStr, &end); // convert to a double
			}
			index+=3; // move three tokens ahead to get to next point	
			p_index++;
		}	
	}

	quicksort(points, 0, numPoints-1, 0);
//	printf("%i points in dataset\n", numPoints);	
//	printf("sorted:\n");
//	for(int i=0; i<numPoints; i++) {
//		double* p = points[i];
//		printf("(");
//		for(int j=0; j<DIM; j++) {
//			if(j == DIM - 1) printf("%.7f)\n", p[j]); 
//			else printf("%.7f, ", p[j]);
//		}
//	}
//
	node* root = (node*) malloc(sizeof(node));	
	root->p = NULL;
	root->visited = 0;
	root->level = 0;
	numNodes++;

	printf("Build tree! los gehts!\n");
	buildTree_r(root, points, 0, 0, numPoints-1);
	printf("Build done. %i points, %i nodes, %i unique, %i duplicates\n", numPoints, numNodes, uniquePoints, duplicates);
	node* queue[numNodes];
	queue[0] = root;	
	int head = 1;
	int tail = 0;
	printTree(queue, &head, &tail, 0);
	
	return root;	
}

void findBetter(node* n, double* p, double** best) {

//	printf("findBetter() --- \n");	
	// calcualate distance to the current node, which can either be a point or a split
	char is_pt = 0; // if this node is a leaf node, we can simply return it if is the best
	double sum = 0;
	double dist = 0;
//	printf("current node: ");
	if(n->p) {
		is_pt = 1;
		for(int i=0; i<DIM; i++) {
			sum += pow((p[i] - n->p[i]), 2);
		}
		dist = sqrt(sum);
//		printPoint(n->p);
	} else {
		// calculate distance to split line
		sum = pow((p[n->axis] - n->split), 2);
		dist = sqrt(sum);
//		printf("split %.3f\n", n->split);	
	}
//	printf("dist from point p to current node = %.3f\n", dist);	

	// calculate the distance of best
	sum = 0;
	for(int i=0; i<DIM; i++) {
		sum += pow((p[i] - (*best)[i]), 2);
	}
	double dist_best = sqrt(sum);
	//printf("current best; ");
//	printPoint(*best); 
//	printf("dist to best = %.3f\n", dist_best);

	if(dist < dist_best) {
		if(is_pt) {
			*best =  n->p; // return the point
			return;
		}
		// explore the other side of the split line
		if(p[n->axis] > n->split) findBetter(n->left, p, best); // explore the other half of tree
		else findBetter(n->right, p, best); 
	}
}

double* findNearestPoint(node* n, double* p) {

//	printf("findNearestNeighbor()\n");		
//	printf("axis %i, level %i ", n->axis, n->level);
//	if(n->p) printPoint(n->p);
//	else printf("split %.3f\n", n->split);

	// base case: leaf node
	if(n->p) return n->p;
	
	int axis = n->axis;
	double* best;
	if(p[n->axis] > n->split) best = findNearestPoint(n->right, p);	
	else best = findNearestPoint(n->left, p);

	findBetter(n, p, &best); // looks for a closer point, if not found, best remains the same; sets best
	
	return best;
}

int main(int argc, char* argv[]) {

	time_t start, end;	
	if(argc < 3) {
		printf("./nf <filename> <dimension n> <d1> ... <dn>\n");
		return 1;
	}
	
	file = argv[1]; // file name	
	printf("filename: %s\n", file);
	DIM = atoi(argv[2]);
	if(DIM > 3 || DIM < 1) {
		printf("Dimension must be between 1-3.\n");
		return 1;
	}
	if(argc != 3 + DIM) {
		printf("Did not provide a %i dimension coordinate.\n", DIM);
		return 1;
	}
	
	ref_pt = malloc(DIM * sizeof(double));
	for(int i=0; i<DIM; i++) {
		ref_pt[i] = (double)rand()/(double)(RAND_MAX/100);
	}
//	// scanf messes this part up.... figure out the issue
//	printf("Finding the nearest neighbor from point (");
//	for(int i=0; i<DIM; i++) {
//	//	p[i] = strtod(argv[3+i], strlen(argv[3+i]));
//		if(i == DIM-1) printf("%f)\n", p[i]);
//		else printf("%f, ", p[i]);
//	}
 
	parseJSON(file);
	printf("Ready to build tree\n");
	start = clock();
	node* kdtree = buildTree();
	end = clock();
	printf("tree build in ");
	elapsed(start, end);

	double a[3];
	for(int i=0; i<DIM; i++) {
		a[i] = atof(argv[3+i]);
	}
	printf("Starting nearest point search\n");
	start = clock();	
	double* b = findNearestPoint(kdtree, a);
	end = clock();
	printf("query in ");
	elapsed(start, end);
	printf("the closest point from point from ");
	printPoint(a);
	printf("is ");
	printPoint(b);

	return 0;		
}
