#include <stdlib.h>
#include <stdio.h>
#include <time.h> 	
#include <string.h>
#include <float.h> // max double value 
#include <math.h>

#include "jsmn.h" // JSON file parser


/*
	Nearest Feature
	Given a dataset of tokenStr and a given point p, find the closest feature from p
*/

#define BUFFER_SIZE 1600000
#define POINTS_BUFFER_SIZE 100000 // the number of points to allocate at once ; arbitrary
#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371 // in km

#define DIM 3 // dimension of cartesian coordinates

int numPoints = 0;
int numNodes = 0;

double** points = NULL;
double** geopoints = NULL;

typedef struct node {
	int axis; // 0=longitude, 1=latitiude
	double split; // the value that determines which nodes go on the left/right
	int geop; // index into the geopoints array (lat, longt) 
	int p; // index into the points array; only contains point if this node has no children // (x,y,z)
	char visited;
	int level;
	struct node* left;
	struct node* right;
} node;

// returns the time elapsed in seconds 
double  elapsed(time_t start, time_t end) {
	double t = ( ((double) (end - start)) / CLOCKS_PER_SEC );	
	printf("%.5f ms elapsed.\n", t*1000);
	return t;
}

double getRadians(double degrees) {
	return (degrees * PI)/180;
}

// calculates the distance between 2 points (Latitude, Longitude) in kilometers
// source: https://www.movable-type.co.uk/scripts/latlong.html
double haversine(double* a, double* b) {
	
	double a_lat_r = getRadians(a[1]);
	double b_lat_r = getRadians(b[1]);
	double delta_lat_r = a_lat_r - b_lat_r;

	double a_long_r = getRadians(a[0]);
	double b_long_r = getRadians(b[0]);
	double delta_long_r = a_long_r - b_long_r;

	double A = pow(sin(delta_lat_r/2), 2) + cos(a_lat_r)*cos(b_lat_r)*pow(sin(delta_long_r/2),2);
	double C = 2 * atan2(sqrt(A), sqrt(1-A));
	return EARTH_RADIUS * C;
	
}

// calculates the distance between cartesian points
 double distance(double* a, double* b) {
	
	double sum = 0;
	for(int i=0; i<DIM; i++) {
		sum += pow((a[i] - b[i]), 2);	
	}
	
	return sqrt(sum);
}

// calculates the closest distance from point p to the split line, which is a line going across Earth's sphere
 double p_to_line(double* p, double split, int axis) {
	// generate the two points that represents the split line
//	return fabs(p[axis]*split) / sqrt(pow(split, 2));
	return fabs(p[axis] - split);
}

void getTokenString(jsmntok_t* tokens, int token_idx,  char* json_str, char* dest) {
	
	jsmntok_t t = tokens[token_idx];
	int ptr = 0;

	for(int i=t.start; i<t.end; i++) {
		dest[ptr] = json_str[i];
		ptr++;
	}
	dest[ptr] = '\0';
}

void printToken(jsmntok_t* tokens, int index, char* jsonString) {
	jsmntok_t t = tokens[index];
	// add a nicer switch statement to print the type of token	
	printf("-----\nTOKEN %i: type %i, start %i, end %i, size %i, ", index, t.type, t.start, t.end, t.size);
	for(int i=t.start; i<t.end; i++) {
		printf("%c", jsonString[i]);
	}
	printf("\n-----\n");
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


void parseJSON(char* file) {

    FILE *f;
    char c;
    int index = 0;
	char* jsonString = (char*) malloc(BUFFER_SIZE);
    f = fopen(file, "rt");
	if(!f) {
		printf("Failed to open file %s\n", file);
		exit(1);
	}
    while((c = fgetc(f)) != EOF){
        jsonString[index] = c;
        index++;
		if(index % BUFFER_SIZE == 0) {
			jsonString = realloc(jsonString, index+BUFFER_SIZE); // allocates more memory
		}    
	}
	fclose(f);
    jsonString[index] = '\0';
	index++; 
	jsonString = realloc(jsonString, index); // adjusts buffer size
	
	jsmn_parser parser;
	jsmn_init(&parser);
	int numTokens = jsmn_parse(&parser, jsonString, strlen(jsonString), NULL, 0); 

	jsmntok_t* tokens = (jsmntok_t*) malloc(numTokens * sizeof(jsmntok_t));
	jsmn_init(&parser); // restart parser	
 	numTokens = jsmn_parse(&parser, jsonString, strlen(jsonString), tokens, numTokens); 
	if(numTokens < 0) {
		printf("Failed to parse JSON file: ");
		switch(numTokens) {
			case JSMN_ERROR_INVAL:
				printf("JSON string is corrupted.\n");
				break;
			case JSMN_ERROR_NOMEM:
				printf("not enough tokens.\n");
				break;
			case JSMN_ERROR_PART:
				printf("JSON string is too short. Expecting more JSON data.\n");
				break;	
			default: printf("error unknown.\n");
				break;
		}
		exit(1);
	}
	// else printf("%i tokens parsed.\n", numTokens);
	
	
	char tokenStr[BUFFER_SIZE];
	int numFeatures = 0;

	// must parse according to geometry type: Point, LineString, Polygon, MultiPoint, MultiLineString, MutliPolygon 
	//https://en.wikipedia.org/wiki/GeoJSON

	// obtain the index into the tokens array where the features start	
	for(int i=1; i<numTokens; i++) { 
		if(tokens[i].type == JSMN_ARRAY && (tokens[i-1].end - tokens[i-1].start) == 8 && jsonString[ tokens[i-1].start ] == 'f') {	
			getTokenString(tokens, i-1, jsonString, tokenStr); 
			if(strcmp(tokenStr, "features") == 0) {
				numFeatures = tokens[i].size;
				index = i;
				break;
			}
		}
	}

	points = (double**) malloc(POINTS_BUFFER_SIZE * sizeof(double*));
	geopoints = (double**) malloc(POINTS_BUFFER_SIZE * sizeof(double*));
	int numPoints = 0; // index into points
	printf("Getting %i features\n", numFeatures);

	for(int i=0; i<numFeatures; i++) {
		while(strcmp(tokenStr, "coordinates") != 0) {
			index++;
			if(tokens[index].type == JSMN_STRING) getTokenString(tokens, index, jsonString, tokenStr);
		}
		
		while(tokens[index].type != JSMN_PRIMITIVE) {
			index++;
		}

		char* end;
		double lat, longt;

		char goto_next_feat = 0;
		while(!goto_next_feat && index != numTokens) { // get all points in this feature	
		
			// extracts the longtitude
			getTokenString(tokens, index, jsonString, tokenStr);
			end = &tokenStr[tokens[index].size - 1]; // get a pointer to the last digit
			longt = strtod(tokenStr, &end); 
//			printf("|longt = %s| is primitive: %i\n", tokenStr, tokens[index].type == JSMN_PRIMITIVE);			
			
			index++;
 
			// extract the latitude
			getTokenString(tokens, index, jsonString, tokenStr);
			end = &tokenStr[tokens[index].size - 1]; // get pointer to the last digit
			lat = strtod(tokenStr, &end); // save latitude value at index 0
//			printf("|lat= %s| is primitive: %i\n", tokenStr, tokens[index].type == JSMN_PRIMITIVE);			

			index++;

//			for(int m=index; m<index+10; m++) {
//				getTokenString(tokens, m, jsonString, tokenStr);
//				printf("%s is primitive? %i\n", tokenStr, tokens[m].type == JSMN_PRIMITIVE);
//			}
	
			// see if this point has already been seen
			char dup = 0;
			for(int k=0; k<numPoints; k++) {
				if(geopoints[k][0] - lat < 1e-32 && geopoints[k][1] - longt < 1e-32) {
	//				printf("%.12f, %.12f is already in the list\n", longt, lat);
					dup = 1;
					break;
				}	
			}

			if(!dup) {
				points[numPoints] = malloc(DIM * sizeof(double));
				geopoints[numPoints] = malloc(2 * sizeof(double));

				// save the geo coordinates	
				geopoints[numPoints][0] = lat;
				geopoints[numPoints][1] = longt;
			
				lat = getRadians(lat);
				longt = getRadians(longt);	
				// convert (lat, long) to (x,y,z)
				points[numPoints][0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
				points[numPoints][1] = EARTH_RADIUS * cos(lat) * sin(longt);
				points[numPoints][2] = EARTH_RADIUS * sin(lat);
			
				numPoints++;
				if(numPoints % POINTS_BUFFER_SIZE == 0) {
					printf("realloc! to hold %i points\n", numPoints + POINTS_BUFFER_SIZE);
					points = realloc(points, (numPoints + POINTS_BUFFER_SIZE) * sizeof(double*));
					geopoints = realloc(points, (numPoints + POINTS_BUFFER_SIZE) * sizeof(double*));
					if(!points || !geopoints) {
						printf("%i %i realloc() failed for buffer size %lu\n", points == 0, geopoints == 0, (numPoints + POINTS_BUFFER_SIZE) * sizeof(double*));
						exit(1);
					}
				}
			}

			if(index >= numTokens) break;

			// move to next point
			while(tokens[index].type != JSMN_PRIMITIVE) {
				index++;
				if(index >= numTokens) break;
				getTokenString(tokens, index, jsonString, tokenStr);
				if(tokens[index].type == JSMN_STRING && strcmp(tokenStr, "Feature") == 0) {
					goto_next_feat = 1;	
					break;
				}
			}
		}
		
		while(tokens[index].type != JSMN_STRING) {
			index++;
			getTokenString(tokens, index, jsonString, tokenStr);
		}	

	}

	free(jsonString);
	free(tokens);

	quicksort(points, 0, numPoints-1, 0);
  	
	printf("%i points in dataset\n", numPoints);			
	printf("sorted:\n");
	for(int i=0; i<numPoints; i++) {
		double* p = points[i];
		printf("%i (", i);
		for(int j=0; j<DIM; j++) {
			if(j == DIM - 1) printf("%.7f)\n", p[j]); 
			else printf("%.7f, ", p[j]);
		}
	}


}
double max(double a, double b) {
	if(a > b) return a;
	else return b;
}

void printPoint(double* p) {
	printf("(");
	for(int i=0; i<DIM; i++) {
		if(i == DIM - 1) printf("%.13f)\n", p[i]);
		else printf("%.13f, ", p[i]);
	}
}

void printTree(node* queue[], int* head, int* tail, int count) {
	
	if(count == numNodes) return;

	// take a node out of the queue
	node* n = queue[ (*tail) ];
	(*tail) = ((*tail) + 1) % numNodes;
	
	printf("axis %i, level %i ", n->axis, n->level);
	if(!n->left && !n->right) printPoint(points[n->p]);
	else printf("split %.12f\n", n->split);

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

void buildTree_r(double* kdarray, int axis, int l, int r, int root_index) { // points = sorted points in x_axis, p_y = sorted points in y_axis

	// leaf node
	if(r == l || l > r || r-l == 0) {
		//printf("level %i leaf: ", root->level);
		//printPoint(points[r]); 
		kdarray[numPoints*0 + root_index] = points[r][0]; // x
		kdarray[numPoints*1 + root_index] = points[r][1]; // y
		kdarray[numPoints*2 + root_index] = points[r][2]; // z 
		return;
	}
	
	int split_index = l + (r-l)/2; // median index

	if(((r-l)+1)%2 == 0) kdarray[root_index] = (points[split_index][axis] + points[split_index+1][axis])/2; // even number of elements
	else kdarray[root_index] = points[split_index][axis]; // odd number of elments; clear median
	
	// left child	
	quicksort(points, l, split_index, (axis+1)%DIM); // sort in the next axis 
	buildTree_r(kdarray, (axis+1)%DIM, l, split_index, 2*root_index+1);	

	// right child
	quicksort(points, split_index+1, r, (axis+1)%DIM); // sort in the next axis 
	buildTree_r(kdarray, (axis+1)%DIM, split_index+1, r, 2*root_index+2);	
	
}

double* buildTree() {
	

	double* kdarray = (double*) malloc(numPoints * DIM * sizeof(double));	
	numNodes++;

	buildTree_r(kdarray, 0, 0, numPoints-1, 0); //buildTree_r(double* kdarray, int axis, int l, int r, int root_index) 
	printf("Build done. %i points, %i nodes\n", numPoints, numNodes);
//	node* queue[numNodes];
//	queue[0] = root;	
//	int head = 1;
//	int tail = 0;
//	printTree(queue, &head, &tail, 0);
	
	return kdarray;	
}

void printNode(node* n) {
	printf("level %i, ", n->level);
	if(n->p == -1) printf("split %.12f\n", n->split);
	else printPoint(points[n->p]);
}

void findBetter(node* n, double* p, int* best, double* dist_best) {

//	printf("best %i, current node: ", *best);
//	printNode(n);
//
	double dist = 0;

	if(!n->left && !n->right) { // is a leaf node	
		for(int i=0; i<DIM; i++) {
			dist += pow((p[i] - points[n->p][i]), 2);	
		}	
		dist = sqrt(dist);

		if(dist < *dist_best) {
			*best = n->p;
			*dist_best = dist;
//			printf("best distance %.12f at point ", dist);
//			printPoint(points[*best]);
			return;

		}
	}
	else { // is a split node
		dist = fabs(p[n->axis] - n->split); // distance to line
		
		if(dist < *dist_best) {
			// explore the other side of the split line
			findBetter(n->left, p, best, dist_best); // explore the other half of tree
			findBetter(n->right, p, best, dist_best); 
		}	
	}	
	
}

void findNearestPoint(node* n, double* p, int* best, double* dist_best) {

	//printf("findNearestNeighbor()\n");		
	//printf("axis %i, level %i ", n->axis, n->level);
	//if(!n->left && !n->right) printPoint(points[n->p]);
	//else printf("split %.3f\n", n->split);

	// base case: leaf node
	if(!n->left && !n->right) {
		//n->colored = 1;
		double dist = distance(p, points[n->p]);	
		if(dist < *dist_best) {
			*dist_best = dist;
			*best = n->p;
		}
		return;
	}
	
	if(p[n->axis] > n->split) {
		findNearestPoint(n->right, p, best, dist_best);	
		findBetter(n->left, p, best, dist_best); 
	} else {
		findNearestPoint(n->left, p, best, dist_best);
		findBetter(n->right, p, best, dist_best);
	}
}

int main(int argc, char* argv[]) {

	time_t start, end;	
	if(argc < 4) {
		printf("./nf <filename> <latitude> <longitude> <max distance>\n");
		return 1;
	}
	
	char* file = argv[1]; // file name		
	parseJSON(file);
	printf("Parsing done\n");
	
	return 0;

	start = clock();
	double* kdtree = buildTree();
	end = clock();
	printf("%i points, %i nodes: tree build in ", numPoints, numNodes);
	elapsed(start, end);

	// test print
	for(int i=0; i<numPoints*DIM; i++) {
		if(i%numPoints == 0) printf("|");
		printf("%.12f, ", kdtree[i]);
	}
	printf("\n");

	return 0;

	double lat, longt;
	lat = atof(argv[2]); // lat
	longt = atof(argv[3]); // longt
//	printf("(lat=%.12f, longt=%.12f) to cartesian: ", lat, longt);
	
	lat = getRadians(lat);
	longt = getRadians(longt);

	double* a = malloc(DIM * sizeof(double));
	a[0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
	a[1] = EARTH_RADIUS * cos(lat) * sin(longt); // y
	a[2] = EARTH_RADIUS * sin(lat); // z
//	printf("(%.12f, %.12f, %.12f)\n", a[0], a[1], a[2]);
	
	double range = DBL_MAX;
	if(argv[4]) range = atof(argv[4]);
	double kdtree_time, bf_time;

	int best;
	double dist_best = range;
	start = clock();	
	//findNearestPoint(kdtree, a, &best, &dist_best);
	end = clock();

	double dist = distance(a, points[best]);
	printf("---kdtree result:\nindex %i, distance: %.12f\n", best, dist);
//	printf("%.12f, %.12f\n", geopoints[b][1], geopoints[b][0]);
	printPoint(points[best]);
	kdtree_time = elapsed(start, end);

	// brute force check
	double min_dist = DBL_MAX;
	int min_index = 0;
	start = clock();	
	for(int i=0; i<numPoints; i++) {
		dist = distance(points[i], a);
		if(dist < min_dist) {
			min_dist  = dist;
			min_index = i;
		}
	}
	end = clock();

	printf("---brute force result:\nindex %i, distance: %.12f\n", min_index, min_dist);
	//printf("%.12f, %.12f\n", geopoints[min_index][1], geopoints[min_index][0]);
	printPoint(points[min_index]);
	bf_time = elapsed(start, end);

	if(best != min_index) {
		printf("kdtree solution is INCORRECT... :(\n");
		return 0;
	}

	if(kdtree_time < bf_time) printf("YES! kdtree is faster!!!\n");
	else printf("brute force is faster by %.2f percent... how can we optimize?\n", fabs(kdtree_time - bf_time)/((kdtree_time+bf_time)/2) * 100);
//	printf("%.12f, %.12f\n", geopoints[min_index][0], geopoints[min_index][1]);
	return 0;		
}
