//#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h> 	
#include <string.h>
#include "jsmn.h" // JSON file parser

/*

	Nearest Feature
	Given a dataset of tokenStr and a given point p, find the closest feature from p

*/

#define BUFFER_SIZE 1600000

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

typedef struct point {
	float x;
	float y;
} point;

typedef enum axis_t {x_axis, y_axis} axis_t;

typedef struct node {
	axis_t axis; // either x_axis or y_axis
	float split_value;
	point p; // only contains point if this node has no children
	struct node* left;
	struct node* right;
} node;

// returns the time elapsed in seconds 
void elapsed(time_t start, time_t end) {
	double t = ( ((double) (end - start)) / CLOCKS_PER_SEC );	
	printf("%.5f s elapsed.\n", t);
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

void buildTree() {
	
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
	for(int i=0; i<tokens[features_index].size; i++) { // for each feature
		while(strcmp(tokenStr, "coordinates") != 0) {
			if(tokens[index].type == JSMN_STRING) getTokenString(index, tokenStr);
			index++;
		}
		index+=3; // goes to the token where the x-coordinate starts
		point p;

		for(int j=0; j<2; j++) { // for each point; 2 points per line
			getTokenString(index, tokenStr);
			p.x = atof(tokenStr);	
			getTokenString(index+1, tokenStr);
			p.y = atof(tokenStr);
			
			points_x[i*2 + j] = p; 	
			index+=3; // move three tokens ahead to get to next point
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
	
	
}

int main(int argc, char* argv[]) {
		
	if(argv[1]) file = argv[1];	
	if(argc == 4) {
		p[0] = atoi(argv[2]); // x-coordinate
		p[1] = atoi(argv[3]); // y coordinate
	}

	printf("Finding the nearest neighbor from point (%i, %i) in dataset %s\n", p[0], p[1], file);
	parseJSON(file);
	buildTree();

	return 0;		
}
