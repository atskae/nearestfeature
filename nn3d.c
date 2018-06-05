#include <stdlib.h>
#include <stdio.h>
#include <time.h> 	
#include <string.h>
#include <float.h> // max double value 
#include <math.h>

#include "jsmn.h" // JSON file parser
#include "kdtree.h"

/*
	Implementation of "GPU-accelerated Nearest Neighbor Search for 3D Registration" Qiu 2009	
	
	Notes;
	- Must convert (lat, longt) to (x, y, z) to implement kdtree
	
*/

#define BUFFER_SIZE 1600000
#define POINTS_BUFFER_SIZE 10000 // the number of points to allocate at once ; arbitrary
#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371 // in km

#define DIM 3 // dimension of cartesian coordinates

double** geopoints = NULL; // saving only for debugging purposes

// returns the time elapsed in seconds 
double elapsed(time_t start, time_t end) {
	double t = ( ((double) (end - start)) / CLOCKS_PER_SEC );	
	printf("%.5f ms elapsed.\n", t*1000);
	return t;
}

double getRadians(double degrees) { return (degrees * PI)/180; }

// calculates the distance between cartesian points
 double distance(double* a, double* b) {
	
	double sum = 0;
	for(int i=0; i<DIM; i++) {
		sum += pow((a[i] - b[i]), 2);	
	}
	
	return sqrt(sum);
}

// calculates the closest distance from point p to the a line
double p_to_line(double* p, double split, int axis) { return fabs(p[axis] - split); }

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

double** parseJSON(char* file, int* num_points) {

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

	// must parse according to geometry type: Point, LineString, Polygon, MultiPoint, MultiLineString, MutliPolygon ; https://en.wikipedia.org/wiki/GeoJSON

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

	double** points = (double**) malloc(POINTS_BUFFER_SIZE * sizeof(double*));
	geopoints = (double**) malloc(POINTS_BUFFER_SIZE * sizeof(double*));
	*num_points = 0; // index into points

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
			
			index++;
 
			// extract the latitude
			getTokenString(tokens, index, jsonString, tokenStr);
			end = &tokenStr[tokens[index].size - 1]; // get pointer to the last digit
			lat = strtod(tokenStr, &end); // save latitude value at index 0

			index++;

			// see if this point has already been seen
			char dup = 0;
			for(int k=0; k<*num_points; k++) {
				if(geopoints[k][0] - lat < 1e-32 && geopoints[k][1] - longt < 1e-32) {
					dup = 1;
					break;
				}	
			}

			if(!dup) {
				int i = *num_points;
				points[i] = malloc(DIM * sizeof(double));
				geopoints[i] = malloc(2 * sizeof(double));

				// save the geo coordinates	
				geopoints[i][0] = lat;
				geopoints[i][1] = longt;

				// convert (lat, long) to (x,y,z)
			
				lat = getRadians(lat);
				longt = getRadians(longt);	
				points[i][0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
				points[i][1] = EARTH_RADIUS * cos(lat) * sin(longt);
				points[i][2] = EARTH_RADIUS * sin(lat);
			
				(*num_points)++;
				if(*num_points % POINTS_BUFFER_SIZE == 0) {
					points = realloc(points, (*num_points + POINTS_BUFFER_SIZE) * sizeof(double*));
					geopoints = realloc(points, (*num_points + POINTS_BUFFER_SIZE) * sizeof(double*));
					if(!points || !geopoints) {
						printf("%i %i realloc() failed for buffer size %lu\n", points == 0, geopoints == 0, (*num_points + POINTS_BUFFER_SIZE) * sizeof(double*));
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
		} // get all the points for this feature ; end
	
	} // for each feature ; end

	free(jsonString);
	free(tokens);

	if((*num_points < 20)) {
		quicksort(points, 0, (*num_points)-1, 0); 
		printf("%i points in dataset; sorted by x-axis\n", *num_points);			
		for(int i=0; i<*num_points; i++) {
			double* p = points[i];
			double* g = geopoints[i];
			printf("(%.12f, %.12f, %.12f) | (%.12f, %.12f)\n", p[0], p[1], p[2], g[1], g[0]);
		}
	}

	return points;
}

//void findBetter(node* n, double* p, int* best, double* dist_best) {
//
////	printf("best %i, current node: ", *best);
////	printNode(n);
////
//	double dist = 0;
//
//	if(!n->left && !n->right) { // is a leaf node	
//		for(int i=0; i<DIM; i++) {
//			dist += pow((p[i] - points[n->p][i]), 2);	
//		}	
//		dist = sqrt(dist);
//
//		if(dist < *dist_best) {
//			*best = n->p;
//			*dist_best = dist;
////			printf("best distance %.12f at point ", dist);
////			printPoint(points[*best]);
//			return;
//
//		}
//	}
//	else { // is a split node
//		dist = fabs(p[n->axis] - n->split); // distance to line
//		
//		if(dist < *dist_best) {
//			// explore the other side of the split line
//			findBetter(n->left, p, best, dist_best); // explore the other half of tree
//			findBetter(n->right, p, best, dist_best); 
//		}	
//	}	
//	
//}
//
//void findNearestPoint(node* n, double* p, int* best, double* dist_best) {
//
//	//printf("findNearestNeighbor()\n");		
//	//printf("axis %i, level %i ", n->axis, n->level);
//	//if(!n->left && !n->right) printPoint(points[n->p]);
//	//else printf("split %.3f\n", n->split);
//
//	// base case: leaf node
//	if(!n->left && !n->right) {
//		//n->colored = 1;
//		double dist = distance(p, points[n->p]);	
//		if(dist < *dist_best) {
//			*dist_best = dist;
//			*best = n->p;
//		}
//		return;
//	}
//	
//	if(p[n->axis] > n->split) {
//		findNearestPoint(n->right, p, best, dist_best);	
//		findBetter(n->left, p, best, dist_best); 
//	} else {
//		findNearestPoint(n->left, p, best, dist_best);
//		findBetter(n->right, p, best, dist_best);
//	}
//}

int main(int argc, char* argv[]) {

	time_t start, end;	
	if(argc < 4) {
		printf("./nf <filename> <latitude> <longitude> <max distance>\n");
		return 1;
	}
	
	char* file = argv[1]; // file name
	int num_points;
	double** points = parseJSON(file, &num_points);

	start = clock();
	kdtree* kdt = kdtree_build(points, num_points); 
	end = clock();
	printf("kdtree build in ");
	elapsed(start, end);

	double lat, longt;
	lat = atof(argv[2]); // lat
	longt = atof(argv[3]); // longt
	
	lat = getRadians(lat);
	longt = getRadians(longt);

	double* a = malloc(DIM * sizeof(double));
	a[0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
	a[1] = EARTH_RADIUS * cos(lat) * sin(longt); // y
	a[2] = EARTH_RADIUS * sin(lat); // z
	
	double range = DBL_MAX;
	if(argv[4]) range = atof(argv[4]);
	double kdtree_time, bf_time;

	/* running the nearest neightbor algorithm */

	// brute force
	start = clock();
	double min = DBL_MAX;
	double dist;
	int best;
	for(int i=0; i<kdt->num_points; i++) {
		dist = distance(a, points[i]);	
		if(dist < min) {
			min = dist; 
			best = i;
		}
	}
	end = clock();
	printf("brute force result %i: ", best);
	elapsed(start, end);

	return 0;		
}
