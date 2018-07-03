#include <cuda_runtime.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h> 	
#include <string.h>
#include <float.h> // max double value 
#include <math.h>

extern "C" {
	#include "jsmn.h" // JSON file parser
	#include "kdtree.h"
}

#include "pq.cuh"

/*
	Implementation of "GPU-accelerated Nearest Neighbor Search for 3D Registration" Qiu 2009	
	
	Notes;
	- Must convert (lat, longt) to (x, y, z) to implement kdtree
	
*/

/*

	Device

*/
// Error handling macro
#define CHECK(function) {															\
	const cudaError_t error = function;												\
	if(error != cudaSuccess) {														\
		printf("ERROR in %s:%d\n", __FILE__, __LINE__); 							\
		printf("error code:%d, reason %s\n", error, cudaGetErrorString(error)); 	\
		exit(1);																	\
	}																				\
}

__device__ __constant__ int ARRAY_LIM; // size of the kdtree array in each dimension
__device__ __constant__ int LEN_NODES;
__device__ __constant__ double INVALID_X_d; // copy of INVALID_X for device read

/*

	Host

*/
#define BUFFER_SIZE 1600000
#define DIM 3 // dimension of cartesian coordinates
extern double INVALID_X; // from kdtree.h

//#define PI 3.14159265358979323846
//#define EARTH_RADIUS 6371 // in km
//#define PQ_SIZE 1 // priority queue size

// returns the time elapsed in miliseconds 
double elapsed(time_t start, time_t end) {
	double t = ( ((double) (end - start)) / CLOCKS_PER_SEC );	
	printf("%.5f ms elapsed.\n", t*1000);
	return t*1000;
}

//double getRadians(double degrees) { return (degrees * PI)/180; }

// calculates the distance between cartesian points
double distance(double* a, double* b) {
	
	double sum = 0;
	for(int i=0; i<DIM; i++) {
		sum += pow((a[i] - b[i]), 2);	
	}
	
	return sqrt(sum);
}

double distance_by_idx(kdtree* kdt, int idx, double* p) {
		
	double sum = 0;
	sum += pow((kdt->x[idx] - p[0]), 2);	
	sum += pow((kdt->y[idx] - p[1]), 2);	
	sum += pow((kdt->z[idx] - p[2]), 2);	

	return sqrt(sum);

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
			jsonString = (char*) realloc(jsonString, index+BUFFER_SIZE); // allocates more memory
		}    
	}
	fclose(f);
    jsonString[index] = '\0';
	index++; 
	jsonString = (char*) realloc(jsonString, index); // adjusts buffer size
	
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
	int points_buffer_size = 10000; // the number of points to allocate at once ; arbitrary

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

	double** points = (double**) malloc(points_buffer_size * sizeof(double*));
	double** geopoints = (double**) malloc(points_buffer_size * sizeof(double*));
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
				points[i] = (double*) malloc(DIM * sizeof(double));
				geopoints[i] = (double*) malloc(2 * sizeof(double));

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
				if(*num_points == points_buffer_size) {					
					points_buffer_size = points_buffer_size * 2;
					points = (double**) realloc(points, points_buffer_size * sizeof(double*));
					geopoints = (double**) realloc(geopoints, points_buffer_size * sizeof(double*));
					if(!points || !geopoints) {
						printf("%i %i realloc() failed for buffer size %lu\n", points == 0, geopoints == 0, (*num_points + points_buffer_size) * sizeof(double*));
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
	free(geopoints);
	points = (double**) realloc(points, (*num_points) * sizeof(double*));

//	quicksort(points, 0, (*num_points)-1, 0); 	
//	printf("%i points in dataset; sorted by x-axis\n", *num_points);			
//	for(int i=25600; i<25601+10; i++) {
//		double* p = points[i];
//		printf("(%.12f, %.12f, %.12f)\n", p[0], p[1], p[2]);
//	}
//
	return points;
}

void findBetter(kdtree* kdt, int idx, int axis, double* p, int* best, double* dist_best) {

//	printf("findBetter() ");
//	print_point(kdt, idx);	
//	printf("best distance so far: %.12f at idx %i\n---\n", *dist_best, *best);

	double dist = 0;
	if(kdt->emptys[idx*2+1] && kdt->emptys[idx*2+2]) { // check if this is a leaf node
		
		double sum = 0;
		sum += pow((kdt->x[idx] - p[0]), 2);	
		sum += pow((kdt->y[idx] - p[1]), 2);	
		sum += pow((kdt->z[idx] - p[2]), 2);
		dist = sqrt(sum);
		
		if(dist < *dist_best) {
			*best = idx;
			*dist_best = dist;
			return;
		}
	}
	else { // is a split node

//		double split = kdt->x[idx] + kdt->y[idx] + kdt->z[idx];
		double split;
	
		switch(axis) {
			case 0:
				split = kdt->x[idx];
				break;	
			case 1:
				split = kdt->x[idx];
				break;	
			case 2:
				split = kdt->x[idx];
				break;	
			default:
				printf("This should never happen.\n");
				break;
		}

		dist = fabs(p[axis] - split); // distance to line
	
		if(dist < *dist_best) {

			if(idx*2+1 > kdt->array_lim || idx*2+2 > kdt->array_lim) return;
			if(kdt->emptys[idx*2+1] || kdt->emptys[idx*2+2]) return;

			// explore the other side of the split line
			findBetter(kdt, (idx*2+1), (axis+1)%3, p, best, dist_best); // left child
			findBetter(kdt, (idx*2+2), (axis+1)%3, p, best, dist_best); // right child
 		}	
		
	}	
	
}

void findNearestPoint_r(kdtree* kdt, int idx, int axis, double* p, int* best, double* dist_best) {

//	printf("findNearestPoint_r() ");
//	print_point(kdt, idx);	
	
	// base case: leaf node
	if(kdt->emptys[idx*2+1] && kdt->emptys[idx*2+2]) { // check if this node is a leaf
//		printf("leaf node: ");
//		print_point(kdt, idx);
		
		double sum = 0;
		sum += pow((kdt->x[idx] - p[0]), 2);	
		sum += pow((kdt->y[idx] - p[1]), 2);	
		sum += pow((kdt->z[idx] - p[2]), 2);
		double dist = sqrt(sum);
		
		if(dist < *dist_best) {
			*dist_best = dist;
			*best = idx;
		}
		return;
	}

	
	
//	double split = kdt->x[idx] + kdt->y[idx] + kdt->z[idx];
	double split;
	switch(axis) {
		case 0:
			split = kdt->x[idx];
			break;	
		case 1:
			split = kdt->x[idx];
			break;	
		case 2:
			split = kdt->x[idx];
			break;	
		default:
			printf("This should never happen.\n");
			break;
	}

	int next_axis = (axis+1)%3;	

	if(idx*2+1 > kdt->array_lim || idx*2+2 > kdt->array_lim) return;
	if(kdt->emptys[idx*2+1] || kdt->emptys[idx*2+2]) return;

	if(p[axis] > split) {
		findNearestPoint_r(kdt, (idx*2 + 2), next_axis, p, best, dist_best); // go to right side	
		findBetter(kdt, (idx*2 + 1), next_axis, p, best, dist_best); 
	} else {
		findNearestPoint_r(kdt, (idx*2 + 1), next_axis, p, best, dist_best); // go to left side
		findBetter(kdt, (idx*2 + 2), next_axis, p, best, dist_best);
	}
}

// returns the index in the kdtree of the nearest point
int findNearestPoint(kdtree* kdt, double* p, double range) {

	int best = -1;
	double dist_best = range;
	
	int root = 0;
	int axis = 0;

	findNearestPoint_r(kdt, root, axis, p, &best, &dist_best);

	return best;	
}
//
//// using priority queue
//int findNearestPoint_pq(kdtree* kdt, double* p, double range) {
//
////	printf("findNearestPoint_pq()\n");
//
//	int best = -1;
//	double dist_best = range;
//
//	// initalize priority queue and place root node
//	pqueue* q = pq_build(PQ_SIZE);
//	pq_insert(q, 0, distance_by_idx(kdt, 0, p)); // distance to root
//
//	int current_node;
//	double left_dist, right_dist;
//	int left_idx, right_idx;
//	double current_dist;
//	while(q->num_elems != 0) { // while the queue is not empty
//	// extract the current minimum
//		current_node = pq_extract();
////		pq_print(q);
//
//		// descend this subtree until we reach a leaf; add the higher distanced' sibling in the queue as we descend the nodes
//		while(current_node*2+1 <= kdt->array_lim) { // descend until we reach a leaf node
//			
//			left_idx = current_node*2 + 1;
//			right_idx = current_node*2 + 2;
//
//			if(kdt->emptys[left_idx] && kdt->emptys[right_idx]) break; // is a leaf node
//	
//			left_dist = DBL_MAX;	
//			right_dist = DBL_MAX;
//				
//			// check which children has a shorter distance from the query point
//			if(!kdt->emptys[left_idx]) left_dist = distance_by_idx(kdt, left_idx, p);
//			if(!kdt->emptys[right_idx]) right_dist = distance_by_idx(kdt, right_idx, p);
//	
//			if(left_dist < right_dist) {
//				if(!kdt->emptys[right_idx]) pq_insert(q, right_idx, right_dist); // insert the other sibling
//				current_node = left_idx; // descend to lower distance'd node
//			} else {
//				if(!kdt->emptys[left_idx]) pq_insert(q, left_idx, left_dist);				
//				current_node = right_idx;
//			}
//
//		} // while not at leaf ; end		
//
//		// reached leaf node ; update best seen so far
//		current_dist = distance_by_idx(kdt, current_node, p);
//
//		if(current_dist < dist_best) {
//			best = current_node;
//			dist_best = current_dist;	
//		}
//	}
//
//	return best;
//}

// initiated from the host
__global__ void device_findNearestPoint(int* results, kdtree* kdt, int num_queries, double* pts_x, double* pts_y, double* pts_z, double range) {

	int thread_idx = blockDim.x * blockIdx.x + threadIdx.x; // unique thread across grid ; blockDim.x = 1024 
	if(thread_idx > num_queries-1) return;

//	printf("thread %i's query point: (%.2f, %.2f, %.2f)\n", thread_idx, pts_x[thread_idx], pts_y[thread_idx], pts_z[thread_idx]);
	
	int best = -1;
	double dist_best = range;

	// initalize priority queue and place root node
	pqueue* q = (pqueue*) malloc(sizeof(pqueue));
	if(!q) {
		printf("Thread %i: malloc pqueue failed\n");
		return;
	}

	q->max_size = 1000;
	q->num_elems = 0;
	q->elems = (int*) malloc(sizeof(int) * q->max_size);
	q->dists = (double*) malloc(sizeof(double) * q->max_size);
	
	memset(q->elems, 0, q->max_size * sizeof(int));
	memset(q->dists, 0, q->max_size * sizeof(double));

	// calculate distance from root
	double dist = pts_x[thread_idx] - kdt->x[0];
	pq_insert(q, 0, dist); // distance to root

	int current_node;
	double current_dist;
//	double left_dist, right_dist;
	int current_axis;
	int left_idx, right_idx;
	while(q->num_elems != 0) { // while the queue is not empty

		// extract the current minimum
		pq_extract(q, &current_node, &current_dist);

		if(current_dist >= dist_best) break;

		// descend this subtree until we reach a leaf; add the higher distanced' sibling in the queue as we descend the nodes
		while(current_node*2+1 <= ARRAY_LIM) { // descend until we reach a leaf node	
			
			left_idx = current_node*2 + 1;
			right_idx = current_node*2 + 2;

			if(kdt->emptys[left_idx] && kdt->emptys[right_idx]) break; // is a leaf node

			current_axis = kdt->axes[current_node];	
				
			// check which side of the kd-tree to visit
			double* array;
			double* query;
			switch(current_axis) {
				case 0:
					array = kdt->x;
					query = pts_x;
					break;
				case 1:
					array = kdt->y;
					query = pts_y;
					break;
				case 2:
					array = kdt->z;
					query = pts_z;
					break;
				default:
					printf("Invalid axis %i\n", current_axis);
					break;
			}			
		
			current_dist = query[thread_idx] - array[current_node]; // distance to split node
			if(current_dist < 0) { // visit the left child ; add the RIGHT child to the queue for later
				if(fabs(current_dist) < dist_best && !kdt->emptys[current_node*2 + 2] && q->num_elems < q->max_size) { // there is potential point that is closer than best on the other side
					pq_insert(q, right_idx, current_dist);
				}
				current_node = current_node*2 + 1; // go to left side	
			} else { // go to right side ; add the LEFT child to the queue later
				if(fabs(current_dist) < dist_best && !kdt->emptys[current_node*2 + 1] && q->num_elems < q->max_size) { // there is potential point that is closer than best on the other side
					pq_insert(q, left_idx, current_dist);
				}
				current_node = current_node*2 + 2; // go to right side		
			}

						

		} // while not at leaf ; end		

		// reached leaf node ; update best seen so far		
		// calculate distance to the current node	
		current_dist = pow((kdt->x[current_node] - pts_x[thread_idx]), 2); // overwrite old current_dist	
		current_dist += pow((kdt->y[current_node] - pts_y[thread_idx]), 2);	
		current_dist += pow((kdt->z[current_node] - pts_z[thread_idx]), 2);	
		current_dist = sqrt(current_dist);		

		if(current_dist < dist_best) {
			best = current_node;
			dist_best = current_dist;	
		}
	}
	
//	printf("thread %i) best %i, dist_best %.2f, point (%.2f, %.2f, %.2f) empty %i\n", thread_idx, best, dist_best, kdt->x[best], kdt->y[best], kdt->z[best], kdt->emptys[best]);
	results[thread_idx] = best;
}

__global__ void device_findNearestPoint2(int* results, double* nodes, char* leaves, double* qpoints, int num_queries, double range) {

	int thread_idx = blockDim.x * blockIdx.x + threadIdx.x; // unique thread across grid ; blockDim.x = 1024 
	if(thread_idx > num_queries-1) return;

	//printf("thread %i's query point: (%.2f, %.2f, %.2f)\n", thread_idx, qpoints[thread_idx], qpoints[num_queries + thread_idx], qpoints[num_queries*2 + thread_idx]);
	
	int best = -1;
	double dist_best = range;

	// initalize priority queue and place root node
	pqueue* q = (pqueue*) malloc(sizeof(pqueue));
	//if(!q) {
	//	printf("Thread %i: malloc pqueue failed\n");
	//	return;
	//}

	q->max_size = 1000;
	q->num_elems = 0;
	q->elems = (int*) malloc(sizeof(int) * q->max_size);
	q->dists = (double*) malloc(sizeof(double) * q->max_size);
	
	//memset(q->elems, 0, q->max_size * sizeof(int));
	//memset(q->dists, 0, q->max_size * sizeof(double));

	// calculate distance from root
	double dist = fabs( qpoints[thread_idx] - nodes[0] ); // distance to split node in the x-dimension
	pq_insert(q, 0, dist); // distance to root

	int current_node;
	double current_dist;
	int current_axis;
	int left_idx, right_idx;
	while(q->num_elems != 0) { // while the queue is not empty

		// extract the current minimum
		pq_extract(q, &current_node, &current_dist);

		if( current_dist >= dist_best) break;

		// descend this subtree until we reach a leaf; add the higher distanced' sibling in the queue as we descend the nodes
		while(current_node*2+1 < LEN_NODES) { // descend until we reach a leaf node	
		
			if( leaves[current_node/8] & 1<<(current_node%8) ) break; // is a leaf node
			
			left_idx = current_node*2 + 1;
			right_idx = current_node*2 + 2;
			current_axis = (int) floor(log2( (float)(current_node + 1) )) % 3; // log base 2						
			current_dist = qpoints[num_queries*current_axis + thread_idx] - nodes[LEN_NODES*current_axis + current_node]; // distance to split node

			if(current_dist < 0) { // visit the left child ; add the RIGHT child to the queue for later
				
				if(fabs(current_dist) < dist_best && q->num_elems < q->max_size && fabs(nodes[current_node] - INVALID_X_d) > 1e-32) { // there is potential point that is closer than best on the other side
					pq_insert(q, right_idx, fabs(current_dist) );
				}
				current_node = current_node*2 + 1; // go to left side	
			} else { // go to right side ; add the LEFT child to the queue later
				if(fabs(current_dist) < dist_best && q->num_elems < q->max_size && fabs(nodes[current_node] - INVALID_X_d) > 1e-32) { // there is potential point that is closer than best on the other side
					pq_insert(q, left_idx, fabs(current_dist) );
				}
				current_node = current_node*2 + 2; // go to right side		
			}	

		} // while not at leaf ; end		

		// reached leaf node ; update best seen so far		
		
		// calculate distance to this leaf node	
		current_dist = pow((nodes[current_node] - qpoints[thread_idx]), 2); // overwrite old current_dist	
		current_dist += pow((nodes[LEN_NODES + current_node] - qpoints[num_queries + thread_idx]), 2);	
		current_dist += pow((nodes[LEN_NODES*2 + current_node] - qpoints[num_queries*2 + thread_idx]), 2);	
		current_dist = sqrt(current_dist);		

		if(current_dist < dist_best) {
			best = current_node;
			dist_best = current_dist;	
		}

	} // while pq not empty ; end
	
	//printf("thread %i) best %i, dist_best %.2f, point (%.2f, %.2f, %.2f) isleaf %i\n", thread_idx, best, dist_best, nodes[best], nodes[LEN_NODES + best], nodes[LEN_NODES*2 + best], leaves[best/8]& 1<<(best%8));
	results[thread_idx] = best;
}


int main(int argc, char* argv[]) {

	time_t start, end;
	srand(time(NULL)); // initialize random number generator
	
	if(argc < 3) {
		printf("./nf <filename> <num query points>\n");
		return 1;
	}
	
	char* file = argv[1]; // file name
	int num_queries = atoi(argv[2]); // number of random query points to generate

	int num_points;
	double** points = parseJSON(file, &num_points);
	printf("Parsing complete.\n");
	
	start = clock();
	kdtree* kdt = kdtree_build(points, num_points); 
	end = clock();
	printf("Build complete. %i points, %i nodes in %.2f ms\n", kdt->num_points, kdt->num_nodes, elapsed(start, end));
	//kdtree_print(kdt);
	
	/* running the nearest neightbor algorithm */ 
// 
//	// stats
//	int num_correct[2] = {0, 0};
//	double avr_pd_pq[2] = {0.0, 0.0}; // average percent difference for priority queue implementation (b/c it is approx)
//	double total_times[3] = {0.0, 0.0, 0.0};
//	
	double range = DBL_MAX;

	double min, dist;
	int best[num_queries];
//
	double lat, longt;	
	double** query_pts = (double**) malloc(sizeof(double*) * num_queries);
	
	// format for coalesced access
	double* query_pts_x = (double*) malloc(sizeof(double) * num_queries);
	double* query_pts_y = (double*) malloc(sizeof(double) * num_queries);
	double* query_pts_z = (double*) malloc(sizeof(double) * num_queries);

	// generate a random query
	for(int i=0; i<num_queries; i++) {
			
		lat = (double) (-90) + ( (rand() * (1.0/(RAND_MAX + 1.0))) * (90 - (-90)));
		longt = (double) (-180) + ( (rand() * (1.0/(RAND_MAX + 1.0))) * (180 - (-180)));

		lat = getRadians(lat);
		longt = getRadians(longt);
	
		query_pts[i] = (double*) malloc(sizeof(double) * 3); // x y z
		query_pts[i][0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
		query_pts[i][1] = EARTH_RADIUS * cos(lat) * sin(longt); // y
		query_pts[i][2] = EARTH_RADIUS * sin(lat); // z

		// gpu-optimized format
		query_pts_x[i] = query_pts[i][0];
		query_pts_y[i] = query_pts[i][1];
		query_pts_z[i] = query_pts[i][2];

//		printf("### Query %i: (%.12f, %.12f, %.12f)\n", i, query_pts[i][0], query_pts[i][1], query_pts[i][2]);
	} // generate random query points end
		
	// brute force
	start = clock();
	for(int i=0; i<num_queries; i++) {
		min = range;
		for(int j=0; j<num_points; j++) {
			dist = distance(query_pts[i], points[j]);	
			if(dist < min) {
				min = dist; 
				best[i] = j;
			}
		}
//		printf("bf %i) result (%.2f, %.2f, %.2f), dist %.2f\n", i, points[best[i]][0], points[best[i]][1], points[best[i]][2], min);

	}
	end = clock();
	printf("Brute force time for %i queries: ", num_queries);
	elapsed(start, end);
//
//	printf("brute force results before ###\n ");
//	for(int i=0; i<num_queries; i++) {
//		printf("bf %i) result (%.2f, %.2f, %.2f)\n", i, points[best[i]][0], points[best[i]][1], points[best[i]][2]);
//	}
//		double bf_best[3];
//		bf_best[0] = points[best][0];
//	 	bf_best[1] = points[best][1];
//		bf_best[2] = points[best][2];
		
////		printf("bf result \t(%.12f, %.12f, %.12f) dist %.12f, idx %i: ", points[best][0], points[best][1], points[best][2], min, best);
//		time_t bf = elapsed(start, end);
//		total_times[0] += bf;
//
//		// cpu ; only to verify that a correct tree has been built
//		start = clock();
//		best = findNearestPoint(kdt, query_pts[i], range);	
//		end = clock();
//	
//		double b[3];
//		if(best != -1) {
//			b[0] = kdt->x[best];
//			b[1] = kdt->y[best];
//			b[2] = kdt->z[best];
//			dist = distance(query_pts[i], b);
////			printf("kd result\t(%.12f, %.12f, %.12f) dist %.12f, idx %i: ", kdt->x[best], kdt->y[best], kdt->z[best], dist, best);		
//		} else printf("kd: no points could be found...\n");
//		time_t kd = elapsed(start, end);
//		total_times[1] += kd;
//	
//		double pd = fabs(min - dist)/((min + dist)/2) * 100;
//		avr_pd_pq[0] += pd; 
//
//	//	printf("distance: %.2f percent different\n", pd);	
//		// check for correct answer
////		if(pd != 0) {
////			
////			char found = 0;
////			double d1, d2, d3;
////			
////			for(int i=0; i<kdt->array_lim; i++) {	
////				d1 = fabs(bf_best[0] - kdt->x[i]);
////				d2 = fabs(bf_best[1] - kdt->y[i]);
////				d3 = fabs(bf_best[2] - kdt->z[i]);
////		
////				if(d1 < 1e-32 && d2 < 1e-32 && d3 < 1e-32) {
//////					printf("i=%i, (%.12f, %.12f, %.12f)\n", i, kdt->x[i], kdt->y[i], kdt->z[i]);
////					found = 1;
////					break;
////				}
////			}
////		
//////			if(found) printf("no excuse... there was a point that had a shorter distance in the kdtree :(\n");
////			else printf("WAAAAAS\n");
//		
//		} else num_correct[0]++;	
//		// nearest neighbor using priority queue
//	//	start = clock();
//	//	best = findNearestPoint_pq(kdt, query_pts[i], range);	
//	//	end = clock();
//	//	if(best != -1) {
//	//		dist = distance_by_idx(kdt, best, query_pts[i]);
////	//		printf("pq result\t(%.12f, %.12f, %.12f) dist %.12f, idx %i: ", kdt->x[best], kdt->y[best], kdt->z[best], dist, best);		
//	//	} else printf("kd pq: no points could be found...\n");
//	//	kd = elapsed(start, end);
//	//	total_times[2] += kd;
//
//	//	pd = fabs(min - dist)/((min + dist)/2) * 100; // percent difference with brute force
//	//	avr_pd_pq[1] += pd; 
//	//	if(pd == 0) num_correct[1]++;
//	//	printf("distance: %.2f percent different (pq and brute force)\n", pd);	
//	//	kdtree_print(kdt);
//	//	printf("###\n\n");		
//

//	printf("kd: %i correct out of %i queries\t(%.2f %% accuracy)\taverage pd: %.2f %%\n", num_correct[0], num_queries, ((float) num_correct[0]/num_queries) * 100, (float)avr_pd_pq[0]/num_queries);
//	printf("pq: %i correct out of %i queries\t(%.2f %% accuracy)\taverage pd: %.2f %%\n", num_correct[1], num_queries, ((float) num_correct[1]/num_queries) * 100, (float)avr_pd_pq[1]/num_queries);
//	printf("bf %.2f ms, kd %.2f ms, pq %.2f ms\n", (double) total_times[0]/num_queries, (double) total_times[1]/num_queries, (double) total_times[2]/num_queries);	

	/* GPU */

	// get device information and checks if GPU is good to run	
	int dev = 0;
	cudaDeviceProp dp;
	CHECK( cudaGetDeviceProperties(&dp, dev) ); // properties, device number
	printf("Using device %i:%s\n", dev, dp.name);
	CHECK( cudaSetDevice(dev) );

//	start = clock();
	// transfer query points to device
//	double *x_d, *y_d, *z_d;
	int numBytes = num_queries * sizeof(double);
//	CHECK( cudaMalloc( (double**)&x_d, numBytes ) );
//	CHECK( cudaMalloc( (double**)&y_d, numBytes ) );
//	CHECK( cudaMalloc( (double**)&z_d, numBytes ) );
//
//	// copy over each dimension of the query points to the device	
//	CHECK( cudaMemcpy(x_d, query_pts_x, numBytes, cudaMemcpyHostToDevice)  );	
//	CHECK( cudaMemcpy(y_d, query_pts_y, numBytes, cudaMemcpyHostToDevice)  );	
//	CHECK( cudaMemcpy(z_d, query_pts_z, numBytes, cudaMemcpyHostToDevice)  );	

//	// transfer kdtree to page-locked device memory
//	kdtree* kdt_d;	
//	CHECK( cudaMallocHost( (kdtree**)&kdt_d, sizeof(kdtree) ) );
//	CHECK( cudaMemcpy(kdt_d, kdt, sizeof(kdtree), cudaMemcpyHostToDevice)  );	
//
//	numBytes = sizeof(double) * (kdt->array_lim+1);
//	CHECK( cudaMallocHost( (double**)&kdt_d->x, numBytes) );
//	CHECK( cudaMemcpy(kdt_d->x, kdt->x, numBytes, cudaMemcpyHostToDevice)  );	
//
//	CHECK( cudaMallocHost( (double**)&kdt_d->y, numBytes) );
//	CHECK( cudaMemcpy(kdt_d->y, kdt->y, numBytes, cudaMemcpyHostToDevice)  );	
//
//	CHECK( cudaMallocHost( (double**)&kdt_d->z, numBytes) );
//	CHECK( cudaMemcpy(kdt_d->z, kdt->z, numBytes, cudaMemcpyHostToDevice)  );
//
//	CHECK( cudaMemcpyToSymbol( ARRAY_LIM, &kdt->array_lim, sizeof(int), 0, cudaMemcpyHostToDevice) );	
//
//	numBytes = sizeof(char) * (kdt->array_lim+1);
//	CHECK( cudaMallocHost( (char**)&kdt_d->emptys, numBytes) );
//	CHECK( cudaMemcpy(kdt_d->emptys, kdt->emptys, numBytes, cudaMemcpyHostToDevice)  );
//
//	numBytes = sizeof(int) * (kdt->array_lim+1);
//	CHECK( cudaMallocHost( (char**)&kdt_d->axes, numBytes) );
//	CHECK( cudaMemcpy(kdt_d->axes, kdt->axes, numBytes, cudaMemcpyHostToDevice)  );
//	
//	// define block and grid size ; max is 1024 per dimension ; 1024 max per block ; so max query points = 1048576 for simple grid dimensions 
//	if(num_queries > 1048576) {
//		// figure out grid calculation...
//		return 0;
//	}
	// assumes num_queries can fit in 1 dimension (x) in the grid
//	printf("Launching %i threads\n", num_queries);
	dim3 block(num_queries); // a strip of threads
	dim3 grid(1024/num_queries + 1); // a long array of a strip of threads....
//	printf("Block dimension (x:%d, y:%d)\n", block.x, block.y);	
//	printf("Grid dimension: (x:%d, y:%d)\n", grid.x, grid.y);

	// results array
	int* results_d;
	CHECK( cudaMalloc( (int**)&results_d, sizeof(int) * num_queries) );
//	end = clock();
//	printf("Tranfer to device time: ");
//	elapsed(start, end);

//	start = clock();
//	// run on GPU
////	printf("Executing kernel\n");
//	device_findNearestPoint <<< grid, block >>> (results_d, kdt_d, num_queries, x_d, y_d, z_d, range);
	cudaError_t err = cudaGetLastError();
//	if (err != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(err));
//	cudaDeviceSynchronize(); // wait for all the threads to finish
//	// results array on host
//	int* results = (int*) malloc(sizeof(int) * num_queries); 
//	CHECK( cudaMemcpy(results, results_d, sizeof(int) * num_queries, cudaMemcpyDeviceToHost) );	
//	end = clock();
//
//	printf("GPU time for %i queries: ", num_queries);
//	elapsed(start, end);
//
//	int num_correct = 0;
//	double avr_pdif = 0.0; // average percent difference
//	for(int i=0; i<num_queries; i++) {
//		if(results[i] == -1) printf("%i) error\n", i);
//		else {
////			printf("%i) result (%.2f, %.2f, %.2f) idx = %i\n", i, kdt->x[results[i]], kdt->y[results[i]], kdt->z[results[i]], results[i]);
//			char correct = 1;
//			if( abs(points[best[i]][0] - kdt->x[results[i]]) > 1e-32) correct = 0;
//			if( abs(points[best[i]][1] - kdt->y[results[i]]) > 1e-32) correct = 0;
//			if( abs(points[best[i]][2] - kdt->z[results[i]]) > 1e-32) correct = 0;
//			
//			if(correct) num_correct++;
//			else {
//				double dist = distance(query_pts[i], points[best[i]]);	
//				double dist2;	
//				dist2 = pow((kdt->x[results[i]] - query_pts[0][i]), 2);	
//				dist2 += pow((kdt->y[results[i]] - query_pts[1][i]), 2);	
//				dist2 += pow((kdt->z[results[i]] - query_pts[2][i]), 2);	
//				dist2 = sqrt(dist2);
//					
//				avr_pdif += (float) (abs(dist - dist2)/((float)(dist + dist2)/2)) * 100;
//			} 
//		}
//	}
//
//	printf("GPU accuracy: %.2f %% (%i correct out of %i queries)\n", (float)num_correct/num_queries * 100, num_correct, num_queries);
//	printf("Average GPU distance error: %.2f %%\n", (double)avr_pdif/num_queries);

	/*

		GPU2		

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
	} kdtree;

	*/	

	// set up query point to be 1D array
	numBytes = num_queries * sizeof(double);
	double* qpoints = (double*) malloc(numBytes * 3); // 3 dimensions x,y,z
	memcpy(qpoints + num_queries*0, query_pts_x, numBytes);
	memcpy(qpoints + num_queries*1, query_pts_y, numBytes);
	memcpy(qpoints + num_queries*2, query_pts_z, numBytes);
	
	double* qpoints_d;
	CHECK( cudaMallocHost( (double**)&qpoints_d, numBytes * 3) );
	CHECK( cudaMemcpy(qpoints_d, qpoints, numBytes * 3, cudaMemcpyHostToDevice)  );	

	start = clock();
	kdtree* kdt2 = kdtree_build_gpu(points, num_points);
	end = clock();
	printf("Build kdt2 complete\n");

	double* nodes_d;
	char* leaves_d;
	numBytes = kdt2->len_nodes * sizeof(double) * 3; // for 3 dimensions x,y,z
	CHECK( cudaMallocHost( (double**)&nodes_d, numBytes) );
	CHECK( cudaMemcpy(nodes_d, kdt2->nodes, numBytes, cudaMemcpyHostToDevice)  );	

	CHECK( cudaMemcpyToSymbol( LEN_NODES, &kdt2->len_nodes, sizeof(int), 0, cudaMemcpyHostToDevice) );	
	CHECK( cudaMemcpyToSymbol( INVALID_X_d, &INVALID_X, sizeof(double), 0, cudaMemcpyHostToDevice) );	

	numBytes = kdt2->len_nodes/8 + 1;
	CHECK( cudaMallocHost( (char**)&leaves_d, numBytes) );
	CHECK( cudaMemcpy(leaves_d, kdt2->leaves, numBytes, cudaMemcpyHostToDevice)  );	

	start = clock();
	device_findNearestPoint2 <<< grid, block >>> (results_d, nodes_d, leaves_d, qpoints_d, num_queries, range); 
	cudaDeviceSynchronize(); // wait for all the threads to finish
	err = cudaGetLastError();
	if (err != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(err));

	// results array on host
	int* results = (int*) malloc(sizeof(int) * num_queries); 
	CHECK( cudaMemcpy(results, results_d, sizeof(int) * num_queries, cudaMemcpyDeviceToHost) );	
	end = clock();

	printf("GPU time for %i queries: ", num_queries);
	elapsed(start, end);

//	printf("brute force results: ");
//	for(int i=0; i<num_queries; i++) {
//		printf("bf %i) result (%.2f, %.2f, %.2f)\n", i, points[best[i]][0], points[best[i]][1], points[best[i]][2]);
//	}

	int num_correct = 0;
	double avr_pdif = 0.0; // average percent difference
	for(int i=0; i<num_queries; i++) {
		if(results[i] == -1) printf("%i) error\n", i);
		else {
//			printf("thread %i result: %.2f, %.2f, %.2f\n", i, kdt2->nodes[results[i]], kdt2->nodes[kdt2->len_nodes + results[i]], kdt2->nodes[kdt2->len_nodes*2 + results[i]]);
//			printf("brute force %i result: %.2f, %.2f, %.2f\n", i, points[best[i]][0], points[best[i]][1], points[best[i]][2]);
			char correct = 1;
			if( abs(points[best[i]][0] - kdt2->nodes[results[i]]) > 1e-32) correct = 0;
			if( abs(points[best[i]][1] - kdt2->nodes[kdt2->len_nodes + results[i]]) > 1e-32) correct = 0;
			if( abs(points[best[i]][2] - kdt2->nodes[kdt2->len_nodes*2 + results[i]]) > 1e-32) correct = 0;
			
			if(correct) num_correct++;
			else {
				double dist = distance(query_pts[i], points[best[i]]);	
				double dist2;	
				dist2 = pow((kdt2->nodes[results[i]] - query_pts_x[i]), 2);	
				dist2 += pow((kdt2->nodes[kdt2->len_nodes + results[i]] - query_pts_y[i]), 2);	
				dist2 += pow((kdt2->nodes[kdt2->len_nodes*2 + results[i]] - query_pts_z[i]), 2);	
				dist2 = sqrt(dist2);
					
				avr_pdif += (float) (abs(dist - dist2)/((float)(dist + dist2)/2)) * 100;
			} 
		}
	}

	printf("GPU2 accuracy: %.2f %% (%i correct out of %i queries)\n", (float)num_correct/num_queries * 100, num_correct, num_queries);
	printf("Average GPU2 distance error: %.2f %%\n", (double)avr_pdif/num_queries);

	// Always call this function when ending program
	cudaDeviceReset();

	return 0;		
}
