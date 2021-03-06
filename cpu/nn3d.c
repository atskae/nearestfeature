#include <stdlib.h>
#include <stdio.h>
#include <time.h> 	
#include <string.h>
#include <float.h> // max double value 
#include <math.h>

#include "jsmn.h" // JSON file parser
#include "kdtree.h"
#include "pq.h"

/*
	Implementation of "GPU-accelerated Nearest Neighbor Search for 3D Registration" Qiu 2009	
	
	Notes;
	- Must convert (lat, longt) to (x, y, z) to implement kdtree
	
*/

#define BUFFER_SIZE 1600000
#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371 // in km
#define PQ_SIZE 20 // priority queue size

#define DIM 3 // dimension of cartesian coordinates

// returns the time elapsed in miliseconds 
double elapsed(time_t start, time_t end) {
	double t = ( ((double) (end - start)) / CLOCKS_PER_SEC ) * 1000;	
//	printf("%.5f ms elapsed.\n", t);
	return t;
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
				points[i] = malloc(DIM * sizeof(double));
				geopoints[i] = malloc(2 * sizeof(double));

				// save the geo coordinates	
				geopoints[i][0] = lat;
				geopoints[i][1] = longt;

//				printf("longt=%.2f, lat=%.2f | ", longt, lat);

				// convert (lat, long) to (x,y,z)
			
				lat = getRadians(lat);
				longt = getRadians(longt);	
				points[i][0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
				points[i][1] = EARTH_RADIUS * cos(lat) * sin(longt);
				points[i][2] = EARTH_RADIUS * sin(lat);
		
//				printf("(%.2f, %.2f, %.2f)\n", points[i][0], points[i][1], points[i][2]);
	
				(*num_points)++;
				if(*num_points == points_buffer_size) {					
					points_buffer_size = points_buffer_size * 2;
					points = realloc(points, points_buffer_size * sizeof(double*));
					geopoints = realloc(geopoints, points_buffer_size * sizeof(double*));
					if(!points || !geopoints) {
//						printf("%i %i realloc() failed for buffer size %lu\n", points == 0, geopoints == 0, (*num_points + points_buffer_size) * sizeof(double*));
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
	points = realloc(points, (*num_points) * sizeof(double*));

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
				split = kdt->y[idx];
				break;	
			case 2:
				split = kdt->z[idx];
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
			split = kdt->y[idx];
			break;	
		case 2:
			split = kdt->z[idx];
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
		if(fabs(p[axis] - split) < *dist_best) findBetter(kdt, (idx*2 + 1), next_axis, p, best, dist_best); 
	} else {
		findNearestPoint_r(kdt, (idx*2 + 1), next_axis, p, best, dist_best); // go to left side
		if(fabs(p[axis] - split) < *dist_best) findBetter(kdt, (idx*2 + 2), next_axis, p, best, dist_best);
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

/*
	CPU 
*/

void findNearestPoint_cpu(double** points, node* n, double* p, int* best, double* dist_best) {

	// base case: leaf node
	if(n->leaf) {
		double dist = distance(p, points[n->p]);	
		if(dist < *dist_best) {
			*dist_best = dist;
			*best = n->p;
		}
		return;
	}

	double dist = fabs(p[n->axis] - n->split);	
	if(p[n->axis] > n->split) { // visit the right side
		if(n->right) findNearestPoint_cpu(points, n->right, p, best, dist_best);
		// check for potentially closer nodes on the other side	
		if(n->left && dist < *dist_best) findNearestPoint_cpu(points, n->left, p, best, dist_best); 
	} else { // visit the left side
		if(n->left) findNearestPoint_cpu(points, n->left, p, best, dist_best);
		// check for potentially closer nodes on the other side
		if(n->right && dist < *dist_best) findNearestPoint_cpu(points, n->right, p, best, dist_best);
	}
}

// using priority queue
int findNearestPoint_pq(kdtree* kdt, double* p, double range) {

//	printf("findNearestPoint_pq()\n");
	int best = -1;
	double dist_best = range;

	// initalize priority queue and place root node
	pqueue* q = pq_build(PQ_SIZE);
	pq_insert(q, 0, fabs(p[0] - kdt->x[0]) ); // distance to root

	int current_node;
//	double left_dist, right_dist;
	int left_idx, right_idx;
	double current_dist;
	int current_axis;
	while(q->num_elems != 0) { // while the queue is not empty
	// extract the current minimum
		current_node = pq_extract(q, &current_dist);	
		if(current_dist >= dist_best) break;

		// descend this subtree until we reach a leaf; add the higher distanced' sibling in the queue as we descend the nodes
		while(current_node*2+1 <= kdt->array_lim) { // descend until we reach a leaf node
			left_idx = current_node*2 + 1;
			right_idx = current_node*2 + 2;

			if(kdt->emptys[left_idx] && kdt->emptys[right_idx]) break; // is a leaf node

			current_axis = kdt->axes[current_node];	
				
			// check which side of the kd-tree to visit
			double* array;
			switch(current_axis) {
				case 0:
					array = kdt->x;
					break;
				case 1:
					array = kdt->y;
					break;
				case 2:
					array = kdt->z;
					break;
				default:
					printf("Invalid axis %i\n", current_axis);
					break;
			}			
		
			current_dist = p[current_axis] - array[current_node]; // distance to split node
			if(current_dist < 0) { // visit the left child ; add the RIGHT child to the queue for later
				if(fabs(current_dist) < dist_best && !kdt->emptys[current_node*2 + 2] && q->num_elems < q->max_size) { // there is potential point that is closer than best on the other side
					pq_insert(q, right_idx, fabs(current_dist) );
				}
				current_node = current_node*2 + 1; // go to left side	
			} else { // go to right side ; add the LEFT child to the queue later
				if(fabs(current_dist) < dist_best && !kdt->emptys[current_node*2 + 1] && q->num_elems < q->max_size) { // there is potential point that is closer than best on the other side
					pq_insert(q, left_idx, fabs(current_dist) );
				}
				current_node = current_node*2 + 2; // go to right side		
			}

		} // while not at leaf ; end		

		// reached leaf node ; update best seen so far
		current_dist = distance_by_idx(kdt, current_node, p);

		if(current_dist < dist_best) {
			best = current_node;
			dist_best = current_dist;	
		}
	}

	return best;
}

/* NEED TO FIX THIS*/
// using priority queue with a GPU tree
int findNearestPoint_pq2(double* nodes, int LEN_NODES, char* leaves, int* axes, double INVALID_X, double* p, double range) {

	int best = -1;
	double dist_best = range;

	// initalize priority queue and place root node
	pqueue* q = pq_build(PQ_SIZE);
	pq_insert(q, 0, fabs(p[0] - nodes[0])); // distance to root ; x-split

	int current_node;
	int left_idx, right_idx;
	double current_dist;
	int current_axis;

	while(q->num_elems != 0) { // while the queue is not empty

		// extract the current minimum
		current_node = pq_extract(q, &current_dist);
		if( current_dist >= dist_best) break;

		// descend this subtree until we reach a leaf; add the higher distanced' sibling in the queue as we descend the nodes
		while(current_node*2 + 1 < LEN_NODES) { // descend until we reach a leaf node	
		
			if(leaves[current_node/8] & 1<<(current_node%8)) break;

			left_idx = current_node*2 + 1;
			right_idx = current_node*2 + 2;
//			current_axis = (int) floor(log2( (float)(current_node + 1) )) % 3; // log base 2						
			current_axis = axes[current_node];
			current_dist = p[current_axis] - nodes[LEN_NODES*current_axis + current_node]; // distance to split node

			if(current_dist < 0) { // visit the left child ; add the RIGHT child to the queue for later
				current_node = current_node*2 + 1; // go to left side	
			
				if(q->num_elems < q->max_size && fabs(current_dist) < dist_best && fabs(nodes[LEN_NODES*current_axis + current_node] - INVALID_X) > 1e-32) {
					pq_insert(q, right_idx, fabs(current_dist) );
				}
				
			} else { // go to right side ; add the LEFT child to the queue later
				current_node = current_node*2 + 2; // go to right side
		
				if(q->num_elems < q->max_size && fabs(current_dist) < dist_best && fabs(nodes[LEN_NODES*current_axis + current_node] - INVALID_X) > 1e-32) {
						pq_insert(q, left_idx, fabs(current_dist) );
				}
			}	

		} // while not at leaf ; end		

		// reached leaf node ; update best seen so far		
		
		// calculate distance to this leaf node	
		current_dist = pow((nodes[current_node] - p[0]), 2); // overwrite old current_dist	
		current_dist += pow((nodes[LEN_NODES + current_node] - p[1]), 2);	
		current_dist += pow((nodes[LEN_NODES*2 + current_node] - p[2]), 2);	
		current_dist = sqrt(current_dist);		

		if(current_dist < dist_best) {
			best = current_node;
			dist_best = current_dist;	
		}
	}

	return best;
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
//	kdtree_print(kdt);
	end = clock();
	printf("Build complete. %i points, %i nodes in %.2f ms\n", kdt->num_points, kdt->num_nodes, elapsed(start, end));

	start = clock();
	node* kdt_cpu = kdtree_build_cpu(points, num_points);
	end = clock();
	printf("CPU kdtree build complete in %.2f ms\n", elapsed(start, end));
	int best_cpu;
	double dist_best_cpu;

	kdtree* kdt_gpu= kdtree_build_gpu(points, num_points);
	
	/* running the nearest neightbor algorithm */ 
	 
	// stats
//	int num_correct[2] = {0, 0};
//	double avr_pd_pq[2] = {0.0, 0.0}; // average percent difference for priority queue implementation (b/c it is approx)
//	double total_times[3] = {0.0, 0.0, 0.0};
//
//	int num_correct_cpu = 0;
//	double total_time_cpu = 0.0;
//	double pd_cpu = 0.0;	
//
//	double min, dist;
//	int best;

	double lat, longt;	
	double query_pt[3];
	double range = DBL_MAX;

	/*
		write data to .csv file
	*/
	FILE* output= fopen("results.csv", "w");
	fprintf(output, "num queries,brute force,kdtree cpu,kdtree cpu accuracy,priority queue,priority queue accuracy,priority queue percent difference\n");

	for(int double_rate = 1; num_queries*double_rate <= 12800; double_rate*=2) {

		// stats
		int num_correct[2] = {0, 0};
		double avr_pd_pq[2] = {0.0, 0.0}; // average percent difference for priority queue implementation (b/c it is approx)
		double total_times[3] = {0.0, 0.0, 0.0};
	
		int num_correct_cpu = 0;
		double total_time_cpu = 0.0;
		double pd_cpu = 0.0;	
	
		double min, dist;
		int best;
	
		for(int i=0; i<num_queries*double_rate; i++) {
	
			// generate a random query
		
			lat = (double) (-90) + ( (rand() * (1.0/(RAND_MAX + 1.0))) * (90 - (-90)));
			longt = (double) (-180) + ( (rand() * (1.0/(RAND_MAX + 1.0))) * (180 - (-180)));
	
			lat = getRadians(lat);
			longt = getRadians(longt);
		
			query_pt[0] = EARTH_RADIUS * cos(lat) * cos(longt); // x
			query_pt[1] = EARTH_RADIUS * cos(lat) * sin(longt); // y
			query_pt[2] = EARTH_RADIUS * sin(lat); // z
		
	//		printf("### Query %i: (%.12f, %.12f, %.12f)\n", i, query_pt[0], query_pt[1], query_pt[2]);
			
			// brute force
			start = clock();
			min = range;
			best = -1;
			for(int j=0; j<kdt->num_points; j++) {
				dist = distance(query_pt, points[j]);	
				if(dist < min) {
					min = dist; 
					best = j;
				}
			}
			end = clock();
			double bf = elapsed(start, end);
			total_times[0] += bf;	
	
			if(best == -1) {
				printf("Something went wrong...\n");
				continue;
			}
	
			double bf_best[3];
			bf_best[0] = points[best][0];
		 	bf_best[1] = points[best][1];
			bf_best[2] = points[best][2];
			
//			printf("bf result \t(%.12f, %.12f, %.12f) dist %.12f, idx %i: ", points[best][0], points[best][1], points[best][2], min, best);
	
			/*
				CPU kdtree (node structs)
			*/
			best_cpu = -1; // index of the points array
			dist_best_cpu = DBL_MAX;
			start = clock();
			findNearestPoint_cpu(points, kdt_cpu, query_pt, &best_cpu, &dist_best_cpu);
			end = clock();
			total_time_cpu += elapsed(start, end);
			//printf("CPU elapsed: %.2f ms\n", elapsed(start, end));;
			if(best != -1) {
				if(fabs(dist_best_cpu - min) <= 1e-32) num_correct_cpu++;
		
			} else if(best != -1) {
				pd_cpu += fabs(min - dist_best_cpu)/((min + dist_best_cpu)/2) * 100;	
			} else printf("CPU error\n");
	
			// cpu but gpu-optimized kdtree; only to verify that a correct tree has been built
			start = clock();
			best = findNearestPoint(kdt, query_pt, range);	
			end = clock();
		
			double b[3];
			if(best != -1) {
				b[0] = kdt->x[best];
				b[1] = kdt->y[best];
				b[2] = kdt->z[best];
				dist = distance(query_pt, b);
	//			printf("kd result\t(%.12f, %.12f, %.12f) dist %.12f, idx %i: ", kdt->x[best], kdt->y[best], kdt->z[best], dist, best);		
			} else printf("kd: no points could be found...\n");
			double kd;
			kd = elapsed(start, end);
			total_times[1] += kd;
			double pd;
			pd = fabs(min - dist)/((min + dist)/2) * 100;
			avr_pd_pq[0] += pd; 
	
			//printf("distance: %.2f percent different\n", pd);	
		  // check for correct answer
			if(pd != 0) {
				
				char found = 0;
				double d1, d2, d3;
				
				for(int i=0; i<kdt->array_lim; i++) {	
					d1 = fabs(bf_best[0] - kdt->x[i]);
					d2 = fabs(bf_best[1] - kdt->y[i]);
					d3 = fabs(bf_best[2] - kdt->z[i]);
			
					if(d1 < 1e-32 && d2 < 1e-32 && d3 < 1e-32) {
	//					printf("i=%i, (%.12f, %.12f, %.12f)\n", i, kdt->x[i], kdt->y[i], kdt->z[i]);
						found = 1;
						break;
					}
				}
			
	//			if(found) printf("no excuse... there was a point that had a shorter distance in the kdtree :(\n");
	//			else printf("WAAAAAS\n");
			
			} else num_correct[0]++;
		
			// nearest neighbor using priority queue
			start = clock();
//			best = findNearestPoint_pq(kdt, query_pt, range);	
			best = findNearestPoint_pq2(kdt_gpu->nodes, kdt_gpu->len_nodes, kdt_gpu->leaves, kdt->axes, kdt->invalid_x, query_pt, range);
			end = clock();
			if(best != -1) {
				dist = distance_by_idx(kdt, best, query_pt);
//				printf("pq result\t(%.12f, %.12f, %.12f) dist %.12f, idx %i\n", kdt->x[best], kdt->y[best], kdt->z[best], dist, best);		
			} else printf("kd pq: no points could be found...\n");
			kd = elapsed(start, end);
			total_times[2] += kd;
	
			pd = fabs(min - dist)/((min + dist)/2) * 100; // percent difference with brute force
			avr_pd_pq[1] += pd; 
			if(pd == 0) num_correct[1]++;
	//		printf("distance: %.2f percent different (pq and brute force)\n", pd);	
			//kdtree_print(kdt);
	//		printf("###\n\n");		
	
		} // each query point ; end
	
	 	fprintf(output, "%i,%.5f,%.5f,%.2f,%.5f,%.2f,%.2f\n", num_queries*double_rate, total_times[0], total_time_cpu, (double)num_correct_cpu/(num_queries*double_rate) * 100, total_times[2], ((double) num_correct[1]/(num_queries*double_rate) * 100), (double)avr_pd_pq[1]/(double_rate*num_queries) );
		
		printf("bf: %.5f, CPU-optimized: %.2f, pq: %.5f\n", total_times[0], total_time_cpu, total_times[2]);
	//	printf("kd: %i correct out of %i queries\t(%.2f %% accuracy)\taverage pd: %.2f %%\n", num_correct[0], num_queries, ((float) num_correct[0]/num_queries) * 100, (float)avr_pd_pq[0]/num_queries);
		printf("cp: %i correct out of %i queries (%.2f%% accuracy)\taverage pd: %.2f %%\n", num_correct_cpu, num_queries*double_rate, (double)num_correct_cpu/(num_queries*double_rate)* 100, (double)pd_cpu/(double_rate*num_queries));
		printf("pq: %i correct out of %i queries\t(%.2f %% accuracy)\taverage pd: %.2f %%\n", num_correct[1], num_queries*double_rate, ((double) num_correct[1]/(num_queries*double_rate)) * 100, (double)avr_pd_pq[1]/(double_rate*num_queries));
	//	printf("bf %.5f ms, kd %.5f ms, pq %.5f ms\n", (double) total_times[0]/num_queries, (double) total_times[1]/num_queries, (double) total_times[2]/num_queries);	
	} // double rate ; end
	fclose(output);
	
	return 0;		
}
