// Simple sorted array for ResponseTime Calculation


#ifndef SIMPLE_HEAP_H
#define SIMPLE_HEAP_H

/*typedef struct {
  int key_priority;
  void* data;
} simple_heap_node;
*/

#define MAX_ARR_SIZE 31
typedef struct {
  int prio;
  void* data;
} array_node;

extern array_node sorted_arr[];

int insert_node_arr(void *data, int prio);
int search_node_arr(int prio);
void init_node_arr(void);
int get_last_index_arr(void);

#endif
