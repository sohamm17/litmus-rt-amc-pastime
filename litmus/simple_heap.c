/*
A simple implementation of heap for task priority ordering in array;
Just for the sake, nothing fancy.
*/

/*
struct bin_heap_node {
  struct bin_heap_node* left;
  struct bin_heap_node* right;
  struct bin_heap node* parent;
  
  void* data;
}*/

#include <litmus/simple_heap.h>
#include <linux/printk.h>
/*#define LEFT_CHILD_INDEX (2 * i + 1)
#define RIGHT_CHILD_INDEX (2 * i + 2)
#define PARENT_INDEX ((i - 1) / 2)

#define CUR_NODE(i) heap_array[i]
#define LEFT_CHILD(i) heap_array[LEFT_CHILD_INDEX(i)]
#define RIGHT_CHILD(i) heap_array[RIGHT_CHILD_INDEX(i)]
#define PARENT(i) heap_array[PARENT_INDEX(i)]
#define MAX_SIMPLE_HEAP_SIZE 20


simple_heap_node heap_array[MAX_SIMPLE_HEAP_SIZE];
*/


array_node sorted_arr[MAX_ARR_SIZE];
int last_sorted_idx = 0;

int get_last_index_arr(void) {
  return last_sorted_idx;
}

int insert_node_arr(void *data, int prio) {
  if(last_sorted_idx < MAX_ARR_SIZE) {
    sorted_arr[last_sorted_idx].data = data;
    sorted_arr[last_sorted_idx].prio = prio;
    //printk("last idx: %d\n", last_sorted_idx);
    last_sorted_idx++;
    return 0;
  }
  return -1;
}

int bin_search_node_arr(int start, int end, int prio) {
  int middle;

  if(start > end || last_sorted_idx <= 0)
    return -1;

  middle = (start + end) / 2;
  /*printk("Middle: %d\n", middle);
  printk("MiddlePrio: %d\n", sorted_arr[middle].prio);
  printk("prio: %d\n", prio);*/
  if (sorted_arr[middle].prio == prio) {
    return middle;
  } else if(sorted_arr[middle].prio < prio) {
    return bin_search_node_arr(middle + 1, end, prio);
  } else if(sorted_arr[middle].prio > prio) {
    return bin_search_node_arr(start, middle - 1, prio);
  }
}

int search_node_arr(int prio) {
  return bin_search_node_arr(0, last_sorted_idx - 1, prio);
}

void init_node_arr(void) {
  int i;
  for(i = 0; i < MAX_ARR_SIZE; i++) {
    sorted_arr[i].data = 0;
    sorted_arr[i].prio = 0;
  }
  last_sorted_idx = 0;
}

/*int last_heap_node = 0;

inline void init_node(int index, void* data, int key_priority) {
  CUR_NODE(index).key_priority = key_priority;
  CUR_NODE(index).data = data;
}

inline int comparator(simple_heap_node* first, simple_heap_node* second) {
  if(first->key_priority == second->key_priority && first->data == second->data) {
    return 0;
  } else if (first->key_priority < second->key_priority)
    return -1;
  else
    return 1;
}

// this function sifts up the node_indexes until it is <= paren
void sift_up(int node_index) {
  if(comparator(&CUR_NODE(node_index), &PARENT(node_index)) > 1) {
    simple_heap_node temp = CUR_NODE(node_index);
    CUR_NODE(node_index) = PARENT(node_index);
    PARENT(node_index) = CUR_NODE(node_index);
    sift_up(PARENT_INDEX(i));
  }
}

int insert_node(void *data, int key_priority) {
  if(last_heap_node <= 0) {
    init_node(last_heap_node, data, key_priority);
    last_heap_node++;
    return 0;
  } else if(last_heap_node >= MAX_SIMPLE_HEAP_SIZE) {
    return -1;
  }
  init_node(last_heap_node, data, key_priority);
  sift_up(last_heap_node);
  last_heap_node++;
}

// search a particular node from a particular start_index
int search_node(int start_index, simple_heap_node* temp) {
  if(comparator(&CUR_NODE(start_index), &temp) == 0) {
    return start_index;
  } else if (comparator(&LEFT_CHILD(start_index), &temp) <= 0) {
    return search_node(LEFT_CHILD_INDEX(start_index), temp);
  } else if(comparator(&RIGHT_CHILD(start_index), &temp) <= 0) {
    return search_node(RIGHT_CHILD_INDEX(start_index), temp);
  } else
    return -1;
}

// get_node_index by priority and data value
int get_node_index(void* data, int key_priority) {
  simple_heap_node temp;
  temp.data = data; temp.key_priority = key_priority;

  return search_node(0, &temp);
}

// subtree size should not grow more than MAX_SIMPLE_HEAP_SIZE
void get_node_subtree(void* data, int key_priority, simple_heap_node** subtree) {
  int subtree_head_index = get_node_index(data, key_priority);

}*/
