#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <sys/time.h>

int maptox(int mapindex, int width);
int maptoy(int mapindex, int width);
int maptox3d(int mapindex, int height, int width);
int maptoy3d(int mapindex, int height, int width);
int maptoz3d(int mapindex, int height, int width);
int mapIndex3D(int r, int c, int z, int nr, int nc, int nz);
void print_timing(FILE *fp, struct timeval start, struct timeval end);

/* ---------------------------------------------------------------------------
   Index worklist

   Every flood fill and front propagation here used to declare its own local
   `struct list { int num_elem; int *elem; }` and write out the push, the pop
   and the overflow check by hand. They all share this one type instead.

   elem[0 .. num_elem-1] holds every index pushed, in push order, so one list
   can serve as the worklist *and* as the record of what it visited: take it
   in FIFO order, pop it in LIFO order, or read elem[] directly. Capacity is
   fixed at init; push reports overflow rather than growing, because callers
   size the list from the image and a full list means a logic error upstream.
   --------------------------------------------------------------------------- */
struct index_list {
  int *elem;
  int num_elem; /* number of indices pushed */
  int first;    /* next index to hand out in FIFO order */
  int max_elem; /* capacity */
};

/* Returns 0 on success, 1 if the allocation failed. */
int list_init(struct index_list *list, int max_elem);
void list_free(struct index_list *list);
/* Forget everything, keeping the allocation. */
void list_clear(struct index_list *list);
/* Returns 0 on success, 1 if the list is full (nothing is stored). */
int list_push(struct index_list *list, int mapindex);
/* LIFO: removes and returns the most recently pushed index, -1 if empty. */
int list_pop(struct index_list *list);
/* FIFO: returns the next unvisited index without erasing it, -1 if empty. */
int list_take(struct index_list *list);
void list_swap(struct index_list *a, struct index_list *b);

/* ---------------------------------------------------------------------------
   Adjacency

   These are the only places that turn a linear index into neighbouring linear
   indices. Every neighbour is derived from the pixel's real (row, column) --
   and (plane, row, column) in 3D -- so a pixel in column 0 never picks up the
   last column of the neighbouring row, which is what indexing by raw offset
   and then bounds-checking with maptox/maptoy used to allow.

   Each fills `neighbors` with the in-bounds neighbours of `mapindex` and
   returns how many were written. Neighbours come back in the same order the
   hand-written loops produced them, so traversal order is unchanged.
   --------------------------------------------------------------------------- */
#define MAX_NEIGHBORS_2D 8
#define MAX_NEIGHBORS_3D 26
#define MAX_NEIGHBORS_3D_FACES 6

/* 8-connected. */
int neighbors2D(int mapindex, int height, int width, int *neighbors);
/* 26-connected. */
int neighbors3D(int mapindex, int height, int width, int depth, int *neighbors);
/* 6-connected: the face neighbours only. */
int neighbors3D_faces(int mapindex, int height, int width, int depth, int *neighbors);

#endif /* UTILS_H */
