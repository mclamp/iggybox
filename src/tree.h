#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

typedef struct _Node {

  struct _Node *left;
  struct _Node *right;
  struct _Node *parent;
  char* name;
  float   branch_length;
  
  int    isRoot;
  int    ycount;
  int    count;
  float  dist_from_root;

} Node;

Node *create_Node();
Node *read_Tree(char *file);
Node *read_TreeStr(char *file);
char *read_TreeFile(char *file);
void  print_Node(Node *node);
Node *find_Leaf(Node *node,char *name,Node *found);
Node *delete_Leaf(Node *root,Node *node);
void  track_Nodes(Node *node);
void count_Leaves(Node *node,int *countptr);
void print_Node_string(Node *node, char **str);
void find_total_height(Node *node, float *height);
void find_max_height(Node *node, float *height);

void fill_Boxes(Node *node, char *boxes,int width, int height,float  dist_per_char,int chars_per_ylevel);

void draw_Tree(Node *topnode,int width,int chars_per_ylevel);
void print_Boxes(char *boxes, int width, int height,float dist_per_char);
