#include "tree.h"

//Trees - creation and pruning


int main(int argc, char *argv[]) {

  Node *node;
  Node *mono;
  Node *topnode;

  char **nodestr;
  char **monostr;

  nodestr = (char **)malloc(sizeof(char *));
  monostr = (char **)malloc(sizeof(char *));

  *nodestr = NULL;
  *monostr = NULL;

  node = read_Tree(argv[1]);

  print_Node_string(node,nodestr);

  printf("String %s\n",*nodestr);
  
  topnode = node;

  int num = 2;
  Node *newroot;

  while (num  < argc) {

    int leaves = 0;

    count_Leaves(topnode,&leaves);

    printf("Number of leaves %d\n",leaves);

    mono = find_Leaf(topnode,argv[num],mono);
    
    printf("Mono %s\n",mono->name);
    printf("Parent %s\n",mono->parent->name);
    
    newroot = delete_Leaf(topnode,mono);
    
    printf("New root %s\n",newroot->name);

    monostr = (char **)malloc(sizeof(char *));
    *monostr = NULL;

    print_Node_string(newroot,monostr);
    
    printf("String %s\n",*monostr);

    topnode = newroot;

    num++;
  }
}

