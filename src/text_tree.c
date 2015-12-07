#include "tree.h"

int main(int argc, char *argv[]) {
  int i = 0;

  char *treestr;
  treestr = NULL;
  treestr = read_TreeFile(argv[1]);

  Node *masternode;
  char **nodestr;

  masternode = NULL;
  masternode = read_TreeStr(treestr);

  nodestr = (char **)malloc(sizeof(char *));
  *nodestr = NULL;
  print_Node_string(masternode,nodestr);

  printf("Master tree string is  %s\n",*nodestr);

  float height = 0;

  find_total_height(masternode,&height);

  float max_height = 0;

  find_max_height(masternode,&max_height);

  printf("total height is %f\n",height);
  printf("max height is   %f\n",max_height);

  draw_Tree(masternode,160,2);
}
