#include "tree.h"

Node *create_Node() {
  Node * node;

  node = (Node *)malloc(sizeof(Node));

  node->left   = NULL;
  node->right  = NULL;
  node->parent = NULL;

  node->isRoot = 0;
  node->branch_length = 0;

  node->name   = (char *)malloc(2*sizeof(char));
  node->name[0] = 'N';
  node->name[1] = '\0';

  node->count  = 0;
  node->ycount = 0;

  node->dist_from_root = 0;

  return node;

}

Node *read_Tree(char *file) {
  char *treestr;

  treestr = read_TreeFile(file);

  return read_TreeStr(treestr);
}

Node *read_TreeStr(char *treestr) {
  Node *currnode;
  Node *rootnode;
  Node *tmpnode;

  rootnode = create_Node();
  
  rootnode->isRoot = 1;
  rootnode->name = "root";

  currnode = rootnode;


  int pos = 0;
  char c;

  int ycount = 0;
  
  while (pos < strlen(treestr)) {

    c = treestr[pos];

    if (c == '(') {

      // Add in a new left node

      tmpnode = create_Node();
      tmpnode->parent = currnode;
      
      currnode->left = tmpnode;
      currnode = tmpnode;

      //      printf("Got new left node %s %s\n",tmpnode->name,tmpnode->parent->name);

      pos++;
      //      print_Node(tmpnode);
      //print_Node(rootnode);

    } else if (c == ',') {

      Node *parent;
      
      parent = currnode->parent;

      // Move up and over to the right node

      tmpnode = create_Node();
      
      tmpnode->parent = parent;

      parent->right = tmpnode;

      currnode = tmpnode;



      //      printf("Got new right node\n");

      pos++;

    } else if (c == ':') {
      
      pos++;

      char *diststr = NULL;

      while (treestr[pos] != ',' &&
	     treestr[pos] != ')') {

	if (diststr == NULL) {

	  diststr = (char *)malloc(2*sizeof(char));
	  diststr[0] = treestr[pos];
	  diststr[1] = '\0';

	  pos++;

	} else {

	  int len = strlen(diststr);

	  diststr = (char *)realloc(diststr,(strlen(diststr)+2)*sizeof(char));

	  diststr[len] = treestr[pos];
	  diststr[len+1] = '\0';

	  pos++;
	}
      }


      float dist = atof(diststr);

      //      printf("Got new distance %f for %s\n",dist,currnode->name);      

      currnode->branch_length = dist;
      currnode->dist_from_root = currnode->parent->dist_from_root + dist;

    } else if (c == ')') {

      // Move up to the parent

      Node *parent;

      parent = currnode->parent;

      currnode = parent;

      //      printf("Going up to the parent %s\n",currnode->name);

      if (currnode->count != 1) {
	currnode->count = currnode->left->count + currnode->right->count;
	currnode->ycount = (currnode->left->ycount + currnode->right->ycount)/2;
      }

      pos++;
    } else if (c == ';') {

      track_Nodes(rootnode);
      return rootnode;

    } else {

      // Read org name

      char *org = NULL;

      while (treestr[pos] != ':' &&
	     treestr[pos] != ';') {

	if (org == NULL) {
	  org = (char *)malloc(2*sizeof(char));

	  org[0] = treestr[pos];
	  org[1] = '\0';

	  pos++;

	} else {
	  int len = strlen(org);

	  org = (char *)realloc(org,(strlen(org)+2)*sizeof(char));

	  org[len] = treestr[pos];
	  org[len+1] = '\0';
	

	  pos++;
	}
      }
      //      printf("Got org %s for current node\n",org);

      ycount++;

      currnode->count = 1;
      currnode->ycount = ycount;
      currnode->name = org;
      currnode->left = NULL;
      currnode->right = NULL;
    }
  }

  rootnode->count = rootnode->left->count + rootnode->right->count;

  rootnode->ycount = (rootnode->left->count + rootnode->right->count)/2;

  track_Nodes(rootnode);

  return rootnode;

}

void track_Nodes(Node *node) {
  
  if (node->parent != NULL) {
    node->dist_from_root = node->parent->dist_from_root + node->branch_length;
  }
  if (node->left != NULL) {
    Node *left;

    left = node->left;
    left->parent = node;

    track_Nodes(left);
  }
  if (node->right != NULL) {
    Node *right;
    
    right = node->right;

    right->parent = node;
    
    track_Nodes(right);
  }
  
}

char *read_TreeFile(char *file) {

    FILE *f;

    f = fopen(file,"r");

    char *str = NULL;
    
    char c;

    while ((c = fgetc(f)) != EOF) {

      if (c != '\n') {

	if (str == NULL) {
	  str = (char *)malloc(2*sizeof(char));

	  str[0] = c;
	  str[1] = '\0';

	} else {
	  int len = strlen(str);

	  str = (char *)realloc(str,(strlen(str)+2)*sizeof(char));

	  str[len] = c;
	  str[len+1] = '\0';
	}
      }
    }

    return str;
}

void print_Node(Node *node) {

  char *name;

  name = node->name;

  if (name == NULL) {
    name = (char *)malloc(2*sizeof(char));
    name[0] = 'N';
    name[1] = '\0';
  }
  printf("Node %f\t%s\t%d\t%d\n",node->branch_length,node->name,node->count,node->ycount);

  if (node->left != NULL) {
    print_Node(node->left);
  }
  if (node->right != NULL) {
    print_Node(node->right);
  }
}

void print_Node_string(Node *node, char **str) {

  char *tmpstr;

  tmpstr = *str;

  if (node->left != NULL &&
      node->right != NULL) {

    if (tmpstr == NULL) {
      tmpstr = (char *)malloc(2*sizeof(char));
      tmpstr[0] = '(';
      tmpstr[1] = '\0';

      *str = tmpstr;
    } else {
      int len = strlen(*str);
      tmpstr = (char *)realloc(*str,(len+2)*sizeof(char));

      *str = tmpstr;

      tmpstr[len] = '(';
      tmpstr[len+1] = '\0';
    }

    print_Node_string(node->left,str);

    int len = strlen(*str);

    tmpstr = (char *)realloc(*str,(len+2)*sizeof(char));

    *str = tmpstr;

    tmpstr[len] = ',';
    tmpstr[len+1] = '\0';

    print_Node_string(node->right,str);

    if (node->branch_length > 0) {
      
      int len = strlen(*str);

      tmpstr = (char *)realloc(*str,(len+2)*sizeof(char));
      *str = tmpstr;
    
      tmpstr[len] = ')';
      tmpstr[len+1] = '\0';
    }
  } else {

    int len1 = strlen(*str);               // original string
    int len2 = strlen(node->name);        // node id

    tmpstr = (char *)realloc(*str,((len1+len2+1)*sizeof(char)));
    
    *str = tmpstr;

    strcat(tmpstr,node->name);
  }

  int len1 = strlen(*str);
  int len3 = 9;                         // branch_length

  if (node->isRoot == 0) {
    
    tmpstr = (char *)realloc(*str,((len1+len3+1)*sizeof(char)));
    *str = tmpstr;
    
    strcat(tmpstr,":");
    
    // Now the branch length
    
    char *dist;
    
    dist = (char *)malloc(9*sizeof(char));
    
    sprintf(dist,"%8.6f",node->branch_length);
    
    dist[8] = '\0';
    
    strcat(tmpstr,dist);

  } else {
    tmpstr = (char *)realloc(*str,(len1+3)*sizeof(char));

    *str = tmpstr;

    tmpstr[len1] = ')';
    tmpstr[len1+1] = ';';
    tmpstr[len1+2] = '\0';
  }
}

Node *delete_Leaf(Node *root,Node *node) {

  Node *parent;
  Node *grand_parent;
  Node *othernode;

  printf("Deleting %s\n",node->name);

  parent = node->parent;

  if (parent->parent != NULL) {
    // We're not at root node

    if (strcmp(node->name,parent->left->name) == 0) {
      othernode = parent->right;
    } else {
      othernode = parent->left;
    }
    

    float len = othernode->branch_length + parent->branch_length;

    //printf("Other node is %s\t%8.6f\n",othernode->name,len);

    othernode->branch_length = len;

    grand_parent = parent->parent;

    if (grand_parent->left == parent) {
      grand_parent->left = othernode;
    } else {
      grand_parent->right = othernode;
    }

    othernode->parent = grand_parent;

    return root;
  } else {

    parent = node->parent;
    //    printf("Root node %s\n",node->name);
    //printf("Parent %s\n",parent->name);

    // We're at the root node so reset root to othernode
    if (strcmp(node->name,parent->left->name) == 0) {
      othernode = parent->right;
    } else {
      othernode = parent->left;
    }

    othernode->branch_length = 0;
    othernode->isRoot = 1;
    othernode->parent = NULL;

    return othernode;
  }
}

Node *find_Leaf(Node *node,char *name,Node *found) {

  if (node->left != NULL &&
      node->right != NULL) {

    found =  find_Leaf(node->left,name,found);
    found = find_Leaf(node->right,name,found);
  } else {

    //    printf("Checking %s %s\n",node->name,name);

    if (strcmp(node->name,name) == 0) {
      //printf("Returning %d %s\n",cmp,node->name);
      found = node;
    }
  }
  return found;
}
void count_Leaves(Node *node,int *countptr) {

  if (node->left != NULL &&
      node->right != NULL) {

    count_Leaves(node->left,countptr);
    count_Leaves(node->right,countptr);

  } else {

    *countptr += 1;

  }
}

 void find_total_height(Node *node, float *height) {
   
   if (node->left != NULL  &&
       node->right != NULL) {
     find_total_height(node->left,height);
     find_total_height(node->right,height);
   }
   
   *height += node->branch_length;
 }


    
void draw_Tree(Node *topnode,int width,int chars_per_ylevel) {
  
  float dist_per_char;
  int   height;
  float tree_height;
  int   leaves;

  tree_height = 0;
  find_max_height(topnode,&tree_height);

  dist_per_char = (tree_height + 0.1*tree_height)/width;

  leaves = 0;
  count_Leaves(topnode,&leaves);

  leaves++;

  char *boxes;

  boxes = (char *)malloc(width*leaves*chars_per_ylevel*sizeof(char));

  int j = 0;

  while (j < leaves*chars_per_ylevel) {
    int i = 0;

    while (i < width) {
      char *newboxes;
      newboxes = boxes;

      newboxes += j*width + i;
      *newboxes = ' ';
      i++;
    }
    j++;
  }

  fill_Boxes(topnode,boxes,width,leaves*chars_per_ylevel,dist_per_char,chars_per_ylevel);

  print_Boxes(boxes,width,leaves*chars_per_ylevel,dist_per_char);
}

void print_Boxes(char *boxes, int width, int height,float dist_per_char) {
  int j = 0;

  while (j < height) {
    int i = 0;

    while (i < width) {
      char *newboxes;
      newboxes = boxes;

      newboxes += j*width + i;
      printf("%c",*newboxes);
      i++;
    }
    printf("\n");

    j++;
  }

  printf("\n%f subs/site per character\n",dist_per_char);
  printf("\n|--------| (10 chars) = %f subs/site\n",10*dist_per_char);

  int nchars = (int)(0.1/dist_per_char);

  int i = 0;
  printf("\n");
  while (i < nchars) {
    if (i == 0 || i == nchars-1) {
      printf("|");
    } else {
      printf("-");
    }
    i++;
  }
  printf(" = 0.1 subs/site\n\n");

}

void fill_Boxes(Node *node, char *boxes,int width, int height,float  dist_per_char,int chars_per_ylevel) {

  float startx_f = node->dist_from_root;
  int   startx   = (int)(startx_f/dist_per_char);

  //printf("Filling for %s\t%f\t%d\n",node->name,node->dist_from_root,node->ycount);

  if (node->left != NULL) {
    
    // First draw the vertical line

    int   starty_f = node->ycount;
    int   starty   = node->ycount*chars_per_ylevel;
    int   endy     = node->left->ycount * chars_per_ylevel;
    
    if (endy < starty) {
      int tmp = endy;
      endy = starty;
      starty = tmp;
    }
    
    int y = starty;
    
    //printf("for left node %s starty endy %d %d x %d\n",node->name,starty,endy,startx);
    while (y <= endy) {
      char *newboxes;
      newboxes = boxes;
      
      newboxes += y*width + startx;
      *newboxes = '|';
      //boxes[startx][i] = '|';
      y++;
    }
   
    //print_Boxes(boxes,width,height);
    // Now the across piece
    
    starty_f = node->ycount;
    starty   = node->ycount*chars_per_ylevel;
    endy     = node->left->ycount * chars_per_ylevel;
    int   endx     = (int)(node->left->dist_from_root/dist_per_char);
    
    int x = startx;
    
    while (x <= endx) {
      char *newboxes;
      newboxes = boxes;
      
      newboxes += endy*width + x;
      *newboxes = '-';
      
      //boxes[i][endy] = "-";
      x++;
    }
    //print_Boxes(boxes,width,height);
    fill_Boxes(node->left,boxes,width,height,dist_per_char,chars_per_ylevel);
  } 
  
  if (node->right != NULL) {
    
    // First draw the vertical line

    int starty_f = node->ycount;
    int starty   = node->ycount*chars_per_ylevel;
    int endy     = node->right->ycount * chars_per_ylevel;
    
    if (endy < starty) {
      int tmp = endy;
      endy = starty;
      starty = tmp;
    }
    
    int y = starty; 
    //printf("for right node %s starty endy %d %d x %d\n",node->name,starty,endy,startx);   
    while (y <= endy) {
      char *newboxes;
      newboxes = boxes;
      
      newboxes += y*width + startx;
      *newboxes = '|';
      
      //boxes[startx][i] = "|";
      y++;
    }
    //    print_Boxes(boxes,width,height);
    // Now the across piece

    starty_f = node->ycount;
    starty   = node->ycount*chars_per_ylevel;
    endy     = node->right->ycount * chars_per_ylevel;
    int   endx     = (int)(node->right->dist_from_root/dist_per_char);

    //printf("for right node %s starty endy %d %d x %d %d\n",node->name,starty,endy,startx,endx);    
    int x = startx;
    
    while (x <= endx) {
      char *newboxes;
      newboxes = boxes;
      
      newboxes += endy*width + x;
      *newboxes = '-';
      
      //boxes[i][endy] = "-";
      x++;
    }
    //    print_Boxes(boxes,width,height);
    fill_Boxes(node->right,boxes,width,height,dist_per_char,chars_per_ylevel);
  } 
  
  if (node->left == NULL &&
      node->right == NULL) {
    
    // print the leaf
    
    int startx = (int)(node->dist_from_root/dist_per_char);
    int starty = node->ycount*chars_per_ylevel;
    
    int x = startx;
    
    while (x-startx < strlen(node->name)) {
      char *newboxes;
      newboxes = boxes;
      
      newboxes += width*starty + x;
      *newboxes = (char)(node->name[x-startx]);
      //boxes[i][starty] = "X";
      x++;
    }

  }
  //print_Boxes(boxes,width,height);

}
void  find_max_height(Node *node, float *height) {
  
   if (node->left != NULL  &&
       node->right != NULL) {
     find_max_height(node->left,height);
     find_max_height(node->right,height);
   }

   if (node->dist_from_root > *height) {
     *height = node->dist_from_root;
   }
 }

