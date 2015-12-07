#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <gtk/gtk.h>

#include "colour.h"

ColorScheme *make_ZappoColorScheme(){
    JalColor **cptr;

    GdkColor *gyy     = make_GdkColor_from_RGB(204,255,0);
    GdkColor *bbb     = make_GdkColor_from_RGB(0,0,255);
    GdkColor *brb     = make_GdkColor_from_RGB(204,0,255);
    GdkColor *rrr     = make_GdkColor_from_RGB(255,0,0);
    GdkColor *yyy     = make_GdkColor_from_RGB(255,255,0);
    GdkColor *rbr     = make_GdkColor_from_RGB(255,0,204);
    GdkColor *brr     = make_GdkColor_from_RGB(255,0,102);
    GdkColor *yry     = make_GdkColor_from_RGB(255,153,0);
    GdkColor *gbb     = make_GdkColor_from_RGB(0,102,255);
    GdkColor *gyg     = make_GdkColor_from_RGB(102,255,0);
    GdkColor *ygg     = make_GdkColor_from_RGB(51,255,0);
    GdkColor *rbb     = make_GdkColor_from_RGB(102,0,255);
    GdkColor *ggg     = make_GdkColor_from_RGB(0,255,0);
    GdkColor *ryy     = make_GdkColor_from_RGB(255,204,0);
    GdkColor *yrr     = make_GdkColor_from_RGB(255,51,0);
    GdkColor *ryr     = make_GdkColor_from_RGB(255,102,0);
    GdkColor *bgg     = make_GdkColor_from_RGB(0,204,255);
    GdkColor *gbg     = make_GdkColor_from_RGB(0,255,204);
    GdkColor *ygy     = make_GdkColor_from_RGB(153,255,0);
    GdkColor *white   = make_GdkColor_from_RGB(255,255,255);
    
    ColorScheme * zappo  = make_ColorScheme("Zappo");
  
    JalColor *colA = make_Color('A',gyy);  
    JalColor *colR = make_Color('R',bbb);
    JalColor *colN = make_Color('N',brb);
    JalColor *colD = make_Color('D',rrr);
    JalColor *colC = make_Color('C',yyy);
    JalColor *colQ = make_Color('Q',rbr);
    JalColor *colE = make_Color('E',brr);
    JalColor *colG = make_Color('G',yry);
    JalColor *colH = make_Color('H',gbb);
    JalColor *colI = make_Color('I',gyg);
    JalColor *colJ = make_Color('J',white);
    JalColor *colL = make_Color('L',ygg);  
    JalColor *colK = make_Color('K',rbb);
    JalColor *colM = make_Color('M',ggg);
    JalColor *colF = make_Color('F',bgg);
    JalColor *colP = make_Color('P',ryy);
    JalColor *colS = make_Color('S',yrr);
    JalColor *colT = make_Color('T',ryr);
    JalColor *colW = make_Color('W',bgg);
    JalColor *colY = make_Color('Y',gbg);
    JalColor *colV = make_Color('V',ygy);
    JalColor *colspace = make_Color('-',white);
    JalColor *colX = make_Color('X',white);

    add_Color(zappo,colA);
    add_Color(zappo,colR);
    add_Color(zappo,colN);
    add_Color(zappo,colD);
    add_Color(zappo,colC);
    add_Color(zappo,colQ);
    add_Color(zappo,colE);
    add_Color(zappo,colG);
    add_Color(zappo,colH);
    add_Color(zappo,colI);
    add_Color(zappo,colJ);
    add_Color(zappo,colL);
    add_Color(zappo,colK);
    add_Color(zappo,colM);
    add_Color(zappo,colF);
    add_Color(zappo,colP);
    add_Color(zappo,colS);
    add_Color(zappo,colT);
    add_Color(zappo,colW);
    add_Color(zappo,colY);
    add_Color(zappo,colV);
    add_Color(zappo,colX);
    add_Color(zappo,colspace);

    cptr = zappo->colors;

    return zappo;

}

GdkColor *get_Zappo_Color_by_Char(ColorScheme *zappo, char name) {
  char  tmp;
  int   i = 0;
  
  GdkColor    *color = NULL;
  JalColor   **cptr  = zappo->colors;
  
  JalColor    *tmpcolor;
  
  while (i < zappo->count) {
    tmpcolor = *cptr;
    tmp      = tmpcolor->name;
    
    if (name == tmp) {
      color = tmpcolor->color;
      i = zappo->count;
      break;
    }
    
    cptr++;
    i++;
  }

  return color;

}

ColorScheme *make_ColorScheme(char *name) {

  ColorScheme *colorScheme = (ColorScheme *)malloc(sizeof(ColorScheme));
  
  if (colorScheme == NULL) {
    printf("Failed to allocate memory for colorScheme %s\n",name);
    exit(0);
  }
  
  colorScheme->name    = name;
  colorScheme->count   = 0;
  
  return colorScheme;
}

GdkColor *make_GdkColor_from_RGB(int red, int green, int blue) {

  GdkColor *gdkcolor = (GdkColor *)malloc(sizeof(GdkColor));
  
  gdkcolor->red   = red   * 65535/255;
  gdkcolor->green = green * 65535/255;
  gdkcolor->blue  = blue  * 65535/255;
  
  gdkcolor->pixel = (gulong)(red*65536 + green*256 + blue);

  return gdkcolor;
}

  
JalColor *make_Color(char name, GdkColor *gdkcolor) {

  JalColor *color = (JalColor *)malloc(sizeof(JalColor));
  
  if (color == NULL) {
    printf("Failed to allocate memory for color %c\n",name);
    exit(0);
  }

  color->name  = name;
  color->color = gdkcolor;

  return color;

}

void add_Color(ColorScheme *scheme,JalColor *color) {
  JalColor **cptr = scheme->colors;
  JalColor **nptr;
  JalColor **tptr = cptr;

  int i     = 0;
  int count = scheme->count;

  nptr = (JalColor **)malloc((count+1)*sizeof(JalColor*));

  if (nptr == NULL) {
    printf("Failed to allocate memory for color scheme %s\n",scheme->name);
    exit(0);
  }

  scheme->colors = nptr;

  tptr = nptr;

  i = 0;

  while (i < count) {
    *tptr = *cptr;

    tptr++;
    cptr++;

    i++;
  }

  *tptr = color;
  count++;
  scheme->count = count;

}
  

