#include <gtk/gtk.h>

typedef struct _JalColor {
  char     name;
  GdkColor *color;
} JalColor,*JalColor_ptr;

typedef struct _ColorScheme {
  char      *name;
  JalColor **colors;
  int        count;
} ColorScheme,*ColorScheme_ptr ;

ColorScheme *make_ZappoColorScheme(void);
GdkColor    *get_Zappo_Color_by_Char(ColorScheme *zappo, char name);
GdkGC       *get_Zappo_GC_by_Char(GtkWidget *widget,ColorScheme *zappo, char name);
ColorScheme *make_ColorScheme(char *name);
GdkColor    *make_GdkColor_from_RGB(int red, int green, int blue);
JalColor    *make_Color(char name, GdkColor *gdkcolor);
void         add_Color(ColorScheme *scheme,JalColor *color);
