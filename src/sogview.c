#include <stdlib.h>
#include <gtk/gtk.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>

#include "datamodel.h"
#include "colour.h"
#include "viewport.h"
#include "menu.h"
#include "seqio.h"

static void draw_align_area(GtkWidget *w);
static void zoom_in(void);
static void zoom_out(void);
static void zoom_in_x(void);
static void zoom_out_x(void);
static void open_file(void);
static void save_file(void);
static void toggle_consensus_color(void);
static void toggle_block_color(void);
static void first_sequence_color(void);
static void toggle_hash_color(void);
static void get_main_menu( GtkWidget  *window,
			   GtkWidget **menubar );
static void draw_region(int row_start, int row_end);
static void draw_id(int row,int position);
static void find_hash(Sequence *head,int n,Viewport *viewport);



extern const signed char fasta_encoding[256];

#define _ (-1)
#define B(n) ((1ULL<<(n*2))-1)

const signed char fasta_encoding[] = {
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, 0, _, 1, _, _, _, 2, _, _, _, _, _, _, _, _,
 _, _, _, _, 3, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
};



static unsigned int map_comb(Comb *comb,unsigned long long i, int print);
GdkColor *get_hash_color(Sequence *seq,int pos,Viewport *viewport);

static GtkItemFactoryEntry menu_items[] = {
  { "/File",         NULL,         NULL, 0, "<Branch>" },
  { "/File/_New",     "<control>N", NULL, 0, NULL },
  { "/File/_Open",    "<control>O", open_file, 0, NULL },
  { "/File/_Save",    "<control>S", NULL, 0, NULL },
  { "/File/Save _As", NULL,         save_file, 0, NULL },
  { "/File/sep1",     NULL,         NULL, 0, "<Separator>" },
  { "/File/Quit",     "<control>Q", gtk_main_quit, 0, NULL },
  { "/View",         NULL,         NULL, 0, "<Branch>" },
  { "/View/Zoom _in",  "<control>I",zoom_in, 0, NULL },
  { "/View/Zoom in X","<control>F",zoom_in_x, 0, NULL },
  { "/View/Zoom _out", "<control>U",zoom_out,0 , NULL },
  { "/View/Zoom out X", "<control>T",zoom_out_x,0 , NULL },
  { "/Color/By _Consensus","<control>C",toggle_consensus_color,0,NULL},
  { "/Color/By _Block","<control>B",toggle_block_color,0,NULL},
  { "/Color/By _First Sequence","<control>A",first_sequence_color,0,NULL},
  { "/Calculate/Hash","<control>X",toggle_hash_color,0,NULL},
  { "/Help",         NULL,         NULL, 0, "<LastBranch>" },
  { "/Help/About",   NULL,         NULL, 0, NULL },
};

void get_main_menu( GtkWidget  *window,
                    GtkWidget **menubar )
{
  GtkItemFactory *item_factory;
  GtkAccelGroup  *accel_group;

  gint nmenu_items = sizeof (menu_items) / sizeof (menu_items[0]);

  accel_group = gtk_accel_group_new ();

  item_factory = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<main>", 
                                       accel_group);


  gtk_item_factory_create_items (item_factory, nmenu_items, menu_items, NULL);

  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);

  if (menubar) {
    *menubar = gtk_item_factory_get_widget (item_factory, "<main>");
  }
}

static GtkWidget     *window;

static GdkPixmap     *align_pixmap    = NULL;
static GdkPixmap     *id_pixmap       = NULL;
static GdkPixmap     *scale_pixmap    = NULL;

static GtkWidget     *align_area;
static GtkWidget     *id_area;
static GtkWidget     *scale_area;
static GtkWidget     *id_scale_area;

static GtkWidget     *status_bar;

static GtkWidget     *hscrollbar;
static GtkWidget     *vscrollbar;

static GtkAdjustment *hadjust;
static GtkAdjustment *vadjust;

static GtkWidget     *file_dialog;

static GdkFont       *fixed_font;
static GdkColor      *color;

static GdkGC         *gc;
static ColorScheme   *scheme;


/* ----------------------------------------- */

static Sequence      *head;
static Viewport      *viewport;

int    align_height;

static char          *residue;



static void draw_residue(GtkWidget *widget,
			 int    x,
			 int    y,
			 int        xposition,
			 int        yposition)
{
  int           offset;
  int           tmppos;

  Sequence      *seq;

  Block         **blocks;
  Block          *block;

  seq = Sequence_get_Sequence_at(head,yposition,align_height);

  blocks = viewport->blocks;

  offset = (int)hadjust->value;

  if (residue == NULL) {
    residue = (char *)malloc(2*(sizeof(char)));
    
    if (residue == NULL) {
      printf("Failed to allocate memory for residue\n");
      exit(0);
    }
    
    residue[1] = '\0';
  }

  if (seq == NULL) {
    gdk_draw_rectangle (align_pixmap,
			widget->style->black_gc,
			TRUE,
			x, y,
			viewport->char_width, viewport->char_height);
    
    return;
  } else if (seq != NULL && seq->length <= offset+xposition) {
    gdk_draw_rectangle (align_pixmap,
			widget->style->white_gc,
			TRUE,
			x, y,
			viewport->char_width, viewport->char_height);
    return;
  }
  
  residue[0] = Sequence_get_char_at(seq,offset+xposition);
  
  if (fixed_font == NULL) {
    fixed_font = gdk_font_load ("-misc-fixed-medium-r-*-*-*-100-*-*-*-*-*-*");
  }
  
  if (color == NULL) {
    color = (GdkColor *)malloc(sizeof(GdkColor));
    
    if (color == NULL) {
      printf("Failed to allocate memory for colour\n");
      exit(0);
    }
  }
  
  if (gc == NULL) {
    gc = gdk_gc_new(widget->window);
  }
  
  if (viewport->first_sequence_color) {
    
    if(residue[0] == head->sequence[offset+xposition]) {
      color = get_Zappo_Color_by_Char(scheme,(char)residue[0]); 
    } else{
      color = get_Zappo_Color_by_Char(scheme,'-');
    }
  } else if (viewport->color_by_hash && viewport->hash != NULL) {
    
    color = get_hash_color(seq,offset+xposition,viewport);
    
  } else if (viewport->color_by_consensus && viewport->consensus != NULL) {
    if(residue[0] == viewport->consensus->sequence[offset+xposition]) {
      color = get_Zappo_Color_by_Char(scheme,residue[0]);
    } else{
      color = get_Zappo_Color_by_Char(scheme,'-');
    }
  } else if (viewport->color_by_blocks) {
    
    if (residue[0] != '-') {
      tmppos = offset + xposition - 100;
      
      if (tmppos < 0) {
        tmppos = 0;
      }
      
      if (tmppos > Sequence_max_sequence_length(head)) {
        tmppos = Sequence_max_sequence_length(head);
      }

      if (blocks[offset + xposition] != NULL) {

        block = blocks[offset + xposition];

        if (block->length == 1) {
          color = get_Zappo_Color_by_Char(scheme,'K');
        } else if (block->length > 20) {
          color = get_Zappo_Color_by_Char(scheme,'C');
        } else if (block->length > 10) {
          color = get_Zappo_Color_by_Char(scheme,'D');
        } else if (block->length > 5) {
          color = get_Zappo_Color_by_Char(scheme,'F');
        } else if (block->length > 3) {
          color = get_Zappo_Color_by_Char(scheme,'I');
        } else if (block->length > 2) {
          color = get_Zappo_Color_by_Char(scheme,'G');
        } else if (block->length > 1) {
          color = get_Zappo_Color_by_Char(scheme,'T');
        }
      } else {
        color = get_Zappo_Color_by_Char(scheme,'-');
      }
    } else {
      color = get_Zappo_Color_by_Char(scheme,'-');
    }
  } else {
    color = get_Zappo_Color_by_Char(scheme,residue[0]);
  }
  
  if (color != NULL) {
    gdk_gc_set_foreground(gc,color);
  } else {
    color = get_Zappo_Color_by_Char(scheme,'-');
    gdk_gc_set_foreground(gc,color);
    
  }
  
    
  gdk_draw_rectangle (align_pixmap,
		      gc,
		      TRUE,
		      x, y,
		      viewport->char_width, viewport->char_height);
  
  
  if (viewport->char_width > 5) {
    gdk_draw_text(align_pixmap,
		  fixed_font,
		  widget->style->black_gc,
		  x + (viewport->char_width - gdk_char_width(fixed_font,'W'))/2,
		  y + (viewport->char_height + gdk_char_height(fixed_font,'W'))/2,
		  residue,
		  1);
  }
  
}


static GdkColor *get_hash_color(Sequence *seq,int pos,Viewport *viewport) {
  int e;
  int i;
  int start;

  int m = viewport->comb->m;

  unsigned long long ei = 0;
  unsigned int ecode;
  int count;

  double av = 1.0*Sequence_max_sequence_length(head)/(1.0*nmer2pow(m));
  
  GdkColor *color;

  if (pos < m/2) {
    return NULL;
  }
  if (pos > seq->length - m/2) {
    return NULL;
  }

  start = pos - m/2;
  i = start;

  while (i < start + m) {
    e = fasta_encoding[(int)seq->sequence[i]];
    
    if (e < 0) {
      return NULL;
    }
    ei = ( ei << 2) | e;
    i++;
  }

  ecode = map_comb(viewport->comb,ei,0);

  if (viewport->hash[ecode] != NULL) {
    count = viewport->hash[ecode]->count;

    if (count > 5*av) {
      color = get_Zappo_Color_by_Char(scheme,'K'); //purple
    } else if (count > 4*av) {
      color = get_Zappo_Color_by_Char(scheme,'N'); //purle
    } else if (count > 3*av) {
      color = get_Zappo_Color_by_Char(scheme,'E'); //reddy pink
    } else if (count > 2*av) {
      color = get_Zappo_Color_by_Char(scheme,'S'); //red
    } else if (count > av) {
      color = get_Zappo_Color_by_Char(scheme,'T'); //orange
    } else if (count > av/2) {
      color = get_Zappo_Color_by_Char(scheme,'P'); //lightorange
    } else if (count > av/3) {
      color = get_Zappo_Color_by_Char(scheme,'C');//yellow
    } else {
      color = get_Zappo_Color_by_Char(scheme,'A'); //light green
    }

  } else {
    color = get_Zappo_Color_by_Char(scheme,'-');
  }

  return color;
}
    


static void id_area_event(GtkWidget *widget) {

  int j;
  int row;

  int vvalue = (int)vadjust->value;
  int yend   = vvalue + viewport->height/viewport->char_height + 1;

  // Draw background
  gdk_draw_rectangle (id_pixmap,
		      widget->style->white_gc,
                      TRUE,
                      0,0,
                      widget->allocation.width, widget->allocation.height);
  

  // Draw all ids
  for (j = vvalue; j < yend; j++) {
    row = j % (align_height +2);
    
    if (row < align_height) {
      draw_id(row,j-vvalue);
    }
  }

  // Copy pixmap to window
  gdk_draw_drawable (widget->window,
		     widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		     id_pixmap,
		     0, 0,
		     0, 0,
		     widget->allocation.width,widget->allocation.height);
  
}

static void scale_area_event(GtkWidget *widget) {

  int start;
  int i;
  int hvalue;
  
  gdk_draw_rectangle (scale_pixmap,
		      widget->style->white_gc,
                      TRUE,
                      0,0,
                      widget->allocation.width, widget->allocation.height);

  hvalue = (int)hadjust->value;

  start = (hvalue + 10)%10 * viewport->char_width;
  
  for (i = start; i < widget->allocation.width/viewport->char_width; i+=viewport->char_width) {
    gdk_draw_text(scale_pixmap,
		  fixed_font,
		  widget->style->black_gc,
		  id_area->allocation.width + 3+ (i-hvalue)*viewport->char_width,
		  15,
		  "|",
		  1);
  }

  gdk_draw_drawable (widget->window,
		     widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		     scale_pixmap,
		     0, 0,
		     0, 0,
		     widget->allocation.width,widget->allocation.height);
  
}
    

int pixel_to_sequence_number(int n) {

  int seqnum;

  seqnum = vadjust->value + (int)(n/viewport->char_height);

  return seqnum;

}

int pixel_to_char_number(int n) {
  int charnum;
  
  charnum = hadjust->value + (int)(n/viewport->char_width);
  
  return charnum;
}

void button_press_event(GtkWidget *widget,
			GdkEventButton *event) {
  int seqnum;
  int charnum;

  Sequence *seq;
  char *status_message;
  char c;

  seqnum  = pixel_to_sequence_number(event->y);
  charnum = pixel_to_char_number(event->x);
  
  seq = Sequence_get_Sequence_at(head,seqnum,align_height);

  if (seq == NULL || charnum >= seq->length) {
    gtk_statusbar_push(GTK_STATUSBAR(status_bar),
		       0,
		       "");
    return;
  }
  
  status_message = (char *)malloc((strlen(seq->id) + 20 + 29)*sizeof(char));

  if (status_message == NULL) {
    printf("Couldn't allocate memory for status message\n");
    return;
  }
  
  c = seq->sequence[charnum];
  
  sprintf(status_message,"Id: %s Position : %d Sequence : %d",seq->id,charnum+1,seqnum+1);
  gtk_statusbar_push(GTK_STATUSBAR(status_bar),
                     0,
                     status_message);
  
  free(status_message);
}

  
void hscroll_event(GtkWidget *widget,
		   GdkEventScroll *event) {

  hadjust->value = (int)hadjust->value;
  
  if (hadjust->value != viewport->old_hvalue) {
    draw_align_area(align_area);
    id_area_event   (id_area);
    scale_area_event(scale_area);

    viewport->old_hvalue = (int)hadjust->value;
  }

}

void vscroll_event(GtkWidget *widget,
	     GdkEventScroll *event) 
{

  vadjust->value = (int)vadjust->value;

  if (vadjust->value != viewport->old_vvalue) {
    draw_align_area(align_area);
    id_area_event(id_area);
    scale_area_event(scale_area);
    viewport->old_vvalue = (int)vadjust->value;
  }

}


/* This is some sort of initialization routine - not
   sure where it is called */
static gboolean align_configure_event( GtkWidget         *widget,
				       GdkEventConfigure *event )
{

  int i;
  int j;
  
  int yend;
  //int align_height = Sequence_count_Sequences(head);

  int value = (int)hadjust->value;
  viewport->old_hvalue = value;

  if (align_pixmap) {
    g_object_unref(align_pixmap);
  }
  
  align_pixmap = gdk_pixmap_new (widget->window,
                           widget->allocation.width,
                           widget->allocation.height,
                           -1);

  viewport->width  = widget->allocation.width;
  viewport->height = widget->allocation.height;

  gdk_draw_rectangle (align_pixmap,
		      widget->style->white_gc,
		      TRUE,
		      0, 0,
		      widget->allocation.width,
		      widget->allocation.height);

  yend = event->height/viewport->char_height;

  if ((yend+(int)vadjust->value) > align_height) {
    yend = align_height-(int)vadjust->value;
  }
  
  for (i = 0; i <= (int)(event->width/viewport->char_width+1); i++) {
    for (j = 0; j <= yend; j++) {
      draw_residue(widget,
		   (gdouble)((i) * viewport->char_width),
		   (gdouble)((j) * viewport->char_height),
		   i,j);
    }
  }
  
  return TRUE;
}
static gboolean scale_configure_event(GtkWidget *widget,
				      GdkEventConfigure *event)
{

  if (scale_pixmap) {
    g_object_unref(scale_pixmap);
  }
  scale_pixmap = gdk_pixmap_new (widget->window,
				 widget->allocation.width,
				 widget->allocation.height,
				 -1);

  gdk_draw_rectangle (scale_pixmap,
		      widget->style->black_gc,
		      TRUE,
		      0, 0,
		      widget->allocation.width,
		      widget->allocation.height);
    
  return TRUE;
}

static gboolean id_configure_event( GtkWidget         *widget,
				    GdkEventConfigure *event )
{

  if (id_pixmap) {
    g_object_unref(id_pixmap);
  }
  
  id_pixmap = gdk_pixmap_new (widget->window,
                           widget->allocation.width,
                           widget->allocation.height,
                           -1);

  gdk_draw_rectangle (id_pixmap,
		      widget->style->black_gc,
		      TRUE,
		      0, 0,
		      widget->allocation.width,
		      widget->allocation.height);
    
  return TRUE;
}

static void draw_align_area(GtkWidget *widget) {

  int nrows = viewport->height/viewport->char_height;
  int vdiff;

  vdiff = (int)vadjust->value - viewport->old_vvalue;

  if (vdiff < 0) {
    vdiff *= -1;

    if (vdiff < nrows) {
      // Copy area down vdiff chars
      gdk_draw_drawable (align_pixmap,
			 widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
			 align_pixmap,
			 0,0,
			 0,viewport->char_height*vdiff,
			 viewport->width,viewport->height-viewport->char_height*vdiff);
      
      draw_region(0,vdiff);
    } else {
      draw_region(0,nrows);
    }
    
  } else if (vdiff > 0) {
    if (vdiff < nrows) {
      // Copy region up vdiff chars
      gdk_draw_drawable (align_pixmap,
			 widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
			 align_pixmap,
			 0,viewport->char_height*vdiff,
			 0,0,
			 viewport->width,viewport->height-viewport->char_height*vdiff);
      
      draw_region(nrows-vdiff,nrows);
    } else {
      draw_region(0,nrows);
    }
  } else {
    // Draw the whole thing
    gdk_draw_rectangle (align_pixmap,
			widget->style->white_gc,
			TRUE,
			0,0,
			widget->allocation.width, widget->allocation.height);
    
    draw_region(0,nrows);
  }

  // Copy the pixmaps to the window
  gdk_draw_drawable (widget->window,
		     widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		     align_pixmap,
		     0, 0,
		     0, 0,
		     viewport->width,viewport->height);
  
  viewport->old_vvalue = (int)vadjust->value;
}

void draw_id(row,position) {
  Sequence *seq;
  char     *id;

  seq = Sequence_get_Sequence_at(head,row,align_height);
  id  = seq->id;

  gdk_draw_text(id_pixmap,
		fixed_font,
		id_area->style->black_gc,
		2,
		position*viewport->char_height +(viewport->char_height+gdk_char_height(fixed_font,'W'))/2,
		id,
		strlen(id));
}

static void draw_region(int row_start, int row_end) {

  int offset = (int)vadjust->value;

  int i;
  int j;

  int height          = (row_end - row_start + 1)*viewport->char_height;

  //int align_height    = Sequence_count_Sequences(head);
  int chars_per_row   = viewport->width/viewport->char_width;
  
  int rows_per_chunk  = align_height + 2;

  int chunkstart      = (offset + row_start)/rows_per_chunk;

  int xstart          = chunkstart*chars_per_row;
  int ystart          = (offset + row_start) % rows_per_chunk;

  int toprag          = (align_height + 2 - ystart)%(align_height+2);
  
  int nchunks         = (height/viewport->char_height - toprag)/rows_per_chunk;

  int endrag          = (row_end - row_start + 1) - rows_per_chunk * nchunks - toprag;

  int row    = 0;
  int xcount = xstart;
  int rows   = row_start;
  int chunk;

  //printf("Toprag %d, nchunks %d, endrag %d ystart %d rowstart %d\n",toprag,nchunks,endrag,ystart,row_start);

  //printf("Xstart %d chunkstart %d chars_per_row %d\n",xstart,chunkstart,chars_per_row);
  if (toprag > 0) {
    row = ystart;
    
    while (row < align_height && rows <= row_end) {
      //printf("Drawing toprag %d %d\n",row,rows);
      for (i = 0; i < chars_per_row; i++) {
	draw_residue(align_area,
		     (gdouble)((i)    * viewport->char_width),
		     (gdouble)((rows) * viewport->char_height),
		     xcount + i,
		     row);
      }
      
//      draw_id(row,rows);
      
      rows++;
      row++;
      
    }
    
    if (row >= align_height) {
      //printf("Drawing rect %d %d\n",row,rows);
      gdk_draw_rectangle (align_pixmap,
			  align_area->style->white_gc,
			  TRUE,
			  0,(rows)*viewport->char_height,
			  align_area->allocation.width, 2*viewport->char_height);
    }
    xcount += chars_per_row;
    rows += 2;
    
  }
  
  for (chunk = 0; chunk < nchunks ; chunk++) {
    //printf("Drawsing chunk %d %d\n",chunk,xcount);
    for (j = 0; j < align_height; j++) {    
      for (i = 0; i < chars_per_row; i++) {
	draw_residue(align_area,
		     (gdouble)((i) * viewport->char_width),
		     (gdouble)((j+rows) * viewport->char_height),
		     i+xcount,j);
      }
//      draw_id(j,rows+j);
    }

    
    xcount += chars_per_row;
    rows   += align_height;

    gdk_draw_rectangle (align_pixmap,
			align_area->style->white_gc,
			TRUE,
			0,(rows)*viewport->char_height,
			align_area->allocation.width, 2*viewport->char_height);

    rows+=2;

  }
  
  j = 0;
    
  while (j < endrag && j < align_height) {
    //printf("Drawing endrag char %d row %d realrow %d\n",xcount,j,rows);    
    for (i = 0; i < chars_per_row; i++) {

      draw_residue(align_area,
		   (gdouble)((i) * viewport->char_width),
		   (gdouble)((rows) * viewport->char_height),
		   i+xcount,
		   j);
      
      
      
    }
//    draw_id(j,rows+j);
    j++;
    rows++;
  }
//  printf("J endrag %d %d\n",j,endrag);

  while (j < endrag) {
    gdk_draw_rectangle (align_pixmap,
			align_area->style->white_gc,
			TRUE,
			0,(offset+rows)*viewport->char_height,
			align_area->allocation.width, viewport->char_height);
    rows++;
    j++;
  }

}


void file_save_sel(GtkWidget *w,
		   GtkFileSelection *fs)
{
  printf("%s\n", gtk_file_selection_get_filename (GTK_FILE_SELECTION (fs)));

  save_fasta((char *)gtk_file_selection_get_filename (GTK_FILE_SELECTION (fs)),head);

  gtk_widget_destroy(file_dialog);
}
void file_open_sel(GtkWidget *w,
		   GtkFileSelection *fs) 
{

  free(head);
  free(viewport->consensus);

  printf("%s\n", gtk_file_selection_get_filename (GTK_FILE_SELECTION (fs)));

  head = read_fasta((char *)gtk_file_selection_get_filename (GTK_FILE_SELECTION (fs)));
		
  
  viewport->consensus = NULL;
  viewport->blocks    = NULL;

  viewport->color_by_consensus = 0;
  viewport->color_by_blocks    = 0;
  viewport->first_sequence_color = 0;

  align_height = Sequence_count_Sequences(head);

  vadjust->upper = (gdouble)(1+Sequence_max_sequence_length(head)/(align_area->allocation.width/viewport->char_width) * (2 + align_height));

  hadjust->upper = (gfloat)Sequence_max_sequence_length(head);

  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);

  gtk_widget_destroy(file_dialog);
}

static void open_file() 
{

  file_dialog = gtk_file_selection_new("Select new fasta file");

  /* Connect the ok_button to file_ok_sel function */
  gtk_signal_connect (GTK_OBJECT (GTK_FILE_SELECTION (file_dialog)->ok_button),
		      "clicked", (GtkSignalFunc) file_open_sel, file_dialog);
    
  /* Connect the cancel_button to destroy the widget */
  gtk_signal_connect_object (GTK_OBJECT (GTK_FILE_SELECTION
					 (file_dialog)->cancel_button),
			     "clicked", (GtkSignalFunc) gtk_widget_destroy,
			     GTK_OBJECT (file_dialog));
    
  /* Lets set the filename, as if this were a save dialog, and we are giving
     a default filename */
  gtk_file_selection_set_filename (GTK_FILE_SELECTION(file_dialog), 
				   "penguin.png");
    
  gtk_widget_show(file_dialog);
}
static void save_file() 
{

  file_dialog = gtk_file_selection_new("Save file");

  /* Connect the ok_button to file_ok_sel function */
  gtk_signal_connect (GTK_OBJECT (GTK_FILE_SELECTION (file_dialog)->ok_button),
		      "clicked", (GtkSignalFunc) file_save_sel, file_dialog);
    
  /* Connect the cancel_button to destroy the widget */
  gtk_signal_connect_object (GTK_OBJECT (GTK_FILE_SELECTION
					 (file_dialog)->cancel_button),
			     "clicked", (GtkSignalFunc) gtk_widget_destroy,
			     GTK_OBJECT (file_dialog));
    
  /* Lets set the filename, as if this were a save dialog, and we are giving
     a default filename */
  gtk_file_selection_set_filename (GTK_FILE_SELECTION(file_dialog), 
				   "penguin.png");
    
  gtk_widget_show(file_dialog);
}
static void first_sequence_color() {
  viewport->first_sequence_color = 1;

  viewport->color_by_blocks = 0;
  viewport->color_by_consensus = 0;
  viewport->color_by_hash        = 0;

  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);
}
static void toggle_block_color() {
  if (viewport->blocks == NULL) {
    viewport->blocks = find_Blocks(head);
  }
  viewport->color_by_blocks = 1;
  viewport->color_by_consensus = 0;
  viewport->first_sequence_color = 0;
  viewport->color_by_hash        = 0;

  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);

  
}

static void toggle_hash_color() {
  if (viewport->hash == NULL) {
    find_hash(head,5,viewport);
  } 
  viewport->color_by_hash   = 1;
  viewport->color_by_blocks = 0;
  viewport->color_by_consensus = 0;
  viewport->first_sequence_color = 0;

  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);

}

static void toggle_consensus_color() {
  viewport->first_sequence_color = 0;
  viewport->color_by_blocks      = 0;
  viewport->color_by_hash        = 0;

  if (viewport->consensus == NULL) {
    viewport->consensus = get_consensus_Sequence(head);
    viewport->consensus->next = NULL;

    //printf("Consensus %s\n",viewport->consensus->sequence);
    Sequence_add_Sequence(head,viewport->consensus);
    align_height++;
  }

  if (viewport->color_by_consensus == 0) {
    viewport->color_by_consensus = 1;
  } else {
    viewport->color_by_consensus = 0;
  }

  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);
}
static void  zoom_in()
{
  viewport->char_width++;
  viewport->char_height++;

  vadjust->upper = (gdouble)(1+Sequence_max_sequence_length(head)/(align_area->allocation.width/viewport->char_width) * (2 + align_height));

  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);

}
static void  zoom_in_x()
{

 viewport->char_width++;

 vadjust->upper = (gdouble)(1+Sequence_max_sequence_length(head)/(align_area->allocation.width/viewport->char_width) * (2 + align_height));
 draw_align_area(align_area);
 id_area_event(id_area);
 scale_area_event(scale_area);

}

static void  zoom_out()
{
  int found = 0;
  if (viewport->char_width > 1) {
    viewport->char_width--;
    found = 1;
  }
  if (viewport->char_height > 1) {
    viewport->char_height--;



    found = 1;
  }

  if (found == 1) {
    vadjust->upper = (gdouble)(1+Sequence_max_sequence_length(head)/(align_area->allocation.width/viewport->char_width) * (2 + align_height));

    draw_align_area(align_area);
    id_area_event(id_area);  
    scale_area_event(scale_area);
  }

}
static void  zoom_out_x()
{

  if (viewport->char_width > 1) {
    viewport->char_width--;

    vadjust->upper = (gdouble)(1+Sequence_max_sequence_length(head)/(align_area->allocation.width/viewport->char_width) * (2 + align_height));

    draw_align_area(align_area);
    id_area_event(id_area);  
    scale_area_event(scale_area);
  }
}

static gboolean expose_event( GtkWidget      *widget,
                              GdkEventExpose *event )
{
  draw_align_area(align_area);
  id_area_event(id_area);
  scale_area_event(scale_area);
  return FALSE;
} 

int nmer2pow(int n) {
  int i=0;
  int val = 1;

  while (i < n*2) {
    val = val * 2;
    i++;
  }
  return val;
}
static int *generate_bitstring(int n, int m) {
  /* Generates a random n of m bitstring */

  int rand_num;
  int i = 0;
  int *bits;

  bits = (int *)malloc(m*sizeof(int));
  
  for (i=0; i < m;i++) {
    bits[i] = 1;
  }

  i=0;

  while (i < (m-n)) {
    rand_num = (int)((double)m*rand()*1.0/(1.0*RAND_MAX));
    
    if (bits[rand_num] != 0) {
      bits[rand_num] = 0;
      i++;
    }
  }

  for(i=0;i< m;i++) {
    printf("%d",bits[i]);
  }

  printf("\n\n");

  return bits;
}

static Comb *bits2comb(int *bits, int n,int m) {
  int i      = 0;
  int length = 0;
  int start  = 0;

  int lshift;

  BitBlock *prev;
  Comb     *comb;

  int numblocks = 0;

  lshift = 0;

  comb = (Comb *)malloc(sizeof(Comb));

  comb->n = n;
  comb->m = m;

  for (i=0; i < m; i++) {
    printf("%d",bits[i]);
  }
  
  printf("\n");
  
  i=0;

  while (i < m) {

    if (i > 0) {
      if (bits[i] == 1 && bits[i-1] == 0) {
	start = i;
	length = 1;
	lshift++;
	numblocks++;
      } else if (bits[i] == 0 && bits[i-1] == 1) {
	if (numblocks == 1) {
	  comb->blocks = (BitBlock *)malloc(sizeof(BitBlock));

	  prev = comb->blocks;
	  prev->next = NULL;
	} else {
	  prev->next = (BitBlock *)malloc(sizeof(BitBlock));

	  prev = prev->next;
	  prev->next = NULL;
	}

	prev->length      = length;

	prev->leftshift   = (m - start - length);
	prev->rightshift  = (prev->leftshift  - n + lshift); 

      } else if (bits[i] == 1 && bits[i-1] == 1) {
	length++;
	lshift++;
      } 
    } else {
      if (bits[i] == 1) {
	numblocks++;
	start = i;
	length = 1;
	lshift++;
      }
    }
      
    i++;
  }

  if (bits[m-1] == 1) {
    if (numblocks == 1) {
      comb->blocks = (BitBlock *)malloc(sizeof(BitBlock));
      
      prev = comb->blocks;
      prev->next = NULL;
    } else {
      prev->next = (BitBlock *)malloc(sizeof(BitBlock));
      
      prev = prev->next;
      prev->next = NULL;
    }
    
    prev->length      = length;
    prev->leftshift   = (m - start - length);
    prev->rightshift  = (prev->leftshift  - n + lshift); 
  }
  
  comb->numblocks = numblocks;

  return comb;
}

static unsigned int map_comb(Comb *comb,unsigned long long i, int print) {

  unsigned int j = 0;
  
  BitBlock *block;
  
  block = comb->blocks;
  
  while (block != NULL) {
    j = j | ((i &(B(block->length) << (2*block->leftshift))) >> (2*block->rightshift));

    block = block->next;
  }
  return j;
}

static void find_hash(Sequence *head,int n,Viewport *viewport) {
  Hash ** hash;
  Hash  * tmphash;
  Comb  * comb;
  int   * bits;
  int     i;

  static unsigned long long ei = 0;
  static unsigned int ecode;

  int e;

  Sequence *seq;

  hash = (Hash **)malloc(nmer2pow(n) * sizeof(Hash *));

  bits = (int *)malloc(sizeof(int)*n);
  bits = generate_bitstring(n,n);

  comb = bits2comb(bits,n,n);

  memset(hash,0,nmer2pow(n)*sizeof(Hash *));

  seq = head;

  while (seq != NULL) {
    i = 0;

    while (i < n-1) {
      
      e = fasta_encoding[(int)seq->sequence[i]];
      ei = (ei << 2) | e;
      i++;
    }
    
    i = n-1;
    
    while (i < seq->length) {
      e = fasta_encoding[(int)seq->sequence[i]];
      
      if (e >=0) {
	ei = (ei << 2) | e;
	ecode = map_comb(comb,ei,0);
	
	tmphash = (Hash *)malloc(sizeof(Hash));
	tmphash->pos = i;
	
	if (hash[ecode] != NULL) {
	  tmphash->count = hash[ecode]->count + 1;
	} else {
	  tmphash->count = 1;
	}
	
	tmphash->next = hash[ecode];
	hash[ecode] = tmphash;
	//print_ll(stdout,ei,m);
	//printf(" Ecode %d ",ecode);
	//print_l(stdout,ecode,n);
	//printf("\n");
      }
      i++;
    }
    seq = seq->next;
  }

  viewport->hash = hash;
  viewport->comb = comb;

  tmphash = hash;

  /*  i=0;
  while (i < nmer2pow(n)) {
    tmphash = hash[i];

    if (tmphash != NULL) {
      printf("HASH %d\n",tmphash->count);
    }
    i++;
    }*/
}
    
void quit() {
  exit(0);
}

int create_window(Sequence *head)
{

  GtkWidget *hbox1;
  GtkWidget *hbox2;
  GtkWidget *vbox;

  GtkWidget *pane;
  GtkWidget *menubar;


  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  
  viewport  = (Viewport *)malloc(sizeof(Viewport));

  viewport->char_width  = 9;
  viewport->char_height = 11;
  viewport->consensus   = NULL;
  viewport->first_sequence_color = 0;
  viewport->color_by_consensus   = 0;

  scheme = make_ZappoColorScheme();

  gtk_widget_set_name (window, "Test Alignment Viewer");
  gtk_widget_set_usize( GTK_WIDGET ( window ) , 800 , 200 );

  g_signal_connect (G_OBJECT (window), "destroy",
                    G_CALLBACK (quit), window);

  get_main_menu (window, &menubar);

  align_area = gtk_drawing_area_new ();
  id_area    = gtk_drawing_area_new ();
  
  scale_area    = gtk_drawing_area_new();
  id_scale_area = gtk_drawing_area_new();

  gtk_drawing_area_size (GTK_DRAWING_AREA (scale_area), 20, 20);

  if (head != NULL) {
    align_height = Sequence_count_Sequences(head);

    hadjust = gtk_adjustment_new((gdouble)0.0,
		                 (gdouble) 0.0,
				 (gdouble)Sequence_max_sequence_length(head),
				 (gdouble)1.0,
				 (gdouble)50.0,
				 (gdouble)(viewport->height*1.0/viewport->char_height));

    vadjust = gtk_adjustment_new((gfloat)0.0,(gfloat)0.0,
				 (gfloat)100*align_height,
				 (gfloat)1.0,
				 (gfloat)viewport->height/viewport->char_height,
				 (gdouble)(viewport->height*1.0/viewport->char_height));
  } else {
    hadjust = gtk_adjustment_new(0.0,0.0,
				 (gfloat)100,
				 1.0,
				 (gfloat)viewport->width/viewport->char_width,
				 (gdouble)(viewport->width*1.0/viewport->char_width));

    vadjust = gtk_adjustment_new((gfloat)0.0,(gfloat)0.0,
				 (gfloat)100,
				 (gfloat)1.0,
				 (gfloat)viewport->height/viewport->char_height,
				 (gfloat)viewport->height/viewport->char_height);
  }
  hscrollbar = gtk_hscrollbar_new(hadjust);
  vscrollbar = gtk_vscrollbar_new(vadjust);

  pane = gtk_hpaned_new ();

  gtk_widget_set_size_request (GTK_WIDGET (id_area), 100, 300);

  gtk_widget_set_double_buffered(align_area,TRUE);
  gtk_widget_set_double_buffered(id_area,TRUE);

  status_bar = gtk_statusbar_new();

  hbox1 = gtk_hbox_new(FALSE,0);
  hbox2 = gtk_hbox_new(FALSE,0);
  vbox  = gtk_vbox_new(FALSE,0);


  gtk_container_add(GTK_CONTAINER (window), vbox);
  
  gtk_box_pack_start(GTK_BOX(hbox1),
		     align_area,
		     TRUE,
		     TRUE,
		     0);

  gtk_box_pack_start(GTK_BOX(hbox1),
		     vscrollbar,
		     FALSE,
		     TRUE,
		     0);

  gtk_paned_add1(GTK_PANED(pane),id_area);
  gtk_paned_add2(GTK_PANED(pane),hbox1);

  gtk_box_pack_start(GTK_BOX(vbox),
		     menubar,
		     FALSE,
		     TRUE,
		     0);

  gtk_box_pack_start(GTK_BOX(vbox),
		     scale_area,
		     FALSE,
		     TRUE,
		     0);

  gtk_box_pack_start(GTK_BOX(vbox),
		     pane,
		     TRUE,
		     TRUE,
		     0);

  gtk_box_pack_start(GTK_BOX(vbox),
		     hscrollbar,
		     FALSE,
		     TRUE,
		     0);

  gtk_box_pack_start(GTK_BOX(vbox),
		     status_bar,
		     FALSE,
		     TRUE,
		     3);

  gtk_widget_show(status_bar);
  gtk_widget_show(align_area);
  gtk_widget_show(id_area);
  gtk_widget_show(scale_area);
  gtk_widget_show(hscrollbar);
  gtk_widget_show(vscrollbar);

  gtk_widget_show(pane);
  gtk_widget_show(hbox1);
  gtk_widget_show(vbox);
  gtk_widget_show(menubar);

  /* Signals used to handle backing pixmap */

  g_signal_connect (G_OBJECT (align_area), "expose_event",
		    G_CALLBACK (expose_event), NULL);
  g_signal_connect (G_OBJECT (align_area),"configure_event",
                    G_CALLBACK (align_configure_event), NULL);
  g_signal_connect (G_OBJECT (id_area), "expose_event",
		    G_CALLBACK (id_area_event), NULL);
  g_signal_connect (G_OBJECT (id_area),"configure_event",
                    G_CALLBACK (id_configure_event), NULL);
  g_signal_connect (G_OBJECT (scale_area),"configure_event",
                    G_CALLBACK (scale_configure_event), NULL);
  g_signal_connect (G_OBJECT (scale_area), "expose_event",
		    G_CALLBACK (scale_area_event), NULL);
  g_signal_connect (G_OBJECT (hadjust), "value_changed",
		    G_CALLBACK (hscroll_event), NULL);
  g_signal_connect (G_OBJECT (vadjust), "value_changed",
		    G_CALLBACK (vscroll_event), NULL);
  g_signal_connect (G_OBJECT (align_area), "button_press_event",
		    G_CALLBACK(button_press_event), NULL);

  gtk_widget_set_events (align_area, GDK_EXPOSURE_MASK
			 | GDK_KEY_RELEASE_MASK
			 | GDK_KEY_PRESS_MASK
                         | GDK_LEAVE_NOTIFY_MASK
                         | GDK_ENTER_NOTIFY_MASK
                         | GDK_BUTTON_PRESS_MASK
                         | GDK_POINTER_MOTION_MASK
                         | GDK_POINTER_MOTION_HINT_MASK);

  gtk_widget_set_events (id_area, GDK_EXPOSURE_MASK
                         | GDK_LEAVE_NOTIFY_MASK
                         | GDK_BUTTON_PRESS_MASK
                         | GDK_POINTER_MOTION_MASK
                         | GDK_POINTER_MOTION_HINT_MASK);

  gtk_widget_show (window);

  vadjust->upper = (gdouble)(1+Sequence_max_sequence_length(head)/(align_area->allocation.width/viewport->char_width) * (2 + align_height));

  vadjust->page_increment = (gdouble)(2 + align_height);
  
  return 0;
}

int main(int argc, char*argv[]) {

  gtk_init (&argc, &argv);

  head = read_fasta(argv[1]);

  create_window(head);

  gtk_main ();

  return 1;
}


