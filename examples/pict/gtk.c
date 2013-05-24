/*
    gtk.c -- gtk 3.0 interface of pict, a graphical example using fgen.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2012, Harm Hanemaaijer

    This file is part of fgen.

    fgen is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fgen is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with fgen.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gtk/gtk.h>

#include "pict.h"

static GtkWidget *window;
static GtkWidget *image1_drawing_area, *image2_drawing_area;
static cairo_surface_t *surface1 = NULL, *surface2 = NULL;

static gboolean
delete_event_cb(GtkWidget * widget, GdkEvent * event, gpointer data)
{
	return FALSE;
}

static void destroy_cb(GtkWidget * widget, gpointer data)
{
	gtk_main_quit();
}

static gboolean key_press_cb(GtkWidget * widget, GdkEventKey * event, gpointer data)
{
	switch (event->keyval) {
	case GDK_KEY_q:
	case GDK_KEY_Q:
		gtk_main_quit();
		exit(0);
		break;
	}

	return TRUE;
}

static gboolean area1_configure_event_cb(GtkWidget * widget, GdkEventConfigure * event, gpointer data)
{
	/* We've handled the configure event, no need for further processing. */
	return TRUE;
}

static gboolean area1_draw_cb(GtkWidget * widget, cairo_t * cr, gpointer data)
{
	if (magnification == 1) {
		gdk_cairo_set_source_pixbuf(cr, image_pixbuf, 0, 0);
		cairo_paint(cr);
		return FALSE;
	}
	GdkPixbuf *magnified = gdk_pixbuf_scale_simple(image_pixbuf, image_width * magnification,
		image_height * magnification, GDK_INTERP_NEAREST);
	gdk_cairo_set_source_pixbuf(cr, magnified, 0, 0);
	cairo_paint(cr);
	g_object_unref(magnified);
	return FALSE;
}

static gboolean area2_configure_event_cb(GtkWidget * widget, GdkEventConfigure * event, gpointer data)
{
	/* We've handled the configure event, no need for further processing. */
	return TRUE;
}

static gboolean area2_draw_cb(GtkWidget * widget, cairo_t * cr, gpointer data)
{
	if (magnification == 1) {
		gdk_cairo_set_source_pixbuf(cr, result_pixbuf, 0, 0);
		cairo_paint(cr);
		return FALSE;
	}
	// Magnfication.
	gdk_pixbuf_scale(result_pixbuf, magnified_pixbuf, 0, 0, image_width * magnification,
		image_height * magnification, 0, 0, magnification, magnification, GDK_INTERP_NEAREST);
	gdk_cairo_set_source_pixbuf(cr, magnified_pixbuf, 0, 0);
	cairo_paint(cr);
	return FALSE;
}

void create_window_layout(int width, int height)
{
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "pict");

	g_signal_connect(G_OBJECT(window), "delete_event",
			 G_CALLBACK(delete_event_cb), NULL);

	g_signal_connect(G_OBJECT(window), "destroy",
			 G_CALLBACK(destroy_cb), NULL);

	g_signal_connect(G_OBJECT(window), "key-press-event",
			 G_CALLBACK(key_press_cb), NULL);

	gtk_container_set_border_width(GTK_CONTAINER(window), 0);

	GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_container_add(GTK_CONTAINER(window), hbox);

	image1_drawing_area = gtk_drawing_area_new();
	image2_drawing_area = gtk_drawing_area_new();
	gtk_box_pack_start(GTK_BOX(hbox), image1_drawing_area, TRUE, TRUE, 0);
	gtk_widget_set_size_request(image1_drawing_area, width, height);
	gtk_box_pack_start(GTK_BOX(hbox), image2_drawing_area, TRUE, TRUE, 0);
	gtk_widget_set_size_request(image2_drawing_area, width, height);

	g_signal_connect(image1_drawing_area, "draw", G_CALLBACK(area1_draw_cb),
			 NULL);
	g_signal_connect(image1_drawing_area, "configure-event",
			 G_CALLBACK(area1_configure_event_cb), NULL);
	g_signal_connect(image2_drawing_area, "draw", G_CALLBACK(area2_draw_cb),
			 NULL);
	g_signal_connect(image2_drawing_area, "configure-event",
			 G_CALLBACK(area2_configure_event_cb), NULL);

	gtk_widget_show_all(window);
}

void clear_result_image() {
	guchar *sub_pixels = gdk_pixbuf_get_pixels(result_pixbuf);
	for (int y = 0; y < image_height; y++)
		memset(sub_pixels + y * image_rowstride, 0, image_width * 3);
	update_area2();
}

void update_area2() {
//    memcpy(gdk_pixbuf_get_pixels(result_pixbuf), gdk_pixbuf_get_pixels(scratch_pixbuf),
//	image_height * image_width * 3);
    gtk_widget_queue_draw(GTK_WIDGET(image2_drawing_area));
//    gdk_window_process_updates(window, TRUE);  
}

