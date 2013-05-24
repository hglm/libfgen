/*
    gtk.c -- gtk interface of tsp, a graphical TSP example using fgen.

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
#define _USE_MATH_DEFINES
#include <math.h>
#include <gtk/gtk.h>
#include "tsp.h"

typedef struct {
    int index;
    GtkWidget *image_drawing_area;
    cairo_surface_t *surface;
    double area_width, area_height;
} Window;

static GtkWidget *gtk_window;
static Window *window;
static int GUI_initialized = 0;
static int window_size = 800;

static gboolean delete_event_cb(GtkWidget *widget, GdkEvent *event, gpointer data) {
    return FALSE;
}

static void destroy_cb(GtkWidget *widget, gpointer data) {
    gtk_main_quit();
    exit(0);
}

static gboolean area_configure_event_cb(GtkWidget *widget, GdkEventConfigure *event, Window *window) {
    if (!GUI_initialized) {
        return TRUE;
    }
    window->area_width = event->width;
    window->area_height = event->height;
    cairo_surface_destroy(window->surface);
    window->surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, window->area_width, window->area_height);
    gui_draw_window();
    /* We've handled the configure event, no need for further processing. */
    return TRUE;
}

static gboolean area_draw_cb(GtkWidget *widget, cairo_t *cr, Window *window) {
    cairo_set_source_surface(cr, window->surface, 0, 0);
    cairo_paint(cr);
    return TRUE;
}

static void gui_update_drawing_area(Window *window) {
    gtk_widget_queue_draw(GTK_WIDGET(window->image_drawing_area));
}

void gui_initialize(int *argc, char ***argv) {
    gtk_init(argc, argv);
}

static void menu_item_new_cities_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    if (ga_busy)
        return;
    initialize_cities();
    gui_draw_and_show_window();
}

static void entry_nu_cities_reset_text(GtkEntry *entry) {
    char s[8];
    sprintf(s, "%d", nu_cities);
    gtk_entry_set_text(entry, s);
}

static void entry_nu_cities_activate_cb(GtkEntry *entry, gpointer data) {
    if (ga_busy)
        return;
    const gchar *s = gtk_entry_get_text(entry);
    int value;
    int n = sscanf(s, "%d", &value);
    if (n == 1 && value >= 3 && value <= 1000) {
        nu_cities = value;
        initialize_cities();
        gui_draw_and_show_window();
        return;
    }
    entry_nu_cities_reset_text(entry);
}

static void entry_mutation_prob_reset_text(GtkEntry *entry) {
    char s[8];
    sprintf(s, "%.3f", mutation_probability);
    gtk_entry_set_text(entry, s);
}

static void entry_mutation_prob_activate_cb(GtkEntry *entry, gpointer data) {
    if (ga_busy)
        return;
    const gchar *s = gtk_entry_get_text(entry);
    float value;
    int n = sscanf(s, "%f", &value);
    if (n == 1 && value >= 0 && value <= 1.0) {
        mutation_probability = value;
        return;
    }
    entry_mutation_prob_reset_text(entry);
}

static void entry_crossover_prob_reset_text(GtkEntry *entry) {
    char s[8];
    sprintf(s, "%.3f", crossover_probability);
    gtk_entry_set_text(entry, s);
}

static void entry_crossover_prob_activate_cb(GtkEntry *entry, gpointer data) {
    if (ga_busy)
        return;
    const gchar *s = gtk_entry_get_text(entry);
    float value;
    int n = sscanf(s, "%f", &value);
    if (n == 1 && value >= 0 && value <= 1.0) {
        crossover_probability = value;
        return;
    }
    entry_crossover_prob_reset_text(entry);
}

GtkWidget *entry_mutation_prob;
GtkWidget *entry_crossover_prob;

static void menu_item_run_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    if (ga_busy)
        return;
    entry_mutation_prob_activate_cb(GTK_ENTRY(entry_mutation_prob), NULL);
    entry_crossover_prob_activate_cb(GTK_ENTRY(entry_crossover_prob), NULL);
    start_ga();
}

static void menu_item_stop_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    if (!ga_busy)
        return;
    stop_ga();
}

static void menu_item_benchmark_current_settings_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    if (ga_busy)
        return;
    entry_mutation_prob_activate_cb(GTK_ENTRY(entry_mutation_prob), NULL);
    entry_crossover_prob_activate_cb(GTK_ENTRY(entry_crossover_prob), NULL);
    benchmark_current_settings();
}

static void menu_item_benchmark_operators_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    if (ga_busy)
        return;
    entry_mutation_prob_activate_cb(GTK_ENTRY(entry_mutation_prob), NULL);
    entry_crossover_prob_activate_cb(GTK_ENTRY(entry_crossover_prob), NULL);
    benchmark_operators();
}

static void menu_item_benchmark_rates_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    if (ga_busy)
        return;
    entry_mutation_prob_activate_cb(GTK_ENTRY(entry_mutation_prob), NULL);
    entry_crossover_prob_activate_cb(GTK_ENTRY(entry_crossover_prob), NULL);
    benchmark_rates();
}

static void menu_item_quit_activate_cb(GtkMenuItem *menu_item, gpointer data) {
    gtk_main_quit();
    exit(0);
}

static GtkWidget *population_size_radio_button1;
static GtkWidget *population_size_radio_button2;
static GtkWidget *population_size_radio_button3;
static GtkWidget *population_size_radio_button4;
static GtkWidget *population_size_radio_button5;
static GtkWidget *population_size_radio_button6;
static GtkWidget *population_size_radio_button7;
static GtkWidget *population_size_radio_button8;

static void population_size_radio_button_toggled_cb(GtkButton *button, gpointer user_data) {
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button1)))
        population_size = 64;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button2)))
        population_size = 128;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button3)))
        population_size = 256;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button4)))
        population_size = 512;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button5)))
        population_size = 1024;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button6)))
        population_size = 2048;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button7)))
        population_size = 4096;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(population_size_radio_button8)))
        population_size = 8192;
}

static GtkWidget *mutation_type_radio_button1;
static GtkWidget *mutation_type_radio_button2;
static GtkWidget *mutation_type_radio_button3;
static GtkWidget *mutation_type_radio_button4;
static GtkWidget *mutation_type_radio_button5;

static void mutation_type_radio_button_toggled_cb(GtkButton *button, gpointer user_data) {
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(mutation_type_radio_button1)))
        mutation_type = MUTATION_TYPE_SWAP;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(mutation_type_radio_button2)))
        mutation_type = MUTATION_TYPE_SWAP_NEIGHBOURS;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(mutation_type_radio_button3)))
        mutation_type = MUTATION_TYPE_INSERT;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(mutation_type_radio_button4)))
        mutation_type = MUTATION_TYPE_INVERT;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(mutation_type_radio_button5)))
        mutation_type = MUTATION_TYPE_NOOP;
}

static GtkWidget *crossover_type_radio_button1;
static GtkWidget *crossover_type_radio_button2;
static GtkWidget *crossover_type_radio_button3;

static void crossover_type_radio_button_toggled_cb(GtkButton *button, gpointer user_data) {
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(crossover_type_radio_button1)))
        crossover_type = CROSSOVER_TYPE_ORDER1;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(crossover_type_radio_button2)))
        crossover_type = CROSSOVER_TYPE_PBX;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(crossover_type_radio_button3)))
        crossover_type = CROSSOVER_TYPE_NOOP;
}

static GtkWidget *reporting_frequency_radio_button1;
static GtkWidget *reporting_frequency_radio_button2;
static GtkWidget *reporting_frequency_radio_button3;
static GtkWidget *reporting_frequency_radio_button4;

static void reporting_frequency_radio_button_toggled_cb(GtkButton *button, gpointer user_data) {
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(reporting_frequency_radio_button1)))
        reporting_frequency = 1;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(reporting_frequency_radio_button2)))
        reporting_frequency = 10;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(reporting_frequency_radio_button3)))
        reporting_frequency = 100;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(reporting_frequency_radio_button4)))
        reporting_frequency = 1000;
}

static void steady_state_check_button_toggled_cb(GtkCheckButton *check_button, gpointer user_data) {
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_button)))
        steady_state = 1;
    else
        steady_state = 0;
}

void gui_create_window_layout() {
    gtk_window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(gtk_window), "Travelling Salesman Problem");

    g_signal_connect(G_OBJECT(gtk_window), "delete_event", G_CALLBACK(delete_event_cb), NULL);

    g_signal_connect(G_OBJECT(gtk_window), "destroy", G_CALLBACK(destroy_cb), NULL);

    gtk_container_set_border_width(GTK_CONTAINER(gtk_window), 0);

    // Create the menu vbox for holding the menu and the rest of the application.
    GtkWidget *menu_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_container_add(GTK_CONTAINER(gtk_window), menu_vbox);
    // Create the menu items.
    GtkWidget *menu_item_new_cities = gtk_menu_item_new_with_label("New cities");
    g_signal_connect(G_OBJECT(menu_item_new_cities), "activate", G_CALLBACK(menu_item_new_cities_activate_cb), NULL);
    GtkWidget *menu_item_run = gtk_menu_item_new_with_label("Run");
    g_signal_connect(G_OBJECT(menu_item_run), "activate", G_CALLBACK(menu_item_run_activate_cb), NULL);
    GtkWidget *menu_item_stop = gtk_menu_item_new_with_label("Stop");
    g_signal_connect(G_OBJECT(menu_item_stop), "activate", G_CALLBACK(menu_item_stop_activate_cb), NULL);
    GtkWidget *menu_item_benchmark_current_settings = gtk_menu_item_new_with_label("Benchmark (current settings)");
    g_signal_connect(G_OBJECT(menu_item_benchmark_current_settings), "activate",
        G_CALLBACK(menu_item_benchmark_current_settings_activate_cb), NULL);
    GtkWidget *menu_item_benchmark_operators = gtk_menu_item_new_with_label("Benchmark (operators)");
    g_signal_connect(G_OBJECT(menu_item_benchmark_operators), "activate",
        G_CALLBACK(menu_item_benchmark_operators_activate_cb), NULL);
    GtkWidget *menu_item_benchmark_rates = gtk_menu_item_new_with_label("Benchmark (rates)");
    g_signal_connect(G_OBJECT(menu_item_benchmark_rates), "activate", G_CALLBACK(menu_item_benchmark_rates_activate_cb), NULL);
    GtkWidget *menu_item_quit = gtk_menu_item_new_with_label("Quit");
    g_signal_connect(G_OBJECT(menu_item_quit), "activate", G_CALLBACK(menu_item_quit_activate_cb), NULL);
    GtkWidget *action_menu = gtk_menu_new();
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_new_cities);
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_run);
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_stop);
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_benchmark_current_settings);
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_benchmark_operators);
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_benchmark_rates);
    gtk_menu_shell_append(GTK_MENU_SHELL(action_menu), menu_item_quit);
    // Create the menu bar.
    GtkWidget *menu_bar = gtk_menu_bar_new();
    // Create the Action menu.
    GtkWidget *action_item = gtk_menu_item_new_with_label("Action");
    gtk_menu_item_set_submenu(GTK_MENU_ITEM(action_item), action_menu);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu_bar), action_item);
    // Add the whole menu.
    gtk_box_pack_start(GTK_BOX(menu_vbox), menu_bar, FALSE, FALSE, 0);

    window = (Window *)malloc(sizeof(Window));
    // Create a hbox.
    GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    // Create a vbox that holds the settings to the left.
    GtkWidget *settings_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    // Create the number of cities label and entry.
    GtkWidget *label_nu_cities = gtk_label_new("Number of cities:");
    GtkWidget *entry_nu_cities = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry_nu_cities), 4);
    entry_nu_cities_reset_text(GTK_ENTRY(entry_nu_cities));
    g_signal_connect(G_OBJECT(entry_nu_cities), "activate", G_CALLBACK(entry_nu_cities_activate_cb), window);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_nu_cities, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), entry_nu_cities, FALSE, FALSE, 0);
    // Create the population size label and radio buttons.
    GtkWidget *label_population_size = gtk_label_new("Population size:");
    population_size_radio_button1 = gtk_radio_button_new_with_label(NULL, "64");
    population_size_radio_button2 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "128");
    population_size_radio_button3 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "256");
    population_size_radio_button4 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "512");
    population_size_radio_button5 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "1024");
    population_size_radio_button6 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "2048");
    population_size_radio_button7 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "4096");
    population_size_radio_button8 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(population_size_radio_button1)), "8192");
    g_signal_connect(G_OBJECT(population_size_radio_button1), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button2), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button3), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button4), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button5), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button6), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button7), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(population_size_radio_button8), "toggled",
        G_CALLBACK(population_size_radio_button_toggled_cb), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(population_size_radio_button3), TRUE);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_population_size, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button1, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button2, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button3, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button4, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button5, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button6, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button7, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), population_size_radio_button8, FALSE, FALSE, 0);
    // Create the mutation probability label and entry.
    GtkWidget *label_mutation_prob = gtk_label_new("Mutation probability:");
    entry_mutation_prob = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry_mutation_prob), 5);
    entry_mutation_prob_reset_text(GTK_ENTRY(entry_mutation_prob));
    g_signal_connect(G_OBJECT(entry_mutation_prob), "activate", G_CALLBACK(entry_mutation_prob_activate_cb), window);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_mutation_prob, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), entry_mutation_prob, FALSE, FALSE, 0);
    // Create the crossover probability label and entry.
    GtkWidget *label_crossover_prob = gtk_label_new("Crossover probability:");
    entry_crossover_prob = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry_crossover_prob), 5);
    entry_crossover_prob_reset_text(GTK_ENTRY(entry_crossover_prob));
    g_signal_connect(G_OBJECT(entry_crossover_prob), "activate", G_CALLBACK(entry_crossover_prob_activate_cb), window);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_crossover_prob, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), entry_crossover_prob, FALSE, FALSE, 0);
    // Create the mutation type label and radio buttons.
    GtkWidget *label_mutation_type = gtk_label_new("Mutation type:");
    mutation_type_radio_button1 = gtk_radio_button_new_with_label(NULL, "Swap");
    mutation_type_radio_button2 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(mutation_type_radio_button1)), "Swap neighbours");
    mutation_type_radio_button3 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(mutation_type_radio_button1)), "Insert");
    mutation_type_radio_button4 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(mutation_type_radio_button1)), "Invert");
    mutation_type_radio_button5 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(mutation_type_radio_button1)), "No operation");
    g_signal_connect(G_OBJECT(mutation_type_radio_button1), "toggled",
        G_CALLBACK(mutation_type_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(mutation_type_radio_button2), "toggled",
        G_CALLBACK(mutation_type_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(mutation_type_radio_button3), "toggled",
        G_CALLBACK(mutation_type_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(mutation_type_radio_button4), "toggled",
        G_CALLBACK(mutation_type_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(mutation_type_radio_button5), "toggled",
        G_CALLBACK(mutation_type_radio_button_toggled_cb), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mutation_type_radio_button4), TRUE);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_mutation_type, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), mutation_type_radio_button1, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), mutation_type_radio_button2, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), mutation_type_radio_button3, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), mutation_type_radio_button4, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), mutation_type_radio_button5, FALSE, FALSE, 0);
    // Create the crossover type label and radio buttons.
    GtkWidget *label_crossover_type = gtk_label_new("Crossover type:");
    crossover_type_radio_button1 = gtk_radio_button_new_with_label(NULL, "Order Based 1");
    crossover_type_radio_button2 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(crossover_type_radio_button1)), "Position Based");
    crossover_type_radio_button3 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(crossover_type_radio_button1)), "No operation");
    g_signal_connect(G_OBJECT(crossover_type_radio_button1), "toggled",
        G_CALLBACK(crossover_type_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(crossover_type_radio_button2), "toggled",
        G_CALLBACK(crossover_type_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(crossover_type_radio_button3), "toggled",
        G_CALLBACK(crossover_type_radio_button_toggled_cb), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(crossover_type_radio_button1), TRUE);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_crossover_type, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), crossover_type_radio_button1, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), crossover_type_radio_button2, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), crossover_type_radio_button3, FALSE, FALSE, 0);
    // Create the reporting frequency label and radio buttons.
    GtkWidget *label_reporting_frequency = gtk_label_new("Reporting frequency:");
    reporting_frequency_radio_button1 = gtk_radio_button_new_with_label(NULL, "Every generation");
    reporting_frequency_radio_button2 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(reporting_frequency_radio_button1)), "Every 10th generation");
    reporting_frequency_radio_button3 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(reporting_frequency_radio_button1)), "Every 100th generation");
    reporting_frequency_radio_button4 = gtk_radio_button_new_with_label(
        gtk_radio_button_get_group(GTK_RADIO_BUTTON(reporting_frequency_radio_button1)), "Every 1000th generation");
    g_signal_connect(G_OBJECT(reporting_frequency_radio_button1), "toggled",
        G_CALLBACK(reporting_frequency_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(reporting_frequency_radio_button2), "toggled",
        G_CALLBACK(reporting_frequency_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(reporting_frequency_radio_button3), "toggled",
        G_CALLBACK(reporting_frequency_radio_button_toggled_cb), NULL);
    g_signal_connect(G_OBJECT(reporting_frequency_radio_button4), "toggled",
        G_CALLBACK(reporting_frequency_radio_button_toggled_cb), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(reporting_frequency_radio_button1), TRUE);
    gtk_box_pack_start(GTK_BOX(settings_vbox), label_reporting_frequency, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), reporting_frequency_radio_button1, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), reporting_frequency_radio_button2, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), reporting_frequency_radio_button3, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(settings_vbox), reporting_frequency_radio_button4, FALSE, FALSE, 0);
    // Create the steady-state check box.
    GtkWidget *steady_state_check_button = gtk_check_button_new_with_label("Use steady-state evolution");
    g_signal_connect(G_OBJECT(steady_state_check_button), "toggled",
        G_CALLBACK(steady_state_check_button_toggled_cb), NULL);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(steady_state_check_button), FALSE);
    gtk_box_pack_start(GTK_BOX(settings_vbox), steady_state_check_button, FALSE, FALSE, 0);

    // Add the settings vbox to the hbox
    gtk_box_pack_start(GTK_BOX(hbox), settings_vbox, FALSE, FALSE, 8);

    // Add the drawing area to the hbox.
    window->image_drawing_area = gtk_drawing_area_new();
    gtk_widget_set_size_request(window->image_drawing_area, window_size, window_size);
    gtk_box_pack_start(GTK_BOX(hbox), window->image_drawing_area, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(menu_vbox), hbox, TRUE, TRUE, 0);

    window->surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, window_size, window_size);
    window->area_width = window_size;
    window->area_height = window_size;
    g_signal_connect(window->image_drawing_area, "draw", G_CALLBACK(area_draw_cb), window);
    g_signal_connect(window->image_drawing_area, "configure-event", G_CALLBACK(area_configure_event_cb), window);

    gtk_widget_show_all(gtk_window);
    GUI_initialized = 1;
}

void gui_draw_window() {
    cairo_t *cr;
    cr = cairo_create(window->surface);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_paint(cr);
    cairo_scale(cr, window->area_width, window->area_height);
    cairo_set_line_width(cr, 1.0 / window->area_width);
    // Draw the cities.
    cairo_set_source_rgb(cr, 0, 1, 0);
    for (int i = 0; i < nu_cities; i++) {
        // Draw the city as a green circle.
        double x = city[i].x / (AREA_SIZE + AREA_SIZE / 5) + 0.1;
        double y = city[i].y / (AREA_SIZE + AREA_SIZE / 5) + 0.1;
        cairo_move_to(cr, x + 0.02, y);
        cairo_arc(cr, x, y, 0.02, 0, 2 * M_PI);
        cairo_stroke(cr);
    }
    // Draw the best route.
    if (best_route_valid) {
        cairo_set_line_width(cr, 2.0 / window->area_width);
        cairo_set_source_rgb(cr, 1, 0, 0);
        double x = city[best_route[0]].x / (AREA_SIZE + AREA_SIZE / 5) + 0.1;
        double y = city[best_route[0]].y / (AREA_SIZE + AREA_SIZE / 5) + 0.1;
        cairo_move_to(cr, x, y);
        for (int i = 1; i < nu_cities; i++) {
            double x = city[best_route[i]].x / (AREA_SIZE + AREA_SIZE / 5) + 0.1;
            double y = city[best_route[i]].y / (AREA_SIZE + AREA_SIZE / 5) + 0.1;
            cairo_line_to(cr, x, y);
        }
        cairo_stroke(cr);
        // Show text.
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 0.025);
        cairo_set_source_rgb(cr, 0, 1, 1);
        char s[32];
        if (ga_busy) {
            cairo_move_to(cr, 0.03, 0.03);
            sprintf(s, "Running (gen = %d)", ga_generation);
            cairo_show_text(cr, s);
        }
        cairo_move_to(cr, 0.03, 0.03 + 0.028);
        sprintf(s, "Best distance: %.3lf", best_distance());
        cairo_show_text(cr, s);
    }
    cairo_destroy(cr);
}

void gui_draw_and_show_window() {
    gui_draw_window();
    gui_update_drawing_area(window);
}

void gui_run() {
    gtk_main();
}

void gui_handle_events() {
    while (gtk_events_pending())
        gtk_main_iteration();
}


