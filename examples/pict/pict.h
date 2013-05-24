
// Defined in main.c:

extern GdkPixbuf *image_pixbuf, *result_pixbuf, *magnified_pixbuf;
extern int image_width, image_height, image_rowstride;
extern int image_extended_width, image_extended_height;
extern guchar *image_pixels;
extern int grid_size;
extern int nu_subdivisions;
extern int use_threading;
extern int block4x4_mode;
extern int magnification;
extern char *output_filename;
extern int texture_output_format;
extern int texture_format;
extern int ktx_no_key_and_value;
extern int nu_islands;
extern int verbose;
extern int modal;

// Defined in gtk.c:

void create_window_layout(int width, int height);
void update_area2();
void clear_result_image();

// Defined ing ga.c:

extern int nu_components;
extern int nu_generations_per_subdivision;
void copy_buffer_to_result_pixbuf(unsigned int *buffer, int sub_image_width, int sub_image_height, int sub_x_offset,
	int sub_y_offset);
void clear_result_image();
void start_ga(int n);
void start_ga_grid(int n);

// Defined in block4x4.c:

void start_ga_block4x4();
void compare_files();


