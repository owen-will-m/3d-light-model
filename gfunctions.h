
#ifndef graphics_functions
#define graphics_functions

//structs
typedef struct _color{double r; double g; double b;} color;
typedef struct _point{double x; double y; double z;} point;
typedef struct _polygon{int numpoints; int intersects[100]; point pts[100]; color c;} polygon;


void clip3(polygon * p, double Z);
//functions
void clear(int wid, int hei);
//clears the screen white given width and height
int eliminate(polygon p);

double getlight(polygon p, point source);
double zmid(polygon p);
void copypolygon(polygon * a, polygon * b);
void readfileintopolygons(polygon *  polygons,int * npolys, char * drawing);
//

void assignpointstopoly(polygon p, point * pts, int * intersects,int nsects);
//

void convertpointstoarrays( point p[],int numpoints, double *X, double *Y);
//converts form array of points to array of x and y values

void Ptopoint(double * P, point p);
//
void printpoint(point p);
//prints out a point

double getorigin(double *x, int n);
// given an array, returns the (greatest value + lowest value) / 2 

double getscale(double *x, int n);
//
void copypts(polygon p, point * pts, int n);
//
void clip(double * x,double * y,int * np, double * X, double * Y, int cp);
//
void clip2(double * x,double * y,int * np, double * X, double * Y, int cp);
//
void printpolygon( polygon p);
void printcolor( color c);
int in_out_pt(point * pts, int n, point p);
int sameside(double x1, double y1, double x2, double y2, double px,double py, double midx, double midy);

int in_out (double *x, double *y, int n, double px, double py);


void printarray(double * x, int l);







#endif
