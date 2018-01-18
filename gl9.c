#include <FPT.h>
#include <D2d_matrix.h>
#include <gfunctions.h>
//GLOBAL VARIABLES
int sheight = 500;int center = 250,l=0;
double zwindow = 2;
double xwindow = 5;

color  colors[10] = {  {.r=.2,.g=.3,.b=.4},
                    {.r=.9,.g=.1,.b=.1},
                    {.r=1,.g=1,.b=0},
                    {.r=0,.g=1,.b=0},
                    {.r=0,.g=0,.b=1}    };
int pcolors[10000];
polygon polygons[10][5000];
polygon temps[50000];



///
//
//
//
//
//
//













void drawpolygon(polygon p, int colorindex,point light){
        double x[p.numpoints],y[p.numpoints],x1[p.numpoints],y1[p.numpoints];
        double id[3][3], useless[3][3], scalefact,xdif,ydifi,lmag;
        int i;


                    clip3(&p,5);




        for(i=0;i<p.numpoints;i++){
                p.pts[i].x = p.pts[i].x / p.pts[i].z;
                p.pts[i].y = p.pts[i].y / p.pts[i].z;
        }
        convertpointstoarrays(p.pts,p.numpoints,x,y);
        D2d_make_identity(id);
        D2d_scale(id,useless,750,750);
        D2d_translate(id,useless,center,center);
        D2d_mat_mult_points(x,y,id,x,y,p.numpoints);
        ///////////////////////////
        lmag = getlight(p,light);
        
        G_rgb(lmag*colors[colorindex].r,lmag*colors[colorindex].g,lmag*colors[colorindex].b);
        
        
        ///////////////////////////
        G_fill_polygon(x,y,p.numpoints);
        //G_rgb(0,0,0);
        //G_polygon(x,y,p.numpoints);
}


//==========main function==========
int main(int z, char * args[]){
        //=====var declaration=====
        int p=1,i=0,j=0,k=0,npolys[10], current=1,tindex;
        double  dz[z],dx[z],dy[z],tx[z],ty[z],tz[z],zrads[z],xrads[z],yrads[z],c=1;
        polygon temp;
        point light = (point){.x=10,.y=0,.z=0};
        G_init_graphics(sheight, sheight);
        for(i=1;i<z;i++){
                readfileintopolygons(polygons[i-1],&npolys[i-1],args[i]);
        }
        z--;
        while(1){
                tindex=0;
                G_rgb(0,0,0);
                G_clear(0,0,0);
                G_rgb(0,1,0);
                if(p>='1'&&p<='9'){
                        current = p - 48;
                        printf("current: %d\n", current);
                }
                for(k=0;k<z;k++){
                        if(current-1==k){
                                if(p=='q')dz[k]-=.2;
                                if(p=='e')dz[k]+=.2;
                                if(p=='z')zrads[k]+=.1;
                                if(p=='x')xrads[k]+=.1;
                                if(p=='c')yrads[k]+=.1;
                                if(p=='w')dy[k]+=c;
                                if(p=='s')dy[k]-=c;
                                if(p=='a')dx[k]-=c;
                                if(p=='d')dx[k]+=c;
                        }
                        for(i=0;i<npolys[k];i++){
                                temp.numpoints = polygons[k][i].numpoints;
                                for(j=0;j<temp.numpoints;j++){
                                        temp.pts[j].x = polygons[k][i].pts[j].x;
                                        temp.pts[j].y = polygons[k][i].pts[j].y;
                                        temp.pts[j].z = polygons[k][i].pts[j].z;
                                        tx[k] = temp.pts[j].x * cos(zrads[k]) - temp.pts[j].y * sin(zrads[k]);
                                        ty[k] = temp.pts[j].x * sin(zrads[k]) + temp.pts[j].y * cos(zrads[k]);
                                        temp.pts[j].x = tx[k];
                                        temp.pts[j].y = ty[k];
                                        ty[k] = temp.pts[j].y * cos(xrads[k]) - temp.pts[j].z * sin(xrads[k]);
                                        tz[k] = temp.pts[j].y * sin(xrads[k]) + temp.pts[j].z * cos(xrads[k]);
                                        temp.pts[j].y = ty[k];
                                        temp.pts[j].z = tz[k];
                                        tz[k] = temp.pts[j].z * cos(yrads[k]) - temp.pts[j].x * sin(yrads[k]);
                                        tx[k] = temp.pts[j].z * sin(yrads[k]) + temp.pts[j].x * cos(yrads[k]);
                                        temp.pts[j].z = tz[k];
                                        temp.pts[j].x = tx[k];
                                        temp.pts[j].x+=dx[k];
                                        temp.pts[j].y+=dy[k];
                                        temp.pts[j].z+=dz[k];
                                }
                                pcolors[tindex] = k;
                                copypolygon(&temp, &temps[tindex]);               
                                tindex++;
                        }
                }
                for(i=0;i<tindex-1;i++){
                        for(j=0;j<tindex-i-1;j++){
                                if(zmid(temps[j])<zmid(temps[j+1])){
                                        l = pcolors[j+1];
                                        pcolors[j+1] = pcolors[j];
                                        pcolors[j] = l;


                                        temp = temps[j+1];
                                        temps[j+1] = temps[j];
                                        temps[j]=temp;
                                }
                        }
                }
                l=0;
                for(i=0;i<tindex;i++){
                    drawpolygon(temps[i], pcolors[i],light);
                }
                p = G_wait_key();
        }
        exit(0);
}
//==========main function==========
