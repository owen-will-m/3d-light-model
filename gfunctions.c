#include <gfunctions.h>
#include <FPT.h> 



void copypolygon(polygon * a, polygon * b){
        //copy a into b
        int i=0;
        b -> numpoints = a -> numpoints;
        for(i=0;i<a->numpoints;i++){
                b->pts[i].x = a->pts[i].x;
                b->pts[i].y = a->pts[i].y;
                b->pts[i].z = a->pts[i].z;
        }

}



double zmid(polygon p){
        int i; double sum=0;
        for(i=0;i<p.numpoints;i++){
                sum+=p.pts[i].z;
        }
        return sum / p.numpoints;
}
double getyorigin(polygon p){
        int i; double sum=0;
        for(i=0;i<p.numpoints;i++){
                sum+=p.pts[i].y;
        }
        return sum / p.numpoints;
}
double getxorigin(polygon p){
        int i; double sum=0;
        for(i=0;i<p.numpoints;i++){
                sum+=p.pts[i].x;
        }
        return sum / p.numpoints;
}
void adjustzorigin(polygon p, double offset){
        int i;double zorigin = zmid(p);
        double change = offset - zorigin;
        //printf("\n%lf\n",change);
        puts("polygon points:");
        for(i=0;i<p.numpoints;i++){
                p.pts[i].z = p.pts[i].z + change;
                printf("\n%lf\n",p.pts[i].z);
        }
}

void rotate(polygon p, double rads){
        //about z axis
        int i=0;
        double tx,ty;
        for(i=0;i<p.numpoints;i++){
                tx = p.pts[i].x * cos(rads) - p.pts[i].y * sin(rads);
                ty = p.pts[i].x * sin(rads) + p.pts[i].y * cos(rads);

                p.pts[i].x = tx;
                p.pts[i].y = ty;
        }


}
double getshapeorigin(polygon * p, int n){

}


void copypts(polygon p, point * pts, int n){
        int i;
        p.numpoints = n;
        for(i=0;i<n;i++){
                p.pts[i] = pts[i];
        }
}


double middle(double *x, int n){
        int i; double highest=x[0],lowest=x[0];
        for(i=1;i<n;i++){
                if(highest<x[i]) highest = x[i];
                if(lowest>x[i]) lowest = x[i];
        }
        return (highest + lowest) / 2;
}


int sameside(double x1, double y1, double x2, double y2, double px,double py, double midx, double midy){
        double dy =   y2 - y1;
        double dx = x2 - x1;
        double num = dy * (midx - x1) - dx * (midy -y1);
        double num2= dy * (px - x1) - dx * (py - y1);
        if((num<=0&&num2<=0)||(num>=0&&num2>=0)){
                return 1;
        }else{
                return 0;
        }
}

int in_out (double *x, double *y, int n, double px, double py)
        // return 1 if point P is inside the convex polygon
        // else return 0
{
        double midx = middle(x,n);
        double midy = middle(y,n);
        int i;
        for(i=0;i<n-1;i++){

                if(!sameside(x[i],y[i],x[i+1],y[i+1],px,py,midx,midy)) return 0;

        }
        if(!sameside(x[0],y[0],x[n-1],y[n-1],px,py,midx,midy)) return 0;


        return 1;
}
int in_out_pt(point * pts, int n, point p)
        // return 1 if point P is inside the convex polygon
        // else return 0
{
        double * x; double * y;
        convertpointstoarrays(pts,n,x,y);

        double midx = middle(x,n);
        double midy = middle(y,n);
        int i;
        for(i=0;i<n-1;i++){

                if(!sameside(x[i],y[i],x[i+1],y[i+1],p.x,p.y,midx,midy)) return 0;

        }
        if(!sameside(x[0],y[0],x[n-1],y[n-1],p.x,p.y,midx,midy)) return 0;


        return 1;
}

void clear(int wid, int hei){
        G_rgb(1,1,1);
        G_fill_rectangle(0,0,wid,hei);
}
void convertpointstoarrays( point p[],int numpoints, double *X, double *Y){
        int i;
        for(i=0;i<numpoints;i++){
                X[i] = p[i].x;
                Y[i] = p[i].y;
        }
}
void printpoint(point p){
        printf("\nx: %lf, y: %lf, z: %lf\n", p.x,p.y,p.z);
}
void printcolor( color c){
        printf("\nColor: red:%lf green:%lf blue:%lf\n",c.r,c.g,c.b);
}
void printpolygon( polygon p){
        printf("\n---Polygon---\nnum pts: %d\npoints: ", p.numpoints);
        int i=0;
        for(i=0;i<p.numpoints;i++){
                printpoint(p.pts[i]);
        }
        printcolor(p.c);
        puts("\n");
}
double getorigin(double *x, int n){
        int i; double highest=x[0],lowest=x[0];
        for(i=1;i<n;i++){
                if(highest<x[i]) highest = x[i];
                if(lowest>x[i]) lowest = x[i];
        }
        return (highest + lowest) / 2;
}
double getscale(double *x, int n){
        int i; double highest=x[0],lowest=x[0];
        for(i=1;i<n;i++){
                if(highest<x[i]) highest = x[i];
                if(lowest>x[i]) lowest = x[i];
        }
        return 250 / (highest - lowest); 

}
/*
   void drawpolygon(polygon p){
   double * x; double * y;
   convertpointstoarrays(p.pts,p.numpoints,x,y);
   G_fill_polygon(x,y,p.numpoints);
   }
   */

void Ptopoint(double * P, point p){
        p.x = P[0]; p.y =P[1];
}
double getM(double x1,double y1, double x2,double y2){
        if(x1==x2) x2+=.0001;
        if(y1==y2) y2+=.0001;
        return (y1 - y2) / (x1 - x2);
}
double getB(double m, double y, double x){
        //y=mx+b so b = y - mx
        return y - m*x;
}

void intersect(double *nx, double *ny, double Ax, double Ay, double Bx, double By,double Cx, double Cy, double Dx, double Dy){
        // (y0 - y1) = m * (x0 - x1) 
        //y = mx + b
        double m1, m2, b1, b2;
        m1 = getM(Ax,Ay,Bx,By); b1 = getB(m1,Ay,Ax);
        m2 = getM(Cx,Cy,Dx,Dy); b2 = getB(m2,Cy,Cx);
        if(m1==m2) m1+=.000001;
        *nx = (b2-b1) / (m1-m2);
        *ny = m1 * *nx + b1;
}


int validintersect(double *nx, double *ny, double * x, double * y, int np, double outx, double outy, double inx, double iny){
        int i; double tempx,tempy;

        for(i=1;i<=np;i++){
                //get intersection point for infinite lines
                if(i==np){
                        intersect(&tempx,&tempy,x[0],y[0],x[np-1],y[np-1],outx,outy,inx,iny);
                }else{
                        intersect(&tempx,&tempy,x[i-1],y[i-1],x[i],y[i],outx,outy,inx,iny);
                }
                if((tempx>=outx&&tempx<=inx)||(tempx<=outx&&tempx>=inx)){//w/in bounds of given line
                        *nx = tempx; *ny = tempy;
                        return 1;
                }
                //check that they are within the bounds of the line inx,iny outx,outy
        }
        return 0;
}


void clip(double * x,double * y,int * np, double * X, double * Y, int cp){
        int i,j,k=0,oneinside,twoinside; double ix, iy, newx[1000], newy[1000];

        //add first value to end of the array--
        //that way it looks at last to first pt 
        x[*np] = x[0];
        y[*np] = y[0];
        *np = *np +1;

        for(j=1;j<*np;j++){
                i = j;
                oneinside =  in_out(X,Y,cp,x[j-1],y[j-1]);//true or false
                twoinside =  in_out(X,Y,cp,x[i],y[i]);//true or false
                //in-in (keep both)
                if(oneinside&&twoinside){
                        newx[k] = x[j-1]; newy[k] = y[j-1];
                        k++;
                        newx[k] = x[i]; newy[k] = y[i];
                        k++;
                }
                //in-out (keep first point, keep intersection w clip line
                if(oneinside&&!twoinside){
                        newx[k] = x[j-1]; newy[k] = y[j-1];
                        k++;
                        validintersect(&newx[k],&newy[k],X,Y,cp,x[j-1],y[j-1],x[i],y[i]);
                        k++;
                }
                //out-out (keep none)
                if(!oneinside&&!twoinside){
                        //there could be multiple line intersections!!!! what do?
                        //jeff didn't address this in his class!!


                }
                //out-in (keep intersection w clip line and second point)
                if(!oneinside&&twoinside){
                        validintersect(&newx[k],&newy[k],X,Y,cp,x[j-1],y[j-1],x[i],y[i]);
                        k++;
                        newx[k] = x[i]; newy[k] = y[i];
                        k++;
                }
        }





        printf("\n k is %d\n",k); 
        *np = k;
        for(i=0;i<k;i++){
                x[i] = newx[i];
                y[i] = newy[i];
        }
}
void copyarray(double * in, double * out, int n){
        int i;
        for(i=0;i<n;i++){
                out[i]=in[i];
        }
}


void clip2(double * x,double * y,int * np, double * X, double * Y, int cp){
        int i,j,k=0,bool1,bool2; double midx,midy,tempx[1000], tempy[1000];

        midx=middle(X,cp); midy=middle(Y,cp);
        x[*np] = x[0]; y[*np] = y[0];
        *np = *np + 1;
        for(j=1;j<cp;j++){
                k=0;
                for(i=1;i<*np;i++){
                        bool1= sameside(X[j-1],Y[j-1],X[j],Y[j],x[i-1],y[i-1],midx,midy);            
                        bool2= sameside(X[j-1],Y[j-1],X[j],Y[j],x[i],y[i],midx,midy);
                        if(bool1&&bool2){//in, in keep first
                                tempx[k] = x[i-1]; tempy[k] = y[i-1]; k++;
                        }else if(!bool1&&bool2){// out, in--get intersection of lines 
                                intersect(&tempx[k],&tempy[k],x[i-1],y[i-1],x[i],y[i],X[j-1],Y[j-1],X[j],Y[j]); k++; 
                        }else if(bool1&&!bool2){//in, out-- keep first, get intersection of lines
                                tempx[k] = x[i-1]; tempy[k] = y[i-1]; k++;
                                intersect(&tempx[k],&tempy[k],x[i-1],y[i-1],x[i],y[i],X[j-1],Y[j-1],X[j],Y[j]); k++; 
                        }
                }
                tempx[k] = tempx[0]; tempy[k] = tempy[0]; k++;
                copyarray(tempx, x, k);
                copyarray(tempy, y, k);
                *np = k;
        }
        printf("k: %d\n",k);


}



void readfileintopolygons(polygon *  polygons,int * npolys, char * drawing){
        int i,np,j; FILE *f; int numpoints, intersects[100];point points[2000];
        f = fopen(drawing,"r") ;
        if (f == NULL) {
                printf("filename: %s\n",drawing);
                exit(0) ;
        }

        fscanf(f,"%d",&numpoints) ;

        for(i=0;i<numpoints;i++){
                fscanf(f,"%lf",&points[i].x);
                fscanf(f,"%lf",&points[i].y);
                fscanf(f,"%lf",&points[i].z);
                //printpoint(points[i]);
        }

        fscanf(f,"%d",npolys);
        for(i=0;i<*npolys;i++){
                fscanf(f,"%d",&np);
                //puts("\n intersects:");printf("%d---",np);
                for(j=0;j<np;j++){
                        fscanf(f,"%d",&intersects[j]);
                }

                for (j=0;j<np;j++){
                        polygons[i].pts[j].x = points[intersects[j]].x;
                        polygons[i].pts[j].y = points[intersects[j]].y;
                        polygons[i].pts[j].z = points[intersects[j]].z;
                }




                polygons[i].numpoints = np;
                //                printpolygon(polygons[i]);

        }
        // printf("\n%d-----NUMPTSd\n",*numpoints);
        fclose(f);
}

void copypoint(point a, point b){//copies point a into b
        b.x=a.x;
        b.y=a.y;
        b.z=a.z;
}
void assignpointstopoly(polygon p, point * pts,int * intersects,int nsects){
        int i;
        for (i=0;i<nsects;i++){
                p.pts[i].x = pts[intersects[i]].x;
                p.pts[i].y = pts[intersects[i]].y;
                p.pts[i].z = pts[intersects[i]].z;
        }
        p.numpoints = nsects;
        printpolygon(p);
}
void printarray(double * x, int l){
        int i;
        puts("array:\n");
        for(i=0;i<l;i++){
                printf("||%lf||",x[i]);
        }
}

double dist(point a, point b){
        return sqrt( (a.z-b.z)*(a.z-b.z) + (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) );
}
double magnitude(point a, point b, char dimension){
        if(dimension=='x'){
                return a.x-b.x;
        }else if(dimension=='y'){
                return a.y-b.y;
        }else if(dimension=='z'){
                return a.z-b.z;
        }else{
                return 0;
        }
}

double anglebetween(point a1, point a2, point b1, point b2){
        //theta = acos( ( a dot b ) / ( |a|*|b| ) )
        double amag = dist(a1, a2), bmag = dist(b1,b2); 
        double adotb =  (a1.x - a2.x) * (b1.x - b2.x) +
                        (a1.y - a2.y) * (b1.y - b2.y) +
                        (a1.z - a2.z) * (b1.z - b2.z);
        return acos(adotb / (amag*bmag));
}


int eliminate(polygon p){
        double amag,bmag,adotb,axb,theta,theta2; point c, or;
        int i=1;
        or.x=0;or.y=0;or.z=0;
                c.x = magnitude(p.pts[i],p.pts[i-1],'y') * magnitude(p.pts[i],p.pts[i+1],'z') -
                      magnitude(p.pts[i],p.pts[i-1],'z') * magnitude(p.pts[i],p.pts[i+1],'y');
                c.y = magnitude(p.pts[i],p.pts[i-1],'z') * magnitude(p.pts[i],p.pts[i+1],'x') -
                      magnitude(p.pts[i],p.pts[i-1],'x') * magnitude(p.pts[i],p.pts[i+1],'z');
                c.z = magnitude(p.pts[i],p.pts[i-1],'x') * magnitude(p.pts[i],p.pts[i+1],'y') -
                      magnitude(p.pts[i],p.pts[i-1],'y') * magnitude(p.pts[i],p.pts[i+1],'x');
                c.x = c.x + p.pts[i].x;
                c.y = c.y + p.pts[i].y;
                c.z = c.z + p.pts[i].z;
                theta2 = anglebetween(p.pts[i],c,p.pts[i],or);
                if(theta2>=(1.57))return 1;
        return 0;
}

////////////////
double getlight(polygon p, point source){
double intensity, ambient = .2, diffusemax = .3, specpow = 20; 

        double amag,bmag,adotb,axb,thetare,thetanl,scalarprod; 
        point c, eye = (point){.x=0,.y=0,.z=0}, l;
        int i=1;
                c.x = magnitude(p.pts[i],p.pts[i-1],'y') * magnitude(p.pts[i],p.pts[i+1],'z') -
                      magnitude(p.pts[i],p.pts[i-1],'z') * magnitude(p.pts[i],p.pts[i+1],'y');
                c.y = magnitude(p.pts[i],p.pts[i-1],'z') * magnitude(p.pts[i],p.pts[i+1],'x') -
                      magnitude(p.pts[i],p.pts[i-1],'x') * magnitude(p.pts[i],p.pts[i+1],'z');
                c.z = magnitude(p.pts[i],p.pts[i-1],'x') * magnitude(p.pts[i],p.pts[i+1],'y') -
                      magnitude(p.pts[i],p.pts[i-1],'y') * magnitude(p.pts[i],p.pts[i+1],'x');
                c.x = c.x + p.pts[i].x;
                c.y = c.y + p.pts[i].y;
                c.z = c.z + p.pts[i].z;

                thetanl = anglebetween(p.pts[i],c,p.pts[i],source);
                scalarprod = cos(thetanl)*2;
                
                c.x = c.x - p.pts[i].x;
                c.y = c.y - p.pts[i].y;
                c.z = c.z - p.pts[i].z;
                
                c.x = scalarprod*c.x -  magnitude(p.pts[i],source,'x') ;
                c.y = scalarprod*c.y -  magnitude(p.pts[i],source,'y') ;
                c.z = scalarprod*c.z -  magnitude(p.pts[i],source,'z') ;

                c.x = c.x + p.pts[i].x;
                c.y = c.y + p.pts[i].y;
                c.z = c.z + p.pts[i].z;

                thetare = anglebetween(p.pts[i],eye,p.pts[i],c);
                printf("\n e dot r: %lf \n",cos(thetare)); 
                intensity = ambient + diffusemax*cos(thetanl) + .5*cos(thetare);


return 3 * intensity;
}
////////////////

void intersect3d(point * new,point A, point B, double Z){
       //line AB at Z
       point AB;
       AB.x = A.x - B.x;
       AB.y = A.y - B.y;
       AB.z = A.z - B.z;

       double t = (Z - A.z) / AB.z;

       new->x = A.x + AB.x * t;
       new->y = A.y + AB.y * t;
       new->z = Z;
}


void clip3(polygon * p, double Z){
        int i,j,k=0,bool1,bool2; 
            point points[100];

        p->pts[p->numpoints].x = p->pts[0].x;
        p->pts[p->numpoints].y = p->pts[0].y;
        p->pts[p->numpoints].z = p->pts[0].z;
        p->numpoints = p->numpoints + 1;
        //put the first point again on the end of the array of points
                
                for(i=1;i<p->numpoints;i++){
                        
                        bool1= (p->pts[i-1].z>Z); //is pts[i-1] farther out than Z            
                        bool2= (p->pts[i].z>Z); //is pts[i] farther out than Z
                        
                        if(bool1&&bool2){//in, in keep first
                                
                            points[k].x = p->pts[i-1].x;
                            points[k].y = p->pts[i-1].y; 
                            points[k].z = p->pts[i-1].z; k++;
                        
                        }else if(!bool1&&bool2){// out, in--get intersection of lines 
                             
                            intersect3d(&points[k],p->pts[i-1],p->pts[i],Z); k++;    
                            //put the intersection of that one line and the z= Z plane on temp arrays          
                        
                        }else if(bool1&&!bool2){//in, out-- keep first, get intersection of lines
                        
                            points[k].x = p->pts[i-1].x;
                            points[k].y = p->pts[i-1].y; 
                            points[k].z = p->pts[i-1].z; k++;
                            
                            intersect3d(&points[k],p->pts[i-1],p->pts[i],Z);    
                            k++;
                        }
                }
                for(i=0;i<k;i++){
                            p->pts[i].x = points[i].x;
                            p->pts[i].y = points[i].y;
                            p->pts[i].z = points[i].z;
                }
                p->numpoints = k;
}










