
#include <D2d_matrix.h>



/*

   ( x')          (x)
   ( y')  =   M * (y)  
   ( 1 )          (1)

   instead of (x',y',1) = (x,y,1) * M  

*/

int D2d_mat_mult_points(double *X, double *Y,double m[3][3], double *x, double *y, int n){
        double ytemp[n] , xtemp[n] ;
        int i;
        for(i=0;i<n;i++){      
                xtemp[i]=  m[0][0] * x[i] + m[0][1] * y[i] + m[0][2];
                ytemp[i]=  m[1][0] * x[i] + m[1][1] * y[i] + m[1][2]; 
        }
        for(i=0;i<n;i++){
            X[i]=xtemp[i];
            Y[i]=ytemp[i];
        }
        
        return 1;
}

int D2d_negate_y(double a[3][3], double b[3][3]){
        
        double temp[3][3];
        int i,j,k;
        for(i=0;i<3;i++){
                for(j=0;j<2;j++){
                temp[i][j] = 0;  
                }
        }
        temp[0][0] = 1;
        temp[1][1] = -1;
        temp[2][2]=1;
        D2d_copy_mat(b,a);
        D2d_mat_mult(a,temp,a);
        return 1;
}



int D2d_negate_x(double a[3][3], double b[3][3]){ 
        double temp[3][3];
        int i,j,k;
        for(i=0;i<3;i++){
                for(j=0;j<2;j++){
                        if(i==j){
                                temp[i][k] = 1;
                        }else{
                                temp[i][k] = 0;
                        }
                }
        }
        temp[0][0] = -1;
        D2d_copy_mat(b,a);
        D2d_mat_mult(a,temp,a);
        return 1;
}



int D2d_scale(double a[3][3], double b[3][3], double sx, double sy){
        double temp[3][3];
        D2d_make_identity(temp);
        temp[0][0] = sx;
        temp[1][1] = sy;
        D2d_mat_mult(a, temp, a);
        temp[0][0] = 1/sx;
        temp[1][1] = 1/sy;
        D2d_mat_mult(b, temp, b);
        return 1;
}
int D2d_print_mat (double a[3][3])
{
        int r,c ;
        for (r = 0 ; r < 3 ; r++ ) {
                for (c = 0 ; c < 3 ; c++ ) {
                        printf(" %12.4lf ",a[r][c]) ;
                }
                printf("\n") ;
        }

        return 1 ;
} 
//WORKS
int D2d_mat_mult(double res[3][3], double a[3][3], double b[3][3]){
        //set value of result to be a * b, and must be able to use same parameter for res and a or b
        double temp[3][3];
        int i,j,k,l;
        double sum=0;
        for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                        for(k=0;k<3;k++){
                                sum += a[i][k] * b[k][j];  
                        }
                        temp[i][j] = sum;
                        sum=0;
                }
        }
        for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                        res[i][j] = temp[i][j];
                }
        }
        return 1;
}

int D2d_copy_mat (double a[3][3], double b[3][3])
        // a = b
{
        int r,c ;
        for (r = 0 ; r < 3 ; r++ ) {
                for (c = 0 ; c < 3 ; c++ ) {
                        a[r][c] = b[r][c] ;
                }
        }

        return 1 ;
} 

int D2d_make_identity (double a[3][3])
        // a = I
{
        int r,c ;
        for (r = 0 ; r < 3 ; r++ ) {
                for (c = 0 ; c < 3 ; c++ ) {
                        if (r == c) a[r][c] = 1.0 ;
                        else    a[r][c] = 0.0 ;
                }
        }
        return 1 ;
} 

int D2d_rotate(double a[3][3], double b[3][3], double radians){
        double temp[3][3];
        double temp2[3][3];
        int h,i,j;
        for(h=0;h<2;h++){
                for(i=0;i<3;i++){
                        for(j=0;j<3;j++){
                                if(i==j&&i!=2) temp[i][j] = cos(radians);
                                if(i==1&&j==0) temp[i][j] = sin(radians);
                                if(i==0&&j==1) temp[i][j] = -1 * sin(radians);
                                if(i==2||j==2) temp[i][j]=0;
                                if(i==2&&j==2) temp[i][j]=1;
                        }
                }
                // D2d_print_mat(temp);

                if(h==0){
                        //puts("matrix a");
                        //D2d_print_mat(a);
                        //puts("rotational matrix");
                        //D2d_print_mat(temp);
                        D2d_mat_mult(temp2,temp,a);
                        radians = -1 * radians;
                }else if(h==1){
                        D2d_mat_mult(b,temp,b);

                }
        }
        D2d_copy_mat(a,temp2);

        return 1;
}   

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


int D2d_translate (double a[3][3], double b[3][3], double dx, double dy)
        // a = translation*a  
        // b = b*translation_inverse  
{
        double t[3][3] ;

        D2d_make_identity(t) ;

        t[0][2] =  dx ;  t[1][2] = dy ;  
        D2d_mat_mult(a,  t,a) ;

        t[0][2] = -dx ;  t[1][2] = -dy ;
        D2d_mat_mult(b,  b,t) ;

        return 1 ;
}

