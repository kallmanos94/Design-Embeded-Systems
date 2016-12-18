#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define IW 352
#define IH 288
#define EIW 356
#define EIH 292
#define N 5
#define S 0.88
#define K 3
#define filename "ak.yuv"
#define finalfilename "akout.yuv"

int i,j;
int image[IH][IW];
double is[IH][IW];
double g[N][N];
int Rep[EIH][EIW];

void read(){
  FILE *frame_c;
  if((frame_c=fopen(filename,"rb"))==NULL)
  {
    printf("current frame doesn't exist\n");
    exit(-1);
  }

  for(i=0;i<IH;i++)
  {
    for(j=0;j<IW;j++)
    {
      image[i][j]=fgetc(frame_c);
    }
  }
  fclose(frame_c);
}

void write(){
  FILE *frame_y;
  frame_y=fopen(finalfilename,"wb");

  for(i=0;i<IH;i++)
  {
    for(j=0;j<IW;j++)
    {
      fputc(image[i][j],frame_y);
    }
  }
  fclose(frame_y);
}

void calfin(){
    for(i=0;i<IH;i++){
        for(j=0;j<IW;j++){
            double temp = image[i][j] + K*( image[i][j] - is[i][j] );
            if(temp>255) image[i][j] = 255;
            else if(temp<0) image[i][j] = 0;
            else image[i][j] = temp;
        }
    }
}

void gauss(){
    double c = exp(N*N/(S*S));

    int a = -N/2;
    int b;

    for(i=0;i<N;i++){
        b = -N/2;
        for(j=0;j<N;j++){
            g[i][j] = c * exp( -(a*a+b*b) / (2*S*S) );
            g[i][j] = g[i][j]/501187233627271,5;
            b++;
        }
        a++;
    }
}

void rot90(){
    double tmp;
    for(i=0; i<N/2; i++){
        for(j=i; j<N-i-1; j++){
            tmp=g[i][j];
            g[i][j]=g[j][N-i-1];
            g[j][N-i-1]=g[N-i-1][N-j-1];
            g[N-i-1][N-j-1]=g[N-j-1][i];
            g[N-j-1][i]=tmp;
        }
    }
    for(i=0; i<N/2; i++){
        for(j=i; j<N-i-1; j++){
            tmp=g[i][j];
            g[i][j]=g[j][N-i-1];
            g[j][N-i-1]=g[N-i-1][N-j-1];
            g[N-i-1][N-j-1]=g[N-j-1][i];
            g[N-j-1][i]=tmp;
        }
    }
}

void conv2(){
    int center = N/2;
    int left = center-1;
    int top = center - 1;
    int x,y;
    rot90();
    rot90();
    
    for(i=top+1;i<=IH+top;i++){
        for(j=1+left;j<=IW+left;j++){
            Rep[i][j]= image[i-top-1][j-left-1];
        }
    }
    for(x=0;x<IH;x++){
        for(y=0;y<IW;y++){
            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    is[x][y]=is[x][y]+( Rep[i+x][j+y] * g[i][j] );
                }
            }
        }
    }
}

int main(){
    read();
    gauss();
    conv2();
    calfin();
    write();

    return(0);
}
