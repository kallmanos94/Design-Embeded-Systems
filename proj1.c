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
