#include <iostream>
#include <fstream>
#include <math.h>
// #include<bits/stdc++.h>

using namespace std;

double minimum(double a, double b, double c)
{
    double mini;
    if( (a<b) && (a<c) )
        mini = a;
    else if(b<c)
        mini = b;
    else
        mini = c;
    return mini;
}
double minimum(double a, double b)
{
    double mini;
    if(a<b)
        mini = a;
    else
        mini = b;
    return mini;
}
int findminindexw(double *check, int w){
    int i;
    double mini = check[0];  
    int mindex=0;  
    for (i=1; i<w; i++) 
    {       
        if(check[i] < mini)  
        {  
        mini = check[i]; 
        mindex = i;  
        }
    }  
    return mindex;
}
int findminindexh(double *check2, int h){
    int i;
    double mini = check2[0];  
    int mindex=0;  
    for (i=1; i<h; i++) 
    {       
        if(check2[i] < mini)  
        {  
            mini = check2[i]; 
            mindex = i;  
        }
    }   
    return mindex;  
}

void decrease_width(double **cost, int *removew, double *check, int ***rgb, int h, int w, int c, int h_, int w_, int c_)
{
    int i,j,k,l;			// i=h, j=w, k=c
    int r,g,b;
    double ans;


    // Making Energy and Cost Matrices Together
    for(i=0; i<h; i++)
    {
        for(j=0; j<w; j++)
        {
            int rw_10 = rgb[i][(w+j-1)%w][0];
            int rw10 = rgb[i][(j+1)%w][0];
            int rh10 = rgb[(i+1)%h][j][0];
            int rh_10 = rgb[(h+i-1)%h][j][0];
            int rw_11 = rgb[i][(w+j-1)%w][1];
            int rw11 = rgb[i][(j+1)%w][1];
            int rh_11 = rgb[(h+i-1)%h][j][1];
            int rh11 = rgb[(i+1)%h][j][1];
            int rw_12 = rgb[i][(w+j-1)%w][2];
            int rw12 = rgb[i][(j+1)%w][2];
            int rh_12 = rgb[(h+i-1)%h][j][2];
            int rh12 = rgb[(i+1)%h][j][2];

            double x = ( (rw_10-rw10)*(rw_10-rw10) + (rw_11-rw11)*(rw_11-rw11) + (rw_12-rw12)*(rw_12-rw12) );
            double y = ( (rh_10-rh10)*(rh_10-rh10) + (rh_11-rh11)*(rw_11-rw11) + (rh_12-rh12)*(rh_12-rh12) );
            cost[i][j] = sqrt(x+y);
            
            if(i==0)
                continue;
            else
            {
                if(j == 0)
                    cost[i][j] = minimum( cost[i-1][j], cost[i-1][j+1] ) + cost[i][j];
                else if(j == w-1 )
                    cost[i][j] = minimum( cost[i-1][j-1], cost[i-1][j] ) + cost[i][j];
                else
                    cost[i][j] = minimum( cost[i-1][j-1], cost[i-1][j], cost[i-1][j+1] ) + cost[i][j];
            }
        }
    }

    // Taking value for h to be removed
    // Copying last row

    for(i=0; i<w; i++)
        check[i]=cost[h-1][i];
    

    // Finding minimum index of neighbour's upper row

    int minindex = findminindexw(check, w);
    removew[h-1] = minindex;
    for(i=h-2; i>=0; i--)
    {
        if(minindex == 0)
            ans = minimum( cost[i][minindex], cost[i][minindex+1] );
        else if(minindex == w-1)
            ans = minimum( cost[i][minindex-1], cost[i][minindex] );
        else
            ans = minimum( cost[i][minindex-1], cost[i][minindex], cost[i][minindex+1] );

        if(ans == cost[i][minindex-1])
            minindex = minindex-1;
        else if(ans == cost[i][minindex])
            minindex = minindex;
        else if(ans == cost[i][minindex+1])
            minindex = minindex+1;
        removew[i] = minindex;
    }

    // Removing the single seam
    for(i=0; i<h; i++)
    {
        l = removew[i];
        for(j=l; j<w-1; j++)
            rgb[i][j] = rgb[i][j+1];
        
    }

}

void decrease_height(double **cost, int *removeh, double *check2, int ***rgb, int h, int w, int c, int h_, int w_, int c_)
{
    int i,j,k,l;			// i=h, j=w, k=c
    int r,g,b;
    double ans;

    // Making Energy and Cost Matrices
    for(j=0; j<w; j++)		// Check this convention during Segmentation Fault
    {
        for(i=0; i<h; i++)
        {           
            int rw_10 = rgb[i][(w+j-1)%w][0];
            int rw10 = rgb[i][(j+1)%w][0];
            int rh_10 = rgb[(h+i-1)%h][j][0];
            int rh10 = rgb[(i+1)%h][j][0];
            int rw_11 = rgb[i][(w+j-1)%w][1];
            int rw11 = rgb[i][(j+1)%w][1];
            int rh_11 = rgb[(h+i-1)%h][j][1];
            int rh11 = rgb[(i+1)%h][j][1];
            int rw_12 = rgb[i][(w+j-1)%w][2];
            int rw12 = rgb[i][(j+1)%w][2];
            int rh_12 = rgb[(h+i-1)%h][j][2];
            int rh12 = rgb[(i+1)%h][j][2];
            double x = ( (rw_10-rw10)*(rw_10-rw10) + (rw_11-rw11)*(rw_11-rw11) + (rw_12-rw12)*(rw_12-rw12) );
            double y = ( (rh_10-rh10)*(rh_10-rh10) + (rh_11-rh11)*(rw_11-rw11) + (rh_12-rh12)*(rh_12-rh12) );
            cost[i][j] = sqrt(x+y);

            if(j==0)
                continue;
            else
            {
                if(i == 0)
                    cost[i][j] = minimum( cost[i][j-1], cost[i+1][j-1] ) + cost[i][j];
                else if(i == h-1 )
                    cost[i][j] = minimum( cost[i-1][j-1], cost[i][j-1] ) + cost[i][j];
                else
                    cost[i][j] = minimum( cost[i-1][j-1], cost[i][j-1], cost[i+1][j-1] ) + cost[i][j];
            }
        }
    }

    // Taking value for w to be removed
    // Copying last column
    for(i=0; i<h; i++)
    {
        check2[i]=cost[i][w-1];
    }

    // Finding minimum index of last column    
    // Finding minimum index of neighbour's left row
    int minindex = findminindexh(check2, h);
    removeh[w-1] = minindex;  

    for(j=w-2; j>=0; j--)
    {
        if(minindex == 0)
            ans = minimum( cost[minindex][j], cost[minindex+1][j] );
        else if(minindex == h-1)
            ans = minimum( cost[minindex-1][j], cost[minindex][j] );
        else
            ans = minimum( cost[minindex-1][j], cost[minindex][j], cost[minindex+1][j] );
            
        if(ans == cost[minindex-1][j])
            minindex = minindex-1;
        else if(ans == cost[minindex][j])
            minindex = minindex;
        else if(ans == cost[minindex+1][j])
            minindex = minindex+1;
        removeh[j] = minindex;
    }

    // Removing the single seam
    for(j=0; j<w; j++)
    {
        l = removeh[j];
        for(i=l; i<h-1; i++)
            rgb[i][j] = rgb[i+1][j];
    }
}

void solve(int ***rgb, int h, int w, int c, int h_, int w_, int c_) {

    double **cost;
    cost = new double *[h];

    for(int i = 0; i < h; ++i) 
        cost[i] = new double [w];
    

    int *removew;
    removew = new int [h];

    double *check;
    check = new double [w];

    int *removeh;
    removeh = new int [w];

    double *check2;
    check2 = new double [h];


// Used for Decreasing Width
    while(w_<w) {   
        decrease_width(cost, removew, check, rgb, h, w, c, h_, w_, c_);
        w--;
    }

// Used for Decreasing Height
    while(h_<h){
        decrease_height(cost, removeh, check2, rgb, h, w, c, h_, w_, c_);
        h--;
    }
}   

int main() {
    string ip_dir = "./data/input/";
    string ip_file = "rgb_in.txt";
    ifstream fin(ip_dir + ip_file);

    int H, W, C;
    fin >> H >> W >> C;

    int ***rgb;
    rgb = new int **[H];
    for(int i = 0; i < H; ++i) {
        rgb[i] = new int *[W];
        for(int j = 0; j < W; ++j) {
            rgb[i][j] = new int[C];
            for(int k = 0; k < C; ++k)
                fin >> rgb[i][j][k];
        }
    }
    fin.close();

    int H_, W_, C_;
    cout << "Enter new value for H (must be less than " << H << "): ";
    cin >> H_;
    cout << "Enter new value for W (must be less than " << W << "): ";
    cin >> W_;
    cout << '\n';
    C_ = C;

    solve(rgb, H, W, C, H_, W_, C_);

    string op_dir = "./data/output/";
    string op_file = "rgb_out.txt";
    ofstream fout(op_dir + op_file);
    
    fout << H_ << " " << W_ << " " << C_ << endl;
    for(int i = 0; i < H_; ++i) {
        for(int j = 0; j < W_; ++j) {
            for(int k = 0; k < C_; ++k) {
                fout << rgb[i][j][k] << " ";
            }
        }
        fout << '\n';
    }
    fout.close();

    return 0;
}

