/*      -----TEST - A sous forme de matrice    

#include <conio.h>
#include <iostream.h>
#include <math.h>


// Création d'une matrice(n,m)

	double **ZeroMatrix(int n, int m)
		{
   		 double **result = (double **) malloc(n * sizeof(double*));
    		for (int row = 0; row < n; row++)
       		 result[row] = (double *) calloc(m , sizeof(double));
    		return result;
		}


// Création de la matrice A (avec les coeffs et conditions aux bords)          taille de A ????

 	double **A(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments)
	{
	std::vector<double> a=m.spot_vector();
        int size = a.size();

	  double **result = ZeroMatrix(size+1, size-1);
    		
		for (int i = 0; i < size+1; i++) {
		for (int j=0; j<size-1 ;j++) { 

       		 if (i == 0 || i == dim - 1) 
		 {
           	 result[i][j] = 1;            // conditions aux bords    - quelle valeur ??
       		 } 
		
		else if(i == j) 
		{
           	 result[i][j] = 1.0+arguments[4]*m.get_mesh_dt()*((pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+rate.function_operator(arguments));        // beta sur la diagonale 
       		} 


		else if(i == (j-1)) 
		{
           	 result[i][j] = ((-1.0 / 2.0) * (arguments[4] / (m.get_mesh_dx()*m.get_mesh_dx()))*pow(vol.function_operator(arguments), 2) + (1.0 / (4.0 * m.get_mesh_dx()))*arguments[4] * (pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments)));    //updiagonale 
       		} 


		else if(i == (j+1)) 
		{
           	 result[i][j] = -0.5*arguments[4]*m.get_mesh_dt()*((pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+((pow(vol.function_operator(arguments),2)-rate.function_operator(arguments))/(2.0*m.get_mesh_dx())));     //subdiagonale 
       		} 


	}
    }
    return result;
}






// Inversion de la matrice           <<>> cin pour paramètres entrés par l utilisateur (cf theta) ?



int scanf(float a[100][100]){
	int i,j,n;
	cout<<"\n Enter Length Of Matrix N*N : ";
	cin>>n;
	cout<<"\n--------------------------\n";
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			cout<<" Matrix["<<i+1<<"]["<<j+1<<"] : ";
			cin>>a[i][j];
		}
	cout<<"\n----------------------------------------------------\n";
return n;
}



//calculate minor of matrix OR build new matrix : k-had = minor
void minor(float b[100][100],float a[100][100],int i,int n){
	int j,l,h=0,k=0;
	for(l=1;l<n;l++)
		for( j=0;j<n;j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
}// end function


//calculate determinte of matrix
float det(float a[100][100],int n){
	int i;
	float b[100][100],sum=0;
	if (n == 1)
return a[0][0];
	else if(n == 2)
return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	else
		for(i=0;i<n;i++){
			minor(b,a,i,n);	// read function
			sum = (float) (sum+a[0][i]*pow(-1,i)*det(b,(n-1)));	// read function	// sum = determinte matrix
		}
return sum;
}// end function





//calculate transpose of matrix
void transpose(float c[100][100],float d[100][100],int n,float det){
	int i,j;
	float b[100][100];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			b[i][j] = c[j][i];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			d[i][j] = b[i][j]/det;	// array d[][] = inverse matrix
}// end function





//calculate cofactor of matrix
void cofactor(float a[100][100],float d[100][100],int n,float determinte){
	float b[100][100],c[100][100];
	int l,h,m,k,i,j;
	for (h=0;h<n;h++)
		for (l=0;l<n;l++){
			m=0;
			k=0;
			for (i=0;i<n;i++)
				for (j=0;j<n;j++)
					if (i != h && j != l){
						b[m][k]=a[i][j];
						if (k<(n-2))
							k++;
						else{
							k=0;
							m++;
						}
					}
			c[h][l] = (float) pow(-1,(h+l))*det(b,(n-1));	// c = cofactor Matrix
		}
	transpose(c,d,n,determinte);	// read function
}// end function





//calculate inverse of matrix
void inverse(float a[100][100],float d[100][100],int n,float det){
	if(det == 0)
		cout<<"\nInverse of Entered Matrix is not possible\n";
	else if(n == 1)
		d[0][0] = 1;
	else
		cofactor(a,d,n,det);	// read function
}// end function




*/