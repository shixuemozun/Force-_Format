#include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
using namespace std;
const int Q=9;          //D2Q9ģ��   //��c++�﷨�У�������������ֱ�д����e[Q][2]��Ϊ����ı�����ֻ����Ϊ���� const ,������ֱ��д�� int Q = 9; 
						 
const int NX=256;       
const int NY=256;       
const double U=0.01;     

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; 
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double F_[NX+1][NY+1][2],G[NX+1][NY+1][Q]; 


double rho[NX+1][NY+1],u[NX+1][NY+1][2],F[NX+1][NY+1][Q],f[NX+1][NY+1][Q],f_[NX+1][NY+1][Q],feq[NX+1][NY+1][Q]; 
	
double P[NX+1][NY+1];	//ѹǿ 

int i,j,k,ip,jp,n; 
double c,cs,Re,dx,dy,Lx,Ly,dt,rho0,T0,P0,tau_f,niu,error;  

double eu,uv;
double rho1,u0,u1; 

double t1=0;

int main()
{
	
	dx=1.0; 
 	dy=1.0;
 	Lx=dx*double(NX); 
 	Ly=dy*double(NY);
 	dt=dx;
 	c=dx/dt;
 
 	rho0=1.0;
 	
 	Re=500;
 	niu=U*Lx/Re;
 	//cout << niu; 
 	tau_f=3.0*niu+0.5;
 	cout <<tau_f;

for(i=0;i<=NX;i++) 
	{
	
	 for(j=0;j<=NY;j++) 
	 {
	 	u[i][j][0]=0;   //��ʼ��������ÿ�����x������ٶ�;
	 	u[i][j][1]=0;   //��ʼ��������ÿ�����y������ٶ�; 
	 	
	 	rho[i][j]=rho0;
	 
	 	u[i][NY][0]=U;  //�Զ��ǵ��ٶȽ��г�ʼ��; 
	 	
	 	
	 									      
	  }
	
}

//�Զ��Ǳ�������ʼ��

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		
		F_[i][j][0] = 0;
		F_[i][j][1] = 0;
		F_[0][j][0] = 0;
		F_[0][j][1] = 0;
		F_[NX][j][0] = 0;
		F_[NX][j][1] = 0;
		
	}

	
 } 


//ѹ����ʼ��

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		P[i][j] = 0;
		
		
	}
	
	
 } 







for( i=0;i<=NX;i++)
	{
		for( j=0;j<=NY;j++)
		{
			
			for( k=0;k<9;k++)
			{
				eu = (e[k][0]*u[i][j][0]+e[k][1]*u[i][j][1]);
				
				uv =(u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1]); 
				
				
				feq[i][j][k] = w[k]*rho[i][j]*(1.0+3.0*eu+4.5*eu*eu-1.5*uv); 
			}	
		 } 
			
	 } 

//�������� G

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		for(k=0;k<9;k++)
		{
		
		G[i][j][k] = (((e[k][0]*F_[i][j][0]+e[k][1]*F_[i][j][1])-(u[i][j][0]*F_[i][j][0]+u[i][j][1]*F_[i][j][1]))*feq[i][j][k])/(rho0*c*c);
	}
		
	}
	
	
	
 } 


	for(i=0;i<=NX;i++) 
	{
	
	 for(j=0;j<=NY;j++) 
	 {

	 	for(k=0;k<Q;k++)
	 	{
	 		f[i][j][k]=feq[i][j][k];
   				
			                                    
		}									      
	  }
	
}

//�µķֲ�����

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		
		for(k=0;k<=Q;k++)

           {
			
			f_[i][j][k] = f[i][j][k] - 0.5*G[i][j][k];
			
			
				   }		
		
	}
	
	
 } 



	for(i=0;i<=NX;i++) 
	{
	
	 for(j=0;j<=NY;j++) 
	 {

	 	for(k=0;k<Q;k++)
	 	{
	 		f[i][j][k]=f_[i][j][k];
   				
			                                    
		}									      
	  }
	
}









for( int t=0;t<100000;t++)
{

cout << t << endl; 
 	
//if(t==0)
//{

//}




/*	
	for(int i=0;i<=NX;i++)
	 {
	 	for(int j=0;j<=NY;j++)
	 	{
	 		for(int k=0;k<9;k++)
	 		{
			 
	 		cout << f[i][j][k]<<"===="; 
	 	}
	 	
		 cout << endl;
		 	
		 }
	 	
	 	
	 }
	*/
	
	for(i=1;i<NX;i++)//�ݻ�
	{
	
	  for(j=1;j<NY;j++)
	  {
	  
	    for(k=0;k<Q;k++)
	    {
		
	    
	    	ip=i-e[k][0];
	    	jp=j-e[k][1];
       
F[i][j][k]=f[ip][jp][k]+(feq[ip][jp][k]-f[ip][jp][k])/0.8+(1-double(1)/(2*0.8))*G[i][j][k];	  //����ȷ���ɳ�ʱ��ᵼ�¼��㷢ɢ 
}
}	
}


for(i=1;i<NX;i++)//��������(�ڲ�����) 
{
	
	
	    for(j=1;j<NY;j++)
		{
	 
            for(k=0;k<Q;k++)
			{
				f[i][j][k]=F[i][j][k];

		}
		
		
		 rho1 = 0;
	   for( k=0;k<9;k++)
	   
	   {
	   rho1 = rho1+f[i][j][k];
	//   cout << f[i][j][k]<<"====";
       }
      rho[i][j]=rho1;
	  
	   u0 = 0;
	   for( k=0;k<9;k++)
	   {
	   
	   u0 += e[k][0]*f[i][j][k]+0.5*F_[i][j][0];
	   //cout << u0<<"====";
	}
	 u[i][j][0] = u0/rho[i][j]; //����x������ٶ�
	  //cout << u[i][j][0];
	  
	  
	   u1 = 0;
	  for( k=0;k<9;k++)
	   {
	   
	   u1 += e[k][1]*f[i][j][k]+0.5*F_[i][j][1];
	  // cout << u1 << "===";
	}
	 u[i][j][1] = u1/rho[i][j]; //����y������ٶ�
	  //cout << u[i][j][1]<<"===";
	  }

}

//���ұ߽� 
 for(j=1;j<NY;j++)
 {
 	for(k=0;k<Q;k++)
		 {
		 	rho[NX][j] = rho[NX-1][j];
		 	f[NX][j][k] = feq[NX][j][k]+(f[NX-1][j][k]-feq[NX-1][j][k]); //�߽翼����ײ   
 	        
			 
			 rho[0][j]=rho[1][j];
 	        f[0][j][k] = feq[0][j][k]+(-1*feq[1][j][k]+f[1][j][k]);
			 
			 
 	
 }
 	
 }

for(i=0;i<=NX;i++)//���±߽�
{
		  for(k=0;k<Q;k++)
		  {
	
	rho[i][0]=rho[i][1];
	f[i][0][k] = feq[i][0][k]+(-1*feq[i][1][k]+f[i][1][k]);
	

	rho[i][NY]=rho[i][NY-1];
	f[i][NY][k] = feq[i][NY][k]+(-1*feq[i][NY-1][k]+f[i][NY-1][k]);
	 
	 u[i][NY][0]=U;
}


}


//����ѹ��

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		P[i][j] = rho[i][j]*c;
		
		
		
	}
	
	
	
	
 } 






	for( i=0;i<=NX;i++)
	{
		for( j=0;j<=NY;j++)
		{
			
			for( k=0;k<9;k++)
			{
				eu = (e[k][0]*u[i][j][0]+e[k][1]*u[i][j][1]);
				
				uv =(u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1]); 
				
				
				feq[i][j][k] = w[k]*rho[i][j]*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);   //����ĩβ������ƽ��̬�ֲ����� 
			}	
		 } 
			
	 } 

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		for(k=0;k<9;k++)
		{
		
		G[i][j][k] = (((e[k][0]*F_[i][j][0]+e[k][1]*F_[i][j][1])-(u[i][j][0]*F_[i][j][0]+u[i][j][1]*F_[i][j][1]))*feq[i][j][k])/(rho0*c*c);
	}
		
	}
	
 } //����   �������� 



//���� �µķֲ�����

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		
		for(k=0;k<=Q;k++)

           {
			
			f_[i][j][k] = f[i][j][k] - 0.5*G[i][j][k];
			
			
				   }		
		
	}
	
	
 } 



	for(i=0;i<=NX;i++) 
	{
	
	 for(j=0;j<=NY;j++) 
	 {

	 	for(k=0;k<Q;k++)
	 	{
	 		f[i][j][k]=f_[i][j][k];
   				
			                                    
		}									      
	  }
	
}





}

 ostringstream name;
	  name<<"cavity_"<<6<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n" << "ZONE T=\"BOX\",I=" << NX+1 << ",J=" << NY+1 << ",F=POINT" << endl;
	  for(int j=0;j<=NY;j++)
	     for(int i=0;i<=NX;i++)
	     {
	     	out<<double(i)/Lx<<" "<<double(j)/Ly<<" "<<u[i][j][0]<<" "<<u[i][j][1]<<" "<<P[i][j] <<endl;
		 }




	return 0;
}
