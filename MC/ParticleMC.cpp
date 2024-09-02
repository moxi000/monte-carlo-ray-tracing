// ParticleMC.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include<omp.h>//OpenMP头文件
#include<random>//随机数头文件
#include<vector>//数据结构，类似数组

#include"Vector.h"//向量头文件
#define PI 3.1415926



using namespace std;



const double AbsorptionCoefficient = 125663.704000000;//谱带衰减系数β=4πk/λ，k为吸收系数，即复折射率里的虚部，λ为波长，（k=0.01，波长为1μm）单位米
const double Ref_Index = 1.5;// 折射率
const int Nthreads = 6;//并行线程数目

#include"Function.h"//函数头文件

vector<Vec3<double> > ParticlePosition;//球的位置
vector<Vec3<double> > Point_range;//xmin ymin zmin xmax ymax zmax 组成的立方体12345678点
vector<double> ParticleRadius;//球的半径
vector<Vec2<double> > F_Position;//光线与入瞳面的交点坐标
vector<double>P_theta;//入瞳光线与入瞳面dA的夹角的cos值
vector<int>s_mesh;//到接受面光线对应的散射角
vector<double>P_length;//入瞳点与粒子散射点的距离的平方

int main()
{
	//Test();
	int ID;//组成大粒子的球球的数量编号
	int N=3;////组成大粒子的球的个数
	ParticlePosition.resize(N);
	Point_range.resize(8);
	double P_L = 0.3;////粒子距离入瞳面P_L=300mm
	double P_R = 0.115;//d=230mm,r=115mm=0.115m 
	double SP[3][5];//组成大粒子的球的个数N=多少，就设置数组SP[N][5]，读入球的参数
	FILE* fp = fopen("shape_sphere.txt", "r");
	//将文档上的值读入数组SP[N][5]，第1列是球的编号，第2-4列是球心坐标，第5列是球的半径，组成大粒子的球的半径可相同，可不同

	if (fp == NULL)
	{
		printf("文件读取错误...");
		system("pause");
	}
    for (int i=0;i<N;i++)
	{
		for (int j=0;j<5;j++)
		{
             fscanf(fp, "%lf", &SP[i][j]);/*每次读取一个数，fscanf函数遇到空格或者换行结束*/
		}
		fscanf(fp, "\n");
	}
		fclose(fp);


	
	for (ID=0;ID<N;ID++)
	{
		ParticlePosition[ID].ReConstruct(SP[ID][1],SP[ID][2],SP[ID][3]);
		ParticleRadius.push_back(SP[ID][4]);
	}
	 //质心
	double c_x = 0.0;
	double c_y = 0.0;
	double c_z = 0.0;
	for (int c_i = 0; c_i < N; c_i++)
	{
		c_x = c_x + SP[c_i][1];
		c_y = c_y + SP[c_i][2];
		c_z = c_z + SP[c_i][3];
	}
	c_x = c_x / N;
	c_y = c_y / N;
	c_z = c_z / N;

	///坐标范围输入
	double Range[9];//定义粒子的坐标范围数组
	FILE* fp1 = fopen("position_range.txt", "r");//读取粒子的坐标范围数组
	if (fp1 == NULL)
	{
		printf("文件读取错误...");
		system("pause");
	}
	//1-9列分别代表 X_min，X_max，Y_min，Y_max，Z_min，Z_max,，X_k（1.1*（X_max-X_min）），Y_k(1.1*（Y_max-Y_min)），Z_k(1.1*（Z_max-Z_min）)
   for (int k=0;k<9;k++)
   {
		
        fscanf(fp1, "%lf",&Range[k]);/*每次读取一个数，fscanf函数遇到空格或者换行结束*/
		fscanf(fp1, "\n");
   }
	fclose(fp1);

	//xmin xmax ymin ymax zmin zmax组成的立方体
	Point_range[0].ReConstruct(Range[0], Range[2], Range[4]);
	Point_range[1].ReConstruct(Range[0], Range[3], Range[4]);
	Point_range[2].ReConstruct(Range[1], Range[2], Range[4]);
	Point_range[3].ReConstruct(Range[1], Range[3], Range[4]);
	Point_range[4].ReConstruct(Range[0], Range[2], Range[5]);
	Point_range[5].ReConstruct(Range[0], Range[3], Range[5]);
	Point_range[6].ReConstruct(Range[1], Range[2], Range[5]);
	Point_range[7].ReConstruct(Range[1], Range[3], Range[5]);

	double inc_angle = 0.5*PI;//照射角度
	double n_inc, p_inc;//平行光入射光线的方向向量
	double y_palmax, y_palmin;//平行光以inc_angle照射角入射后的ymin ymax对应的最大以及最小值
	//求平行光入射到粒子上的入射光线的方向向量，以及ymin ymax
	if ((inc_angle>1e-15)&&(inc_angle < PI)&& (inc_angle!=0.5*PI ))//照射角<180，且不等于90
	{
		
		double t_k[8];
		double x_pal[8];
		double y_pal[8];
		
		if (inc_angle > 0.5 * PI)//照射角>90度，方向向量为负
		{
			n_inc = -sin(PI - inc_angle);//圆周角设为90度，即m_inc=0
			p_inc = -cos(PI - inc_angle);
			//求平行光投影在Z=Zmax=Range[5]的xy范围
			for (int p_i = 0; p_i < 8; p_i++)
			{
				t_k[p_i] = (Range[5] - Point_range[p_i].z) / p_inc;
				x_pal[p_i] = Point_range[p_i].x;
				y_pal[p_i] = t_k[p_i] * n_inc + Point_range[p_i].y;
			}
		}
		else //照射角<90度方向向量y为负，z为正
		{
			n_inc = -sin(inc_angle);//圆周角设为90度，即m_inc=0
			p_inc = cos(inc_angle);
			//求平行光投影在Z=Zmin=Range[4]的xy范围
			for (int p_i = 0; p_i < 8; p_i++)
			{
				t_k[p_i] = (Range[4] - Point_range[p_i].z) / p_inc;
				x_pal[p_i] = Point_range[p_i].x;
				y_pal[p_i] = t_k[p_i] * n_inc + Point_range[p_i].y;
			}
		}

		//求平行光以inc_angle照射角入射后的ymin ymax对应的最大以及最小值
		 y_palmax = y_pal[0];
		 y_palmin = y_pal[0];
		for (int m_i = 1; m_i < 8; m_i++)
		{
			if (y_palmax < y_pal[m_i])
			{
				y_palmax = y_pal[m_i];
			}
			else
			{
				y_palmax = y_palmax;
			}
			if (y_palmin > y_pal[m_i])
			{
				y_palmin = y_pal[m_i];
			}
			else
			{
				y_palmin = y_palmin;
			}

		}

	}
	else if (abs(inc_angle - PI) < 1e-15)//照射角=180度，则从z负轴入射 0，0，-1
	{
		n_inc = 0.0;
		p_inc = -1.0;
	}
	else if (abs(inc_angle - 0.5*PI) < 1e-15)//照射角=90 从y负轴入射 0，-1,0
	{
		n_inc = -1.0;
		p_inc = 0.0;
	}
	else
	{
		n_inc = 0.0;
		p_inc = 1.0;
	}
	cout << n_inc << endl;
	cout << p_inc << endl;




	const int RaysNumber = 2E9;//发射的光线总数
	int eachCycle_number = 1E7;
	
	if (RaysNumber % eachCycle_number != 0)
	{
		printf("错误：总光线数不能被每次循环光线数整除！\n");
		system("pause");
	}

	int Number_Non_valid = 0, Number_absorbed = 0, Number_scattered = 0, Number_reflected = 0, Number_debug = 0;
	int i_r = 0;
	            //int mesh_i=0;
				int Phase_fun[181]={0};//对应散射角的光线数
				int Phase_refr[181]={0};//折射出的光线对应天顶角的光线数
				int Phase_refle[181]={0};//反射出的光线对应天顶角的光线数
				int F_Angle[181][361] = {0};//对应天顶角 圆周角的光线数
				double f_phase_fun[181] = {0.0};//散射相函数
				double I_phase[181] = { 0.0 };//散射强度
				double Func_Angle[181][361] = {0.0};//对应天顶角 圆周角的辐射能

	for (int ii = 0; ii < RaysNumber; ii = ii + eachCycle_number)
	{
		printf("%d of %d w\n", ii / 10000, RaysNumber / 10000);

#pragma omp parallel num_threads(Nthreads)
		{
			int thread_id = omp_get_thread_num();//获取当前所在线程的线程编号
			//default_random_engine rand_engine((time(NULL) % 10001 + thread_id * 777));
			default_random_engine rand_engine(thread_id * 777);

			uniform_real_distribution<double> rand_dist(0, 1);
			//对每个线程分别定义随机数引擎与随机数分布

#pragma omp for
			for (int i = ii; i < ii + eachCycle_number; i++)
			{

				//光线进入介质中时的初始位置
				Vec3<double> R;//发射点
				Vec3<double> R1;//落点
			
				if ((inc_angle<0.5*PI)&& (inc_angle > 1e-15))
				{
					R.x = Range[0] + Range[6] * rand_dist(rand_engine);
					R.y = y_palmin + (y_palmax - y_palmin) * rand_dist(rand_engine);
					R.z = Range[4] - 0.2 * abs(Range[4]);
				}
				else if ((inc_angle > 0.5 * PI) && (inc_angle < PI))
				{
					R.x = Range[0] + Range[6] * rand_dist(rand_engine);
					R.y = y_palmin + (y_palmax - y_palmin) * rand_dist(rand_engine);
					R.z = Range[5] + 0.2 * abs(Range[5]);
				}
				else if (abs(inc_angle - PI) < 1e-15)
				{
					R.x = Range[0] + Range[6] * rand_dist(rand_engine);
		            R.y = Range[2] + Range[7] * rand_dist(rand_engine);
	                R.z = Range[5] + 0.1 * abs(Range[5]);
				}
				else if (abs(inc_angle - 0.5 * PI) < 1e-15)
				{
					R.x = Range[0] + Range[6] * rand_dist(rand_engine);
					R.y = Range[3] + 0.2 * abs(Range[3]);
					R.z = Range[4] + Range[8] * rand_dist(rand_engine);
				}
				//if (inc_angle < 1e-15)
				else
				{
					R.x = Range[0] + Range[6] * rand_dist(rand_engine);
					R.y = Range[2] + Range[7] * rand_dist(rand_engine);
					R.z = Range[4]-0.2 * abs(Range[4]);

				}
				
				
				
				Vec3<double> Direction;//光线传播方向
				Direction.ReConstruct(0.0, n_inc, p_inc);// 入射光线的方向向量
				Vec3<double> Direction1;//光线继续传播方向
				Vec3<double> Emission0;//光线发射方向
				Emission0.ReConstruct(0.0, n_inc, p_inc);// 入射光线的方向向量
				Vec3<double> P_direction;//入瞳点与出射点连线的方向向量
				Vec3<double> Emission(0.0, 0.0, 1.0);//Z轴正方向
				
				int flag = -1;//光线有效标志,-1无效，0散射，1吸收
				int flag_refl = -1;//反射光线有效标志,-1无效，0折射，1反射
				
				double absorb_length = -1.0 / AbsorptionCoefficient * log(1.0 - rand_dist(rand_engine));//光子超出此长度将会被吸收	

				MonteCarlo(R, R1, Direction, Direction1, absorb_length, flag,flag_refl,ParticlePosition, ParticleRadius, rand_engine, rand_dist);
				//调用蒙特卡洛函数，光线追踪
				if (flag == 1)//光线有效标志1，光线被吸收
				{
#pragma omp critical//临界区，防止线程间内存竞争
					{
						Number_absorbed++;//被吸收的光线数+1
					}
				}
				else if (flag == 0)//光线有效标志0，光线被散射
				{
#pragma omp critical//临界区，防止线程间内存竞争
					{
						Number_scattered++;//被散射的光线数+1
						int mesh_i=sca_mesh(Direction1, Emission0);//散射光线的散射角计算
					Phase_fun[mesh_i]=Phase_fun[mesh_i]+1;//散射光线对应散射角的光线数+1
					int phi_i=phi_mesh(Direction1,R1);	//散射光线的天顶角和圆周角计算					
					F_Angle[mesh_i][phi_i] = F_Angle[mesh_i][phi_i] + 1;//散射光线对应的天顶角和圆周角光线数+1
					Direction1 = Direction1 / Direction1.length();//归一化
					///计算辐射亮度
					double Z_L = c_z - P_L;//入瞳面所在位置的纵坐标
					double cos_theta_y = Vec3_Dot(Direction1, Emission);//出射光线与y轴正方向夹角的cos(theta)

					if (cos_theta_y < 1e-15)//散射角大于90°的才能入瞳
					{
						double t_r = ( Z_L-R1.z) / Direction1.z;//直线的参数方程参数t的值
						double X_L = t_r * Direction1.x + R1.x;//X=t*m+x0
						double Y_L = t_r * Direction1.y + R1.y;//Y=t*n+y0
						if ((X_L - c_x) * (X_L - c_x) + (Y_L - c_y) * (Y_L - c_y) <= P_R * P_R)//d=230mm,r=115mm=0.115m   在圆内以及圆圈上
						{
							F_Position.resize(i_r+1);
							double PR_m = R1.x - X_L;
							double PR_n = R1.y - Y_L;
							double PR_k = R1.z - Z_L; 
							double PR_length = pow(PR_m, 2) + pow(PR_n, 2) + pow(PR_k, 2);
							P_length.push_back(PR_length);
						    P_direction.ReConstruct(PR_m, PR_n, PR_k);//入瞳点和出射点连线的方向向量
							P_direction = P_direction / P_direction.length();
							double cos_Ptheta = Vec3_Dot(P_direction, Emission);//入瞳点与出射点连线与Z轴间夹角的cos值
							P_theta.push_back(cos_Ptheta);
							F_Position[i_r].x = X_L;
							F_Position[i_r].y = Y_L;
							s_mesh.push_back(mesh_i);//到接受面光线对应的散射角
						
						
							i_r = i_r + 1;
						}
					}
				
					}
				}
				else
				{
#pragma omp critical//临界区，防止线程间内存竞争
					{
						Number_Non_valid++;
					}
				}

				if (flag_refl == 1)//反射光线标志1，光线被反射
				{
#pragma omp critical//临界区，防止线程间内存竞争
					{
						Number_reflected++;//反射光线数+1
						int mesh_i=sca_mesh(Direction1, Emission0);//反射光线对应的散射角计算
					Phase_refle[mesh_i]=Phase_refle[mesh_i]+1;//反射光线对应散射角的光线数+1
					}
				}
			}
		}
	}

	double TotalNumber = RaysNumber - Number_Non_valid;//成功到达介质的总光线数
	double absorption = Number_absorbed / TotalNumber;//吸收率
	double scattered = Number_scattered / TotalNumber;//散射率
	double fanshe = Number_reflected / TotalNumber;//反射率

	double w_w = (scattered + 1) / (scattered + 1 + absorption);//散射反照率  散射+1=Qsca
	double K_I = w_w / (4 * PI);
	double I0 = 1.0; //入射辐射强度
	//屏幕输出总的有效光线数，吸收率，散射率和反射率，以及散射光线数，反射光线数
	printf("%f\n%f\n%f\n%f\n", TotalNumber, absorption, scattered, fanshe);
	cout<<Number_scattered<<'\t'<<Number_reflected<<endl;
	ofstream out("result.txt");
	//打开文件（result.txt），输出结果 总的有效光线数，吸收光线数，散射光线数，反射光线数，吸收率，散射率和反射率，以及散射光线散射角对应的光线数
	out<<"TotalNumber "<<TotalNumber<<'\t'<<" Number_absorbed "<<Number_absorbed<<'\t'<< " Number_scattered "<<Number_scattered<<'\t'<<" Number_reflected "<<Number_reflected<<endl;
	out<<"absorption "<<absorption<<'\t'<<" scattered "<<scattered<<'\t'<<"fanshe "<<fanshe<<endl;
	//相函数光线
	for ( int mesh = 0; mesh <181; mesh++)
	{
		out<< mesh <<'\t'<<Phase_fun[mesh]<<endl;
	}
	
	out.close();//关闭文件，保存文件。
	
	//求相函数
	f_phase_fun[0] = 4 * Phase_fun[0] / (sin(0.25 * PI / 180) * (PI / 180) * Number_scattered);
	f_phase_fun[180] = 4 * Phase_fun[180] / (sin(179.75 * PI / 180) * (PI / 180) * Number_scattered);
	for (int mesh = 1; mesh < 180; mesh++)
	{
		f_phase_fun[mesh] = 2 * Phase_fun[mesh] / (sin(mesh * PI / 180) * (PI / 180) * Number_scattered);
	}
	//散射方向单位立体角的散射强度
	for (int I_mesh = 0; I_mesh < 181; I_mesh++)
	{
		I_phase[I_mesh] = K_I * I0 * f_phase_fun[I_mesh];
	}

	//输出相函数和对应辐射能
	out.open("f_function.txt");

	for (int mesh = 0; mesh < 181; mesh++)
	{
		out << mesh << '\t' << f_phase_fun[mesh]<<'\t' <<I_phase[mesh]<< endl;
	}
	out.close();//关闭文件，保存文件。
	
	//输出 天顶角 圆周角 对应光线数

	out.open("F_angle.txt");
	for (int i_a = 0; i_a < 181; i_a++)
	{
		for (int i_b = 0; i_b < 361; i_b++)
		{
			out << F_Angle[i_a][i_b] << '\t';

		}
		out << '\n';
	}
	out.close();

	// 计算天顶角 圆周角对应的辐射能

	for (int i_b = 0; i_b < 361; i_b++)
	{
		Func_Angle[0][i_b] = K_I * I0* 4 * F_Angle[0][i_b] / (sin(0.25 * PI / 180) * (PI / 180) * Number_scattered);
	}

	for (int i_b = 0; i_b < 361; i_b++)
	{
		Func_Angle[180][i_b] = K_I * I0* 4 * F_Angle[180][i_b] / (sin(0.25 * PI / 180) * (PI / 180) * Number_scattered);
	}

	for (int i_a = 1; i_a < 180; i_a++)
	{
		for (int i_b = 0; i_b < 361; i_b++)
		{
			Func_Angle[i_a][i_b] = K_I * I0 * 2 * F_Angle[i_a][i_b] / (sin(i_a * PI / 180) * (PI / 180) * Number_scattered);
		}
	}
	//输出 天顶角 圆周角 对应辐射能

	out.open("q_angle.txt");
	for (int i_a = 0; i_a < 181; i_a++)
	{
		for (int i_b = 0; i_b < 361; i_b++)
		{
			out << Func_Angle[i_a][i_b] << '\t';

		}
		out << '\n';
	}
	out.close();


	//到达入瞳面的位置
	out.open("f_position.txt");
	//输出f_position.txt，包含到达入瞳面的位置(x,y),到接受面光线对应的散射角以及入瞳点与出射点连线与A轴间夹角的cos值
	for (int j_r = 0; j_r < F_Position.size(); j_r++)
	{
		out << j_r << '\t' << F_Position[j_r].x << '\t' << F_Position[j_r].y << '\t' << s_mesh[j_r] <<'\t'<< P_theta[j_r] << '\t'<< P_length[j_r] << endl;
	}
	out.close();

	double d_l = 0.01;//接收面格子的间距d=10mm=10/1000m；
	const int N_d = 24;//分成网格子的边长格子数ceil(2*P_R/d_l)+1
	double L_I[N_d][N_d] = { 0.0 };//格子分布矩阵（矩形的）
	for (int j_r = 0; j_r < F_Position.size(); j_r++)
	{
		int m_l = ceil((F_Position[j_r].x + c_x + P_R) / (d_l));//x值属于第几列格子，m_l代表列
		int n_l = ceil((F_Position[j_r].y - c_y - P_R) / (-d_l));//y值属于第几行格子，n_l代表行
		
		if (s_mesh[j_r] < 1e-15)
		{
			L_I[n_l][m_l] = L_I[n_l][m_l] + 2 * w_w * I0 * P_theta[j_r] * (1 / (2 * PI * sin(0.25 * PI / 180) * (PI / 180) * Number_scattered)) / P_length[j_r];
		}
		else if ((s_mesh[j_r] - 180) < 1e-15)
		{
			L_I[n_l][m_l] = L_I[n_l][m_l] + 2 * w_w * I0 * P_theta[j_r] * (1 / (2 * PI * sin(179.25 * PI / 180) * (PI / 180) * Number_scattered)) / P_length[j_r];
		}
		else
		{
			L_I[n_l][m_l] = L_I[n_l][m_l] + w_w * I0 * P_theta[j_r] * (1 / (2 * PI * sin(s_mesh[j_r] * PI / 180) * (PI / 180) * Number_scattered)) / P_length[j_r];
		}
		
		
	}
	out.open("L_I.txt");
	//输出L_I.txt，将入瞳面划分为N_d*N_d个网格，单位面积单位立体角所接受的辐射能
	for (int i_a = 0; i_a < N_d; i_a++)
	{
		for (int i_b = 0; i_b < N_d; i_b++)
		{
			out << L_I[i_a][i_b] << '\t';
		}
		out << '\n';
	}
	out.close();




	
    cout << "Hello World!\n";

	system("pause");
}

