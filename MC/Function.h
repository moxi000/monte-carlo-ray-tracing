#pragma once

using namespace std;
bool SphereRay(Vec3<double> Position, Vec3<double> Direction, Vec3<double> O, double R, double& t);
void Refraction(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, double n1, double n2, vector<Vec3<double> > ParticlePosition, int SphereIndex);
void Reflection(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, vector<Vec3<double> > ParticlePosition, int SphereIndex);
double reflectivity(double Angle_i, double Angle_t);
int sca_mesh(Vec3<double>Direction1, Vec3<double>Emission0);
int phi_mesh(Vec3<double>Direction1, Vec3<double>R1);



//cos转换成sin函数
inline double CosToSin(double cos)
{
	return 	sqrt(1 - cos * cos);
}
//优化反cos函数
inline double acos_optimized(double value)
{
	if (value < -1.0)value = -1.0;
	else if (value > 1.0)value = 1.0;
	value = acos(value);
	return value;
}



bool InOtherSphere(vector<Vec3<double> > ParticlePosition, vector<double> ParticleRadius, Vec3<double> Position, int sphereIndex)
{
	//除当前所在球之外，是否在其他球之中
	for (int i = 0; i < ParticlePosition.size(); i++)
	{
		if (i != sphereIndex)
		{
			if ((Position - ParticlePosition[i]).length() <= ParticleRadius[i] + 1e-15) return true;
		}
	}
	return false;
}

////判断光线是否与粒子相交，并确定交点在那个单元球

bool Intersection(Vec3<double> Position, Vec3<double> Direction, double& t, int& sphereIndex, vector<Vec3<double> > ParticlePosition, vector<double> ParticleRadius)
{
	t = 999999;
	double t_temp;
	for (int i = 0; i < ParticlePosition.size(); i++)
	{
		if (SphereRay(Position, Direction, ParticlePosition[i], ParticleRadius[i], t_temp))
		{
			if (t_temp < t && t_temp>1e-15)
			{
				t = t_temp;
				sphereIndex = i;//交点球编号
			}
		}
	}
	if (t > 99999) return false;
	else return true;
}
//求光线与球的交点，得到参数t
bool SphereRay(Vec3<double> Position, Vec3<double> Direction, Vec3<double> O, double R, double& t)
{
	Vec3<double> temp = Position - O;//A,B,C


	double a = Vec3_Dot(Direction, Direction);
	double b = 2.0 * Vec3_Dot(temp, Direction);
	double c = Vec3_Dot(temp, temp) - R * R;

	double delta = b * b - 4.0 * a * c;

	

	if (delta < -1e-16) return false;
	else if (abs(delta) < 1e-16)
	{
		t = -b / 2.0 / a;
	}
	else
	{
		double t_temp = (-b - sqrt(delta)) / 2.0 / a;

		t = (-b + sqrt(delta)) / 2.0 / a;

		if (t_temp < 0 && t < 0) return false;//反向无效  t<0 & t_temp<0
		else if (t < 1e-15) t = t_temp;//t<0 & t_temp>0
		else
		{
			if (t_temp > 1e-15)//t>0 & t_temp>0
			{
				if (t_temp < t) t = t_temp;//二者都大于0，取更小的一个
			}
		}
	}
	//InterPosition = Position + Direction * t;
	return true;
}
//调用蒙特卡洛函数，光线追踪
void MonteCarlo(Vec3<double> Position, Vec3<double>& R1, Vec3<double> Direction, Vec3<double>& Direction1, double absorb_length, int& flag, int& flag_refl, vector<Vec3<double> > ParticlePosition, vector<double> ParticleRadius, default_random_engine rand_engine, uniform_real_distribution<double> rand_dist)
{
	double PathLenth = 0;//传输长度
	double t;
	int sphereIndex;//交点球编号
	
	bool inSphere = false;
	double n1, n2;
	
	while (Intersection(Position, Direction, t, sphereIndex, ParticlePosition, ParticleRadius))//光线与粒子有交点的情况下，计算
	{
		flag = 0;//散射
		R1 = Position + Direction * t;
		if (!InOtherSphere(ParticlePosition, ParticleRadius, R1, sphereIndex))//R1在粒子表面
		{
			//printf("R1在粒子表面——");

			

			Vec3<double> Normal = R1 - ParticlePosition[sphereIndex];//边界外法向量
			Normal = Normal / Normal.length();//归一化

			double cos_theta_i = Vec3_Dot(Normal, Direction);
			double sin_theta_i = CosToSin(cos_theta_i);
			
			if (inSphere)//从内部撞击边界
			{
				//printf("从内部撞击边界\n");
				PathLenth += (R1 - Position).length();//出粒子的时候加上在粒子中传播的长度

				n1 = Ref_Index;//粒子介质折射率
				n2 = 1.0;//空气折射率

				if (sin_theta_i >= n2 / n1)
				{//全反射
					Reflection(R1, Direction, Direction1, ParticlePosition, sphereIndex);//调用反射函数
				}
				else {
					double Angle_i = acos_optimized(cos_theta_i);//入射角
					double Angle_t = asin(n1 * sin_theta_i / n2);

					if (rand_dist(rand_engine) > reflectivity(Angle_i, Angle_t))
					{//发生折射
						Refraction(R1, Direction, Direction1, n1, n2, ParticlePosition, sphereIndex);//调用折射函数
						inSphere = false;
						
				
					}
					else
					{//发生反射
						Reflection(R1, Direction, Direction1, ParticlePosition, sphereIndex);//调用反射函数
					}
				}
			}
			else//从外部撞击边界
			{
				//printf("从外部撞击边界, ");
				n1 = 1.0;
				n2 = Ref_Index;

				double Angle_i = acos_optimized(-cos_theta_i);//入射角
				double Angle_t = asin(n1 * sin_theta_i / n2);
				
				//printf("反射率：%f\n", reflectivity(Angle_i, Angle_t));

				if (rand_dist(rand_engine) > reflectivity(Angle_i, Angle_t))
				{//发生折射

					if (flag_refl == -1) flag_refl = 0;
					Refraction(R1, Direction, Direction1, n1, n2, ParticlePosition, sphereIndex);//调用折射函数
					inSphere = true;
				}
				else
				{//发生反射
					if (flag_refl == -1) flag_refl = 1;
					Reflection(R1, Direction, Direction1, ParticlePosition, sphereIndex);//调用反射函数
			
					
					
				
				}
			}

			Direction = Direction1;//发射点方向改变
		}
		else//R1处于粒子内部（两球重合区域）
		{
			//printf("两球重合区域\n");
			PathLenth += (R1 - Position).length();//传输长度大于概率光程，光线被吸收
		}

		if (PathLenth >= absorb_length)
		{
			flag = 1;
			return;
		}

		Position = R1;//发射点置于此处
	}
	//cout<<N_sca<<endl;
	
	return;
	
}



//折射函数，求折射方向（满足Snell定律）
void Refraction(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, double n1, double n2, vector<Vec3<double> > ParticlePosition, int SphereIndex)
{
	Vec3<double> Normal = Position - ParticlePosition[SphereIndex];//折射面法向量
	Normal = Normal / Normal.length();//归一化

	double cos_theta1 = Vec3_Dot(Normal, Direction);
	if (cos_theta1 < 0)
	{
		Normal = Normal*(-1.0);//取法向为与入射、折射方向呈锐角的方向
	}
	double sin_theta1 = sqrt(1.0 - cos_theta1 * cos_theta1);
	double sin_theta2 = n1 * sin_theta1 / n2;
	double cos_theta2 = sqrt(1.0 - sin_theta2 * sin_theta2);

	Vec3<double> N = Vec3_Cross(Direction, Normal);

	double A = Normal.y * N.x - Normal.x * N.y;

	double y2_a = N.x * cos_theta2 / A;
	double y2_b = (Normal.x * N.z - Normal.z * N.x) / A;

	double x2_a = -N.y * cos_theta2 / A;
	double x2_b = -N.z / N.x - (Normal.x * N.z - Normal.z * N.x) * N.y / A / N.x;

	double z_c = x2_a * x2_a + y2_a * y2_a - 1.0;
	double z_b = 2.0 * (x2_a * x2_b + y2_a * y2_b);
	double z_a = x2_b * x2_b + y2_b * y2_b + 1.0;

	double delta = z_b * z_b - 4.0 * z_a * z_c;//判据
	if (delta < 0) printf("计算折射时，折射方向计算报错！\n");
	else if (abs(delta) < 1e-16)
	{
		double z2 = -z_b / 2.0 / z_a;
		double x2 = x2_a + x2_b * z2;
		double y2 = y2_a + y2_b * z2;
		Direction1.ReConstruct(x2, y2, z2);
	}
	else
	{
		double z2 = (-z_b + sqrt(delta)) / 2.0 / z_a;
		double x2 = x2_a + x2_b * z2;
		double y2 = y2_a + y2_b * z2;
		Direction1.ReConstruct(x2, y2, z2);

		Vec3<double> Direction1_x_Normal = Vec3_Cross(Direction1, Normal);
		Vec3<double> Direction_x_Normal = Vec3_Cross(Direction, Normal);

		if (Vec3_Dot(Direction1_x_Normal, Direction_x_Normal) < 0)
		{
			z2 = (-z_b - sqrt(delta)) / 2.0 / z_a;
			x2 = x2_a + x2_b * z2;
			y2 = y2_a + y2_b * z2;
			Direction1.ReConstruct(x2, y2, z2);
		}
	}

}
//反射函数，求反射方向（Fresnel反射定律)
void Reflection(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, vector<Vec3<double> > ParticlePosition, int SphereIndex)
{
	Vec3<double> Normal = Position - ParticlePosition[SphereIndex];//反射面法向量
	Normal = Normal / Normal.length();//归一化
	double cos_theta1 = Vec3_Dot(Normal, Direction);
	if (cos_theta1 > 0)
	{
		Normal = Normal * (-1.0);
	}
	else cos_theta1 = -cos_theta1;

	Direction1 = Direction + Normal * 2.0 * cos_theta1;
}
//反射率的计算
double reflectivity(double Angle_i, double Angle_t)
{
	return 0.5 * (pow(tan(Angle_i - Angle_t) / tan(Angle_i + Angle_t), 2) + pow(sin(Angle_i - Angle_t) / sin(Angle_i + Angle_t), 2));
}

//求天顶角，并分到网格中
int sca_mesh(Vec3<double>Direction1, Vec3<double> Emission0)
{
	Direction1=Direction1/Direction1.length();//归一化
	double cos_theta_I=Vec3_Dot(Emission0,Direction1);//cos(theta)
	double Angle = acos_optimized(cos_theta_I)*180/PI;//散射角
	return round(Angle);//整数值
	 
}
//求圆周角，并分布到网格中
int phi_mesh(Vec3<double>Direction1, Vec3<double>R1)
{
	Vec3<double>X_axis(1.0, 0.0, 0.0);//与x轴的夹角
	Direction1 = Direction1 / Direction1.length();//归一化
	double cos_phi = Vec3_Dot(X_axis, Direction1);//cos
	double phi_angle = acos_optimized(cos_phi);
	//cout<<phi_angle<<endl;
	if (R1.y < 1e-16) { phi_angle = (2 * PI - phi_angle) * 180 / PI; }//圆周角
	else phi_angle = phi_angle * 180 / PI;
	//cout<<phi_angle<<endl;
	return round(phi_angle);//整数值
}


