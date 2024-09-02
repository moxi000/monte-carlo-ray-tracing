#pragma once

using namespace std;
bool SphereRay(Vec3<double> Position, Vec3<double> Direction, Vec3<double> O, double R, double& t);
void Refraction(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, double n1, double n2, vector<Vec3<double> > ParticlePosition, int SphereIndex);
void Reflection(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, vector<Vec3<double> > ParticlePosition, int SphereIndex);
double reflectivity(double Angle_i, double Angle_t);
int sca_mesh(Vec3<double>Direction1, Vec3<double>Emission0);
int phi_mesh(Vec3<double>Direction1, Vec3<double>R1);



//cosת����sin����
inline double CosToSin(double cos)
{
	return 	sqrt(1 - cos * cos);
}
//�Ż���cos����
inline double acos_optimized(double value)
{
	if (value < -1.0)value = -1.0;
	else if (value > 1.0)value = 1.0;
	value = acos(value);
	return value;
}



bool InOtherSphere(vector<Vec3<double> > ParticlePosition, vector<double> ParticleRadius, Vec3<double> Position, int sphereIndex)
{
	//����ǰ������֮�⣬�Ƿ���������֮��
	for (int i = 0; i < ParticlePosition.size(); i++)
	{
		if (i != sphereIndex)
		{
			if ((Position - ParticlePosition[i]).length() <= ParticleRadius[i] + 1e-15) return true;
		}
	}
	return false;
}

////�жϹ����Ƿ��������ཻ����ȷ���������Ǹ���Ԫ��

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
				sphereIndex = i;//��������
			}
		}
	}
	if (t > 99999) return false;
	else return true;
}
//���������Ľ��㣬�õ�����t
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

		if (t_temp < 0 && t < 0) return false;//������Ч  t<0 & t_temp<0
		else if (t < 1e-15) t = t_temp;//t<0 & t_temp>0
		else
		{
			if (t_temp > 1e-15)//t>0 & t_temp>0
			{
				if (t_temp < t) t = t_temp;//���߶�����0��ȡ��С��һ��
			}
		}
	}
	//InterPosition = Position + Direction * t;
	return true;
}
//�������ؿ��庯��������׷��
void MonteCarlo(Vec3<double> Position, Vec3<double>& R1, Vec3<double> Direction, Vec3<double>& Direction1, double absorb_length, int& flag, int& flag_refl, vector<Vec3<double> > ParticlePosition, vector<double> ParticleRadius, default_random_engine rand_engine, uniform_real_distribution<double> rand_dist)
{
	double PathLenth = 0;//���䳤��
	double t;
	int sphereIndex;//��������
	
	bool inSphere = false;
	double n1, n2;
	
	while (Intersection(Position, Direction, t, sphereIndex, ParticlePosition, ParticleRadius))//�����������н��������£�����
	{
		flag = 0;//ɢ��
		R1 = Position + Direction * t;
		if (!InOtherSphere(ParticlePosition, ParticleRadius, R1, sphereIndex))//R1�����ӱ���
		{
			//printf("R1�����ӱ��桪��");

			

			Vec3<double> Normal = R1 - ParticlePosition[sphereIndex];//�߽��ⷨ����
			Normal = Normal / Normal.length();//��һ��

			double cos_theta_i = Vec3_Dot(Normal, Direction);
			double sin_theta_i = CosToSin(cos_theta_i);
			
			if (inSphere)//���ڲ�ײ���߽�
			{
				//printf("���ڲ�ײ���߽�\n");
				PathLenth += (R1 - Position).length();//�����ӵ�ʱ������������д����ĳ���

				n1 = Ref_Index;//���ӽ���������
				n2 = 1.0;//����������

				if (sin_theta_i >= n2 / n1)
				{//ȫ����
					Reflection(R1, Direction, Direction1, ParticlePosition, sphereIndex);//���÷��亯��
				}
				else {
					double Angle_i = acos_optimized(cos_theta_i);//�����
					double Angle_t = asin(n1 * sin_theta_i / n2);

					if (rand_dist(rand_engine) > reflectivity(Angle_i, Angle_t))
					{//��������
						Refraction(R1, Direction, Direction1, n1, n2, ParticlePosition, sphereIndex);//�������亯��
						inSphere = false;
						
				
					}
					else
					{//��������
						Reflection(R1, Direction, Direction1, ParticlePosition, sphereIndex);//���÷��亯��
					}
				}
			}
			else//���ⲿײ���߽�
			{
				//printf("���ⲿײ���߽�, ");
				n1 = 1.0;
				n2 = Ref_Index;

				double Angle_i = acos_optimized(-cos_theta_i);//�����
				double Angle_t = asin(n1 * sin_theta_i / n2);
				
				//printf("�����ʣ�%f\n", reflectivity(Angle_i, Angle_t));

				if (rand_dist(rand_engine) > reflectivity(Angle_i, Angle_t))
				{//��������

					if (flag_refl == -1) flag_refl = 0;
					Refraction(R1, Direction, Direction1, n1, n2, ParticlePosition, sphereIndex);//�������亯��
					inSphere = true;
				}
				else
				{//��������
					if (flag_refl == -1) flag_refl = 1;
					Reflection(R1, Direction, Direction1, ParticlePosition, sphereIndex);//���÷��亯��
			
					
					
				
				}
			}

			Direction = Direction1;//����㷽��ı�
		}
		else//R1���������ڲ��������غ�����
		{
			//printf("�����غ�����\n");
			PathLenth += (R1 - Position).length();//���䳤�ȴ��ڸ��ʹ�̣����߱�����
		}

		if (PathLenth >= absorb_length)
		{
			flag = 1;
			return;
		}

		Position = R1;//��������ڴ˴�
	}
	//cout<<N_sca<<endl;
	
	return;
	
}



//���亯���������䷽������Snell���ɣ�
void Refraction(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, double n1, double n2, vector<Vec3<double> > ParticlePosition, int SphereIndex)
{
	Vec3<double> Normal = Position - ParticlePosition[SphereIndex];//�����淨����
	Normal = Normal / Normal.length();//��һ��

	double cos_theta1 = Vec3_Dot(Normal, Direction);
	if (cos_theta1 < 0)
	{
		Normal = Normal*(-1.0);//ȡ����Ϊ�����䡢���䷽�����ǵķ���
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

	double delta = z_b * z_b - 4.0 * z_a * z_c;//�о�
	if (delta < 0) printf("��������ʱ�����䷽����㱨��\n");
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
//���亯�������䷽��Fresnel���䶨��)
void Reflection(Vec3<double> Position, Vec3<double> Direction, Vec3<double>& Direction1, vector<Vec3<double> > ParticlePosition, int SphereIndex)
{
	Vec3<double> Normal = Position - ParticlePosition[SphereIndex];//�����淨����
	Normal = Normal / Normal.length();//��һ��
	double cos_theta1 = Vec3_Dot(Normal, Direction);
	if (cos_theta1 > 0)
	{
		Normal = Normal * (-1.0);
	}
	else cos_theta1 = -cos_theta1;

	Direction1 = Direction + Normal * 2.0 * cos_theta1;
}
//�����ʵļ���
double reflectivity(double Angle_i, double Angle_t)
{
	return 0.5 * (pow(tan(Angle_i - Angle_t) / tan(Angle_i + Angle_t), 2) + pow(sin(Angle_i - Angle_t) / sin(Angle_i + Angle_t), 2));
}

//���춥�ǣ����ֵ�������
int sca_mesh(Vec3<double>Direction1, Vec3<double> Emission0)
{
	Direction1=Direction1/Direction1.length();//��һ��
	double cos_theta_I=Vec3_Dot(Emission0,Direction1);//cos(theta)
	double Angle = acos_optimized(cos_theta_I)*180/PI;//ɢ���
	return round(Angle);//����ֵ
	 
}
//��Բ�ܽǣ����ֲ���������
int phi_mesh(Vec3<double>Direction1, Vec3<double>R1)
{
	Vec3<double>X_axis(1.0, 0.0, 0.0);//��x��ļн�
	Direction1 = Direction1 / Direction1.length();//��һ��
	double cos_phi = Vec3_Dot(X_axis, Direction1);//cos
	double phi_angle = acos_optimized(cos_phi);
	//cout<<phi_angle<<endl;
	if (R1.y < 1e-16) { phi_angle = (2 * PI - phi_angle) * 180 / PI; }//Բ�ܽ�
	else phi_angle = phi_angle * 180 / PI;
	//cout<<phi_angle<<endl;
	return round(phi_angle);//����ֵ
}


