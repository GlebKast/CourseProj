#include "CommonHeader.h"
#include "Point.h"
//---------------------------------------------------------------------------
// ������� ������:
const double	E = 10e-5;						// �������� 
const double	PI = 3.14159265358979323846;	// ����� ��	
const int		N_ksi = 30;						// ����� ��������� �� ��� ���
const int		N_eta = 40;						// ����� ��������� �� ��� ���
double			h_ksi = 1. / N_ksi;				// ��� �� ��� ���
double			h_eta = 1. / N_eta;				// ��� �� ��� ���
const double	R1_1 = 1.;						// ������� ������ �� ���
const double	R1_2 = 3.;						// ������� ������ �� ���
const double	R2_1 = 3.;						// ������� ������ �� ���
const double	R2_2 = 6.;						// ������� ������ �� ���
double			hRksi_1 = (R1_2 - R1_1) / N_ksi;		// ��� �� ������� ����� ���
double			hRksi_2 = (R2_2 - R1_2) / N_eta;		// ��� �� ������� ����� ���

//---------------------------------------------------------------------------
// ����� ��� ������������� ��������:
Point netK[N_ksi + 1][N_eta + 1];				// ����� �� k-� ����
Point netK1[N_ksi + 1][N_eta + 1];				// ����� �� (k+1)-� ����
//---------------------------------------------------------------------------
// ������������� ��������� ��������:
void initBorderValues()
{
	// ��������� ����� �����:
	for (int i = 0; i <= N_ksi; i++)
	{
		netK[i][0].x = (R1_1 + hRksi_1 * i) * cos(PI / (2 * N_eta) * 0.);
		netK[i][N_eta].x = (R1_1 + hRksi_1 * i) * cos(PI / (2 * N_eta) * N_eta);
		netK[i][0].y = (R2_1 + hRksi_2 * i) * sin(PI / (2 * N_eta) * 0.);
		netK[i][N_eta].y = (R2_1 + hRksi_2 * i) * sin(PI / (2 * N_eta) * N_eta);

		netK1[i][0].x = netK[i][0].x;
		netK1[i][N_eta].x = netK[i][N_eta].x;
		netK1[i][0].y = netK[i][0].y;
		netK1[i][N_eta].y = netK[i][N_eta].y;
	}
	for (int j = 0; j <= N_eta; j++)
	{
		netK[0][j].x = (R1_1 + hRksi_1 * 0.) * cos(PI / (2 * N_eta) * j);
		netK[N_ksi][j].x = (R1_1 + hRksi_1 * N_ksi) * cos(PI / (2 * N_eta) * j);
		netK[0][j].y = (R2_1 + hRksi_2 * 0.) * sin(PI / (2 * N_eta) * j);
		netK[N_ksi][j].y = (R2_1 + hRksi_2 * N_ksi) * sin(PI / (2 * N_eta) * j);
		
		netK1[0][j].y = netK[0][j].y;
		netK1[N_ksi][j].y = netK[N_ksi][j].y;
		netK1[0][j].x = netK[0][j].x;
		netK1[N_ksi][j].x = netK[N_ksi][j].x;
	}

	// �� �� �������:
	for (int i = 1; i <= N_ksi - 1; i++)
	{
		if (netK[i + 1][0].x == netK[i - 1][0].x)
		{
			netK[i][0].fi = 0.;
		}
		else
		{
			netK[i][0].fi = 2 / h_ksi * (2 * netK[i][0].x 
				- netK[i + 1][0].x - netK[i - 1][0].x)
				/ (netK[i + 1][0].x - netK[i - 1][0].x);
		}
		if (netK[i + 1][N_eta].x == netK[i - 1][N_eta].x)
		{
			netK[i][N_eta].fi = 0.;
		}
		else
		{
			netK[i][N_eta].fi = 2 / h_ksi * (2 * netK[i][N_eta].x 
				- netK[i + 1][N_eta].x - netK[i - 1][N_eta].x) 
				/ (netK[i + 1][N_eta].x - netK[i - 1][N_eta].x);
		}		
		netK1[i][0].fi = netK[i][0].fi;
		netK1[i][N_eta].fi = netK[i][N_eta].fi;
	}

	// ��� �� �������:
	for (int j = 1; j <= N_eta - 1; j++)
	{
		if (netK[0][j + 1].y == netK[0][j - 1].y)
		{
			netK[0][j].psi = 0.;
		}
		else
		{
			netK[0][j].psi = 2 / h_eta * (2 * netK[0][j].y 
				- netK[0][j + 1].y - netK[0][j - 1].y) 
				/ (netK[0][j + 1].y - netK[0][j - 1].y);
		}
		if (netK[N_ksi][j + 1].y == netK[N_ksi][j - 1].y)
		{
			netK[N_ksi][j].psi = 0.;
		}
		else
		{
			netK[N_ksi][j].psi = 2 / h_eta * (2 * netK[N_ksi][j].y - 
				netK[N_ksi][j + 1].y - netK[N_ksi][j - 1].y) 
				/ (netK[N_ksi][j + 1].y - netK[N_ksi][j - 1].y);
		}
		netK1[0][j].psi = netK[0][j].psi;
		netK1[N_ksi][j].psi = netK[N_ksi][j].psi;
	}
}
//---------------------------------------------------------------------------
// ������� ���������� �������������:
void GetInnerValues()
{
	for (int i = 1; i <= N_ksi - 1; i++)
	{
		for (int j = 1; j <= N_eta - 1; j++)
		{
			// �����, ����, �����:
			netK[i][j].alpha = (pow(netK[i][j + 1].x - netK[i][j - 1].x, 2) 
				+ pow(netK[i][j + 1].y - netK[i][j - 1].y, 2))
				/ (4 * h_eta * h_eta);
			netK[i][j].beta = ((netK[i][j + 1].x - netK[i][j - 1].x) 
				* (netK[i + 1][j].x - netK[i - 1][j].x) 
				+ (netK[i][j + 1].y - netK[i][j - 1].y) 
				* (netK[i + 1][j].y - netK[i - 1][j].y)) 
				/ (4 * h_ksi * h_eta);
			netK[i][j].gamma = (pow(netK[i + 1][j].x - netK[i - 1][j].x, 2) 
				+ pow(netK[i + 1][j].y - netK[i - 1][j].y, 2)) 
				/ (4 * h_ksi * h_ksi);

			// ��:
			if (netK[N_ksi][j].x == netK[0][j].x)
			{
				netK[i][j].fi = 0.;
			}
			else
			{
				netK[i][j].fi = ((netK[i][j].x - netK[N_ksi][j].x) 
					* (netK[i][N_eta].fi - netK[i][0].fi))
					/ (netK[N_ksi][j].x - netK[0][j].x) + netK[i][N_eta].fi;
			}
			
			// ���:
			if (netK[i][N_eta].y == netK[i][0].y)
			{
				netK[i][j].psi = 0.;
			}
			else
			{
				netK[i][j].psi = ((netK[i][j].y - netK[i][N_eta].y) 
					* (netK[N_ksi][j].psi - netK[0][j].psi)) 
					/ (netK[i][N_eta].y - netK[i][0].y) + netK[N_ksi][j].psi;
			}
		}
	}
}
//---------------------------------------------------------------------------
// ����� ����� ����� �� ������:
void printNet()
{
	ofstream out("Net.txt");
	out.precision(6);
	out.setf(ios::fixed);
	for (int i = 0; i <= N_ksi; i++)
	{
		for (int j = 0; j <= N_eta; j++)
		{
			out << "(" << netK[i][j].x << ", " << netK[i][j].y << ")  ";
		}
		out << "\n";
	}
	out.close();
}
//---------------------------------------------------------------------------
// ����� ����� ��� ������������ �������������:
void printNetForGraph()
{
	ofstream out("NetGraph.txt");
	out.precision(6);
	out.setf(ios::fixed);
	out << "{";
	for (int j = 0; j <= N_eta; j++)
	{
		out << "{";
		for(int i = 0; i < N_ksi; i++)
			out << "{" << netK[i][j].x << ", " << netK[i][j].y << "}, ";
		out << "{" << netK[N_ksi][j].x << ", " << netK[N_ksi][j].y << "} ";
		out << "}, ";
	}
	for (int i = 0; i < N_ksi; i++)
	{
		out << "{";
		for (int j = 0; j < N_eta; j++)
			out << "{" << netK[i][j].x << ", " << netK[i][j].y << "}, ";
		out << "{" << netK[i][N_eta].x << ", " << netK[i][N_eta].y << "} ";
		out << "}, ";
	}
	out << "{";
	for(int j = 0; j < N_eta; j++)
		out << "{" << netK[N_ksi][j].x << ", " << netK[N_ksi][j].y << "}, ";
	out << "{" << netK[N_ksi][N_eta].x << ", " << netK[N_ksi][N_eta].y << "}";
	out << "} ";
	out << "} ";
	out.close();
}
//---------------------------------------------------------------------------
// ������� ��������� ������������� ��������:
bool StopCondition()
{
	for (int i = 0; i <= N_ksi; i++)
	{
		for (int j = 0; j <= N_eta; j++)
		{
			if (fabs(netK1[i][j].x - netK[i][j].x) >= E 
				|| fabs(netK1[i][j].y - netK[i][j].y) >= E)
			{
				return true;
			}
		}
	}
	return false;
}
//---------------------------------------------------------------------------
// ����������� ������������� ��������:
int nOfIter = 0;
void DoIter()
{
	double K = 0.;
	do
	{
		for (int i = 1; i <= N_ksi - 1; i++)
		{
			for (int j = 1; j <= N_eta - 1; j++)
			{
				netK[i][j].x = netK1[i][j].x;
				netK[i][j].y = netK1[i][j].y;
			}
		}

		GetInnerValues();

		for (int i = 1; i <= N_ksi - 1; i++)
		{
			for (int j = 1; j <= N_eta - 1; j++)
			{
				if (netK[i][j].alpha * h_eta * h_eta + netK[i][j].gamma * h_ksi * h_ksi == 0)
				{
					K = 0.;
				}
				else
				{
					K = 1 / (4 * (netK[i][j].alpha * h_eta * h_eta + netK[i][j].gamma * h_ksi * h_ksi));
				}
				netK1[i][j].x = K * (netK[i][j].alpha * h_eta * h_eta * ((2 + netK[i][j].fi * h_ksi) * netK[i + 1][j].x + (2 - netK[i][j].fi * h_ksi) * netK1[i - 1][j].x)
					+ netK[i][j].gamma * h_ksi * h_ksi * ((2 + netK[i][j].psi * h_eta) * netK[i][j + 1].x + (2 - netK[i][j].psi * h_eta) * netK1[i][j - 1].x) 
					- netK[i][j].beta * h_ksi * h_eta * (netK[i + 1][j + 1].x - netK1[i + 1][j - 1].x - netK[i - 1][j + 1].x + netK1[i - 1][j - 1].x));
				netK1[i][j].y = K * (netK[i][j].alpha * h_eta * h_eta * ((2 + netK[i][j].fi * h_ksi) * netK[i + 1][j].y + (2 - netK[i][j].fi * h_ksi) * netK1[i - 1][j].y)
					+ netK[i][j].gamma * h_ksi * h_ksi * ((2 + netK[i][j].psi * h_eta) * netK[i][j + 1].y + (2 - netK[i][j].psi * h_eta) * netK1[i][j - 1].y)
					- netK[i][j].beta * h_ksi * h_eta * (netK[i + 1][j + 1].y - netK1[i + 1][j - 1].y - netK[i - 1][j + 1].y + netK1[i - 1][j - 1].y));
			}
		}
		nOfIter++;

	} while (StopCondition());
}
//---------------------------------------------------------------------------
// ������� ����� ���������:
int main()
{
	initBorderValues();

	DoIter();

	printNet();
	printNetForGraph();

	cout << nOfIter << endl;

	return 0;
}