#include "CommonHeader.h"
#include "Point.h"

//---------------------------------------------------------------------------
// входные данные:
const double	E = 10e-5;						// точность 
const double	PI = 3.14159265358979323846;		
const int		N_ksi = 10;						// число разбиений по оси кси
const int		N_eta = 10;						// число разбиений по оси эта
double			h_ksi = 1. / N_ksi;				// шаг по оси кси
double			h_eta = 1. / N_eta;				// шаг по оси кси
const double	R1 = 1.;						// меньший радиус
const double	R2 = 3.;						// больший радиус
double			hRksi = (R2 - R1) / N_ksi;		// шаг по радиусу вдоль кси
double			hReta = (R2 - R1) / N_eta;		// шаг по радиусу вдоль эта


//const double k0 = 0.2;
//const double kN = 0.4;
//const double x0 = 0.;
//const double xN = 1.;
//double xLow(double y)
//{
//	return (0.025 * sin(5 * PI * y));
//}
//double xUp(double y)
//{
//	return (1 + 0.05 * cos(5 * PI / 2 * y));
//}
//double yLow(double x)
//{
//	return (k0 + cos(PI * x));
//}
//double yUp(double x)
//{
//	return (1. + kN * sin(PI * x));
//}
//double l0 = (yUp(x0) - yLow(x0)) / N_eta;
//double lN = (yUp(xN) - yLow(xN)) / N_eta;
//double h0 = (xN - x0) / N_ksi;
//double hN = h0;




//---------------------------------------------------------------------------
// сетки для итерационного процесса:
Point netK[N_ksi + 1][N_eta + 1];				// сетка на k-м шаге
Point netK1[N_ksi + 1][N_eta + 1];				// сетка на (k+1)-м шаге
//---------------------------------------------------------------------------
// инициализация граничных значений:
void initBorderValues()
{
	// граничные точки сетки:
	for (int i = 0; i <= N_ksi; i++)
	{
		netK[i][0].x = (R1 + hReta * i) * cos(PI / (2 * N_eta) * 0.);
		netK[i][N_eta].x = (R1 + hReta * i) * cos(PI / (2 * N_eta) * N_eta);
		netK[i][0].y = (R1 + hReta * i) * sin(PI / (2 * N_eta) * 0.);
		netK[i][N_eta].y = (R1 + hReta * i) * sin(PI / (2 * N_eta) * N_eta);

		/*netK[i][0].x = i * h0;
		netK[i][N_eta].x = i * hN;
		netK[i][0].y = yLow(netK[i][0].x);
		netK[i][N_eta].y = yUp(netK[i][N_eta].x);*/

		netK1[i][0].x = netK[i][0].x;
		netK1[i][N_eta].x = netK[i][N_eta].x;
		netK1[i][0].y = netK[i][0].y;
		netK1[i][N_eta].y = netK[i][N_eta].y;
	}
	for (int j = 0; j <= N_eta; j++)
	{
		netK[0][j].x = (R1 + hRksi * 0.) * cos(PI / (2 * N_ksi) * j);
		netK[N_ksi][j].x = (R1 + hRksi * N_ksi) * cos(PI / (2 * N_ksi) * j);
		netK[0][j].y = (R1 + hRksi * 0.) * sin(PI / (2 * N_ksi) * j);
		netK[N_ksi][j].y = (R1 + hRksi * N_ksi) * sin(PI / (2 * N_ksi) * j);

		/*netK[0][j].y = yLow(x0) + j * l0;
		netK[N_ksi][j].y = yLow(xN) + j * lN;

		netK[0][j].x = xLow(netK[0][j].y); 
		netK[N_ksi][j].x = xUp(netK[N_ksi][j].y);*/
		

		netK1[0][j].y = netK[0][j].y;
		netK1[N_ksi][j].y = netK[N_ksi][j].y;
		netK1[0][j].x = netK[0][j].x;
		netK1[N_ksi][j].x = netK[N_ksi][j].x;
	}

	// фи на границе:
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

	// пси на границе:
	for (int j = 0; j <= N_eta - 1; j++)
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
// подсчёт внутренних коэффициентов:
void GetInnerValues()
{
	for (int i = 1; i <= N_ksi - 1; i++)
	{
		for (int j = 1; j <= N_eta - 1; j++)
		{
			// альфа, бета, гамма:
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

			// фи:
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
			//---
			if (netK[N_ksi][j].y == netK[0][j].y)
			{
				netK[i][j].fi = 0;
			}
			else
			{
				netK[i][j].fi = ((netK[i][j].y - netK[N_ksi][j].y) * (netK[i][N_eta].fi - netK[i][0].fi)) / (netK[N_ksi][j].y - netK[0][j].y) + netK[i][N_eta].fi;
			}
			//---
			// пси:
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
			//---
			if (netK[i][N_eta].x == netK[i][0].x)
			{
				netK[i][j].psi = 0;
			}
			else
			{
				netK[i][j].psi = ((netK[i][j].x - netK[i][N_eta].x) * (netK[N_ksi][j].psi - netK[0][j].psi)) / (netK[i][N_eta].x - netK[i][0].x) + netK[N_ksi][j].psi;
			}
			//---
		}
	}
}
//---------------------------------------------------------------------------
// вывод узлов сетки на печать:
void printNet()
{
	ofstream out("Net.txt");
	out.precision(6);
	out.setf(ios::fixed);
	for (int i = 0; i <= N_ksi; i++)
	{
		for (int j = 0; j <= N_eta; j++)
		{
			out << "{" << netK[i][j].x << ", " << netK[i][j].y << "}, ";
		}
		//out << "\n";
	}
	out.close();
}
//---------------------------------------------------------------------------
// вывод сетки для графического представления:
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
	for (int i = 0; i <= N_ksi; i++)
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
// условие остановки итерационного процесса:
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
// организация итерационного процесса:
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
				if (netK[i][j].alpha == 0 && netK[i][j].gamma == 0)
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
// входная точка программы:
int main()
{
	initBorderValues();

	DoIter();

	printNet();
	printNetForGraph();

	cout << nOfIter << endl;


	return 0;
}