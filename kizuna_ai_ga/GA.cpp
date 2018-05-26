#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<int> _varMax, std::vector<int> _varMin, std::vector<_TCHAR> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num)
{
	//もらった変数をクラス内変数に格納
	model = _model;
	varMax = _varMax;
	varMin = _varMin;

	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//遺伝子の初期設定
		}
	}
	prev_data = data;
	calcResult();

	displayValues(false);
}

bool GA::selection()
{
	int max_num = 0;//最も評価の良い個体の番号
	bool ret = isChanged;//最も評価の良い個体の変化の監視(デバッグ用)
#ifndef __ENABLE_LIGHT_WAIGHT_MODE__
	isChanged = false;
#endif 
	eliteData = searchRank(0);//最も評価の良い個体を保持
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		double selector = random(0.0, 1.0);//乱数を生成
		double needle = 0;//ルーレットの針を生成
		int j = 0;
		for (; ; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//ルーレットの針を乱数の値まで進める
			if (needle > selector)
				break;
			if (j == (data.size() - 1))
				break;
		}
		data[i] = prev_data[j];
	}
	return ret;
}

void GA::blxAlphaCrossover()
{
	prev_data = data;

	for (int i = 0; i < data.size(); i += 2)//2個ずつ交叉
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			double ave = (data[i].x[j] + data[i + 1].x[j]) / 2;
			double length = std::abs((data[i].x[j] - data[i + 1].x[j]));

			data[i].x[j] = (int)random(ave - length * (1 + alpha * 2) / 2, ave + length * (1 + alpha * 2) / 2);
			data[i + 1].x[j] = (int)random(ave - length * (1 + alpha * 2) / 2, ave + length * (1 + alpha * 2) / 2);
		}
	}
}

void GA::mutation()
{
	for (int i = 0; i < data.size(); i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//個体突然変異率の計算
		{
#ifdef __ENABLE_SINGLE_POINT_MUTATION__		
			int pos = random(0, (int)data[i].x.size() - 1);
			data[i].x[pos] += (int)random(varMin[pos] * random(0, genomMutationRate), varMax[pos] * random(0, genomMutationRate));
#else
			for (int j = 0; j < data[i].x.size(); j++)
			{
				data[i].x[j] = random(varMin[j], varMax[j]);
			}
#endif
		}
	}
}

void GA::calc(bool enableDisplay, bool enableOnleLine)
{
	int minNum = 0;
	calcResult();
	Data maxData = searchRank(0);
	for (int i = 0; i < data.size(); i++)//評価関数が最小の奴と最大のやつを検索
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	//評価関数が最もいいやつを保存
	data[minNum] = eliteData;
#ifndef  __ENABLE_LIGHT_WAIGHT_MODE__
	if (searchRank(0).functionValue - eliteData.functionValue > 0)
		isChanged = true;
#endif
	calcResult();

	if (enableDisplay)
		displayValues(enableOnleLine);
}

void GA::calcResult(bool enableSort)
{
	int maxNum = 0;
	for (int i = 0; i < data.size(); i++)
	{
		data[i].functionValue = 0;
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].functionValue += std::pow(data[i].x[j] - model[j], 2.0);
		}
		if (data[maxNum].functionValue < data[i].functionValue)//座標の中で最も関数が大きいやつを検索
			maxNum = i;

	}
	resultSumValue = 0;
	double coefficient = 0.1 / data[0].x.size();//評価関数用の定数

	for (int i = 0; i < data.size(); i++)
	{
		bool flag = true;
		for (int j = 0; j < data[i].x.size(); j++)
		{
#ifndef	__ENABLE_LIGHT_WAIGHT_MODE__
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//座標が場外にいるやつの処理
			{
				flag = false;
				if (data[i].x[j] > varMax[j])
					data[i].x_str[j] = (_TCHAR)varMax[j];
				else
					data[i].x_str[j] = (_TCHAR)varMin[j];
			}
#else
			if (data[i].x[j] > varMax[j])
			{
				data[i].x_str[j] = (_TCHAR)varMax[j];
				flag = false;
			}
			else if (data[i].x[j] < varMin[j])
			{
				data[i].x_str[j] = (_TCHAR)varMin[j];
				flag = false;
			}
#endif
			else
				data[i].x_str[j] = (_TCHAR)data[i].x[j];

		}
		data[i].result = data[i].functionValue == 0 ? std::pow(10 / coefficient, 2.0) : std::pow(1 / (data[i].functionValue*coefficient), 2.0);
		//data[i].result = std::pow((data[i].functionValue - seg),2.0);//与えられた関数の値から切片で設定した値を引いて2乗する→与えられた関数の値が小さいやつが強くなる

		if (!flag)//場外に出たやつの処理
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
#ifndef __ENABLE_LIGHT_WAIGHT_MODE__
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
#endif 
}

int GA::random(int min, int max)
{
	//乱数の設定
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(engine);
}

double GA::random(int min, double max) { return random((double)min, max); }
double GA::random(double min, int max) { return random(min, (double)max); }
double GA::random(double min, double max)
{
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(engine);
}

void GA::displayValues(bool enableOneLine)
{
	setlocale(LC_ALL, "Japanese");
	std::wcout.imbue(std::locale("Japanese", std::locale::ctype));

	std::vector<Data> data_temp;
#ifndef __ENABLE_LIGHT_WAIGHT_MODE__
	if (enableOneLine)
	{
		data_temp.push_back(searchRank(0));
	}
	else
#endif
	{
		data_temp = data;
		std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
	}

	for (int i = 0; i < data_temp.size(); i++)
	{
		for (int j = 0; j < data_temp[i].x.size(); j++)
		{
			wprintf_s(L"%4x,", data_temp[i].x[j]);
		}
		wprintf_s(L"%ls\t f(x,y)=%4.0lf \tResult=%4.0lf\n", &data_temp[i].x_str[0], data_temp[i].functionValue, data_temp[i].result);
		//wprintf_s(L"\n");
	}
}

GA::Data GA::searchRank(int num)
{
	std::vector<Data> data_temp = data;
	std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue < y.functionValue; });
	return data_temp[num];
}

GA::~GA()
{

}
