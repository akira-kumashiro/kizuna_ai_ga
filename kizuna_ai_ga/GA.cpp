#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<int> _varMax, std::vector<int> _varMin, std::vector<char> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num)/*,
	prevElite(_var_num)*/
{
	//もらった変数をクラス内変数に格納
	//max_genom_list = _max_genom_list;
	//var_num = _var_num;
	model = _model;
	varMax = _varMax;
	varMin = _varMin;
	//prevElite = eliteData;

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
	isChanged = false;

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
	/*	int max_num = 0;// , prev_max_num = 0;//最も評価の良い個体の番号
		//bool ret = false;
		bool ret = isChanged;
		isChanged = false;

		//calc(false);

		for (int i = 0; i < data.size(); i++)//ルーレット選択用に評価関数の合計と一番評価の良い番号を取得
		{
			if (data[i].result > data[max_num].result)
				max_num = i;
			//		if (prev_data[i].result > prev_data[prev_max_num].result)
			//			prev_max_num = i;
		}
		//prevElite = eliteData;
		eliteData = data[max_num];//最も評価の良い個体を保持
		//if (eliteData.functionValue - prevElite.functionValue != 0)//最も評価の良い個体の変化の監視(デバッグ用)
		//	ret = true;
		//if (data[maxNum].functionValue - prev_data[prev_max_num].functionValue > 10 * var_num)
		//	ret = true;
		//std::copy(data.begin(), data.end(), prev_data.begin());
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
			//data[i].prev_pos = j;
		}
		//displayValues();
		return ret;*/
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
		/*if (data[i].result > data[maxNum].result)
			maxNum = i;*/
	}
	//評価関数が最もいいやつを保存
	data[minNum] = eliteData;

	if (searchRank(0).functionValue - eliteData.functionValue > 0)
		isChanged = true;

	calcResult();

	if (enableDisplay)
		displayValues(enableOnleLine);
}

void GA::calcResult(bool enableSort)
{
	int maxNum = 0;
	//double seg;
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
	//double seg = data[maxNum].functionValue;//評価関数の切片を与えられた関数が最も大きいやつにセット
	resultSumValue = 0;
	double coefficient = 0.1 / data[0].x.size();//評価関数用の定数

	for (int i = 0; i < data.size(); i++)
	{
		bool flag = true;
		//double coefficient = 0.1 / data[i].x.size();//評価関数用の定数
		for (int j = 0; j < data[i].x.size(); j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//座標が場外にいるやつの処理
			{
				flag = false;
				if (data[i].x[j] > varMax[j])
					data[i].x_str[j] = (char)varMax[j];
				else
					data[i].x_str[j] = (char)varMin[j];
			}
			else
				data[i].x_str[j] = (char)data[i].x[j];

		}
		data[i].result = data[i].functionValue == 0 ? 10 / coefficient : 1 / (data[i].functionValue*coefficient);
		//data[i].result = std::pow((data[i].functionValue - seg),2.0);//与えられた関数の値から切片で設定した値を引いて2乗する→与えられた関数の値が小さいやつが強くなる

		if (!flag)//場外に出たやつの処理
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
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
	std::vector<Data> data_temp;
	if (enableOneLine)
	{
		data_temp.push_back(searchRank(0));
	}
	else
	{
		data_temp = data;
		std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
	}

	for (int i = 0; i < data_temp.size(); i++)
	{
		for (int j = 0; j < data_temp[i].x.size(); j++)
		{
			printf_s("%4d,", data_temp[i].x[j]);
		}
		std::cout << "\t" << std::string(data_temp[i].x_str.begin(), data_temp[i].x_str.end()) << "\t";
		printf_s(" \t f(x,y)=%6.2lf\t Result=%6.2lf\n", data_temp[i].functionValue, data_temp[i].result);
	}
	/*if (enableOneLine)
	{
		for (int j = 0; j < data[0].x.size(); j++)
		{
			//printf_s("%6.2lf,", data[i].x[j]);//デバッグ用
			printf_s("%4d,", eliteData.x[j]);
		}
		std::cout << "\t" << std::string(eliteData.x_str.begin(), eliteData.x_str.end()) << "\t";
		printf_s(" \t f(x,y)=%6.2lf\t Result=%6.2lf\n", eliteData.functionValue, eliteData.result);
	}
	else
	{
		std::vector<Data> data_temp = data;
		std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });

		for (int i = 0; i < data.size(); i++)
		{
			for (int j = 0; j < data[i].x.size(); j++)
			{
				//printf_s("%6.2lf,", data[i].x[j]);//デバッグ用
				printf_s("%4d,", data_temp[i].x[j]);
			}
			std::cout << "\t" << std::string(data_temp[i].x_str.begin(), data_temp[i].x_str.end()) << "\t";
			printf_s(" \t f(x,y)=%6.2lf\t Result=%6.2lf\n", data_temp[i].functionValue, data_temp[i].result);
			//printf_s(" \t f(x,y)=%d\t Result=%d\n", data[i].functionValue, data[i].result);
		}
	}*/
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
