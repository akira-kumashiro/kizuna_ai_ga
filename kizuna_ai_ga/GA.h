#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <locale.h>

#define __ENABLE_SINGLE_POINT_MUTATION__
#define __ENABLE_LIGHT_WAIGHT_MODE__

class GA
{
private:
	double individualMutationRate = 0.3;//個体突然変異率
	double genomMutationRate = 0.1;
	double alpha = 1.0;
	std::vector<int> varMax, varMin;
	std::vector<_TCHAR> model;
	bool isChanged = false;
public:
	double resultSumValue;//評価関数の合計

	class Data//データ格納用クラス
	{
	public:
		std::vector<int> x;//座標
		//std::vector<_TCHAR> x_str;
		std::wstring x_str;
		double functionValue;//与えられた関数の値
		double result;

		Data(int _var_num)//コンストラクタ
		{
			x.resize(_var_num);//isIncludedの配列の長さの設定
			x_str.resize(_var_num);
		}
	};

	std::vector<Data> data, prev_data;//操作前後で値を保持するために2個
	Data eliteData;
	GA(int _max_genom_list, int _var_num, std::vector<int> _varMax, std::vector<int> _varMin, std::vector<_TCHAR> _model);	//コンストラクタ
	bool selection();//選択
	void blxAlphaCrossover();
	void mutation();//突然変異
	void calc(bool enableDisplay, bool enableOneLine = false);//評価関数の計算
private:
	void calcResult(bool enableSort = false);
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	void displayValues(bool enableOneLine);
	Data searchRank(int num);
public:
	~GA();//デコンストラクタ
};


