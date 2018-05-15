#include "stdafx.h"
#include "GA.h"

#define MAX_GENERATION 30000
#define MAX_GENOM_LIST 50
#define VAR_NUM 1

int main()
{
	char model[] = "キズナアイ";

	//配列をstd::vectorへ変換
	std::vector<char> mdl(model, std::end(model));
	std::vector<int> vMax(mdl.size(), 127);
	std::vector<int> vMin(mdl.size(), -127);

	for (int i = 0; i < mdl.size(); i++)
	{
		printf_s("%d,", mdl[i]);
	}
	std::cout << std::endl;
	GA ga(MAX_GENOM_LIST, mdl.size(), vMax, vMin, mdl);//遺伝的アルゴリズム諸関数をまとめたクラスの宣言

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = ga.selection();//選択

		ga.blxAlphaCrossover();//交叉

		ga.mutation();//突然変異

		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			ga.calc(true, i != MAX_GENERATION);//評価関数の計算
		}
		else
		{
			ga.calc(false);//評価関数の計算
		}
	}

	while (1)
	{
		if (_kbhit() && _getch() == 27)
			break;
	}
	return 0;
}

