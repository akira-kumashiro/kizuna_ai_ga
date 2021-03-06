#include "stdafx.h"
#include "GA.h"

#define MAX_GENERATION 30000
#define MAX_GENOM_LIST 50
#define __ENABLE_MUTATION__

int main()
{
	_TCHAR model[] = L"キズナアイ";

	setlocale(LC_ALL, "");
	std::wcout.imbue(std::locale("", std::locale::ctype));

	//配列をstd::vectorへ変換
	std::vector<_TCHAR> mdl(model, std::end(model));

	for (int i = 0; i < mdl.size(); i++)
	{
		printf_s("%d,", mdl[i]);
	}
	std::wcout << model << std::endl;

	std::cout << std::endl;
	GA ga(MAX_GENOM_LIST, mdl.size(), std::vector<int>(mdl.size(), std::numeric_limits<_TCHAR>::max()), std::vector<int>(mdl.size(), std::numeric_limits<_TCHAR>::min()), mdl);//遺伝的アルゴリズム諸関数をまとめたクラスの宣言

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = ga.selection();//選択

		ga.blxAlphaCrossover();//交叉
#ifdef __ENABLE_MUTATION__
		ga.mutation();//突然変異
#endif
		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			ga.calc(true, change);//評価関数の計算
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

