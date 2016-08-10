// #include <cstdlib>

// template<typename T>
// T**** make4DTable(int n1, int n2, int n3, int n4)
// {
// //allocate memory for rnaWins
// 	T **** rnaWins;	//rnaWins[1st rna][2nd rna][w1][w2];
// 	T INF_1 = 9999;
// 	rnaWins = (T ****) malloc(sizeof (T ***) * (n1 + 1));
// 	int i,j,k;
// 	for (i = 0; i <= n1; i++)
// 	{
// 		rnaWins[i] = (T ***) malloc(sizeof (T **) * (n2 + 1));
// 		for (j = 0; j <= n2; j++)
// 		{
// 			rnaWins[i][j] = (T **) malloc(sizeof (T *) * (n3+1));
// 			for (k = 0; k <= 25; k++)
// 			{
// 				rnaWins[i][j][k] = (T *) malloc(sizeof (T) * (n4+1));
// 				int k2 = 0;
// 				for (k2 = 0; k2 <= 25; k2++)
// 				{
// 					rnaWins[i][j][k][k2] = INF_1;
// 				}
// 			}
// 		}
// 	}
// 	return rnaWins;
// }


// template<typename T>
// string vecToString(vector<T> v)
// {
// 	string s;
// 	for(auto e : v)
// 		s += e + " ";
// 	return s; 
// }


// // int main()
// // {
// // 	double ****p = make4DTable<double>(10,10,10,10);
// // }