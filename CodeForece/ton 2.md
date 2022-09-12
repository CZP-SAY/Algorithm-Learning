## D



### 思路

​	结论：经过一次操作1,数组的前缀和数组的总和是不会变的,经过一次操作2,数组的前缀和数组值和-1. 通过它我们就能求出谁用了操作2,那个数组算出的前缀和数组总和最小,那个就用了2,差多少,就用了多少次.

​	前缀和数组:sum[i],可以用sum[r]-sum[l-1]算出一个数组[l,r]的区域之和.它的和:

sum[0]+..+sum[i]+sum[i+1]+sum[i+2] ..

​	(写题的时候,要观察两个操作对数组产生的变化的区别,而对区间操作的,可以联想到前缀和/差分)

​	产生这个结论的原因是,



### 代码

```c++
#include <bits/stdc++.h>
using namespace std;
int n, m, t;
long long a[233333], x;
void slove()
{
     cin >> n >> m;
     vector<long long> sum(n);
     for (int id = 0; id < n; id++)
     {
          int pre = 0;
          for (int i = 1; i <= m; i++)
          {
               int x;
               cin >> x;
               pre += x;
               sum[id] += pre;
          }
     }
     int mx = *max_element(sum.begin(), sum.end());
     int mi = *min_element(sum.begin(), sum.end());
     cout << min_element(sum.begin(), sum.end()) - sum.begin() + 1 << " ";
     cout << mx - mi << endl;
}
signed main()
{
     cin >> t;
     while (t--)
     {
          slove();
     }
}
```







## E



### 思路

