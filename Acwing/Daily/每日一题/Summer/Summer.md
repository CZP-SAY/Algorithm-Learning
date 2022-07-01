### 活动链接:

https://www.acwing.com/activity/content/punch_the_clock/1934/







### 4268:性感素数



#### 题目

“[性感素数](http://mathworld.wolfram.com/SexyPrimes.html) ”是指形如 (p,p+6)这样的一对素数。

之所以叫这个名字，是因为拉丁语管“六”叫“sex”（即英语的“性感”）。

现给定一个整数，请你判断其是否为一个性感素数。

输入格式

输入在一行中给出一个正整数 N。

输出格式

若 N是一个性感素数，则在一行中输出 `Yes`，并在第二行输出与 N配对的另一个性感素数（若这样的数不唯一，输出较小的那个）。

若 N不是性感素数，则在一行中输出 `No`，然后在第二行输出大于 N 的最小性感素数。

数据范围

1≤N≤10^8^

输入样例1：

```
47
```

输出样例1：

```
Yes
41
```

输入样例2：

```
21
```

输出样例2：

```
No
23
```



#### 思路

判断素数

但注意细节：hack的数据：99993431。答案是99993433

```cpp
  for (int i = 2; i <= sqrt(n); i++)
     {
          if (n % i == 0)
               return 0;
     }
```

不优化一下的话,会tle



#### 代码

```c++

#include <bits/stdc++.h>
using namespace std;

bool check(int n)
{
     if (n < 2)
          return 0;
     for (int i = 2; i <= sqrt(n); i++)
     {
          if (n % i == 0)
               return 0;
     }
     return 1;
}

int main()
{

     int n;
     cin >> n;

     if (check(n) && (check(n + 6) || check(n - 6)))
     {
          cout << "Yes" << endl;
          if (check(n - 6))
               cout << n - 6 << endl;
          else
               cout << n + 6 << endl;
     }
     else
     {
          cout << "No" << endl;
          for (int i = n + 1; i <= 1e8; i++)
          {
               if (check(i) && check(i - 6) && i - 6 < n)
               //错在这个小细节了,要再判断一下-6的情况,要是减去6且减去后它是<n的,那么这个i才是n后的最小性感素数
               {
                    cout << i << endl;
                    break;
               }
               if (check(i) && (check(i + 6)))
               {
                    cout << i << endl;
                    break;
               }
          }
     }
}
```









### 4269:校庆



#### 题目



2019 年浙江大学将要庆祝成立 122 周年。

为了准备校庆，校友会收集了所有校友的身份证号。

现在需要请你编写程序，根据来参加校庆的所有人士的身份证号，统计来了多少校友。



输入格式

输入在第一行给出正整数 N。

随后 N 行，每行给出一位校友的身份证号（1818 位由数字和大写字母 XX 组成的字符串）。题目保证身份证号不重复。

随后给出前来参加校庆的所有人士的信息：

首先是一个正整数 M。

随后 M行，每行给出一位人士的身份证号。题目保证身份证号不重复。

输出格式

首先在第一行输出参加校庆的校友的人数。

然后在第二行输出最年长的校友的身份证号 —— 注意身份证第 7−14 位给出的是 `yyyymmdd` 格式的生日。

如果没有校友来，则在第二行输出最年长的来宾的身份证号。题目保证这样的校友或来宾必是唯一的。

数据范围

1≤N,M≤10^5^



输入样例：

```
5
372928196906118710
610481197806202213
440684198612150417
13072819571002001X
150702193604190912
6
530125197901260019
150702193604190912
220221196701020034
610481197806202213
440684198612150417
370205198709275042
```



输出样例：

```
3
150702193604190912
```





#### 思路

一道简单的字符串题目





#### 代码

```c++
#include <bits/stdc++.h>
using namespace std;
const int N = 1e5 + 10;
map<string, int> m;

int main()
{
     int n;
     cin >> n;
     for (int i = 0; i < n; i++)
     {
          string s;
          cin >> s;
          m[s] = 1;
     }

     int n1;
     cin >> n1;
     int cnt = 0;
     string s1 = "99999999";
     string s[N], ans;
     for (int i = 0; i < n1; i++)
     {
          cin >> s[i];
          if (m[s[i]] == 1)
          {
               if (s[i].substr(6, 8) < s1)
               {
                    s1 = s[i].substr(6, 8);
                    ans = s[i];
               }
               cnt++;
          }
     }
     if (cnt > 0)
     {
          cout << cnt << endl;
          cout << ans;
     }
     else
     {

          cout << cnt << endl;
          for (int i = 0; i < n1; i++)
          {

               if (s[i].substr(6, 8) < s1)
               {
                    s1 = s[i].substr(6, 8);
                    ans = s[i];
               }
          }
          cout << ans << endl;
     }
}

```











### 4273:链表合并



#### 题目



给定两个单链表` L1=a1→a2→…→an−1→an` 和 `L2=b1→b2→…→bm−1→bm`，满足：n≥2m。

你的任务是将较短的那个链表逆序，然后将之并入较长的链表，得到形如 `a1→a2→bm→a3→a4→bm−1` 的结果。

例如给定两个链表分别为 6→7和 1→2→3→4→5，你应该输出 1→2→7→3→4→6→5。

**补充**
本题中可能包含不在两个单链表中的节点，这些节点无需考虑。



输入格式

输入首先在第一行中给出两个链表 L1 和 L2 的头结点的地址，以及正整数 N，即给定的结点总数。

一个结点的地址是一个 5 位数的非负整数（可能包含前导 0），空地址 `NULL` 用 −1 表示。

随后 N 行，每行按以下格式给出一个结点的信息：

```
Address Data Next
```

其中 `Address` 是结点的地址，`Data` 是不超过 10^5^ 的正整数，`Next` 是下一个结点的地址。

题目保证没有空链表，并且较长的链表至少是较短链表的两倍长。

输出格式

按顺序输出结果链表，每个结点占一行，格式与输入相同。

数据范围

1≤N≤10^5^



输入样例：

```
00100 01000 7
02233 2 34891
00100 6 00001
34891 3 10086
01000 1 02233
00033 5 -1
10086 4 00033
00001 7 -1
```

输出样例：

```
01000 1 02233
02233 2 00001
00001 7 34891
34891 3 10086
10086 4 00100
00100 6 00033
00033 5 -1
```





#### 思路

orz

本题建立了三个vector来存链表,两个用来存取题目数据,第三个就是答案,把合成后的链表放在这里面.



```c++
typedef pair<int, int> PII;
vector<PII> a,b,c;
```



本题用vector方便读取输出数据(没想到这点,hhh,用vector会明了很多)







#### 代码

```c++
//类似手写链表
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#define x first
#define y second

using namespace std;

typedef pair<int, int> PII;

const int N = 100010;

int h1, h2, n;
int v[N], ne[N];

int main()
{
    scanf("%d%d%d", &h1, &h2, &n);
    while (n -- )
    {
        int addr, val, next;
        scanf("%d%d%d", &addr, &val, &next);
        v[addr] = val, ne[addr] = next;
    }

    vector<PII> a, b;
    for (int i = h1; i != -1; i = ne[i]) a.push_back({i, v[i]});
    for (int i = h2; i != -1; i = ne[i]) b.push_back({i, v[i]});

    if (a.size() < b.size()) swap(a, b);

    vector<PII> c;
    for (int i = 0, j = b.size() - 1; i < a.size(); i += 2, j -- )
    {
        c.push_back(a[i]);
        if (i + 1 < a.size()) c.push_back(a[i + 1]);
        if (j >= 0) c.push_back(b[j]);
    }

    for (int i = 0; i < c.size(); i ++ )
    {
        printf("%05d %d ", c[i].x, c[i].y);
        if (i + 1 < c.size()) printf("%05d\n", c[i + 1].x);
        else puts("-1");
    }

    return 0;
}

```







### 后缀表达式





#### 题目

给定一个二叉表达式树，请你输出相应的后缀表达式，要求使用括号反映运算符的优先级。



输入格式

第一行包含整数 N，表示节点数量。节点编号 1∼N。

接下来 N 行，每行给出一个节点的信息（第 ii 行对应第 ii 个节点），格式为：

```
data left_child right_child
```

其中，`data` 是一个不超过 1010 个字符的字符串，`left_child` 和 `right_child` 分别是该节点的左右子节点的编号。

没有子节点（即 NULL），则用 −1−1 表示。

下面两图分别对应给出的两个样例。

![4d1c4a98-33cc-45ff-820f-c548845681ba.JPG](https://cdn.acwing.com/media/article/image/2022/01/10/19_ec0f334e71-4d1c4a98-33cc-45ff-820f-c548845681ba.JPG) 

![b5a3c36e-91ad-494a-8853-b46e1e8b60cc.JPG](https://cdn.acwing.com/media/article/image/2022/01/10/19_effd4b2471-b5a3c36e-91ad-494a-8853-b46e1e8b60cc.JPG)

输出格式

在一行中输出答案，表达式符号之间不得有空格。

数据范围

1≤N≤20



输入样例1：

```
8
* 8 7
a -1 -1
* 4 1
+ 2 5
b -1 -1
d -1 -1
- -1 6
c -1 -1
```

输出样例1：

```
(((a)(b)+)((c)(-(d))*)*)
```

输入样例2：

```
8
2.35 -1 -1
* 6 1
- -1 4
% 7 8
+ 2 3
a -1 -1
str -1 -1
871 -1 -1
```

输出样例2：

```
(((a)(2.35)*)(-((str)(871)%))+)
```



#### 题意

**dfs**

本题考查后缀表达式.给你一颗二叉树,输出它的后缀表达式,也就是输出它的后序遍历罢了.而一道算数表达式的后续表达式,就是:

对于一个算术表达式我们的一般写法是这样的

> (3 + 4) × 5 - 6

这种写法是**中序表达式**
而**后序表达式**则是将运算符放在操作数的后面,如

> 3 4 + 5 × 6 -

这种写法转换可以用栈来实现:

[后缀表达式](https://blog.csdn.net/u012507347/article/details/52245233)



本题给定了二叉树,我们就用递归输出它的后序遍历.(还要用括号表现优先级)



满足输出后序表达式的二叉树:其中的的节点有两种情况：

+ 有左右儿子；
+ 只有右儿子。

只有右儿子的情况,这个右儿子就是负号,代表这个右儿子的整个子树表达的是个负数.

如题目所给出的两个例子



```c++
string dfs(int t)
{
     if (rs[t] == -1 && ls[t] == -1)
          return '(' + w[t] + ')';
     if (ls[t] == -1)
          return '(' + w[t] + '(' + dfs(rs[t]) + ')' + ')';
     return '(' + dfs(ls[t]) + dfs(rs[t]) + w[t] + ')';
}
```

`return '(' + w[t] + ')';`用括号把每一个操作数都给括上

`return '(' + dfs(ls[t]) + dfs(rs[t]) + w[t] + ')';`操作符的话在这条语句输出

要是删除最后一行这两个小括号`return dfs(ls[t]) + dfs(rs[t]) + w[t]; `,就会输出:`(a)(b)+(c)(-(d))**`

原样:`(((a)(b)+)((c)(-(d))*)*)`





#### 代码

```c++
#include <bits/stdc++.h>
using namespace std;
const int N = 20 + 10;
int ls[N], rs[N];
string w[N];
int n;
int in[N];

string dfs(int t)
{
     if (rs[t] == -1 && ls[t] == -1)
          return '(' + w[t] + ')';
     if (ls[t] == -1)
          return '(' + w[t] + '(' + dfs(rs[t]) + ')' + ')'; //碰到负号情况
     return '(' + dfs(ls[t]) + dfs(rs[t]) + w[t] + ')';
}

int main()
{
     int n;
     cin >> n;
     for (int i = 1; i <= n; i++)
     {
          string a;
          int b, c;
          cin >> a >> b >> c;
          w[i] = a;
          ls[i] = b, rs[i] = c;
          in[b]++, in[c]++;
     }
     int root;
     for (int i = 1; i <= n; i++)	//找根
     {
          if (in[i] == 0)
               root = i;
     }

     cout << dfs(root);
}

```

