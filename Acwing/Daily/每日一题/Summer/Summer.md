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







### 4274:后缀表达式





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









### 4275:dijkstra序列



#### 题目

Dijkstra 算法是非常著名的贪心算法之一。

它用于解决单源最短路径问题，即指定一个特定源顶点，求该顶点到给定图的所有其他顶点的最短路径。

它由计算机科学家 Edsger W. Dijkstra 于 1956 年构思并在三年后出版。

在该算法中，我们需要不断维护一个包含最短路径树中顶点的集合。

在每一步中，我们找到一个尚未在集合内且与源顶点距离最小的顶点，并将其收于集合中。

因此，通过 Dijkstra 算法，我们可以逐步生成一个有序的顶点序列，我们称之为 Dijkstra 序列。

对于一个给定的图，可能有多个 Dijkstra 序列。

例如，{5,1,3,4,2} 和 {5,3,1,2,4}都是给定图的 Dijkstra 序列。

注意，序列中的第一个顶点即为指定的特定源顶点。

你的任务是检查给定的序列是否是 Dijkstra 序列。



输入格式

第一行包含两个整数 N 和 M，表示图中点和边的数量。

点的编号 1∼N。

接下来 MM 行，每行包含三个整数 a,b,c，表示点 a 和点 b 之间存在一条**无向**边，长度为 c。

再一行包含整数 K，表示需要判断的序列个数。

接下来 K 行，每行包含一个 1∼N 的排列，表示一个给定序列。

输出格式

共 K 行，第 ii 行输出第 K 个序列的判断，如果序列是 Dijkstra 序列则输出 `Yes`，否则输出 `No`。



数据范围

1≤N≤1000,1≤M≤10^5^,1≤a,b≤N,1≤c≤100,1≤K≤100,

保证给定无向图是连通图，

保证无重边和自环（官网没提，但是经实测，官网数据符合此条件）。



输入样例：

```
5 7
1 2 2
1 5 1
2 3 1
2 4 1
2 5 2
3 5 1
3 4 1
4
5 1 3 4 2
5 3 1 2 4
2 3 4 5 1
3 2 1 5 4
```

输出样例：

```
Yes
Yes
Yes
No
```



#### 思路

orz

一道稍微变形的dij问题.(稍微变通一下,但核心代码都与朴素版dij差不多)

题目给你了一个已知序列(一条路径),让你判断这是不是dij路径

复杂度是O(m * n * n)



#### 代码

```c++
#include <bits/stdc++.h>
using namespace std;
const int N = 1e3 + 10;
int c[N][N];
int n, m;
int a[N], d[N], st[N];

int dij(int x)
{
     memset(d, 127, sizeof d);
     memset(st, 0, sizeof st);
     d[x] = 0;
     st[x] = 1;
     for (int i = 1; i <= n; i++)
     {
          int t = a[i];
          st[t] = 1;

          for (int j = 1; j <= n; j++)
          {
               if (!st[j] && d[j] < d[t])//发现这个点不是最短距离的点,直接退出
                    return 0;
               
          }
          for (int j = 1; j <= n; j++)
          {
               if (c[t][j] != 0 && d[j] > d[t] + c[t][j])
                    d[j] = d[t] + c[t][j];
          }
     }
     return 1;
}

int main()
{
     cin >> n >> m;
     for (int i = 0; i < m; i++)
     {
          int x, y, z;
          cin >> x >> y >> z;
          c[x][y] = z;
          c[y][x] = z;
     }

     int k;
     cin >> k;
     for (int i = 0; i < k; i++)
     {
          for (int j = 1; j <= n; j++)
               cin >> a[j];
          if (dij(a[1]))
               cout << "Yes" << endl;
          else
               cout << "No" << endl;
     }
}
```









### 4276:删除C



#### 题目



​	当你被面试官要求用 CC 写一个 `Hello World` 时，有本事像下图显示的那样写一个出来吗？

![ba3b8678-061d-4fc6-a87e-ce08e1434410.jpg](https://cdn.acwing.com/media/article/image/2022/01/10/19_0b2aaf6771-ba3b8678-061d-4fc6-a87e-ce08e1434410.jpg)



输入格式

​	输入首先给出 2626 个英文大写字母 A−ZA−Z，每个字母用一个 7×57×5 的、由 `C` 和 `.` 组成的矩阵构成。

​	最后在一行中给出一个句子，以回车结束。句子是由若干个单词（每个包含不超过 10 个连续的大写英文字母）组成的，单词间以任何非大写英文字母分隔。

​	题目保证至少给出一个单词。

输出格式

​	对每个单词，将其每个字母用矩阵形式在一行中输出，字母间有一列空格分隔。单词的首尾不得有多余空格。

​	相邻的两个单词间必须有一空行分隔。输出的首尾不得有多余空行。



数据范围

​	最后一行句子的总长度范围 \[1,5000]，
​	给出的单词数量范围 \[1,300]。



输入样例：(前面输入的就是ABCD...)

```
..C..
.C.C.
C...C
CCCCC
C...C
C...C
C...C
CCCC.
C...C
C...C
CCCC.
C...C
C...C
CCCC.
.CCC.
C...C
C....
C....
C....
C...C
.CCC.
CCCC.
C...C
C...C
C...C
C...C
C...C
CCCC.
CCCCC
C....
C....
CCCC.
C....
C....
CCCCC
CCCCC
C....
C....
CCCC.
C....
C....
C....
CCCC.
C...C
C....
C.CCC
C...C
C...C
CCCC.
C...C
C...C
C...C
CCCCC
C...C
C...C
C...C
CCCCC
..C..
..C..
..C..
..C..
..C..
CCCCC
CCCCC
....C
....C
....C
....C
C...C
.CCC.
C...C
C..C.
C.C..
CC...
C.C..
C..C.
C...C
C....
C....
C....
C....
C....
C....
CCCCC
C...C
C...C
CC.CC
C.C.C
C...C
C...C
C...C
C...C
C...C
CC..C
C.C.C
C..CC
C...C
C...C
.CCC.
C...C
C...C
C...C
C...C
C...C
.CCC.
CCCC.
C...C
C...C
CCCC.
C....
C....
C....
.CCC.
C...C
C...C
C...C
C.C.C
C..CC
.CCC.
CCCC.
C...C
CCCC.
CC...
C.C..
C..C.
C...C
.CCC.
C...C
C....
.CCC.
....C
C...C
.CCC.
CCCCC
..C..
..C..
..C..
..C..
..C..
..C..
C...C
C...C
C...C
C...C
C...C
C...C
.CCC.
C...C
C...C
C...C
C...C
C...C
.C.C.
..C..
C...C
C...C
C...C
C.C.C
CC.CC
C...C
C...C
C...C
C...C
.C.C.
..C..
.C.C.
C...C
C...C
C...C
C...C
.C.C.
..C..
..C..
..C..
..C..
CCCCC
....C
...C.
..C..
.C...
C....
CCCCC
HELLO~WORLD!
```

输出样例：

```
C...C CCCCC C.... C.... .CCC.
C...C C.... C.... C.... C...C
C...C C.... C.... C.... C...C
CCCCC CCCC. C.... C.... C...C
C...C C.... C.... C.... C...C
C...C C.... C.... C.... C...C
C...C CCCCC CCCCC CCCCC .CCC.

C...C .CCC. CCCC. C.... CCCC.
C...C C...C C...C C.... C...C
C...C C...C CCCC. C.... C...C
C.C.C C...C CC... C.... C...C
CC.CC C...C C.C.. C.... C...C
C...C C...C C..C. C.... C...C
C...C .CCC. C...C CCCCC CCCC.
```





#### 思路

一道模拟题~(真ex,有许多小细节

[对应不同题意的输入格式](https://blog.csdn.net/qq_46009744/article/details/121705608?spm=1001.2014.3001.5502)

题目说是输出出现的单词,然后换行,再输出下一个单词.所以把当出现一个完整单词的时候,我们就把它存储进字符串,然后一个一个字母输出.再输出下一个单词.

前提是先储存好题目给你的每一个字母的矩阵表达。g\[t]\[x][y]:t表示哪一个字母,1代表A,2代表B x和y就用来存储这个字母的矩阵形式.

​	要用getchar读入字符串,防止cin过滤掉空格(getchar也可以读入换行)

```c++
    while ((c = getchar()) != EOF) //也可以写成(c = getchar()) != EOF && c!='\n'),不能用cin和scanf,会过滤空格 .当然!=EOF和!=-1是同一个意思
    
        //利用while不断读入字符,但为了防止它无休止的读下去,加上!=EOF.它是文本结束符,在自己终端这用ctrl+z表示
     {
          if (c >= 'A' && c <= 'Z')
               word += c;
          else
          {
               output(word);
               word = "";
          }
     }
```

[粗略讲解EOF用法](https://blog.csdn.net/henu1710252658/article/details83040281)



#### 代码

```c++
#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

char g[26][7][6];
bool is_first = true;
string word; //存储每个单词,当出现一整个单词的时候,就把他输出
void output(string word)
{
     if (word.empty())
          return;

     if (is_first)
          is_first = false;
     else
          cout << endl;

     char str[7][60] = {0};
     for (int i = 0; i < word.size(); i++)
          for (int j = 0; j < 7; j++)
               for (int k = 0; k < 5; k++)
                    str[j][i * 6 + k] = g[word[i] - 'A'][j][k];

     for (int i = 1; i < word.size(); i++) //中间夹着一行空格
          for (int j = 0; j < 7; j++)
               str[j][i * 6 - 1] = ' ';

     for (int i = 0; i < 7; i++)
          cout << str[i] << endl;
}

int main()
{
     for (int i = 0; i < 26; i++) //把每个字符的格式存储起来
          for (int j = 0; j < 7; j++)
               cin >> g[i][j];

     string word;
     char c;
     // getchar(); 如果这里加上getchar,会清除掉上面cin遗留在缓冲区的换行
     while ((c = getchar()) != EOF) //也可以写成(c = getchar()) != EOF && c!='\n'),不能用cin和scanf,会过滤空格
     //利用while不断读入字符,但为了防止它无休止的读下去,加上!=EOF.它是文本结束符.而在自己终端要表示出EOF可以键入ctrl+z,再换行
     {
          if (c >= 'A' && c <= 'Z')
               word += c;
          else
          {
               output(word);
               word = "";
          }
     }

     output(word);

     return 0;
}
```



### 4278:峰会

#### 题目



峰会是国家元首或政府首脑的会议。

为峰会安排休息区可不是一件简单的工作。

一共有 N个首脑参加峰会，编号 1∼N。

这些首脑之间存在 MM 对两两之间的直接朋友关系。

在划分区域时，我们希望被安排在同一休息区域的首脑们满足，任意两人之间都是直接朋友关系。

现在，给定 K 个关于划分休息区域的安排，请你依次判断每个安排是否合理。



输入格式

第一行包含两个整数 N 和 M。

接下来 M 行，每行包含两个整数 a,b，表示首脑 a 和首脑 b 之间存在直接朋友关系。

再一行包含整数 K。

接下来 K 行，每行描述一个区域安排，首先包含一个整数 L，表示该安排打算将 L 个首脑安排在同一区域休息，然后包含 L 个整数，表示这些首脑的编号。

输出格式

共 K 行，第 i 行输出对第 i 个安排的判断，具体格式为

- 如果安排满足其中的任意两人之间都是直接朋友关系并且不存在额外的人与被安排的所有人都是直接朋友关系（即无法安排更多的人在这一区域休息），则输出 `Area X is OK.`
- 如果安排满足其中的任意两人之间都是直接朋友关系并且存在额外的人与被安排的所有人都是直接朋友关系（即可以安排更多的人在这一区域休息），则输出 `Area X may invite more people, such as H.`，其中 `H` 是额外可被安排的人的编号（如果不唯一，则输出最小的那个）。
- 如果安排无法满足其中的任意两人之间都是直接朋友关系，则输出 `Area X needs help.`。

`X` 表示组别编号，从 1 到 K。



数据范围

1≤N≤200,
1≤M≤(N(N−1))/2,
1≤a,b≤N,
a≠b,
1≤K≤100,
1≤L≤N,
同一对直接朋友关系不会在输入中重复出现。



输入样例：

```
8 10
5 6
7 8
6 4
3 6
4 5
2 3
8 2
2 7
5 3
3 4
6
4 5 4 3 6
3 2 8 7
2 2 3
1 1
2 4 6
3 3 2 1
```

输出样例：

```
Area 1 is OK.
Area 2 is OK.
Area 3 is OK.
Area 4 is OK.
Area 5 may invite more people, such as 3.
Area 6 needs help.
```



#### 思路

Orz

看到是小数据,就采用暴力枚举

```c++
任意两个的关系
/*
for (int i = 1; i <= n; i++)
     for (int j = i + 1; j <= n; j++)
*/
    
/*
for (int i = 1; i <= n; i++)
     for (int j = 1; j <= n; j++)
*/
    
```

时间复杂度:O(n^2^k)

但有一些小细节(要增强写代码能力)



#### 代码

```c++
//任意两个的关系
/*
for (int i = 1; i <= n; i++)
     for (int j = i + 1; j <= n; j++)
*/
//数据小:暴力枚举

#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

const int N = 210;

int n, m;
bool g[N][N];
bool st[N];

int main()
{
     scanf("%d%d", &n, &m);
     while (m--)
     {
          int a, b;
          scanf("%d%d", &a, &b);
          g[a][b] = g[b][a] = true;
     }

     scanf("%d", &m);
     for (int T = 1; T <= m; T++)
     {
          int cnt;
          scanf("%d", &cnt);
          memset(st, 0, sizeof st);
          while (cnt--)
          {
               int x;
               scanf("%d", &x);
               st[x] = true;
          }

          bool is_clique = true;
          for (int i = 1; i <= n; i++)
               for (int j = i + 1; j <= n; j++)
                    if (st[i] && st[j] && !g[i][j])
                         is_clique = false;
          if (!is_clique)
               printf("Area %d needs help.\n", T);
          else
          {
               int id = 0;
               for (int i = 1; i <= n; i++)
                    if (!st[i]) // i要不存在
                    {
                         bool all = true;
                         for (int j = 1; j <= n; j++) //这里就不是j=i+1了
                              if (st[j] && !g[i][j])
                              {
                                   all = false;
                                   break;
                              }
                         if (all)
                         {
                              id = i;
                              break;
                         }
                    }
               if (id)
                    printf("Area %d may invite more people, such as %d.\n", T, id);
               else
                    printf("Area %d is OK.\n", T);
          }
     }

     return 0;
}

```









### 4279:笛卡尔树

#### 题目

[笛卡尔树](https://baike.baidu.com/item/笛卡尔树/7579802?fr=aladdin) 是由一系列不同数字构成的二叉树。

树满足堆的性质，中序遍历返回原始序列。
最小笛卡尔树表示满足小根堆性质的笛卡尔树。

例如，给定序列 {8,15,3,4,1,5,12,10,18,6}则生成的最小堆笛卡尔树如图所示。

![6a99f68a-6578-46e0-9232-fbf0adf3691f.jpg](https://cdn.acwing.com/media/article/image/2022/01/11/19_e5d5fcaf72-6a99f68a-6578-46e0-9232-fbf0adf3691f.jpg)

现在，给定一个长度为 N 的原始序列，请你生成最小堆笛卡尔树，并输出其层序遍历序列。



输入格式

第一行包含整数 N。

第二行包含 N 个两两不同的整数，表示原始序列。

输出格式

共一行，输出最小堆笛卡尔树的层序遍历序列。



数据范围

1≤N≤30,
原始序列中元素的取值范围 \[−2147483648,2147483647]。



输入样例：

```
10
8 15 3 4 1 5 12 10 18 6
```

输出样例：

```
1 3 5 8 4 6 15 10 12 18
```



#### 思路

​	只已知中序遍历,我们是无法确定一棵树的,因为可以选任何一个点当root.(要知道中序+后序,中序+前序才能构建),而本题还给了你一个要求,构建最小笛卡尔树(小根堆性质:root小于子节点),这样我们才能确定一棵树:

​	找到序列最小的数字，作为分割点，它左边的作为左儿子，它右边的作为右儿子,递归下去.

​	一般要建树,可能要构造数据结构(采用指针):

```c++
typedef struct TreeNode {
    int data;
    struct TreeNode *left;
    struct TreeNode *right;
    struct TreeNode *parent;
} TreeNode;

void middle_order(TreeNode *Node) {
    if(Node != NULL) {
        middle_order(Node->left);
        printf("%d ", Node->data);
        middle_order(Node->right);
    }
}
```

​	但其实也可以采用数组构建,写的方便一点:

```c++

struct Node {
    int l;
    int r;
    int data;
};
Node tree[N];
```




先用 dfs 建树，然后用 bfs 输出答案

时间复杂度:
无论 dfs 还是 bfs 每个点都只遍历一次，所以时间复杂度 O(n)。

空间复杂度:
dfs 的空间复杂度是树的深度，最多为 O(n)，bfs 的空间复杂度最多是节点数也就是 O(n)，总的空间复杂度就是 O(n)。

​	

#### 代码

```c++
//模仿y总
#include <bits/stdc++.h>
using namespace std;
const int INF = 2147483647;
const int N = 15;
int w[N];
vector<int> v[N];
int f = 0;
void dfs(int l, int r, int d)//还有一个层数参数
{
     if (r < l) //!注意
          return;
     if (d > f)
          f = d;
     int minv = INF;
     int idx;
     for (int i = l; i <= r; i++)
     {
          if (w[i] < minv)
          {
               minv = w[i];
               idx = i;
          }
     }
     v[d].push_back(minv);
     dfs(l, idx - 1, d + 1);
     dfs(idx + 1, r, d + 1);
}

int main()
{
     int n;
     cin >> n;
     for (int i = 1; i <= n; i++)
          cin >> w[i];
     dfs(1, n, 1);
     for (int i = 1; i <= f; i++)
          for (auto x : v[i])
               cout << x << ' ';
}
```



```c++
//dfs+bfs,wp

#include <iostream>
#include <queue> // bfs要用
using namespace std;

const int N = 35;

int n;
int s[N]; // 原始序列

struct Node
{
     int l, r, dat;
} tree[N]; // 笛卡尔树
int cnt;   // 树的长度，便于加入元素

int dfs(int l, int r) // 建s[l] ~ s[r]的笛卡尔树，返回根节点
{
     if (l >= r)
          return -1; // 如果区间不合法返回-1

     int minv = 2147483647, idx;     // minv: s[l] ~ s[r]中最小的数，idx: 它的编号
     for (int i = l; i < r; i++)     // 遍历s[l] ~ s[r]
          if (minv >= s[i])          // 如果能更新
               minv = s[i], idx = i; // 更新

     int tmp = cnt; // 根节点编号需要保存一下，因为往深层递归的时候可能会改动cnt

     tree[cnt++] = {dfs(l, idx), dfs(idx + 1, r), minv}; // 往深层递归

     return tmp; // 返回根节点编号
}

void bfs() // 输出最小堆笛卡尔树的层序遍历序列
{
     queue<int> q; // bfs要用的队列
     q.push(0);    // 把根节点加入队列（因为根节点是最先遍历的，所以它的编号总是0

     while (!q.empty()) // 还能扩展
     {
          int t = q.front();          // 要扩展的点
          cout << tree[t].dat << ' '; // 输出
          q.pop();                    // 出队

          if (~tree[t].l)
               q.push(tree[t].l); // 如果左儿子非空则加入队列
          if (~tree[t].r)
               q.push(tree[t].r); // 同理
     }
}

int main()
{
     cin >> n;
     for (int i = 0; i < n; i++)
          cin >> s[i]; // 输入原始序列

     dfs(0, n); // dfs
     bfs();     // bfs

     return 0;
}
```







### 691:立方体IV





#### 思路

​	忘了,Orz

​	从题干可以看出,走过的格子绝对不会再走第二遍,可以想到记忆化搜索

​	时间复杂度n^2^.

​	本题类似题目:滑雪,dfs类似返还:子树大小,



​	读入每个点的数值后，dfs每一个点最多能到达的房间数量。
​	dfs时，如果f值非空（即dfs前面的点时，已经遍历过该点），可直接返回该点的f值，以免重复搜索该点。
​	否则，先赋初值1（即该房间）。然后依次遍历4个方向，如果到达的点的g值比出发点大1，则dfs该点，同时该点取dfs返回值+1与原先的点的最大值
最后根据要求依次遍历每个点即可找到最大值。

​	详细的记忆化搜索讲解可以查看 [901. 滑雪 - AcWing题库](https://www.acwing.com/problem/content/903/)



#### 代码

```c++
//dfs+记忆化搜索
//f[i][j]:存储以这个点为起点,最长走过几个房间
#include <algorithm>
#include <bitset>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
typedef pair<int, int> pii;
typedef vector<int> vi;
const int MAXN = 1005;
const int dx[4] = {-1, 0, 1, 0}, dy[4] = {0, 1, 0, -1};

int n;
int g[MAXN][MAXN], f[MAXN][MAXN];

bool check(int x, int y)
{
     if (x < 0 || x >= n || y < 0 || y >= n)
          return false;
     return true;
}

int dfs(int x, int y) // int,类似返还子树大小
{
     if (f[x][y])
          return f[x][y];
     f[x][y] = 1; //初始化
     for (int i = 0; i < 4; i++)
     {
          int xx = x + dx[i], yy = y + dy[i];
          if (check(xx, yy) && g[xx][yy] == g[x][y] + 1)
               f[x][y] = f[x][y] + dfs(xx, yy); //子树大小
         		//f[x][y]=max(f[x][y],dfs(xx,yy)+1)
     }
     return f[x][y];
}

/*

int dfs(int x, int y, int sum)
{
    for (int i = 0; i < 4; i++)
    {
        int nx = x + dx[i], ny = y + dy[i];
        if (nx < 1 || nx > n || ny < 1 || ny > n) continue;
        if (g[nx][ny] != g[x][y] - 1) continue;

        return dfs(nx, ny, sum + 1);
    }

    if (maxn <= sum)
    {
        ans = g[x][y];
        maxn = sum;
    }

    return sum;
}
*/

int main()
{
     int T;
     scanf("%d", &T);
     for (int t = 1; t <= T; t++)
     {
          scanf("%d", &n);
          for (int i = 0; i < n; i++)
               for (int j = 0; j < n; j++)
                    scanf("%d", &g[i][j]);
          memset(f, 0, sizeof f);
          for (int i = 0; i < n; i++)
               for (int j = 0; j < n; j++)
                    dfs(i, j);
          int res = 0, id = n * n;
          for (int i = 0; i < n; i++)
               for (int j = 0; j < n; j++)
                    if (f[i][j] > res || f[i][j] == res && g[i][j] < id)
                    {
                         res = f[i][j];
                         id = g[i][j];
                    }
          printf("Case #%d: %d %d\n", t, id, res);
     }
     return 0;
}

```





### 3311:最长算数





#### 思路

简单的双指针题目,O(n)时间复杂度

双指针大多用到for+while

刚开始还是看了一眼标程,感觉要是自己写,结果可能就变成了判断a[i]等不等于a[j]了hhhh,j动一下i动一下,会被绕进去,而不是直接清爽的写出`a[j] - a[j - 1] == a[j - 1] - a[j - 2]`



#### 代码

```c++

#include <bits/stdc++.h>
using namespace std;
const int N = 2e5 + 10;
int a[N];
int main()
{
     int t;
     cin >> t;
     for (int e = 1; e <= t; e++)
     {
          int n;
          cin >> n;
          for (int i = 1; i <= n; i++)
               cin >> a[i];
          int maxx = 0;
          for (int i = 1; i <= n; i++)
          {
               int j = i + 2; //此时初始最长算数长度为:2(j-i)
               while (j <= n)
               {
                    if (a[j] - a[j - 1] == a[j - 1] - a[j - 2])//直接两段一起比
                         j++;
                    else
                         break;
               }
               if (j - i >= maxx)
               {
                    maxx = j - i;
               }
               i = j - 2;
          }
          printf("Case #%d: %d\n", e, maxx);
     }
}
```







### 行程排序



#### 思路

​	先建立(哈希表)每一张机票的起点和终点的映射关系,然后找到最初的起点,最后直接一条连着一条输出就可以了

​	小模拟



#### 代码

```c++
#include <bits/stdc++.h>
using namespace std;
const int N = 1e4 + 10;

string b[N];
int main()
{

     int t;
     cin >> t;
     for (int x = 0; x < t; x++)
     {
          map<string, int> s, e;
          map<string, string> a;
          int n;
          cin >> n;
          for (int i = 0; i < n; i++)
          {
               string s1, s2;
               cin >> s1 >> s2;
               b[i] = s1;     
               a[s1] = s2; //把一张机票的起点和终点映射起来
               e[s2]++;    //找最初的起点
          }
          printf("Case #%d: ", x + 1);
          string start1;
          int j = 0;
          for (int i = 0; i < n; i++)
          {
               if (e[b[i]] == 0)
               {
                    start1 = b[i];
                    while (j < n) //找到起点后直接输出答案
                    {
                         cout << start1 << '-' << a[start1] << ' ';
                         start1 = a[start1];
                         j++;
                    }
                    break;
               }
          }
          cout << endl;
     }
}
```

