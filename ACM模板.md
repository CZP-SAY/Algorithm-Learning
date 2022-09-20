# ACM模板

## CF代码模板

```cpp
#pragma GCC target("avx")
#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")

#include <bits/stdc++.h>
#define between(x, l, r) (x >= l && x <= r)
#define point(a, b, m) ((a - 1) * m + b)
#define IDX(x, m) (PII) { (x - 1) / m + 1, (x - 1) % m + 1 }
#define all(x) (x).begin(), (x).end()
#define np next_permutation
#define lb lower_bound
#define ub upper_bound
#define eb emplace_back
#define pb push_back
#define mk make_pair
#define y second
#define x first

using namespace std;

typedef long long LL;
typedef pair<long double, long double> PDD;
typedef pair<int, int> PII;
mt19937 mrand(random_device{}());
int rnd(int x) { return mrand() % x;}

constexpr int N = 1e6 + 10, M = 4 * N, mod = 1e9 + 7;
constexpr LL INF = 1e18;
constexpr long double eps = 1e-9;

void solve()
{

}

int main()
{
    // freopen("C:\\AKICPC\\test5.in", "r", stdin), freopen("C:\\AKICPC\\test5.out", "w", stdout);
    ios::sync_with_stdio(0), cin.tie(0);
    cout.precision(2), cout.setf(ios::fixed);
    // auto now = clock();

    int t = 1;
    cin >> t;
    while (t--) solve();

    // cout << "Program run for " << (clock() - now) / (long double)CLOCKS_PER_SEC * 1000 << " ms." << endl;
    return 0;
}
```



## 数学

### 矩阵

**用法:**  具体可以参考CF718C。

```cpp
constexpr int mod = 1e9 + 7;

#define T (*this)
struct matrix
{
    static constexpr int n = 2;
    int a[n][n];

    matrix()
    {
        for (auto &i : a)
            for (auto &j : i)
                j = 0;
    }

    int *operator[](const int i) { return a[i]; }

    const int *operator[](const int i) const { return a[i]; }

    friend matrix operator*(const matrix &a, const matrix &b)
    {
        matrix ans;
        for (int i = 0; i < n; i ++ )
            for (int k = 0; k < n; k ++ )
            {
                if(a[i][k] == 0) continue;
                for (int j = 0; j < n; j++)
                    ans[i][j] = (ans[i][j] + (LL)a[i][k] * b[k][j]) % mod;
            }
        return ans;
    }

    friend matrix operator+(const matrix &a, const matrix &b)
    {
        matrix ans;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                ans[i][j] = a[i][j] + b[i][j];
                if (ans[i][j] >= mod) ans[i][j] -= mod;
            }
        return ans;
    }

    bool operator!=(const matrix &a)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (T[i][j] != a[i][j])
                    return 1;
        return 0;
    }

    matrix fpow(LL r) const;
};

const matrix E = []
{
    matrix E;
    for (int i = 0; i < matrix::n; i ++ ) E[i][i] = 1;
    return E;
}();

matrix matrix::fpow(LL r) const
{
    matrix ans = E, x = T;
    while (r)
    {
        if (r & 1) ans = ans * x;
        r >>= 1;
        x = x * x;
    }
    return ans;
}

const matrix F = [] // 转移
{
    matrix F;
    F[0][0] = F[0][1] = F[1][0] = 1;
    return F;
}();

matrix c[N << 2], lazy[N << 2];
#undef T

```



### 快速幂

```
LL qmi(LL a, LL k, LL p)
{
    LL res = 1;
    while (k)
    {
        if (k & 1)res = res * a % p;
        k >>= 1;
        a = a * a % p;
    }
    return res;
}
```



### 组合数预处理

```c++
const int mod = 1e9 + 7;

LL inv[N], fac[N], ifac[N];

void init(int n)
{
    inv[1] = fac[0] = fac[1] = ifac[0] = ifac[1] = 1;
    for (int i = 2; i <= n; i ++ )
    {
        inv[i] = (mod - mod / i) * inv[mod % i] % mod;
        fac[i] = fac[i - 1] * i % mod;
        ifac[i] = ifac[i - 1] * inv[i] % mod;
    }
}

LL C(int m, int n) { if(m < n) return 0; return fac[m] * ifac[n] % mod * ifac[m - n] % mod; }

```



### μ函数线性筛&杜教筛

```
const int N = 5e6 + 10;

int n, p;
bool st[N];
int mu[N], primes[N], cnt;
map<LL, LL> mp;

void get_mu(int n)
{
    mu[1] = 1;
    for(int i = 2; i <= n; i ++ )
    {
        if(!st[i]) primes[++ cnt] = i, mu[i] = -1;
        
        for(int j = 1; j <= cnt && primes[j] * i <= n; j ++ )
        {
            st[primes[j] * i] = 1;
            if(i % primes[j] == 0) break;
            else mu[i * primes[j]] = -mu[i];
        }
    }
    
    for (int i = 1; i < n; i ++ ) mu[i] += mu[i - 1]; // mu的前缀和
}

LL Mu(LL n) // mu的前缀和
{
    if(n < N) return mu[n];
    if(mp[n])return mp[n];
    LL s = 0;
    for (LL l = 1, r; l <= n; l = r + 1)
    {
    	r = n / (n / l);
    	s += (r - l + 1) * Mu(n / x);
    }
    mp[n] = 1 - s;
    return mp[n];
}
```



### 线性筛与因式分解

```cpp
int vis[M], primes[M], minp[M], cnt;

void get_primes(int n)  // 线性筛质数
{
    minp[1] = 1;
    for (int i = 2; i <= n; i ++ )
    {
        if (!vis[i]) primes[cnt ++ ] = i, minp[i] = i;
        for (int j = 0; i * primes[j] <= n; j ++ )
        {
            vis[primes[j] * i] = true, minp[i * primes[j]] = primes[j];
            if (i % primes[j] == 0) break;
        }
    }
}

vector<PII> get(int x)
{
    vector<PII> res;
    while(x > 1)
    {
        int cnt = 0, cur = minp[x];
        while(cur == minp[x]) cnt ++, x /= minp[x];
        res.push_back({cur, cnt});
    }
    return res;
}

```



### 容斥求互质数量

```cpp
int vis[M], primes[M], minp[M], cnt;

void get_primes(int n)  // 线性筛质数
{
    minp[1] = 1;
    for (int i = 2; i <= n; i ++ )
    {
        if (!vis[i]) primes[cnt ++ ] = i, minp[i] = i;
        for (int j = 0; i * primes[j] <= n; j ++ )
        {
            vis[primes[j] * i] = true, minp[i * primes[j]] = primes[j];
            if (i % primes[j] == 0) break;
        }
    }
}

vector<int> get(int x) //获得x的质因子
{
    vector<int> res;
    while(x > 1)
    {
        int cur = minp[x];
        while(cur == minp[x]) x /= minp[x];
        res.push_back(cur);
    }
    return res;
}


int get_coprime(int x, int n) // [1, n]中与x互质的数
{
    vector<int> p = get(x); // 得到x的质因子
    int k = p.size();
    int res = 0;
    for (int i = 1; i < (1 << k); i ++ )
    {
        int cnt = 0, t = 1; //被选取的质因子的数量 以及 所选取的质因子的乘积

        for (int j = 0; j < k; j ++ )
            if(i >> j & 1) // i的二进制表示下第j位为1
            {
                cnt ++ ;
                t *= p[j];
            }
        // (n / t) 是[1, n]内被t整除的数 的数量
        if(cnt & 1) res += n / t;
        else res -= n / t;
    }
    return n - res;
}
```



### exgcd

```cpp
void exgcd(LL a, LL b, LL &d, LL &x, LL &y)
{
    if (!b)
    {
        d = a;
        x = 1;
        y = 0;
    }
    else
    {
        exgcd(b, a % b, d, y, x);
        y -= x * (a / b);
    }
}

LL inv(LL a, LL p) // -1 代表没有逆元，同 gcd(a,p) == 1
{
    LL d, x, y;
    exgcd(a, p, d, x, y);
    return d == 1 ? (x + p) % p : -1;
}

LL MinPair(LL a, LL b, LL n, LL &x, LL& y) //最小非负整数x满足 ax + by = n
{
    LL d;
    exgcd(a, b, d, x, y);
    x *= n / d;
    LL t = b / d;
    x = (x % t + t) % t;
    y = (n - a * x) / b;
    return d;
}

```



### 多项式

#### FFT

```c++
using db = double;
struct cp {
    db x, y;
    cp(db real = 0, db imag = 0) : x(real), y(imag){};
    cp operator+(cp b) const { return {x + b.x, y + b.y}; }
    cp operator-(cp b) const { return {x - b.x, y - b.y}; }
    cp operator*(cp b) const { return {x * b.x - y * b.y, x * b.y + y * b.x}; }
};
using vcp = vector<cp>;
using Poly = vector<LL>;
namespace FFT {
    const db pi = acos(-1);
    vcp Omega(int L) {
        vcp w(L); w[1] = 1;
        for (int i = 2; i < L; i <<= 1) {
            auto w0 = w.begin() + i / 2, w1 = w.begin() + i;
            cp wn(cos(pi / i), sin(pi / i));
            for (int j = 0; j < i; j += 2)
                w1[j] = w0[j >> 1], w1[j + 1] = w1[j] * wn;
        }
        return w;
    }
    auto W = Omega(1 << 20); // NOLINT
    void DIF(cp *a, int n) {
        cp x, y;
        for (int k = n >> 1; k; k >>= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    x = a[i + j], y = a[i + j + k],
                    a[i + j + k] = (a[i + j] - y) *  W[k + j], a[i + j] = x + y;
    }
    void IDIT(cp *a, int n) {
        cp x, y;
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    x = a[i + j], y = a[i + j + k] * W[k + j],
                    a[i + j + k] = x - y, a[i + j] = x + y;
        const db Inv = 1. / n;
        for(int i = 0;i <= n - 1; i ++ ) a[i].x *= Inv, a[i].y *= Inv;
        reverse(a + 1, a + n);
    }
}
namespace Polynomial {
    // basic operator
    void DFT(vcp &a) { FFT::DIF(a.data(), a.size()); }
    void IDFT(vcp &a) { FFT::IDIT(a.data(), a.size()); }
    int norm(int n) { return 1 << (__lg(n - 1) + 1); }

    // Poly mul
    vcp &dot(vcp &a, vcp &b) { for(int i = 0; i <= a.size() - 1; i ++ ) a[i] = a[i] * b[i]; return a; }
    Poly operator*(Poly &a, Poly &b) {
        int n = a.size() + b.size() - 1;
        vcp c(norm(n));
        for(int i = 0;i <= a.size() - 1; i ++ ) c[i].x = a[i];
        for(int i = 0;i <= b.size() - 1; i ++ ) c[i].y = b[i];
        DFT(c), dot(c, c), IDFT(c), a.resize(n);
        for(int i = 0;i <= n - 1; i ++ ) a[i] = LL(c[i].y * .5 - .5);
        return a;
    }
}
using namespace Polynomial;

```



#### NTT全家桶

```c++
namespace polybase {//范围为1e9需要先取模,别忘记改模数和原根, 以及数组大小
    constexpr LL mod = 998244353, g = 3, L = 1 << 22;

    LL fpow(LL x, LL r) // 快速幂
    {
        LL result = 1;
        while (r)
        {
            if (r & 1) result = result * x % mod;
            r >>= 1;
            x = x * x % mod;
        }
        return result;
    }

    int w[L], _ = []
    {
        LL x = fpow(g, (mod - 1) / L); w[L / 2] = 1;
        for (int i = L / 2 + 1; i < L; i ++ ) w[i] = w[i - 1] * x % mod;
        for (int i = L / 2 - 1; i >= 0; i -- ) w[i] = w[i << 1];
        return 0;
    }();

    inline int norm(int n) { return 1 << __lg(n * 2 - 1); }

    void dft(LL *a, int n) // 多项式 -> 点
    {
        for (int k = n >> 1; k; k >>= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; j ++ )
                {
                    LL &x = a[i + j], y = a[i + j + k];
                    a[i + j + k] = (x - y + mod) * w[k + j] % mod;
                    x += y;
                    if (x >= mod) x -= mod;
                }
    }

    void idft(LL *a, int n) // 点 -> 多项式
    {
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; j ++ )
                {
                    LL x = a[i + j], y = a[i + j + k] * w[k + j] % mod;
                    a[i + j + k] = x - y < 0 ? x - y + mod : x - y;
                    a[i + j] += y;
                    if (a[i + j] >= mod) a[i + j] -= mod;
                }

        for (int i = 0, inv = mod - (mod - 1) / n; i < n; i ++ ) a[i] = a[i] * inv % mod;
        reverse(a + 1, a + n);
    }

    struct poly : public vector<LL>
    {
        using vector<LL>::vector;
        #define T (*this)

        poly modxk(int k) const {
            k = min(k, (int)size());
            return poly(begin(), begin() + k);
        }

        poly rev() const { return poly(rbegin(), rend()); }
        friend void dft(poly &a) { dft(a.data(), a.size()); }
        friend void idft(poly &a) { idft(a.data(), a.size()); }

        friend poly operator*(const poly &x, const poly &y)
        {
            if (x.empty() || y.empty()) return poly();

            poly a(x), b(y);

            for (auto &it : a) if(it < 0) it += mod;
            for (auto &it : b) if(it < 0) it += mod;

            int len = a.size() + b.size() - 1;
            int n = norm(len);
            a.resize(n), b.resize(n);
            dft(a), dft(b);
            for (int i = 0; i < n; i ++ ) a[i] = a[i] * b[i] % mod;
            idft(a);
            a.resize(len);
            return a;
        }

        poly operator+(const poly &b)
        {
            poly a(T);
            if (a.size() < b.size()) a.resize(b.size());
            for (int i = 0; i < b.size(); i ++ ) //用b.size()防止越界
            {
                a[i] += b[i];
                if (a[i] >= mod) a[i] -= mod;
            }
            return a;
        }

        poly operator-(const poly &b)
        {
            poly a(T);
            if (a.size() < b.size()) a.resize(b.size());
            for (int i = 0; i < b.size(); i ++ )
            {
                a[i] -= b[i];
                if (a[i] < 0) a[i] += mod;
            }
            return a;
        }

        poly operator*(const LL p)
        {
            poly a(T);
            for (auto &x : a) x = x * p % mod;
            return a;
        }

        poly &operator<<=(int r) { return insert(begin(), r, 0), T; }//注意逗号,F(x)*(x^r)
        poly operator<<(int r) const { return poly(T) <<= r; }
        poly operator>>(int r) const { return r >= size() ? poly() : poly(begin() + r, end()); }
        poly &operator>>=(int r) { return T = T >> r; }//F[x]/(x^r)

        poly deriv() //求导
        {
            if (empty()) return T;
            poly a(size() - 1);
            for (int i = 1; i < size(); i ++ ) //注意是size()
                a[i - 1] = T[i] * i % mod;
            return a;
        }

        poly integ()//积分
        {
            poly a(size() + 1);
            for (int i = 1; i < a.size(); i ++ ) //注意是a.size()
                a[i] = T[i - 1] * fpow(i, mod - 2) % mod;
            return a;
        }

        poly inv(int n) //求逆
        {
            poly a{fpow(T[0], mod - 2)};
            int k = 1;
            while (k < n)
            {
                k <<= 1;
                a = (a * 2 - modxk(k) * a * a).modxk(k);
            }

            return a.modxk(n);
        }

        poly sqrt(int n) //f[0]必须等于1
        {
            poly a{1};
            int k = 1;
            const LL inv2 = fpow(2, mod - 2);
            while (k < n)
            {
                k <<= 1;
                a = ((modxk(k) * a.inv(k)).modxk(k) + a) * inv2;
            }
            return a.modxk(n);
        }

        poly ln(int n) //需要保证f[0]=1
        {
            return (deriv() * inv(n)).integ().modxk(n);
        }

        poly exp(int n)//需要保证f[0]=0
        {
            poly a{1};
            int k = 1;
            while (k < n)
            {
                k <<= 1;
                a = (a * (poly{1} - a.ln(k) + modxk(k))).modxk(k);
            }
            return a.modxk(n);
        }

        poly power(LL k, int n)//若a[0]!=1,只能做k未对mod取模情况,因为(1)处fpow的k要对phi(mod)取模,此时用power2
        {
            if (k == 0)return poly(n) + poly{1};
            int d = 0;
            poly a(T);
            a.resize(n);
            while (d < n && !a[d])d++;
            if (d * k >= n)return poly(n);
            a = poly(a.begin() + d, a.end() - d * (k - 1));
            LL inv = fpow(a[0], mod - 2);
            a = ((a * inv).ln(n) * k).exp(n) * fpow(inv, mod - 1 - k % (mod - 1));//(1)
            a.insert(a.begin(), d * k, 0);
            a.resize(n);
            return a;
        }

        poly power2(LL k, LL k2, int n)//k对mod取模,k2对phi(mod)取模,若k在取模前已经大于等于n并且F[0]==0,则返回全0,见P5273 main函数
        {
            if (k == 0)return poly(n) + poly{fpow(T[0], k2)};
            int d = 0;
            poly a(T);
            a.resize(n);
            while (d < n && !a[d])d++;
            if (d * k >= n)return poly(n);
            a = poly(a.begin() + d, a.end() - d * (k - 1));
            LL inv = fpow(a[0], mod - 2);
            a = ((a * inv).ln(n) * k).exp(n) * fpow(inv, mod - 1 - k2 % (mod - 1));//(1)
            a.insert(a.begin(), d * k, 0);
            a.resize(n);
            return a;
        }

        poly inversion() // 二项式反演, i = k -> n
        {
            if(!size()) return poly();

            //处理阶乘
            vector<LL> fac(size()), ifac(size());

            fac[0] = 1;
            for (int i = 1; i < size(); i ++ ) fac[i] = fac[i - 1] * i % mod;
            ifac[size() - 1] = fpow(fac[size() - 1], mod - 2);
            for (int i = size() - 1; i; i -- ) ifac[i - 1] = ifac[i] * i % mod;

            // 计算答案
            poly A(size()), B(size());
            for (int i = 0; i < size(); i ++ )
            {
                A[i] = T[i] * fac[i] % mod;
                B[i] = (i & 1) ? mod - ifac[i] : ifac[i];
            }

            poly G = A.rev() * B;
            G.resize(size());
            reverse(all(G));

            for (int i = 0; i < size(); i ++ ) G[i] = G[i] * ifac[i] % mod;
            return G;
        }
#undef T
    };

    poly operator/(const poly &x, const poly &y)
    {
        int n = x.size(), m = y.size();
        if (n < m)return {0};
        poly a(x.rev()), b(y.rev());
        int k = norm(n - m + 1);
        a.resize(k);
        return (a * b.inv(k)).modxk(n - m + 1).rev();
    }

    pair<poly, poly> div(const poly &a, const poly &b)
    {
        int m = b.size();
        poly c = a / b;
        return {c, a.modxk(m - 1) - (b * c).modxk(m - 1)};
    }

    poly operator%(const poly &a, const poly &b) { return div(a, b).second; }

    struct SegTree
    {
        vector<poly> p;
        int n, raw_n;

        SegTree(const poly &a)
        {
            n = norm(raw_n = a.size());
            p.resize(n << 1);
            for (int i = 0; i < n; i ++ )
                p[i + n] = poly{1, i < raw_n ? mod - a[i] : 0};
            for (int i = n - 1; i; i--)
            {
                int l = i << 1, r = l | 1, k = (p[l].size() - 1) << 1;
                p[l].resize(k), dft(p[l]);
                p[r].resize(k), dft(p[r]);
                p[i].resize(k);
                for (int j = 0; j < k; j ++ )
                    p[i][j] = p[l][j] * p[r][j] % mod;
                idft(p[i]);
                p[i].emplace_back((p[i][0] - 1 + mod) % mod);
                p[i][0] = 1;
            }
        }

        poly eval(const poly &f)//多点求值
        {
            int m = f.size();
            if (m == 1)return poly(raw_n, f[0]);//raw_n个f[0]
            poly q = f.rev() * p[1].inv(m);
            q.resize(m);
            if (m > n)q >>= m - n;
            else q <<= n - m;
            for (int k = n, o = 1; k > 1; k >>= 1)
                for (int i = 0; i < n; i += k, o++)
                {
                    if (i >= raw_n)continue;
                    LL *a = &q[i], *l = p[o << 1].data(), *r = p[o << 1 | 1].data();
                    dft(a, k);
                    poly x(k), y(k);
                    for (int j = 0; j < k; j ++ )x[j] = a[j] * r[j] % mod;
                    for (int j = 0; j < k; j ++ )y[j] = a[j] * l[j] % mod;
                    idft(x), idft(y);
                    for (int j = k / 2; j < k; j ++ )*a++ = x[j];
                    for (int j = k / 2; j < k; j ++ )*a++ = y[j];
                }
            return q.modxk(raw_n);
        }

        poly interpolate(const poly &b)//快速插值
        {
            assert(b.size() == raw_n);
            poly q = eval(p[1].modxk(raw_n + 1).rev().deriv());
            for (int i = 0; i < raw_n; i ++ )q[i] = fpow(q[i], mod - 2) * b[i] % mod;
            q.resize(n);
            for (int k = 1, h = n >> 1; k < n; k <<= 1, h >>= 1)
                for (int i = 0, o = h; i < n; i += k << 1, o++)
                {
                    if (i >= raw_n)continue;
                    LL *a = &q[i], *b = &q[i + k], *l = p[o << 1].data(), *r = p[o << 1 | 1].data();
                    poly x(k << 1), y(k << 1);
                    for (int j = 0; j < k; j ++ )
                        x[j] = a[j], y[j] = b[j];
                    dft(x), dft(y);
                    for (int j = 0; j < k * 2; j ++ )
                        x[j] = (x[j] * r[j] + y[j] * l[j]) % mod;
                    idft(x);
                    for (int j = 0; j < k * 2; j ++ )a[j] = x[j];
                }
            q.resize(raw_n);
            return q.rev();
        }
    };

    poly cdq(const poly &G, int f0, int n)//F[0...n-1]
    {
        poly F(n);
        F[0] = f0;
        function<void(int, int)> solve = [&](int l, int r)//左闭右开
        {
            if (l + 1 >= r)return;
            int m = (l + r) >> 1;
            solve(l, m);
            poly A = poly(F.begin() + l, F.begin() + m) * poly(G.begin() + 1, G.begin() + r - l);
            for (int i = m; i < r; i ++ )
            {
                F[i] += A[i - l - 1];
                if (F[i] >= mod)F[i] -= mod;
            }
            solve(m, r);
        };
        solve(0, n);//别忘记
        return F;
    }
}

using namespace polybase;
```



#### FWT

```
#include <bits/stdc++.h>
#define cpy(a, b, n) memcpy(a, b, sizeof(LL) * n);
#define maxx(a, b, c) max(max(a, b), c)
#define y second
#define x first

using namespace std;

typedef long long LL;
typedef pair<int, int> PII;

const int N = 2e5 + 10, mod = 998244353, inv2 = 499122177;

LL
Cor[2][2]  ={{1, 0},{1, 1}},
Cand[2][2] ={{1, 1},{0, 1}},
Cxor[2][2] ={{1, 1},{1, mod - 1}},
ICor[2][2] ={{1, 0},{mod - 1, 1}},
ICand[2][2]={{1, mod - 1},{0, 1}},
ICxor[2][2]={{inv2, inv2},{inv2, mod - inv2}};

void FWT(LL *F, LL c[2][2], int n)
{
    for (int len = 1; len < n; len <<= 1)
        for (int p = 0; p < n; p += len << 1)
            for (int i = p; i < p + len; i ++ )
            {
                LL sav0 = F[i];
                F[i] = (c[0][0] * F[i] + c[0][1] * F[i + len]) % mod;
                F[i + len] = (c[1][0] * sav0 + c[1][1] * F[i + len]) % mod;
            }
}

void bitmul(LL *F,LL *G,LL C[2][2],LL IC[2][2], int n)
{
    FWT(F, C, n), FWT(G, C, n);
    for (int i = 0; i < n; i ++ ) F[i] = F[i] * G[i] % mod;
    FWT(F, IC, n);
}

int n;
LL f[N], g[N], a[N], b[N];

void print(LL arr[], int n)
{
    for (int i = 0; i < n; i ++ ) cout << arr[i] << ' ';
    cout << '\n';
}

int main()
{
    ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
    cin >> n;
    n = 1 << n;
    for (int i = 0; i < n; i ++ ) cin >> f[i];
    for (int i = 0; i < n; i ++ ) cin >> g[i];
    cpy(a, f, n); cpy(b, g, n);
    bitmul(a, b, Cor, ICor, n);
    print(a, n);
    cpy(a, f, n); cpy(b, g, n);
    bitmul(a, b, Cand, ICand, n);
    print(a, n);
    cpy(a, f, n); cpy(b, g, n);
    bitmul(a, b, Cxor, ICxor, n);
    print(a, n);
    return 0;
}
```



## 数据结构

### ST表

```cpp
int a[N];
int f[N][20];

LL get_gcd(int l, int r) //注意l<r
{
    int k = __lg(r - l + 1);
    return __gcd(f[l][k], f[r - (1 << k) + 1][k]);
}

void Init_ST(int n)
{
    for (int i = 1; i <= n; i ++ ) f[i][0] = a[i];
    for (int j = 1; j < 20; j ++ )
        for (int i = 1; i + (1 << j) - 1 <= n; i ++ )
            f[i][j] = __gcd(f[i][j - 1], f[i + (1 << (j - 1))][j - 1]);
}
```



### DSU

```cpp
struct DSU {
    vector<int> p, siz;
    DSU(int n) : p(n + 1), siz(n + 1, 1) {iota(all(p), 0);}
    int find(int x) {
        while(p[x] != x) x = p[x] = p[p[x]];
        return p[x];
    }
    int operator[](const int x) { return find(x); }
    bool same(int x, int y) {return find(x) == find(y);}
    bool merge(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) return false;
        siz[y] += siz[x];
        p[x] = y;
        return true;
    }
    int size(int x) { return siz[find(x)];}
};
```



### ODT

```
struct Node{
	int l, r, val;
	Node(int a = -1, int b = -1, int c = 0){
		l = a, r = b, val = c;
	}
	bool operator < (const Node &a) const {
		return l < a.l;
	}
};
 
set<Node> st;
set<Node>::iterator split(int pos) // [l, r]分成[l, pos), [pos, r], 返回[pos, r]
{
	set<Node>::iterator it = st.lower_bound(Node(pos));
	if (it != st.end() && it->l == pos) return it;
	--it; Node tmp = *it; st.erase(it);
	st.insert(Node(tmp.l, pos - 1, tmp.val));
	return st.insert(Node(pos, tmp.r, tmp.val)).first; //first return iterator
}
 
void assign(int l, int r, int val)
{
	set<Node>::iterator itr = split(r + 1), itl = split(l);
    for (auto it = itl; it != itr; it = st.erase(it))
    {
        auto[L, R, Val] = (*it);
        add(R + 1, -d[Val]);
        add(L, d[Val]);
    }
 
    st.insert((Node){l, r, val});
    add(r + 1, d[val]);
    add(l, -d[val]);
}
```



### BIT

```cpp
struct BIT // 注意范围
{
    int N;
    vector<int> tr;
    BIT(int n) { N = n; tr.resize(n + 1); }
    int lowbit(int x) { return x & -x;}
    void add(int x, int c) { for (int i = x; i <= N; i += lowbit(i)) tr[i] += c; }
    int sum(int x)  // 返回前x个数的和
    {
        int res = 0;
        for (int i = x; i; i -= lowbit(i)) res += tr[i];
        return res;
    }
};
```



### 线段树

```cpp
struct Node
{
    int l, r;
    // TODO: 需要维护的信息和懒标记
    #define ls u << 1
    #define rs u << 1 | 1
}tr[N << 2];

void pushup(int u) 
{ 
    // TODO: 利用左右儿子信息维护当前节点的信息
}

void pushdown(int u)
{
    if(!tr[u].lazy) return;
    // TODO: 将懒标记下传
}

void build(int u, int l, int r)
{
    if (l == r) tr[u] = {l, r}; // TODO 赋值
    else
    {
        tr[u] = {l, r};
        int mid = l + r >> 1;
        build(u << 1, l, mid), build(u << 1 | 1, mid + 1, r);
        pushup(u);
    }
}

void update(int u, int l, int r, int d)
{
    if(tr[u].l > r || tr[u].r < l) return;
    if (tr[u].l >= l && tr[u].r <= r)
    {
        //TODO: 影响与懒标记
        return;
    }

    pushdown(u);
    update(u << 1, l, r, d), update(u << 1 | 1, l, r, d);
    pushup(u);
}

int query(int u, int l, int r)
{
    if(tr[u].l > r || tr[u].r < l) return 0;
    if (tr[u].l >= l && tr[u].r <= r) return ;  // TODO: 需要补充返回值

    pushdown(u);
    
    return query(u << 1, l, r) + query(u << 1 | 1, l, r); // 返回两个区间合并为1个区间的操作
}
```



### 主席树

```cpp
#include<bits/stdc++.h>

using namespace std;

const int N = 1e5 + 10;

int n, m;
int rt[N], idx;

struct node
{
    int l, r;
    int sum;
} tr[N << 5];

int ins(int u, int l, int r, int x, int k) // 二分值域，在值等于x处加k
{
    tr[ ++ idx] = tr[u];
    u = idx;
    
    tr[u].sum += k;
    if(l == r) return u;
    
    int mid = l + r >> 1;
    if(x <= mid) tr[u].l = ins(tr[u].l, l, mid, x, k);
    else tr[u].r = ins(tr[u].r, mid + 1, r, x, k);
    return u;
}

int query(int u, int v, int l ,int r, int k)
{
    if(l == r) return l;
    int l1 = tr[u].l, l2 = tr[v].l;
    int r1 = tr[u].r, r2 = tr[v].r;
    int mid = l + r >> 1;
    if(tr[l2].sum - tr[l1].sum >= k) return query(l1, l2, l, mid, k);
    return query(r1, r2, mid + 1, r, k - (tr[l2].sum - tr[l1].sum));
}

int main()
{
    cin >> n >> m;
    for (int i = 1; i <= n; i ++ )
    {
        int x;
        cin >> x;
        rt[i] = ins(rt[i - 1], -1e9, 1e9, x, 1);
    }
    
    while (m -- )
    {
        int l, r, k;
        cin >> l >> r >> k;
        cout << query(rt[l - 1], rt[r], -1e9, 1e9, k) << endl;
    }
    
    return 0;
}
```



### 可持久化01Trie树

```cpp
#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>

using namespace std;

const int N = 600010, M = N * 25;

int n, m;
int s[N];
int tr[M][2], max_id[M];
int root[N], idx;

void insert(int i, int k, int p, int q)
{
    if (k < 0)
    {
        max_id[q] = i;
        return;
    }
    int v = s[i] >> k & 1;
    if (p) tr[q][v ^ 1] = tr[p][v ^ 1];
    tr[q][v] = ++ idx;
    insert(i, k - 1, tr[p][v], tr[q][v]);
    max_id[q] = max(max_id[tr[q][0]], max_id[tr[q][1]]);
}

int query(int root, int C, int L)
{
    int p = root;
    for (int i = 23; i >= 0; i -- )
    {
        int v = C >> i & 1;
        if (max_id[tr[p][v ^ 1]] >= L) p = tr[p][v ^ 1];
        else p = tr[p][v];
    }

    return C ^ s[max_id[p]];
}

int main()
{
    scanf("%d%d", &n, &m);

    max_id[0] = -1;
    root[0] = ++ idx;
    insert(0, 23, 0, root[0]);

    for (int i = 1; i <= n; i ++ )
    {
        int x;
        scanf("%d", &x);
        s[i] = s[i - 1] ^ x;
        root[i] = ++ idx;
        insert(i, 23, root[i - 1], root[i]);
    }

    char op[2];
    int l, r, x;
    while (m -- )
    {
        scanf("%s", op);
        if (*op == 'A')
        {
            scanf("%d", &x);
            n ++ ;
            s[n] = s[n - 1] ^ x;
            root[n] = ++ idx;
            insert(n, 23, root[n - 1], root[n]);
        }
        else
        {
            scanf("%d%d%d", &l, &r, &x);
            printf("%d\n", query(root[r - 1], s[n] ^ x, l - 1));
        }
    }

    return 0;
}
```



### 线性基

```cpp
const int N = 1e5 + 10, M = N * 2, S = 62;

void insert(LL x)
{
    for (int i = S - 1; i >= 0 && x; i -- )
        if(x >> i & 1)
        {
            if(!p[i])
            {
                p[i] = x;
                return;
            }
            x ^= p[i];
        }
}

bool check(LL x)
{
    for (int i = S - 1; i >= 0; i -- )
    {
        if((x >> i & 1) == 0) continue;
        if(!p[i]) return 0;
        x ^= p[i];
    }
    return 1;
}

LL get_max(LL x)
{
    LL res = x;
    for (int i = S - 1; i >= 0; i -- ) res = max(res, res ^ p[i]);
    return res;
}
```



### 无旋Treap

```cpp

struct node
{
    int l, r;
    int val, sz, key;
    #define ls tr[u].l
    #define rs tr[u].r
} tr[N];

int rt, idx;

void pushup(int u) { tr[u].sz = tr[ls].sz + tr[rs].sz + 1; }

int add(int x) {
    tr[++ idx] = {0, 0, x, 1, rnd(1e9)};
    return idx;
}

int merge(int x, int y)
{
    if(!x || !y) return x + y;
    if(tr[x].key < tr[y].key) {
        tr[x].r = merge(tr[x].r, y);
        pushup(x);
        return x;
    } else {
        tr[y].l = merge(x, tr[y].l);
        pushup(y);
        return y;
    }
}

void split(int u, int val, int &x, int &y) // x为小于等于val的树, y为大于val的树
{
    if(!u) {
        x = y = 0;
        return;
    }

    // 该节点的val <= val, 因为右树的值都是大于等于tr[u].val的, 所以x要选tr[u].r
    if(tr[u].val <= val) x = u, split(rs, val, rs, y);
    else y = u, split(ls, val, x, ls);
    pushup(u);
}

int val2rk(int val) {
    int x, y;
    split(rt, val - 1, x, y);
    int res = tr[x].sz + 1;
    rt = merge(x, y);
    return res;
}

int rk2val(int rk){
    int u = rt;
    while (u)
    {
        if(tr[ls].sz + 1 == rk) return tr[u].val;
        else if(tr[ls].sz >= rk) u = ls;
        else rk -= tr[ls].sz + 1, u = rs;
    }
    return -1;
}

void insert(int val)
{
    int x, y;
    split(rt, val, x, y);
    rt = merge(merge(x, add(val)), y);
}

void remove(int val){
    int x, y, z;
    split(rt, val, x, z);
    split(x, val - 1, x, y);
    y = merge(tr[y].l, tr[y].r); // 去掉中间点
    rt = merge(merge(x, y), z);
}

int prv(int val) {
    int x, y;
    split(rt, val - 1, x, y);
    // 往右走
    int u = x;
    while (rs) u = rs;

    int res = tr[u].val;
    rt = merge(x, y);
    return res;
}

int nxt(int val) {
    int x, y;
    split(rt, val, x, y);
    int u = y;
    while (ls) u = ls;
    int res = tr[u].val;
    rt = merge(x, y);
    return res;
}
```



### 莫队

#### 普通莫队

```cpp
#include <iostream>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cmath>

using namespace std;

const int N = 50010, M = 200010, S = 1000010;

int n, m, len;
int w[N], ans[M];
struct Query
{
    int id, l, r;
}q[M];
int cnt[S];

int get(int x)
{
    return x / len;
}

bool cmp(const Query& a, const Query& b)
{
    int i = get(a.l), j = get(b.l);
    if (i ^ j) return a.l < b.l;
    else if(i & 1) return a.r < b.r;
    return a.r > b.r;
}

void add(int x, int& res)
{
    if (!cnt[x]) res ++ ;
    cnt[x] ++ ;
}

void del(int x, int& res)
{
    cnt[x] -- ;
    if (!cnt[x]) res -- ;
}

int main()
{
    scanf("%d", &n);
    for (int i = 1; i <= n; i ++ ) scanf("%d", &w[i]);
    scanf("%d", &m);
    len = sqrt((double)n * n / m);

    for (int i = 0; i < m; i ++ )
    {
        int l, r;
        scanf("%d%d", &l, &r);
        q[i] = {i, l, r};
    }
    sort(q, q + m, cmp);

    for (int k = 0, i = 0, j = 1, res = 0; k < m; k ++ )
    {
        int id = q[k].id, l = q[k].l, r = q[k].r;
        while (i < r) add(w[ ++ i], res);
        while (i > r) del(w[i -- ], res);
        while (j < l) del(w[j ++ ], res);
        while (j > l) add(w[ -- j], res);
        ans[id] = res;
    }

    for (int i = 0; i < m; i ++ ) printf("%d\n", ans[i]);
    return 0;
}
```



#### 回滚莫队

```cpp
#include <iostream>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

typedef long long LL;
const int N = 100010;

int n, m, len;
int w[N], cnt[N];
int ans[N];
struct Query
{
    int id, l, r;
}q[N];

vector<int> nums;

int get(int x)
{
    return x / len;
}

bool cmp(const Query& a, const Query& b)
{
    int i = get(a.l), j = get(b.l);
    if (i != j) return i < j;
    return a.r < b.r;
}

void add(int x, int& res)
{
    cnt[x] ++ ;
    res = max(res, cnt[x]);
}

int main()
{
    scanf("%d%d", &n, &m);
    len = sqrt(n);
    for (int i = 1; i <= n; i ++ ) scanf("%d", &w[i]), nums.push_back(w[i]);
    sort(nums.begin(), nums.end());
    nums.erase(unique(nums.begin(), nums.end()), nums.end());
    for (int i = 1; i <= n; i ++ ) w[i] = lower_bound(nums.begin(), nums.end(), w[i]) - nums.begin();

    for (int i = 0; i < m; i ++ )
    {
        int l, r;
        scanf("%d%d", &l, &r);
        q[i] = {i, l, r};
    }
    
    sort(q, q + m, cmp);
    
    for (int x = 0; x < m;)
    {
        int y = x;
        while (y < m && get(q[y].l) == get(q[x].l)) y ++ ;
        int right = get(q[x].l) * len + len - 1;

        // 暴力求块内的询问
        while (x < y && q[x].r <= right)
        {
            int res = 0;
            int id = q[x].id, l = q[x].l, r = q[x].r;
            for (int k = l; k <= r; k ++ ) add(w[k], res);
            ans[id] = res;
            for (int k = l; k <= r; k ++ ) cnt[w[k]] -- ;
            x ++ ;
        }

        // 求块外的询问
        int res = 0;
        int i = right, j = right + 1;
        while (x < y)
        {
            int id = q[x].id, l = q[x].l, r = q[x].r;
            while (i < r) add(w[ ++ i], res);
            int backup = res;
            while (j > l) add(w[ -- j], res);
            ans[id] = res;
            while (j < right + 1) cnt[w[j ++ ]] -- ;
            res = backup;
            x ++ ;
        }
        
        memset(cnt, 0, sizeof cnt);
    }

    for (int i = 0; i < m; i ++ ) printf("%d\n", -1 * ans[i]);
    return 0;
}

```



#### 带修莫队

```cpp
#include <iostream>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cmath>

using namespace std;

const int N = 10010, S = 1000010;

int n, m, mq, mc, len;
int w[N], cnt[S], ans[N];
struct Query
{
    int id, l, r, t;
}q[N];
struct Modify
{
    int p, c;
}c[N];

int get(int x)
{
    return x / len;
}

bool cmp(const Query& a, const Query& b)
{
    int al = get(a.l), ar = get(a.r);
    int bl = get(b.l), br = get(b.r);
    if (al != bl) return al < bl;
    if (ar != br) return ar < br;
    return a.t < b.t;
}

void add(int x, int& res)
{
    if (!cnt[x]) res ++ ;
    cnt[x] ++ ;
}

void del(int x, int& res)
{
    cnt[x] -- ;
    if (!cnt[x]) res -- ;
}

int main()
{
    scanf("%d%d", &n, &m);
    for (int i = 1; i <= n; i ++ ) scanf("%d", &w[i]);
    for (int i = 0; i < m; i ++ )
    {
        char op[2];
        int a, b;
        scanf("%s%d%d", op, &a, &b);
        if (*op == 'Q') mq ++, q[mq] = {mq, a, b, mc};
        else c[ ++ mc] = {a, b};
    }

    len = cbrt((double)n * mc) + 1;
    sort(q + 1, q + mq + 1, cmp);

    for (int i = 0, j = 1, t = 0, k = 1, res = 0; k <= mq; k ++ )
    {
        int id = q[k].id, l = q[k].l, r = q[k].r, tm = q[k].t;
        while (i < r) add(w[ ++ i], res);
        while (i > r) del(w[i -- ], res);
        while (j < l) del(w[j ++ ], res);
        while (j > l) add(w[ -- j], res);
        while (t < tm)
        {
            t ++ ;
            if (c[t].p >= j && c[t].p <= i)
            {
                del(w[c[t].p], res);
                add(c[t].c, res);
            }
            swap(w[c[t].p], c[t].c);
        }
        while (t > tm)
        {
            if (c[t].p >= j && c[t].p <= i)
            {
                del(w[c[t].p], res);
                add(c[t].c, res);
            }
            swap(w[c[t].p], c[t].c);
            t -- ;
        }
        ans[id] = res;
    }

    for (int i = 1; i <= mq; i ++ ) printf("%d\n", ans[i]);
    return 0;
}
```





## 图论

### 树上问题

#### 树的重心

```cpp
vector<int> g[N];

int siz[N];   // 这个节点的“大小”（所有子树上节点数 + 该节点）
int weight[N]; // 这个节点的“重量”
int centroid[2];  // 用于记录树的重心（存的是节点编号）

void GetCentroid(int u, int last, int sum) // u 表示当前节点
{
    if(last == 0) centroid[0] = centroid[1] = 0;

    siz[u] = 1;
    weight[u] = 0;
    for (auto &v : g[u])
    {
        if(v == last) continue;
        // v 表示这条有向边所通向的节点。
        GetCentroid(v, u, sum);
        siz[u] += siz[v];
        weight[u] = max(weight[u], siz[v]);
    }

    weight[u] = max(weight[u], sum - siz[u]);
    if (weight[u] <= sum / 2) centroid[centroid[0] != 0] = u;
}
```



#### 树链剖分

```cpp

int siz[N], dep[N], id[N], cnt;
int son[N], top[N], fa[N];
vector<int> g[N];

void dfs1(int u, int last)
{
    dep[u] = dep[last] + 1, fa[u] = last, siz[u] = 1;
    for (auto &v : g[u])
    {
        if(v == last) continue;
        dfs1(v, u);
        siz[u] += siz[v];
        if(siz[v] > siz[son[u]]) son[u] = v;
    }
}

void dfs2(int u, int topf) // 将树转化为序列
{
    id[u] = ++ cnt, top[u] = topf;
    if(!son[u]) return;
    dfs2(son[u], topf);
    for (auto &v : g[u])
    {
        if(v == fa[u] || v == son[u]) continue;
        dfs2(v, v);
    }
}

void change(int opt, int x, int y)
{
    if(opt == 1) // 以x为根的子树
    {
        update(1, id[x], id[x] + siz[x] - 1); // 区间
    }
    else if(opt == 2) // x, y的路径上的询问
    {
        while (top[x] != top[y])
        {
            if(depth[top[x]] < depth[top[y]]) swap(x, y);
            query(1, id[top[x]], id[x]); // 一条链
            x = fa[top[x]];
        }

        if(id[x] > id[y]) swap(x, y);
        query(1, id[x], id[y]);
    }
}
```



#### LCA & Kruskal重构树

```c++
int n, m;
struct edge
{
    int a, b, w;
}ed[M];

vector<int> e[N];//即为重构树

int d[N];
int p[N];
int q[N];
int depth[N];
int fa[N][22];

int find(int x) {
    return p[x] == x ? p[x] : p[x] = find(p[x]);
}

int kruskal()
{
    for (int i = 1; i <= 2 * n; i ++ ) p[i] = i;

    int cnt = n;
    for(int i = 1; i <= m; i ++ )
    {
        int a = find(ed[i].a), b = find(ed[i].b);
        if (a ^ b)
        {
            d[++ cnt] = ed[i].w;
            p[a] = p[b] = cnt;
            e[cnt].push_back(a);
            e[cnt].push_back(b);
        }
    }
    return cnt;
}

void bfs(int root)  // 预处理倍增数组
{
    memset(depth, 0x3f, sizeof depth);
    depth[0] = 0, depth[root] = 1;  // depth存储节点所在层数
    int hh = 0, tt = 0;
    q[0] = root;
    while (hh <= tt)
    {
        int t = q[hh ++ ];
        for (auto &j : e[t])
        {
            if (depth[j] > depth[t] + 1)
            {
                depth[j] = depth[t] + 1;
                q[ ++ tt] = j;
                fa[j][0] = t;  // j的第二次幂个父节点
                for (int k = 1; k <= 20; k ++ ) fa[j][k] = fa[fa[j][k - 1]][k - 1];
            }
        }
    }
}

int lca(int a, int b)  // 返回a和b的最近公共祖先
{
    if (depth[a] < depth[b]) swap(a, b);
    for (int k = 20; k >= 0; k -- )
        if (depth[fa[a][k]] >= depth[b]) a = fa[a][k];

    if (a == b) return a;
    for (int k = 20; k >= 0; k -- )
        if (fa[a][k] != fa[b][k])
        {
            a = fa[a][k];
            b = fa[b][k];
        }

    return fa[a][0];
}
```



#### 虚树

```cpp
// Virtual Tree
//
// Comprime uma arvore dado um conjunto S de vertices, de forma que
// o conjunto de vertices da arvore comprimida contenha S e seja
// minimal e fechado sobre a operacao de LCA
// Se |S| = k, a arvore comprimida tem O(k) vertices
//
// O(k log(k))

template<typename T> struct rmq {
	vector<T> v;
	int n; static const int b = 30;
	vector<int> mask, t;

	int op(int x, int y) { return v[x] < v[y] ? x : y; }
	int msb(int x) { return __builtin_clz(1)-__builtin_clz(x); }
	rmq() {}
	rmq(const vector<T>& v_) : v(v_), n(v.size()), mask(n), t(n) {
		for (int i = 0, at = 0; i < n; mask[i++] = at |= 1) {
			at = (at<<1)&((1<<b)-1);
			while (at and op(i, i-msb(at&-at)) == i) at ^= at&-at;
		}
		for (int i = 0; i < n/b; i++) t[i] = b*i+b-1-msb(mask[b*i+b-1]);
		for (int j = 1; (1<<j) <= n/b; j++) for (int i = 0; i+(1<<j) <= n/b; i++)
			t[n/b*j+i] = op(t[n/b*(j-1)+i], t[n/b*(j-1)+i+(1<<(j-1))]);
	}
	int small(int r, int sz = b) { return r-msb(mask[r]&((1<<sz)-1)); }
	T query(int l, int r) {
		if (r-l+1 <= b) return small(r, r-l+1);
		int ans = op(small(l+b-1), small(r));
		int x = l/b+1, y = r/b-1;
		if (x <= y) {
			int j = msb(y-x+1);
			ans = op(ans, op(t[n/b*j+x], t[n/b*j+y-(1<<j)+1]));
		}
		return ans;
	}
};

namespace lca {
	vector<int> g[N];
	int v[2*N], pos[N], dep[2*N];
	int t;
	rmq<int> RMQ;

	void dfs(int i, int d = 0, int p = -1) {
		v[t] = i, pos[i] = t, dep[t++] = d;
		for (int j : g[i]) if (j != p) {
			dfs(j, d+1, i);
			v[t] = i, dep[t++] = d;
		}
	}
	void build(int n, int root) {
		t = 0;
		dfs(root);
		RMQ = rmq<int>(vector<int>(dep, dep+2*n-1));
	}
	int lca(int a, int b) {
		a = pos[a], b = pos[b];
		return v[RMQ.query(min(a, b), max(a, b))];
	}
	int dist(int a, int b) {
		return dep[pos[a]] + dep[pos[b]] - 2*dep[pos[lca(a, b)]];
	}
}

vector<int> virt[N];

//warning lembrar de buildar o LCA antes
int build_virt(vector<int> v) {
	auto cmp = [&](int i, int j) { return lca::pos[i] < lca::pos[j]; };
	sort(v.begin(), v.end(), cmp);
	for (int i = v.size()-1; i; i--) v.push_back(lca::lca(v[i], v[i-1]));
	sort(v.begin(), v.end(), cmp);
	v.erase(unique(v.begin(), v.end()), v.end());
	for (int i : v) virt[i].clear();
	for (int i = 1; i < v.size(); i++) {
    //warning soh to colocando aresta descendo
		virt[lca::lca(v[i-1], v[i])].push_back(v[i]);
	}
	return v[0];
}

```



#### DSU ON TREE

```cpp
#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

typedef long long LL;

const int N = 2e5 + 10, M = 4e5 + 10;

int n, m;
int h[N], e[M], ne[M], idx;
int siz[N], son[N], cnt[N], w[N], del[N];
LL res[N], tmp, cn, mx;

void add(int a, int b)  // 添加一条边a->b
{
    e[idx] = b, ne[idx] = h[a], h[a] = idx ++ ;
}

void dfs1(int u, int last)
{
    siz[u] = 1;
    for (int i = h[u]; ~i; i = ne[i])
    {
        int j = e[i];
        if(j == last) continue;
        dfs1(j, u);
        siz[u] += siz[j];
        if(siz[son[u]] < siz[j]) son[u] = j;
    }
}

void sum(int u, int last, int Son)
{
    if(!cnt[w[u]]) del[++ cn] = w[u];

    cnt[w[u]] ++ ;
    if(cnt[w[u]] > mx)
    {
        mx = cnt[w[u]];
        tmp = w[u];
    }
    else if(cnt[w[u]] == mx) tmp += w[u];

    for (int i = h[u]; ~i; i = ne[i])
    {
        int j = e[i];
        if(j == last || j == Son) continue;
        sum(j, u, Son);
    }
}

void dfs2(int u, int last, int type) // type = 0不存
{
    for (int i = h[u]; ~i; i = ne[i]) // 第一次遍历轻儿子, 计算轻儿子的答案
    {
        int j = e[i];
        if(j == last || j == son[u]) continue;
        dfs2(j, u, 0);
    }

    if(son[u]) dfs2(son[u], u, 1); // 统计重儿子对 u 的答案的影响

    // 统计轻儿子的贡献
    sum(u, last, son[u]);
    res[u] = tmp;

    // 清空
    if(type == 0)
    {
        for (int i = 1; i <= cn; i ++ ) cnt[del[i]] = 0;
        cn = tmp = mx = 0;
    }
}

int main()
{
    ios::sync_with_stdio(0), cin.tie(0);
    memset(h, -1, sizeof h);
    cin >> n;
    for (int i = 1; i <= n; i ++ ) cin >> w[i];
    for (int i = 1; i < n; i ++ )
    {
        int a, b;
        cin >> a >> b;
        add(a, b), add(b, a);
    }

    dfs1(1, 0), dfs2(1, 0, 0);

    for (int i = 1; i <= n; i ++ ) cout << res[i] << ' ';
    return 0;
}
```



#### 树hash

**用法:** 通过两个重心为根跑出来的$hash$值+该树的大小来判断。

```cpp
int vis[M], primes[M], cnt;

void get_primes(int n)  // 线性筛质数
{
    for (int i = 2; i <= n; i ++ )
    {
        if (!vis[i]) primes[cnt ++ ] = i;
        for (int j = 0; i * primes[j] <= n; j ++ )
        {
            vis[primes[j] * i] = true;
            if (i % primes[j] == 0) break;
        }
    }
}

vector<int> g[N];

int n;
int siz[N];   // 这个节点的“大小”（所有子树上节点数 + 该节点）
int weight[N]; // 这个节点的“重量”
int centroid[2];  // 用于记录树的重心（存的是节点编号）

void GetCentroid(int u, int last) // u 表示当前节点
{
    if(last == 0) centroid[0] = centroid[1] = 0;

    siz[u] = 1;
    weight[u] = 0;
    for (auto &v : g[u])
    {
        if(v == last) continue;
        // v 表示这条有向边所通向的节点。
        GetCentroid(v, u);
        siz[u] += siz[v];
        weight[u] = max(weight[u], siz[v]);
    }

    weight[u] = max(weight[u], n - siz[u]);
    if (weight[u] <= n / 2) centroid[centroid[0] != 0] = u;
}

LL f[N];
void dfs(int u, int last)
{
    f[u] = 1, siz[u] = 1;
    for (auto &it : g[u])
    {
        if(it == last) continue;
        dfs(it, u);
        siz[u] += siz[it];
        f[u] = (f[u] + f[it] * primes[siz[it]]) % mod;
    }
}

void solve()
{
    get_primes(1000);

    int m;
    cin >> m;
    map<array<int, 3>, int> mp;
    for (int t = 1; t <= m; t ++ )
    {
        cin >> n;
        for (int i = 1; i <= n; i ++ ) g[i].clear();

        for (int i = 1; i <= n; i ++ )
        {
            int fa;
            cin >> fa;
            if(fa == 0) continue;
            g[fa].push_back(i);
            g[i].push_back(fa);
        }

        GetCentroid(1, 0);

        array<int, 3> tmp = {0, 0, n};
        dfs(centroid[0], 0);
        tmp[0] = f[centroid[0]];
        if(centroid[1])
        {
            dfs(centroid[1], 0);
            tmp[1] = f[centroid[1]];
        }

        if(tmp[0] > tmp[1]) swap(tmp[0], tmp[1]);

        if(!mp[tmp]) mp[tmp] = t;
        cout << mp[tmp] << '\n';
    }
}
```





### 网络流

#### Dinic

```cpp
#include<bits/stdc++.h>

using namespace std;

const int N = 10010, M = 200010, INF = 1e8; // M为两倍边的数量

int n, m, S, T;
int h[N], e[M], f[M], ne[M], idx;
int q[N], d[N], cur[N];

void add(int a, int b, int c)
{
    e[idx] = b, f[idx] = c, ne[idx] = h[a], h[a] = idx ++ ;
    e[idx] = a, f[idx] = 0, ne[idx] = h[b], h[b] = idx ++ ; //无向图可将0改成c
}

bool bfs()
{
    int hh = 0, tt = 0;
    memset(d, -1, sizeof d);
    q[0] = S, d[S] = 0, cur[S] = h[S];
    while (hh <= tt)
    {
        int t = q[hh ++ ];
        for (int i = h[t]; ~i; i = ne[i])
        {
            int ver = e[i];
            if (d[ver] == -1 && f[i])
            {
                d[ver] = d[t] + 1;
                cur[ver] = h[ver];
                if (ver == T)  return true;
                q[ ++ tt] = ver;
            }
        }
    }
    return false;
}

int find(int u, int limit)
{
    if (u == T) return limit;
    int flow = 0;
    for (int i = cur[u]; ~i && flow < limit; i = ne[i])
    {
        cur[u] = i;  // 当前弧优化
        int ver = e[i];
        if (d[ver] == d[u] + 1 && f[i])
        {
            int t = find(ver, min(f[i], limit - flow));
            if (!t) d[ver] = -1;
            f[i] -= t, f[i ^ 1] += t, flow += t;
        }
    }
    return flow;
}

int dinic()
{
    int r = 0, flow;
    while (bfs()) while (flow = find(S, INF)) r += flow;
    return r;
}

int main()
{
    scanf("%d%d%d%d", &n, &m, &S, &T);
    memset(h, -1, sizeof h);
    while (m -- )
    {
        int a, b, c;
        scanf("%d%d%d", &a, &b, &c);
        add(a, b, c);
    }

    printf("%d\n", dinic());

    return 0;
}
```

#### 费用流

```
/*
我们知道，费用流的的一种常见解法是用 SPFA 找到从 S 到 T 的以费用为权的最短路（最小费用），将这条路上的流量乘以最小费用累加到答案中，直到 SPFA 最后找不到这个最短路。
如果流量比较多的话可以考虑一边加边一边跑
因为spfa一次增广出一条路，其时间复杂度和边数有关
*/
#include <bits/stdc++.h>
#define maxx(a, b, c) max(max(a, b), c)
#define minn(a, b, c) min(min(a, b), c)
#define between(x, l, r) (x >= l && x <= r)
#define point(a, b, m) ((a - 1) * m + b)
#define y second
#define x first
 
using namespace std;
 
const int N = 10010, M = 1e6, INF = 1e9; //M为两倍边的上限
 
int n, m, S, T;
int h[N], e[M], f[M], w[M], ne[M], idx;
int q[N], d[N], pre[N], incf[N];
bool st[N];
 
void add(int a, int b, int c, int d)
{
    e[idx] = b, f[idx] = c, w[idx] = d, ne[idx] = h[a], h[a] = idx ++ ;
    e[idx] = a, f[idx] = 0, w[idx] = -d, ne[idx] = h[b], h[b] = idx ++ ;
}
 
bool spfa()
{
    int hh = 0, tt = 1;
    memset(d, 0x3f, sizeof d);
    memset(incf, 0, sizeof incf);
    q[0] = S, d[S] = 0, incf[S] = INF;
    while (hh != tt)
    {
        int t = q[hh ++ ];
        if (hh == N) hh = 0;
        st[t] = false;
 
        for (int i = h[t]; ~i; i = ne[i])
        {
            int ver = e[i];
            if (f[i] && d[ver] > d[t] + w[i])
            {
                d[ver] = d[t] + w[i];
                pre[ver] = i;
                incf[ver] = min(f[i], incf[t]);
                if (!st[ver])
                {
                    q[tt ++ ] = ver;
                    if (tt == N) tt = 0;
                    st[ver] = true;
                }
            }
        }
    }

    return incf[T] > 0;
}
 
void EK(int& flow, int& cost)
{
    flow = cost = 0;
    while (spfa())
    {
        int t = incf[T];
        flow += t, cost += t * d[T];
        for (int i = T; i != S; i = e[pre[i] ^ 1])
        {
            f[pre[i]] -= t;
            f[pre[i] ^ 1] += t;
        }
    }
}

void print() // 输出有流量的边
{
    for (int i = 0; i < idx; i += 2)
        if(f[i ^ 1]) cout << e[i ^ 1] << ' ' << e[i] << '\n';
}

int main()
{
    ios::sync_with_stdio(0), cin.tie(0);
    // cout.precision(2), cout.setf(ios::fixed); 控制精度
    
    S = 0, T = N - 1;
    memset(h, -1, sizeof h);

    // 建边
    
    //输出答案
    int flow, cost;
    EK(flow, cost);
    
    cout << cost << '\n';
    return 0;
}


------------atcoder板子--------------------------
#include <algorithm>
#include <cassert>
#include <limits>
#include <queue>
#include <vector>
#include <cstdio>

namespace atcoder {

template <class Cap, class Cost> struct mcf_graph {
  public:
    mcf_graph() {}
    mcf_graph(int n) : _n(n), g(n) {}

    int add_edge(int from, int to, Cap cap, Cost cost) {
        assert(0 <= from && from < _n);
        assert(0 <= to && to < _n);
        int m = int(pos.size());
        pos.push_back({from, int(g[from].size())});
        int from_id = int(g[from].size());
        int to_id = int(g[to].size());
        if (from == to) to_id++;
        g[from].push_back(_edge{to, to_id, cap, cost});
        g[to].push_back(_edge{from, from_id, 0, -cost});
        return m;
    }

    struct edge {
        int from, to;
        Cap cap, flow;
        Cost cost;
    };

    edge get_edge(int i) {
        int m = int(pos.size());
        assert(0 <= i && i < m);
        auto _e = g[pos[i].first][pos[i].second];
        auto _re = g[_e.to][_e.rev];
        return edge{
            pos[i].first, _e.to, _e.cap + _re.cap, _re.cap, _e.cost,
        };
    }
    std::vector<edge> edges() {
        int m = int(pos.size());
        std::vector<edge> result(m);
        for (int i = 0; i < m; i++) {
            result[i] = get_edge(i);
        }
        return result;
    }

    std::pair<Cap, Cost> flow(int s, int t) {
        return flow(s, t, std::numeric_limits<Cap>::max());
    }
    std::pair<Cap, Cost> flow(int s, int t, Cap flow_limit) {
        return slope(s, t, flow_limit).back();
    }
    std::vector<std::pair<Cap, Cost>> slope(int s, int t) {
        return slope(s, t, std::numeric_limits<Cap>::max());
    }
    std::vector<std::pair<Cap, Cost>> slope(int s, int t, Cap flow_limit) {
        assert(0 <= s && s < _n);
        assert(0 <= t && t < _n);
        assert(s != t);
        std::vector<Cost> dual(_n, 0), dist(_n);
        std::vector<int> pv(_n), pe(_n);
        std::vector<bool> vis(_n);
        auto dual_ref = [&]() {
            std::fill(dist.begin(), dist.end(),
                      std::numeric_limits<Cost>::max());
            std::fill(pv.begin(), pv.end(), -1);
            std::fill(pe.begin(), pe.end(), -1);
            std::fill(vis.begin(), vis.end(), false);
            struct Q {
                Cost key;
                int to;
                bool operator<(Q r) const { return key > r.key; }
            };
            std::priority_queue<Q> que;
            dist[s] = 0;
            que.push(Q{0, s});
            while (!que.empty()) {
                int v = que.top().to;
                que.pop();
                if (vis[v]) continue;
                vis[v] = true;
                if (v == t) break;
                for (int i = 0; i < int(g[v].size()); i++) {
                    auto e = g[v][i];
                    if (vis[e.to] || !e.cap) continue;
                    Cost cost = e.cost - dual[e.to] + dual[v];
                    if (dist[e.to] - dist[v] > cost) {
                        dist[e.to] = dist[v] + cost;
                        pv[e.to] = v;
                        pe[e.to] = i;
                        que.push(Q{dist[e.to], e.to});
                    }
                }
            }
            if (!vis[t]) {
                return false;
            }

            for (int v = 0; v < _n; v++) {
                if (!vis[v]) continue;
                dual[v] -= dist[t] - dist[v];
            }
            return true;
        };
        Cap flow = 0;
        Cost cost = 0, prev_cost_per_flow = -1;
        std::vector<std::pair<Cap, Cost>> result;
        result.push_back({flow, cost});
        while (flow < flow_limit) {
            if (!dual_ref()) break;
            Cap c = flow_limit - flow;
            for (int v = t; v != s; v = pv[v]) {
                c = std::min(c, g[pv[v]][pe[v]].cap);
            }
            for (int v = t; v != s; v = pv[v]) {
                auto& e = g[pv[v]][pe[v]];
                e.cap -= c;
                g[v][e.rev].cap += c;
            }
            Cost d = -dual[s];
            flow += c;
            cost += c * d;
            if (prev_cost_per_flow == d) {
                result.pop_back();
            }
            result.push_back({flow, cost});
            prev_cost_per_flow = d;
        }
        return result;
    }

  private:
    int _n;

    struct _edge {
        int to, rev;
        Cap cap;
        Cost cost;
    };

    std::vector<std::pair<int, int>> pos;
    std::vector<std::vector<_edge>> g;
};

}

int main(){
  int n,m,s,t;
  scanf("%d%d%d%d",&n,&m,&s,&t);
  s--,t--;
  atcoder::mcf_graph<long long,long long> g(n);
  for(int i=0;i<m;++i){
    int x,y,cap,cost;
    scanf("%d%d%d%d",&x,&y,&cap,&cost);
    x--,y--;
    g.add_edge(x,y,cap,cost);
  }
  auto an=g.flow(s,t);
  printf("%lld %lld\n",an.first,an.second);
  return 0;
}
```



### dijkstra

```cpp
vector<PII> g[N];
int dist[N], st[N];

void dij()  // 求1号点到n号点的最短路距离
{
    memset(st, 0, sizeof st);
    memset(dist, 0x3f, sizeof dist);
    
    dist[1] = 0;
    priority_queue<PII, vector<PII>, greater<PII>> heap;
    heap.push({0, 1});

    while (heap.size())
    {
        auto t = heap.top(); heap.pop();
        
        int u = t.second;
        if (st[u]) continue;
        st[u] = true;

        for (auto &[v, w] : g[u])
        {
            if (dist[v] > dist[u] + w)
            {
                dist[v] = dist[u] + w;
                heap.push({dist[v], v});
            }
        }
    }
}
```



### topsort

技巧 : $i → p[i]$ 连一条边的图可以使用 $topsort$ ，剩下的都是环。

```cpp
int n, m;
vector<int> g[N];
int d[N], q[N];

bool topsort()
{
    int hh = 0, tt = -1;

    // d[i] 存储点i的入度
    for (int i = 1; i <= n; i ++ )
        if (!d[i]) q[ ++ tt] = i;

    while (hh <= tt)
    {
        int t = q[hh ++ ];

        for (auto &v : g[t])
        {
            if (-- d[v] == 0) q[ ++ tt] = v;
        }
    }
    
    return tt == n - 1;
}
```



### 最大团

```cpp
/*
最大团 = 补图G的最大独立集数
———>最大独立集数 = 补图G'最大团
*/

//最大团模板
bool g[N][N];//a为图的邻接表(从1开始)
int ans, cnt[N], group[N], n, m, vis[N]; //ans表示最大团, cnt[N]表示当前最大团的节点数, group[N]用以寻找一个最大团集合
bool dfs(int u, int pos) //u为当从前顶点开始深搜，pos为深搜深度（即当前深搜树所在第几层的位置）
{
    int i, j;
    for(i = u + 1; i <= n; i ++ ) //按递增顺序枚举顶点
	{
        if(cnt[i] + pos <= ans) return 0; //剪枝
        if(g[u][i])
		{
            // 与目前团中元素比较, 取 Non-N(i)
            for(j = 0; j < pos; j ++ )
                if(!g[i][vis[j]]) break;

            if(j == pos)
			{     // 若为空，则皆与 i 相邻，则此时将i加入到 最大团中
                vis[pos] = i;//深搜层次也就是最大团的顶点数目，vis[pos] = i表示当前第pos小的最大团元素为i（因为是按增顺序枚举顶点 ）
                if(dfs(i, pos + 1)) return 1;
            }
        }
    }

    if(pos > ans)
	{
        for(i = 0; i < pos; i++ ) group[i] = vis[i]; // 更新最大团元素
        ans = pos;
        return 1;
    }

    return 0;
}

void maxclique()//求最大团
{
    ans = -1;
    for(int i = n; i > 0; i -- )
    {
        vis[0] = i;
        dfs(i, 1);
        cnt[i] = ans;
    }
}

```



### Prufer 序列（有标号无根树计数）

考虑这样一种序列 : 长度为 $n-2$ ，值域在 $[1,n]$ 以内。

奇妙的是，这类序列和带标号无根树之间存在一一对应(双射)。这对计数大有帮助。

- 树 $\rightarrow$ 序列 :

  每次移去所有叶子节点中标号最小的顶点和相连的边，并把与它相邻的点的编号加入 `prufer` 序列中。重复以上步骤直到原图仅剩 $2$ 个顶点。

  我们能够发现，节点 $u$ 会在序列中出现 $(\text{度数}-1)$ 次，这是个重要的性质。

- 序列 $\rightarrow$ 树 :

  建立集合 $G$ 含有节点 $\{1...n\}$。找出集合中最小的，未在 `prufer` 序列中出现过的数，将该点与 `prufer` 序列中首项连一条边，并将该点和序列首项删除。

  重复操作 $n-2$ 次，最后将集合中剩余的两个点之间连边。

如何构造并不是重点，我们只需要记住**双射**的性质，和有关**出现次数**和**度数**的性质就好了。这样，就能把有标号树计数问题转化成了序列计数。

```cpp
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

const int N = 100010;

int n, m;
int f[N], d[N], p[N];

void tree2prufer()
{
    for (int i = 1; i < n; i ++ )
    {
        scanf("%d", &f[i]);
        d[f[i]] ++ ;
    }

    for (int i = 0, j = 1; i < n - 2; j ++ )
    {
        while (d[j]) j ++ ;
        p[i ++ ] = f[j];
        while (i < n - 2 && -- d[p[i - 1]] == 0 && p[i - 1] < j) p[i ++ ] = f[p[i - 1]];
    }

    for (int i = 0; i < n - 2; i ++ ) printf("%d ", p[i]);
}

void prufer2tree()
{
    for (int i = 1; i <= n - 2; i ++ )
    {
        scanf("%d", &p[i]);
        d[p[i]] ++ ;
    }
    p[n - 1] = n;

    for (int i = 1, j = 1; i < n; i ++, j ++ )
    {
        while (d[j]) j ++ ;
        f[j] = p[i];
        while (i < n - 1 && -- d[p[i]] == 0 && p[i] < j) f[p[i]] = p[i + 1], i ++ ;
    }

    for (int i = 1; i <= n - 1; i ++ ) printf("%d ", f[i]);
}

int main()
{
    scanf("%d%d", &n, &m);
    if (m == 1) tree2prufer();
    else prufer2tree();
    return 0;
}
```



### 求补图连通块

核心：找出原图上的度数最小的点 ( $d[u] ≤ \frac{m}{n}$ )，然后找出补图上所有与这个点有边相邻的点，将它们合并成一个连通块。剩下 $d[u]$ 个点暴力即可。

时间复杂度 : $O(n * \frac{m}{n})$

```c++
struct DSU {
    vector<int> p, siz;
    DSU(int n) : p(n + 1), siz(n + 1, 1) {iota(all(p), 0);}
    int find(int x) {
        while(p[x] != x) x = p[x] = p[p[x]];
        return p[x];
    }
    bool same(int x, int y) {return find(x) == find(y);}
    bool merge(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) return false;
        siz[y] += siz[x];
        p[x] = y;
        return true;
    }
};

void solve()
{
    int n, m;
    cin >> n >> m;

    LL xr = 0;
	vector<vector<int>> g(n + 1);
    for (int i = 1; i <= m; i ++ )
    {
        int u, v;
        cin >> u >> v;
		g[u].eb(v);
		g[v].eb(u);
    }

    DSU d(n);
	int r = 1;
	for (int i = 1; i <= n; i ++ )
		if(g[i].size() < g[r].size())
			r = i;

	vector<int> vis(n + 1);
	for (auto &it : g[r]) vis[it] = 1;
	for (int i = 1; i <= n; i ++ )
		if(!vis[i]) d.merge(i, r);

	for (auto &u : g[r])
	{
		vis.assign(n + 1, 0);
		for (auto &v : g[u]) vis[v] = 1;

		for (int i = 1; i <= n; i ++ )
			if(!vis[i]) d.merge(u, i); // 补图中没有
	}

	vis.assign(n + 1, 0);

	vector<int> res;
	for(int i = 1; i <= n; i ++ )
	{
		int p = d.find(i);
		if(vis[p]) continue;
		vis[p] = 1;
		res.eb(d.siz[p]);
	}

	sort(all(res));
	cout << res.size() << '\n';
	for (auto &it : res) cout << it << ' '; // 每个连通块的大小
	cout << '\n';
}

```



## 字符串

### 字符串hash

**用法: ** 先init()随机模数，然后再操作。

```c++
namespace Hashbase
{
    int mod[] = {1000000009, 1004535809, 469762049, 917120411, 515880193};
    int p1[N], p2[N];
    int ha1, ha2, bas1, bas2;
    void init(int a = -1, int b = -1, int c = -1, int d = -1)
    {
        int z1 = rnd(5);
        int z2 = (z1 + rnd(4) + 1) % 5;
        if (a < 0) a = mod[z1];
        if (b < 0) b = mod[z2];
        if (c < 0) c = rnd(114514) + 23333;
        if (d < 0) d = rnd(1919810) + 23333;
        ha1 = a, ha2 = b, bas1 = c, bas2 = d;

        for (int i = 1; i < N; i ++ )
        {
            p1[0] = p2[0] = 1;
            p1[i] = (LL)p1[i - 1] * bas1 % ha1;
            p2[i] = (LL)p2[i - 1] * bas2 % ha2;
        }
    }

    struct StringDoubleHash
    {
        vector<int> h1, h2;
        template<class Tp>
        StringDoubleHash(int n, Tp &s) : h1(n + 1), h2(n + 1)
        {
            for (int i = 1; i <= n; i ++ )
            {
                h1[i] = ((LL)h1[i - 1] * bas1 + s[i]) % ha1;
                h2[i] = ((LL)h2[i - 1] * bas2 + s[i]) % ha2;
            }
        }

        PII get(int l, int r)
        {
            int res1 = (h1[r] - (LL)h1[l - 1] * p1[r - l + 1] % ha1 + ha1) % ha1;
            int res2 = (h2[r] - (LL)h2[l - 1] * p2[r - l + 1] % ha2 + ha2) % ha2;
            return mk(res1, res2);
        }

        template<class Tp>
        PII calc_Hash(int n, Tp &s) // s[0, n)
        {
            int res1 = 0, res2 = 0;
            for (int i = 0; i < n; i ++ )
            {
                res1 = ((LL)res1 * bas1 + s[i]) % ha1;
                res2 = ((LL)res2 * bas2 + s[i]) % ha2;
            }
            return mk(res1, res2);
        }

        PII merge(PII hs1, PII hs2, int len) //s1+s2的哈希, len为hs2的长度, s1+s2最好小于
        {
            int res1 = hs1.x, res2 = hs1.y;
            res1 = ((LL)res1 * p1[len] + hs2.x) % ha1;
            res2 = ((LL)res2 * p2[len] + hs2.y) % ha2;
            return mk(res1, res2);
        }
    };
};

using namespace Hashbase;
```



### SAM

```cpp
struct SuffixAutomation
{
	int last, cnt;
	// map<int, int> ch[N << 1];
	//如果是字符串, 可以空间换时间用以下注释代码
    int ch[N << 1][26];
    int fa[N << 1], len[N << 1], pos[N << 1];
    int sz[N << 1], a[N << 1], c[N << 1];
    LL f[N << 1];

    // len[x]代表集合right[x]最长的公共部分
    // topsort: a为top序，c为入度
    // sz[x]代表集合right[x]的个数
    // fa[x]表示可以包含节点x的right集合，且集合元素最小的点
    // len越大, right的集合越小, 所以一直向上, 可以匹配的越宽泛

	void init()
	{
	    last = cnt = 1;
        // ch[1].clear();
        memset(ch[1], 0, sizeof ch[1]);
	    fa[1] = len[1] = 0;
        a[1] = c[1] = 0;
    }

    int inline newnode(int idx)
    {
        ++cnt;
        // ch[cnt].clear();
        memset(ch[cnt], 0, sizeof ch[cnt]);
        fa[cnt] = len[cnt] = sz[cnt] = a[cnt] = c[cnt] = 0;
        pos[cnt] = idx;
        return cnt;
    }

	void ins(int c)
	{
		int p = last, np = newnode(pos[last] + 1);
		last = np, len[np] = len[p] + 1;
		for(; p && !ch[p][c]; p = fa[p]) ch[p][c] = np;
		if(!p) fa[np] = 1;
		else
		{
			int q = ch[p][c];
			if(len[p] + 1 == len[q]) fa[np] = q;
			else
			{
				int nq = newnode(pos[p] + 1);
				len[nq] = len[p] + 1;
                // ch[nq] = ch[q];
                memcpy(ch[nq], ch[q], sizeof ch[q]);
				fa[nq] = fa[q], fa[q] = fa[np] = nq;
				for(; ch[p][c] == q; p = fa[p]) ch[p][c] = nq;
			}
		}
        sz[np] = 1;
	}

    void Toposort(int T) // T=0则不同位置的相同子串算作一个, 否则算多个
    {
        long long ans = 0;
        for(int i = 1; i <= cnt; i ++ ) c[len[i]] ++;
        for(int i = 1; i <= cnt; i ++ ) c[i] += c[i - 1];
        for(int i = 1; i <= cnt; i ++ ) a[c[len[i]] --] = i;
        for(int i = cnt; i; i -- )
        {
            int p = a[i];
            sz[fa[p]] += sz[p];
        }

        for (int i = 1; i <= cnt; i ++ )
            f[i] = (T == 0 ? sz[i] = 1 : sz[i]);

        sz[1] = f[1] = 0;
        for (int i = cnt; i; i -- )
        {
            int p = a[i];
            for (int j = 0; j < 26; j ++ )
            {
                if(!ch[p][j]) continue;
                f[p] += f[ch[p][j]];
            }
        }
    }

    void print(int u, int k) // 第k小的子串
    {
        if(k <= sz[u]) return;
        k -= sz[u];

        for (int j = 0; j < 26; j ++ )
        {
            int ne = ch[u][j];
            if(!ne) continue;

            if(k > f[ne]) k -= f[ne];
            else
            {
                cout << char(j + 'a');
                print(ne, k);
                return;
            }
        }
    }
} sam;

```



## 计算几何

### 二维全家桶

```cpp
using _T= double;  //全局数据类型，可修改为 long long 等

constexpr _T eps=1e-8;
constexpr long double PI=3.1415926535897932384l;

// 点与向量
template<typename T> struct point
{
    T x,y;
    bool operator==(const point &a) const {return (abs(x-a.x)<=eps && abs(y-a.y)<=eps);}
    bool operator<(const point &a) const {if (abs(x-a.x)<=eps) return y<a.y-eps; return x<a.x-eps;}
    bool operator>(const point &a) const {return !(*this<a || *this==a);}
    point operator+(const point &a) const {return {x+a.x,y+a.y};}
    point operator-(const point &a) const {return {x-a.x,y-a.y};}
    point operator-() const {return {-x,-y};}
    point operator*(const T k) const {return {k*x,k*y};}
    point operator/(const T k) const {return {x/k,y/k};}
    T operator*(const point &a) const {return x*a.x+y*a.y;}  // 点积
    T operator^(const point &a) const {return x*a.y-y*a.x;}  // 叉积，注意优先级
    int toleft(const point &a) const {const auto t=(*this)^a; return (t>eps)-(t<-eps);}  // to-left 测试
    T len2() const {return (*this)*(*this);}  // 向量长度的平方
    T dis2(const point &a) const {return (a-(*this)).len2();}  // 两点距离的平方

    // 涉及浮点数
    long double len() const {return sqrtl(len2());}  // 向量长度
    long double dis(const point &a) const {return sqrtl(dis2(a));}  // 两点距离
    long double ang(const point &a) const {return acosl(max(-1.0l,min(1.0l,((*this)*a)/(len()*a.len()))));}  // 向量夹角
    point rot(const long double rad) const {return {x*cos(rad)-y*sin(rad),x*sin(rad)+y*cos(rad)};}  // 逆时针旋转（给定角度）
    point rot(const long double cosr,const long double sinr) const {return {x*cosr-y*sinr,x*sinr+y*cosr};}  // 逆时针旋转（给定角度的正弦与余弦）
};

using Point=point<_T>;

// 极角排序
struct argcmp
{
    bool operator()(const Point &a,const Point &b) const
    {
        const auto quad=[](const Point &a)
        {
            if (a.y<-eps) return 1;
            if (a.y>eps) return 4;
            if (a.x<-eps) return 5;
            if (a.x>eps) return 3;
            return 2;
        };
        const int qa=quad(a),qb=quad(b);
        if (qa!=qb) return qa<qb;
        const auto t=a^b;
        // if (abs(t)<=eps) return a*a<b*b-eps;  // 不同长度的向量需要分开
        return t>eps;
    }
};

// 直线
template<typename T> struct line
{
    point<T> p,v;  // p 为直线上一点，v 为方向向量

    bool operator==(const line &a) const {return v.toleft(a.v)==0 && v.toleft(p-a.p)==0;}
    int toleft(const point<T> &a) const {return v.toleft(a-p);}  // to-left 测试
    bool operator<(const line &a) const  // 半平面交算法定义的排序
    {
        if (abs(v^a.v)<=eps && v*a.v>=-eps) return toleft(a.p)==-1;
        return argcmp()(v,a.v);
    }

    // 涉及浮点数
    point<T> inter(const line &a) const {return p+v*((a.v^(p-a.p))/(v^a.v));}  // 直线交点
    long double dis(const point<T> &a) const {return abs(v^(a-p))/v.len();}  // 点到直线距离
    point<T> proj(const point<T> &a) const {return p+v*((v*(a-p))/(v*v));}  // 点在直线上的投影
};

using Line=line<_T>;

//线段
template<typename T> struct segment
{
    point<T> a,b;

    bool operator< (const segment &s) const {return make_pair(a,b)<make_pair(s.a,s.b);}

    // 判定性函数建议在整数域使用

    // 判断点是否在线段上
    // -1 点在线段端点 | 0 点不在线段上 | 1 点严格在线段上
    int is_on(const point<T> &p) const
    {
        if (p==a || p==b) return -1;
        return (p-a).toleft(p-b)==0 && (p-a)*(p-b)<-eps;
    }

    // 判断线段直线是否相交
    // -1 直线经过线段端点 | 0 线段和直线不相交 | 1 线段和直线严格相交
    int is_inter(const line<T> &l) const
    {
        if (l.toleft(a)==0 || l.toleft(b)==0) return -1;
        return l.toleft(a)!=l.toleft(b);
    }

    // 判断两线段是否相交
    // -1 在某一线段端点处相交 | 0 两线段不相交 | 1 两线段严格相交
    int is_inter(const segment<T> &s) const
    {
        if (is_on(s.a) || is_on(s.b) || s.is_on(a) || s.is_on(b)) return -1;
        const line<T> l{a,b-a},ls{s.a,s.b-s.a};
        return l.toleft(s.a)*l.toleft(s.b)==-1 && ls.toleft(a)*ls.toleft(b)==-1;
    }

    // 点到线段距离
    long double dis(const point<T> &p) const
    {
        if ((p-a)*(b-a)<-eps || (p-b)*(a-b)<-eps) return min(p.dis(a),p.dis(b));
        const line<T> l{a,b-a};
        return l.dis(p);
    }

    // 两线段间距离
    long double dis(const segment<T> &s) const
    {
        if (is_inter(s)) return 0;
        return min({dis(s.a),dis(s.b),s.dis(a),s.dis(b)});
    }
};

using Segment=segment<_T>;

// 多边形
template<typename T> struct polygon
{
    vector<point<T>> p;  // 以逆时针顺序存储

    size_t nxt(const size_t i) const {return i==p.size()-1?0:i+1;}
    size_t pre(const size_t i) const {return i==0?p.size()-1:i-1;}

    // 回转数
    // 返回值第一项表示点是否在多边形边上
    // 对于狭义多边形，回转数为 0 表示点在多边形外，否则点在多边形内
    pair<bool,int> winding(const point<T> &a) const
    {
        int cnt=0;
        for (size_t i=0;i<p.size();i++)
        {
            const point<T> u=p[i],v=p[nxt(i)];
            if (abs((a-u)^(a-v))<=eps && (a-u)*(a-v)<=eps) return {true,0};
            if (abs(u.y-v.y)<=eps) continue;
            const Line uv={u,v-u};
            if (u.y<v.y-eps && uv.toleft(a)<=0) continue;
            if (u.y>v.y+eps && uv.toleft(a)>=0) continue;
            if (u.y<a.y-eps && v.y>=a.y-eps) cnt++;
            if (u.y>=a.y-eps && v.y<a.y-eps) cnt--;
        }
        return {false,cnt};
    }

    // 多边形面积的两倍
    // 可用于判断点的存储顺序是顺时针或逆时针
    T area() const
    {
        T sum=0;
        for (size_t i=0;i<p.size();i++) sum+=p[i]^p[nxt(i)];
        return sum;
    }

    // 多边形的周长
    long double circ() const
    {
        long double sum=0;
        for (size_t i=0;i<p.size();i++) sum+=p[i].dis(p[nxt(i)]);
        return sum;
    }
};

using Polygon=polygon<_T>;

//凸多边形
template<typename T> struct convex: polygon<T>
{
    // 闵可夫斯基和
    convex operator+(const convex &c) const
    {
        const auto &p=this->p;
        vector<Segment> e1(p.size()),e2(c.p.size()),edge(p.size()+c.p.size());
        vector<point<T>> res; res.reserve(p.size()+c.p.size());
        const auto cmp=[](const Segment &u,const Segment &v) {return argcmp()(u.b-u.a,v.b-v.a);};
        for (size_t i=0;i<p.size();i++) e1[i]={p[i],p[this->nxt(i)]};
        for (size_t i=0;i<c.p.size();i++) e2[i]={c.p[i],c.p[c.nxt(i)]};
        rotate(e1.begin(),min_element(e1.begin(),e1.end(),cmp),e1.end());
        rotate(e2.begin(),min_element(e2.begin(),e2.end(),cmp),e2.end());
        merge(e1.begin(),e1.end(),e2.begin(),e2.end(),edge.begin(),cmp);
        const auto check=[](const vector<point<T>> &res,const point<T> &u)
        {
            const auto back1=res.back(),back2=*prev(res.end(),2);
            return (back1-back2).toleft(u-back1)==0 && (back1-back2)*(u-back1)>=-eps;
        };
        auto u=e1[0].a+e2[0].a;
        for (const auto &v:edge)
        {
            while (res.size()>1 && check(res,u)) res.pop_back();
            res.push_back(u);
            u=u+v.b-v.a;
        }
        if (res.size()>1 && check(res,res[0])) res.pop_back();
        return {res};
    }

    // 旋转卡壳
    // func 为更新答案的函数，可以根据题目调整位置
    template<typename F> void rotcaliper(const F &func) const
    {
        const auto &p=this->p;
        const auto area=[](const point<T> &u,const point<T> &v,const point<T> &w){return (w-u)^(w-v);};
        for (size_t i=0,j=1;i<p.size();i++)
        {
            const auto nxti=this->nxt(i);
            func(p[i],p[nxti],p[j]);
            while (area(p[this->nxt(j)],p[i],p[nxti])>=area(p[j],p[i],p[nxti]))
            {
                j=this->nxt(j);
                func(p[i],p[nxti],p[j]);
            }
        }
    }

    // 凸多边形的直径的平方
    T diameter2() const
    {
        const auto &p=this->p;
        if (p.size()==1) return 0;
        if (p.size()==2) return p[0].dis2(p[1]);
        T ans=0;
        auto func=[&](const point<T> &u,const point<T> &v,const point<T> &w){ans=max({ans,w.dis2(u),w.dis2(v)});};
        rotcaliper(func);
        return ans;
    }

    // 判断点是否在凸多边形内
    // 复杂度 O(logn)
    // -1 点在多边形边上 | 0 点在多边形外 | 1 点在多边形内
    int is_in(const point<T> &a) const
    {
        const auto &p=this->p;
        if (p.size()==1) return a==p[0]?-1:0;
        if (p.size()==2) return segment<T>{p[0],p[1]}.is_on(a)?-1:0;
        if (a==p[0]) return -1;
        if ((p[1]-p[0]).toleft(a-p[0])==-1 || (p.back()-p[0]).toleft(a-p[0])==1) return 0;
        const auto cmp=[&](const Point &u,const Point &v){return (u-p[0]).toleft(v-p[0])==1;};
        const size_t i=lower_bound(p.begin()+1,p.end(),a,cmp)-p.begin();
        if (i==1) return segment<T>{p[0],p[i]}.is_on(a)?-1:0;
        if (i==p.size()-1 && segment<T>{p[0],p[i]}.is_on(a)) return -1;
        if (segment<T>{p[i-1],p[i]}.is_on(a)) return -1;
        return (p[i]-p[i-1]).toleft(a-p[i-1])>0;
    }

    // 凸多边形关于某一方向的极点
    // 复杂度 O(logn)
    // 参考资料：https://codeforces.com/blog/entry/48868
    template<typename F> size_t extreme(const F &dir) const
    {
        const auto &p=this->p;
        const auto check=[&](const size_t i){return dir(p[i]).toleft(p[this->nxt(i)]-p[i])>=0;};
        const auto dir0=dir(p[0]); const auto check0=check(0);
        if (!check0 && check(p.size()-1)) return 0;
        const auto cmp=[&](const Point &v)
        {
            const size_t vi=&v-p.data();
            if (vi==0) return 1;
            const auto checkv=check(vi);
            const auto t=dir0.toleft(v-p[0]);
            if (vi==1 && checkv==check0 && t==0) return 1;
            return checkv^(checkv==check0 && t<=0);
        };
        return partition_point(p.begin(),p.end(),cmp)-p.begin();
    }

    // 过凸多边形外一点求凸多边形的切线，返回切点下标
    // 复杂度 O(logn)
    // 必须保证点在多边形外
    pair<size_t,size_t> tangent(const point<T> &a) const
    {
        const size_t i=extreme([&](const point<T> &u){return u-a;});
        const size_t j=extreme([&](const point<T> &u){return a-u;});
        return {i,j};
    }

    // 求平行于给定直线的凸多边形的切线，返回切点下标
    // 复杂度 O(logn)
    pair<size_t,size_t> tangent(const line<T> &a) const
    {
        const size_t i=extreme([&](...){return a.v;});
        const size_t j=extreme([&](...){return -a.v;});
        return {i,j};
    }
};

using Convex=convex<_T>;

// 圆
struct Circle
{
    Point c;
    long double r;

    bool operator==(const Circle &a) const {return c==a.c && abs(r-a.r)<=eps;}
    bool operator<(const Circle &a) const {return make_pair(c,r)<make_pair(a.c,a.r);}
    long double circ() const {return 2*PI*r;}
    long double area() const {return PI*r*r;}
    int is_in(const Point &p) const {const long double d=p.dis(c); return abs(d-r)<=eps?-1:d<r-eps;}

    // 直线与圆关系
    // 0 相离 | 1 相切 | 2 相交
    int relation(const Line &l) const
    {
        const long double d=l.dis(c);
        if (d>r+eps) return 0;
        if (abs(d-r)<=eps) return 1;
        return 2;
    }

    // 圆与圆关系
    // -1 相同 | 0 相离 | 1 外切 | 2 相交 | 3 内切 | 4 内含
    int relation(const Circle &a) const
    {
        if (*this==a) return -1;
        const long double d=c.dis(a.c);
        if (d>r+a.r+eps) return 0;
        if (abs(d-r-a.r)<=eps) return 1;
        if (abs(d-abs(r-a.r))<=eps) return 3;
        if (d<abs(r-a.r)-eps) return 4;
        return 2;
    }

    // 直线与圆的交点
    vector<Point> inter(const Line &l) const
    {
        const long double d=l.dis(c);
        const Point p=l.proj(c);
        const int t=relation(l);
        if (t==0) return vector<Point>();
        if (t==1) return vector<Point>{p};
        const long double k=sqrt(r*r-d*d);
        return vector<Point>{p-(l.v/l.v.len())*k,p+(l.v/l.v.len())*k};
    }

    // 圆与圆交点
    vector<Point> inter(const Circle &a) const
    {
        if ((*this)==a) return vector<Point>();
        const long double d=c.dis(a.c);
        const int t=relation(a);
        Point e=a.c-c; e=e/e.len();
        if (t==0 || t==4) return vector<Point>();
        if (t==1) return vector<Point>{c+e*r};
        if (t==3) return vector<Point>{c-e*r};
        const long double costh=(r*r+d*d-a.r*a.r)/(2*r*d),sinth=sqrt(1-costh*costh);
        return vector<Point>{c+e.rot(costh,-sinth)*r,c+e.rot(costh,sinth)*r};
    }

    // 圆与圆交面积
    long double inter_area(const Circle &a) const
    {
        if ((*this)==a) return area();
        const long double d=c.dis(a.c);
        const int t=relation(a);
        if (t<2) return 0;
        if (t>2) return min(area(),a.area());
        const long double costh1=(r*r+d*d-a.r*a.r)/(2*r*d),costh2=(a.r*a.r+d*d-r*r)/(2*a.r*d);
        const long double sinth1=sqrt(1-costh1*costh1),sinth2=sqrt(1-costh2*costh2);
        const long double th1=acos(costh1),th2=acos(costh2);
        return area()*th1/PI-r*r*costh1*sinth1+a.area()*th2/PI-a.r*a.r*costh2*sinth2;
    }
};

// 点集的凸包
// Andrew 算法，复杂度 O(nlogn)
Convex convexhull(vector<Point> p)
{
    vector<Point> st;
    sort(p.begin(),p.end());
    const auto check=[](const vector<Point> &st,const Point &u)
    {
        const auto back1=st.back(),back2=*prev(st.end(),2);
        return (back1-back2).toleft(u-back2)<=0;
    };
    for (const Point &u:p)
    {
        while (st.size()>1 && check(st,u)) st.pop_back();
        st.push_back(u);
    }
    size_t k=st.size();
    p.pop_back(); reverse(p.begin(),p.end());
    for (const Point &u:p)
    {
        while (st.size()>k && check(st,u)) st.pop_back();
        st.push_back(u);
    }
    st.pop_back();
    return {st};
}

// 半平面交
// 排序增量法，复杂度 O(nlogn)
// 输入与返回值都是用直线表示的半平面集合
vector<Line> halfinter(vector<Line> l, const _T lim=1e9)
{
    const auto check=[](const Line &a,const Line &b,const Line &c){return a.toleft(b.inter(c))<0;};
    // 无精度误差的方法，但注意取值范围会扩大到三次方
    /*const auto check=[](const Line &a,const Line &b,const Line &c)
    {
        const Point p=a.v*(b.v^c.v),q=b.p*(b.v^c.v)+b.v*(c.v^(b.p-c.p))-a.p*(b.v^c.v);
        return p.toleft(q)<0;
    };*/
    l.push_back({{-lim,0},{0,-1}}); l.push_back({{0,-lim},{1,0}});
    l.push_back({{lim,0},{0,1}}); l.push_back({{0,lim},{-1,0}});
    sort(l.begin(),l.end());
    deque<Line> q;
    for (size_t i=0;i<l.size();i++)
    {
        if (i>0 && l[i-1].v.toleft(l[i].v)==0 && l[i-1].v*l[i].v>eps) continue;
        while (q.size()>1 && check(l[i],q.back(),q[q.size()-2])) q.pop_back();
        while (q.size()>1 && check(l[i],q[0],q[1])) q.pop_front();
        if (!q.empty() && q.back().v.toleft(l[i].v)<=0) return vector<Line>();
        q.push_back(l[i]);
    }
    while (q.size()>1 && check(q[0],q.back(),q[q.size()-2])) q.pop_back();
    while (q.size()>1 && check(q.back(),q[0],q[1])) q.pop_front();
    return vector<Line>(q.begin(),q.end());
}

// 点集形成的最小最大三角形
// 极角序扫描线，复杂度 O(n^2logn)
// 最大三角形问题可以使用凸包与旋转卡壳做到 O(n^2)
pair<_T,_T> minmax_triangle(const vector<Point> &vec)
{
    if (vec.size()<=2) return {0,0};
    vector<pair<int,int>> evt;
    evt.reserve(vec.size()*vec.size());
    _T maxans=0,minans=numeric_limits<_T>::max();
    for (size_t i=0;i<vec.size();i++)
    {
        for (size_t j=0;j<vec.size();j++)
        {
            if (i==j) continue;
            if (vec[i]==vec[j]) minans=0;
            else evt.push_back({i,j});
        }
    }
    sort(evt.begin(),evt.end(),[&](const pair<int,int> &u,const pair<int,int> &v)
    {
        const Point du=vec[u.second]-vec[u.first],dv=vec[v.second]-vec[v.first];
        return argcmp()({du.y,-du.x},{dv.y,-dv.x});
    });
    vector<size_t> vx(vec.size()),pos(vec.size());
    for (size_t i=0;i<vec.size();i++) vx[i]=i;
    sort(vx.begin(),vx.end(),[&](int x,int y){return vec[x]<vec[y];});
    for (size_t i=0;i<vx.size();i++) pos[vx[i]]=i;
    for (auto [u,v]:evt)
    {
        const size_t i=pos[u],j=pos[v];
        const size_t l=min(i,j),r=max(i,j);
        const Point vecu=vec[u],vecv=vec[v];
        if (l>0) minans=min(minans,abs((vec[vx[l-1]]-vecu)^(vec[vx[l-1]]-vecv)));
        if (r<vx.size()-1) minans=min(minans,abs((vec[vx[r+1]]-vecu)^(vec[vx[r+1]]-vecv)));
        maxans=max({maxans,abs((vec[vx[0]]-vecu)^(vec[vx[0]]-vecv)),abs((vec[vx.back()]-vecu)^(vec[vx.back()]-vecv))});
        if (i<j) swap(vx[i],vx[j]),pos[u]=j,pos[v]=i;
    }
    return {minans,maxans};
}

// 判断多条线段是否有交点
// 扫描线，复杂度 O(nlogn)
bool segs_inter(const vector<Segment> &segs)
{
    if (segs.empty()) return false;
    using seq_t=tuple<_T,int,Segment>;
    const auto seqcmp=[](const seq_t &u, const seq_t &v)
    {
        const auto [u0,u1,u2]=u;
        const auto [v0,v1,v2]=v;
        if (abs(u0-v0) <= eps) return make_pair(u1,u2) < make_pair(v1,v2);
        return u0<v0-eps;
    };
    vector<seq_t> seq;
    for (auto seg:segs)
    {
        if (seg.a.x>seg.b.x+eps) swap(seg.a,seg.b);
        seq.push_back({seg.a.x,0,seg});
        seq.push_back({seg.b.x,1,seg});
    }
    sort(seq.begin(),seq.end(),seqcmp);
    _T x_now;
    auto cmp=[&](const Segment &u, const Segment &v)
    {
        if (abs(u.a.x-u.b.x)<=eps || abs(v.a.x-v.b.x)<=eps) return u.a.y<v.a.y-eps;
        return ((x_now-u.a.x)*(u.b.y-u.a.y)+u.a.y*(u.b.x-u.a.x))*(v.b.x-v.a.x)<((x_now-v.a.x)*(v.b.y-v.a.y)+v.a.y*(v.b.x-v.a.x))*(u.b.x-u.a.x)-eps;
    };
    multiset<Segment,decltype(cmp)> s{cmp};
    for (const auto [x,o,seg]:seq)
    {
        x_now=x;
        const auto it=s.lower_bound(seg);
        if (o==0)
        {
            if (it!=s.end() && seg.is_inter(*it)) return true;
            if (it!=s.begin() && seg.is_inter(*prev(it))) return true;
            s.insert(seg);
        }
        else
        {
            if (next(it)!=s.end() && it!=s.begin() && (*prev(it)).is_inter(*next(it))) return true;
            s.erase(it);
        }
    }
    return false;
}
```



## 杂项

### int128读入输出

```cpp
typedef __int128 LL;
inline __int128 read()
{
    __int128 x=0,f=1;
    char ch=getchar();
    while(ch<'0'||ch>'9')
    {
        if(ch=='-')
            f=-1;
        ch=getchar();
    }
    while(ch>='0'&&ch<='9')
    {
        x=x*10+ch-'0';
        ch=getchar();
    }
    return x*f;
}

inline void print(__int128 x)
{
    if(x < 0) putchar('-'), x = -x;
    if(x > 9) print(x / 10);
    putchar(x % 10 + '0');
}
```



### 快读

```c++
char buf[1<<23],*p1=buf,*p2=buf,obuf[1<<23],*O=obuf;

inline LL read(){
#define getchar() (p1==p2&&(p2=(p1=buf)+fread(buf,1,1<<21,stdin),p1==p2)?EOF:*p1++)
    LL sign = 1, x = 0;char s = getchar();
    while(s > '9' || s < '0' ){if(s == '-')sign = -1;s = getchar();}
    while(s >= '0' && s <= '9'){x = (x << 3) + (x << 1) + s - '0';s = getchar();}
    return x * sign;
#undef getchar
}//快读
void print(LL x) {
    if(x/ 10) print(x / 10);
    *O++=x % 10+'0';
}
void write(LL x) {
    if(x < 0)putchar('-'),x = -x;
    print(x);
    fwrite(obuf,O-obuf,1,stdout);
}
```



### modint

```c++
constexpr int mod = 998244353;

// v∈[-mod, 2 * mod)
int Norm(int v) { if(v < 0) v += mod; if(v >= mod) v -= mod; return v; }
// int/LL 的ksm别使用该函数, 此函数不带自动取模
template<class T> constexpr T power(T a, int b, const int &p = mod, T c = 1)
{
    for (a %= p;b;a *= a, b >>= 1) if (b & 1) c *= a;
    return c;
}

struct Z // 初始化的时候在int范围内
{
    int v;
    Z(int _x = 0) : v(Norm(_x)) {}
    // C++20 可用
    // auto operator<=>(const Z &) const = default;
    Z operator-() const { return Z(Norm(mod - v)); }
    Z inv() const { return power(*this, mod - 2, mod); }
    Z &operator*=(const Z &R) { return v = (long long)v * R.v % mod, *this; }
    Z &operator+=(const Z &R) { return v = Norm(v + R.v), *this; }
    Z &operator-=(const Z &R) { return v = Norm(v - R.v), *this; }
    Z &operator/=(const Z &R) { return *this *= R.inv(); }
    Z &operator%=(const int &R) { return v %= R, *this; }
    friend Z operator*(Z L, const Z &R) { return L *= R; }
    friend Z operator+(Z L, const Z &R) { return L += R; }
    friend Z operator-(Z L, const Z &R) { return L -= R; }
    friend Z operator/(Z L, const Z &R) { return L /= R; }
    friend Z operator%(Z L, const int &R) { return L %= R; }
    // C++ 14
    friend auto &operator>>(istream &i, Z &z) { return i >> z.v; }
    friend auto &operator<<(ostream &o, const Z &z) { return o << z.v; }
};
```



### 杜教BM

```cpp
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cassert>
#include<bits/stdc++.h>
using namespace std;
#define rep(i,a,n) for (int i=a;i<n;i++)
#define per(i,a,n) for (int i=n-1;i>=a;i--)
#define pb push_back
#define mp make_pair
#define all(x) (x).begin(),(x).end()
#define fi first
#define se second
#define SZ(x) ((int)(x).size())
typedef vector<int> VI;
typedef long long ll;
typedef pair<int,int> PII;
const ll mod= 998244353;
ll powmod(ll a,ll b) {ll res=1;a%=mod; assert(b>=0); for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}
// head
int _,n;
namespace linear_seq {
    const int N=10010;
    ll res[N],base[N],_c[N],_md[N];
 
    vector<int> Md;
    void mul(ll *a,ll *b,int k) {
        rep(i,0,k+k) _c[i]=0;
        rep(i,0,k) if (a[i]) rep(j,0,k) _c[i+j]=(_c[i+j]+a[i]*b[j])%mod;
        for (int i=k+k-1;i>=k;i--) if (_c[i])
            rep(j,0,SZ(Md)) _c[i-k+Md[j]]=(_c[i-k+Md[j]]-_c[i]*_md[Md[j]])%mod;
        rep(i,0,k) a[i]=_c[i];
    }
    int solve(ll n,VI a,VI b) { // a 系数 b 初值 b[n+1]=a[0]*b[n]+...
        ll ans=0,pnt=0;
        int k=SZ(a);
        assert(SZ(a)==SZ(b));
        rep(i,0,k) _md[k-1-i]=-a[i];_md[k]=1;
        Md.clear();
        rep(i,0,k) if (_md[i]!=0) Md.push_back(i);
        rep(i,0,k) res[i]=base[i]=0;
        res[0]=1;
        while ((1ll<<pnt)<=n) pnt++;
        for (int p=pnt;p>=0;p--) {
            mul(res,res,k);
            if ((n>>p)&1) {
                for (int i=k-1;i>=0;i--) res[i+1]=res[i];res[0]=0;
                rep(j,0,SZ(Md)) res[Md[j]]=(res[Md[j]]-res[k]*_md[Md[j]])%mod;
            }
        }
        rep(i,0,k) ans=(ans+res[i]*b[i])%mod;
        if (ans<0) ans+=mod;
        return ans;
    }
    VI BM(VI s) {
        VI C(1,1),B(1,1);
        int L=0,m=1,b=1;
        rep(n,0,SZ(s)) {
            ll d=0;
            rep(i,0,L+1) d=(d+(ll)C[i]*s[n-i])%mod;
            if (d==0) ++m;
            else if (2*L<=n) {
                VI T=C;
                ll c=mod-d*powmod(b,mod-2)%mod;
                while (SZ(C)<SZ(B)+m) C.pb(0);
                rep(i,0,SZ(B)) C[i+m]=(C[i+m]+c*B[i])%mod;
                L=n+1-L; B=T; b=d; m=1;
            } else {
                ll c=mod-d*powmod(b,mod-2)%mod;
                while (SZ(C)<SZ(B)+m) C.pb(0);
                rep(i,0,SZ(B)) C[i+m]=(C[i+m]+c*B[i])%mod;
                ++m;
            }
        }
        return C;
    }
    int gao(VI a,ll n) {
        VI c=BM(a);
        c.erase(c.begin());
        rep(i,0,SZ(c)) c[i]=(mod-c[i])%mod;
        return solve(n,c,VI(a.begin(),a.begin()+SZ(c)));
    }
};

int main() 
{
    while (~scanf("%d",&n)) 
	{
        vector<int>v;
        v.push_back(4);
        v.push_back(12);
        v.push_back(33);
        v.push_back(88);
        v.push_back(232);
        v.push_back(609);
        printf("%d\n",linear_seq::gao(v,n-1));
    }
}

```



### poLLard_rho

```cpp
//sqrt(sqrt(n))的时间复杂度分解数的质因子

#include <bits/stdc++.h>

using namespace std;

typedef long long LL;

const int N = 1e5 + 10;

LL factor[100], total = 0; // 计算p的约数

LL gcd(LL a, LL b) { return b ? gcd(b, a % b) : a; }

// 快速幂
LL qmi(LL a, LL k, LL p)  // 求a^k mod p
{
    LL res = 1 % p;
    while (k)
    {
        if (k & 1) res = (__int128)res * a % p;
        a = (__int128)a * a % p;
        k >>= 1;
    }
    return res;
}

int a[] = {2, 3, 5, 7, 11};
// 判断一个数是否为质数

bool miLLer_rabin(LL n)
{
    if (n == 1) return 0;
    for (int i = 0; i < 5; ++i)
    {
        if (n == a[i]) return 1;
        if (!(n % a[i])) return 0;
    }
    
    LL d = n - 1;
    while (!(d & 1)) d >>= 1;
    for (int i = 1; i <= 5; ++i)
    {
        LL a = rand() % (n - 2) + 2;
        bool tmp_k = false;
        LL tmp_d = d;
        LL tmp = qmi(a, tmp_d, n);
        if (tmp == 1) tmp_k = true;
        // 先把所有的2给提出来，然后在一个一个乘上去，这样比直接做快
        while (tmp != 1 && tmp != n - 1 && tmp_d != n - 1)
        {
            tmp = (__int128) tmp * tmp % n;
            tmp_d <<= 1;
        } // 二次探测
        
        if (tmp == n - 1) tmp_k = true;
        if (!tmp_k) return 0;
    }
    return 1;
}

// 找出n的一个因数
const int M = (1 << 7) - 1; // 一个玄学值

LL poLLard_rho(LL n)
{
    for (int i = 0; i < 5; ++i)
        if (n % a[i] == 0) return a[i];
        
    LL x = rand(), y = x, t = 1, a = rand() % (n - 2) + 2;
    for (int k = 2;; k <<= 1, y = x)
    {
        LL q = 1;
        // k用于倍增
        for (int i = 1; i <= k; ++i)
        {
            x = ((__int128)x * x % n + a) % n; // f(x) = x^2+a
            q = (__int128)q * abs(x - y) % n;  // 每一次的|x-y|累加
            // 如果次数超过M
            if (!(i & M))
            {
                t = gcd(q, n);
                if (t > 1) break;
            }
        }
        
        if (t > 1 || (t = gcd(q, n)) > 1) break;
    }
    return t;
}

// 找出x的所有因数
void find(LL x)
{
    if (x == 1 || x <= 0) return;
    if (miLLer_rabin(x))
    {
        factor[++total] = x;
        return;
    }
    
    LL p = x;
    while (p == x) p = poLLard_rho(x);
    while (x % p == 0) x /= p;
    find(p); // 递归遍历所有的因数
    find(x);
}

LL get_phi(LL x)
{
    LL res = x;
    total = 0;
    find(x);
    for (int i = 1; i <= total; i++)
    {
        res = res / factor[i] * (factor[i] - 1);
    }
    
    return res;
}

int main()
{
    cout << get_phi(6516156) << endl;
}
```



### PII的hash值

```cpp
template <typename T>
inline void hash_combine(std::size_t &seed, const T &val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
// auxiliary generic functions to create a hash value using a seed
template <typename T> inline void hash_val(std::size_t &seed, const T &val) {
    hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t &seed, const T &val, const Types &... args) {
    hash_combine(seed, val);
    hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        return hash_val(p.first, p.second);
    }
};

```



### 手写hash

```cpp
template <typename A, typename B>
struct unmap
{
    struct node
    {
        A u;
        B v;
        int nxt;
    };
    static const unsigned uN = 2000003;
    vector<node> e;
    int head[uN], rem[uN], cnt;
    unmap() { memset(head, -1, sizeof(head)); }
    bool count(A u)
    {
        int hs = u % uN;
        for (int i = head[hs]; ~i; i = e[i].nxt)
            if (e[i].u == u) return 1;
        return 0;
    }

    B &operator[](A u)
    {
        int hs = u % uN;
        for (int i = head[hs]; ~i; i = e[i].nxt)
            if (e[i].u == u) return e[i].v;

        e.push_back({u, B(), head[hs]});
        head[hs] = e.size() - 1;
        rem[++ cnt] = hs;
        return e.back().v;
    }

    void clear() {
        e.clear();
        for (int i = 1; i <= cnt; i ++ ) head[rem[i]] = -1;
        cnt = 0;
    }
};
```

### 德州扑克板子

```cpp
struct card{
  char suit;
  int rank;
  card(){
  }
  bool operator <(card C){
    return rank < C.rank || rank == C.rank && suit < C.suit;
  }
};
istream& operator >>(istream& is, card& C){
  string S;
  is >> S;
  if (S[0] == 'A'){
    C.rank = 14;
  } else if (S[0] == 'K'){
    C.rank = 13;
  } else if (S[0] == 'Q'){
    C.rank = 12;
  } else if (S[0] == 'J'){
    C.rank = 11;
  } else if (S[0] == 'T'){
    C.rank = 10;
  } else {
    C.rank = S[0] - '0';
  }
  C.suit = S[1];
  return is;
}
vector<int> hand(vector<card> C){
  sort(C.begin(), C.end());
  set<char> suits;
  for (int i = 0; i < 5; i++){
    suits.insert(C[i].suit);
  }
  if (suits.size() == 1 && C[4].rank - C[0].rank == 4){
    if (C[4].rank == 14){
      return vector<int>{9};
    } else {
      return vector<int>{8, C[4].rank};
    }
  }
  if (suits.size() == 1 && C[3].rank == 5 && C[4].rank == 14){
    return vector<int>{8, 5};
  }
  if (C[0].rank == C[3].rank){
    return vector<int>{7, C[0].rank, C[4].rank};
  }
  if (C[1].rank == C[4].rank){
    return vector<int>{7, C[1].rank, C[0].rank};
  }
  if (C[0].rank == C[2].rank && C[3].rank == C[4].rank){
    return vector<int>{6, C[0].rank, C[3].rank};
  }
  if (C[2].rank == C[4].rank && C[0].rank == C[1].rank){
    return vector<int>{6, C[2].rank, C[0].rank};
  }
  if (suits.size() == 1){
    return vector<int>{5, C[4].rank, C[3].rank, C[2].rank, C[1].rank, C[0].rank};
  }
  if (C[1].rank - C[0].rank == 1 && C[2].rank - C[1].rank == 1 && C[3].rank - C[2].rank == 1 && C[4].rank - C[3].rank == 1){
    return vector<int>{4, C[4].rank};
  }
  if (C[0].rank == 2 && C[1].rank == 3 && C[2].rank == 4 && C[3].rank == 5 && C[4].rank == 14){
    return vector<int>{4, 5};
  }
  if (C[0].rank == C[2].rank){
    return vector<int>{3, C[0].rank, C[4].rank, C[3].rank};
  }
  if (C[1].rank == C[3].rank){
    return vector<int>{3, C[1].rank, C[4].rank, C[0].rank};
  }
  if (C[2].rank == C[4].rank){
    return vector<int>{3, C[2].rank, C[1].rank, C[0].rank};
  }
  if (C[0].rank == C[1].rank && C[2].rank == C[3].rank){
    return vector<int>{2, C[2].rank, C[0].rank, C[4].rank};
  }
  if (C[0].rank == C[1].rank && C[3].rank == C[4].rank){
    return vector<int>{2, C[3].rank, C[0].rank, C[2].rank};
  }
  if (C[1].rank == C[2].rank && C[3].rank == C[4].rank){
    return vector<int>{2, C[3].rank, C[1].rank, C[0].rank};
  }
  if (C[0].rank == C[1].rank){
    return vector<int>{1, C[0].rank, C[4].rank, C[3].rank, C[2].rank};
  }
  if (C[1].rank == C[2].rank){
    return vector<int>{1, C[1].rank, C[4].rank, C[3].rank, C[0].rank};
  }
  if (C[2].rank == C[3].rank){
    return vector<int>{1, C[2].rank, C[4].rank, C[1].rank, C[0].rank};
  }
  if (C[3].rank == C[4].rank){
    return vector<int>{1, C[3].rank, C[2].rank, C[1].rank, C[0].rank};
  }
  return vector<int>{0, C[4].rank, C[3].rank, C[2].rank, C[1].rank, C[0].rank};
}
```

