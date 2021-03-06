#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <utility>
#include <algorithm>
#include <vector>
#include <chrono>

const double eps = 1e-9;
std::pair<long long, long long>
convex_hull(const long long &N, const long long &x1, const long long &y1,
            const long long &x2, const std::function<double(long long)> &f,
            const std::function<double(long long)> &df,
            const std::function<bool(long long, long long)> &inside,
            std::vector<std::pair<long long, long long>> &ret) {
  // x1: from (inclusive)  x2: to (exclusive)
  // Note that for concave function, the inside function means the points on f
  // should be included
  // but for convex function, the inside function means the outside
  // return:
  //   ret: the vertices of the convex hull under the function f
  //   return value: the slope of the convex hull (x,y) at x1

  // Assumptions:
  // f(x) is positive on x1<=x<x2
  // df(x) is monotonic
  // df(x) has the same sign on x1<=x<x2
  // (x1, y1) should be on the convex hull (I don't check it)

  assert(x1 < x2);
  const bool convex = df(x1) < df(x1 + 1);
  if (!convex && inside(x1, y1) || convex && inside(x1, y1))
    ret.emplace_back(std::make_pair(x1, y1));
  if (x1 + 1 == x2)
    return std::make_pair(0, 0);

  const int sign = (df(x1) < 0) ? -1 : 1;
  std::vector<std::pair<long long, long long>> stac; // (x, y)
  std::pair<long long, long long> first_slope;
  if (!convex) {
    long long x = x1, y = y1;
    if (sign == -1) {
      stac.emplace_back(0, 1);
      if ((long long)f(x1) == (long long)f(x1 + 1)) {
        first_slope = std::make_pair(1, 0);
        stac.emplace_back(1, 0);
      } else {
        std::pair<long long, long long> left = std::make_pair(1, 0),
                                        right = std::make_pair(0, 1);
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              inside(x + mid.first, y + sign * mid.second)) {
            right = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first >= right.second)
              break;
            left = mid;
          }
        }
        /*printf("my: ");
        for (const auto &i:stac)
            printf("%lld/%lld, ",i.second,i.first);
        printf("\n");*/
        first_slope = stac.back();
      }
      while (true) {
        std::pair<long long, long long> left, right;
        right = stac.back();
        stac.pop_back();
        bool assigned = false;
        std::pair<long long, long long> last;
        while (x + right.first < x2 &&
               inside(x + right.first, y + sign * right.second)) {
          x += right.first;
          y += sign * right.second;
          if (x < x2) {
            assigned = true;
            last = std::make_pair(x, y);
          }
        }
        if (assigned && last.second >= 0)
          ret.emplace_back(last);
        left = right;
        while (!stac.empty()) {
          right = stac.back();
          if (x + right.first >= x2)
            break;
          if (inside(x + right.first, y + sign * right.second))
            break;
          stac.pop_back();
          left = right;
        }
        if (stac.empty())
          break;
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              inside(x + mid.first, y + sign * mid.second)) {
            right = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first >= right.second)
              break;
            left = mid;
          }
        }
      }
      return first_slope;
    } else {
      stac.emplace_back(1, 0);
      if ((long long)f(x1) == (long long)f(x1 + 1)) {
        first_slope = std::make_pair(1, 0);
      } else {
        std::pair<long long, long long> left = std::make_pair(1, 0),
                                        right = std::make_pair(0, 1);
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              inside(x + mid.first, y + sign * mid.second)) {
            left = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first >= right.second)
              break;
            right = mid;
          }
        }
        /*printf("my: ");
        for (const auto &i:stac)
            printf("%lld/%lld, ",i.second,i.first);
        printf("\n");*/
        first_slope = stac.back();
      }
      while (true) {
        std::pair<long long, long long> left, right;
        left = stac.back();
        stac.pop_back();
        bool assigned = false;
        std::pair<long long, long long> last;
        while (x + left.first < x2 &&
               inside(x + left.first, y + sign * left.second)) {
          x += left.first;
          y += sign * left.second;
          if (x < x2) {
            assigned = true;
            last = std::make_pair(x, y);
          }
        }
        if (assigned && last.second >= 0)
          ret.emplace_back(last);
        right = left;
        while (!stac.empty()) {
          left = stac.back();
          if (x + left.first >= x2)
            break;
          if (inside(x + left.first, y + sign * left.second))
            break;
          stac.pop_back();
          right = left;
        }
        if (stac.empty())
          break;
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              inside(x + mid.first, y + sign * mid.second)) {
            left = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            double v = sign * df(x + mid.first);
            if (sign * df(x + mid.first) * right.first >= right.second)
              break;
            right = mid;
          }
        }
      }
      return first_slope;
    }
  } else {
    const auto &outside = inside;
    if (sign == -1) {
      long long x = x1, y = y1;
      while (!outside(x, y))
        ++y;
      stac.emplace_back(1, 0);
      long long y2 = f(x1 + 1);
      while (!outside(x1 + 1, y2))
        ++y2;
      if (y == y2) {
        first_slope = std::make_pair(1, 0);
        stac.emplace_back(0, 1);
      } else {
        std::pair<long long, long long> left = std::make_pair(0, 1),
                                        right = std::make_pair(1, 0);
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              outside(x + mid.first,
                      y + sign * mid.second)) // take care! It should be
                                              // equivalent to
          // f(x+mid.first)<=y+sign*mid.second
          {
            right = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first <= right.second)
              break;
            left = mid;
          }
        }
        /*printf("my: ");
        for (const auto &i : stac)
          printf("%lld/%lld, ", i.second, i.first);
        printf("\n");*/
        first_slope = stac.back();
      }
      while (true) {
        std::pair<long long, long long> left, right;
        right = stac.back();
        stac.pop_back();
        bool assigned = false;
        std::pair<long long, long long> last;
        while (x + right.first < x2 &&
               outside(x + right.first, y + sign * right.second)) {
          x += right.first;
          y += sign * right.second;
          if (x < x2) {
            assigned = true;
            last = std::make_pair(x, y);
            //if (y >= 0) ret.emplace_back(x, y);
          } else
            return first_slope;
        }
        if (assigned && last.second >= 0)
          ret.emplace_back(last);
        left = right;
        while (!stac.empty()) {
          right = stac.back();
          if (x + right.first >= x2)
            break;
          if (outside(x + right.first, y + sign * right.second))
            break;
          left = right;
          stac.pop_back();
        }
        if (stac.empty())
          break;
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              outside(x + mid.first, y + sign * mid.second)) {
            right = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first <= right.second)
              break;
            left = mid;
          }
        }
      }
      return first_slope;
    } else {
      long long x = x1, y = y1;
      while (!outside(x, y))
        ++y;
      stac.emplace_back(0, 1);
      long long y2 = f(x1 + 1);
      while (!outside(x1 + 1, y2))
        y2++;
      if (y == y2) {
        first_slope = std::make_pair(1, 0);
        stac.emplace_back(1, 0);
      } else {
        std::pair<long long, long long> left = std::make_pair(0, 1),
                                        right = std::make_pair(1, 0);
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              outside(
                  x + mid.first,
                  y + sign * mid.second)) // take care! It should be equivalent
                                          // to
                                          // f(x+mid.first)<=y+sign*mid.second
          {
            left = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first <= right.second)
              break;
            right = mid;
          }
        }
        /*printf("my: ");
        for (const auto &i : stac)
          printf("%lld/%lld, ", i.second, i.first);
        printf("\n");*/
        first_slope = stac.back();
      }
      while (true) {
        std::pair<long long, long long> left, right;
        left = stac.back();
        if (!left.first)
          break;
        stac.pop_back();
        bool assigned = false;
        std::pair<long long, long long> last;
        while (x + left.first < x2 &&
               outside(x + left.first, y + sign * left.second)) {
          x += left.first;
          y += sign * left.second;
          if (x < x2) {
            assigned = true;
            last = std::make_pair(x, y - 1);
          } else
            break;
        }
        if (assigned && last.second >= 0)
          ret.emplace_back(last);
        right = left;
        while (!stac.empty()) {
          left = stac.back();
          if (x + left.first >= x2)
            break;
          if (outside(x + left.first, y + sign * left.second))
            break;
          right = left;
          stac.pop_back();
        }
        if (stac.empty())
          break;
        while (true) {
          std::pair<long long, long long> mid = std::make_pair(
              left.first + right.first, left.second + right.second);
          if (x + mid.first < x2 &&
              outside(
                  x + mid.first,
                  y + sign * mid.second)) // take care! It should be equivalent
                                          // to
                                          // f(x+mid.first)<=y+sign*mid.second
          {
            left = mid;
            stac.emplace_back(mid);
          } else {
            if (x + mid.first >= x2)
              break;
            if (sign * df(x + mid.first) * right.first <= right.second)
              break;
            right = mid;
          }
        }
      }
      return first_slope;
    }
  }
}

inline __int128 S1(__int128 x){return (x&1)?((x+1)/2*x):(x/2*(x+1));}
inline __int128 S2(__int128 x){__int128 a=S1(x),b=2*x+1;return (a%3)?(b/3*a):(a/3*b);}
inline __int128 Sqr(__int128 x){return x*x;}
struct node{
	__int128_t f,g,h;
	node(){f=0,g=0,h=0;}
};
node calc(__int128_t a,__int128_t b,__int128_t c,__int128_t n){
	node ans,res;
	if(!n){ans.f=b/c;ans.g=Sqr(b/c);return ans;}
  __int128_t m,t1,t2,s1,s2;
	if(!a){
		t1=b/c;
		ans.f=(n+1)*t1;
		ans.g=(n+1)*Sqr(t1);
		ans.h=S1(n)*t1;
		return ans;
	}
	if (c==1) {
    ans.f=a*S1(n)+b*(n+1);
    ans.g=a*a*S2(n)+2*a*b*S1(n)+b*b*(n+1);
    ans.h=a*S2(n)+b*S1(n);
    return ans;
  }
	if(a>=c||b>=c){
		t1=a/c;t2=b/c;
		res=calc(a%c,b%c,c,n);
		s1=S1(n);s2=S2(n);
		ans.f=(((s1*t1)+(n+1)*t2)+res.f);
		ans.g=(((Sqr(t1)*s2+(n+1)*Sqr(t2)))+((t1*t2)*2*s1+(t1*2*res.h))+(res.g+t2*2*res.f));
		ans.h=((s2*t1+s1*t2)+res.h);
		return ans;
	}
	m=(n*a+b)/c-1;
	res=calc(c,c-b-1,a,m);
	__int128_t w1=n*(m+1),w2=n*(n+1),w3=m+1;
	ans.f=(w1-res.f);
	ans.g=((w1*w3)-((res.h*2+res.f)));
	ans.h=((w2*w3)-(res.f+res.g))/2;
	return ans;
}

__int128_t T(long long a, __int128_t b, long long c, long long n) {
  // calculate \sum_{x=0}^{n} floor((ax+b)/c)*x where a is negative
  if (n < 10) {
    __int128_t correct = 0;
    for (long long x = 1; x <= n; ++x)
      correct = correct + (a * x + b) / c * x;
    return correct;
  }
  node _node = calc(-a, b + a * n, c, n);
  __int128_t my = _node.f * n - _node.h;
  /*if (my != correct) {
    _node = calc(-a, b + a * n, c, n);
    assert(false);
  }*/
  return my;
}

long long gcd(long long x, long long y) {
  return y ? gcd(y, x % y) : x;
}
__int128_t S(long long N) {
  if (N < 12) {
    long long ret = 0;
    for (long long i = 1; i <= N; ++i)
      ret += N / i * i;
    ret -= N * (N + 1) / 2;
    return ret;
  }
  std::function<double(long long)> f = [&](long long x) { return N * 1.0 / x; };
  std::function<double(long long)> df = [&](long long x) {
    return -N * 1.0 / x / x;
  };
  std::function<double(long long, long long)> outside =
      [&](long long x, long long y) { return (__int128_t)x * y > (__int128_t)N; };

  std::vector<std::pair<long long, long long>> hull;
  long long x0 = sqrtl(N + 0.5), x1 = cbrtl(N) * cbrtl(N);
  convex_hull(N, x0 + 1, N / (x0 + 1) + 1, x1, f, df, outside, hull);
  //printf("size = %d\n", (int)hull.size());

  const std::pair<long long, long long> last = hull.back();
  //std::vector<std::pair<long long, long long>> hull_inv;
  //for (long long x = 2; x < hull.back().second; ++x)
  //  hull_inv.emplace_back(x, N / x + 1);
  //for (auto iter = hull.rbegin(); iter != hull.rend(); ++iter) {
  //  std::pair<long long, long long> p = std::make_pair(iter->second, iter->first);
  //  hull_inv.emplace_back(p);
  //}
  //for (long long y = hull.back().second - 1; y >= 2; --y)
  //  hull.emplace_back(N / y + 1, y);
  //hull.insert(hull.begin(), hull_inv.begin(), hull_inv.end());

  /*for (const auto &i : hull)
    printf("(%lld, %lld), ", i.first, i.second);
  printf("\n");
  for (const auto &i : hull_inv)
    printf("(%lld, %lld), ", i.first, i.second);
  printf("\n");*/
  __int128_t ret = 0;
  __int128_t k = N / 2;
  if (N & 1) ret = ret + ((k & 1) ? ((k + 1) / 2 * (3 * k + 2)) : ((3 * k + 2) / 2 * (k + 1)));
  else ret = ret + ((k & 1) ? ((3 * k + 1) / 2 * k) : (k / 2 * (3 * k + 1)));
  for (long long x = 1; x < last.second; ++x)
    ret += N / x * x;
  auto section = [&](const std::pair<long long, long long> &p1, const std::pair<long long, long long> &p2, const bool &flip=false) -> __int128_t {
    if (!flip) {
      __int128_t a = p2.first - p1.first, b = p1.second - p2.second;
      if (a == 0) return 0;
      //__int128_t c = (__int128_t)p1.first * b + (__int128_t)p1.second * a;
      __int128_t g = gcd(a, b);
      //__int128_t val_true = T(-b, c, a, p2.first - 1) - T(-b, c, a, p1.first - 1) - p1.first * g -  ((g & 1) ? ((g - 1) / 2 * g) : (g / 2 * (g - 1))) * (a / g);
      __int128_t points = (__int128_t)(p2.first - p1.first) + g + p2.second + p1.second;
      __int128_t internal = ((__int128_t)(p2.second + p1.second) * (p2.first - p1.first) - points) / 2 + 1;
      __int128_t val2 = internal + p1.second - 1 + g;
      __int128_t val = T(-b, (__int128_t)p1.second * a, a, a - 1) + val2 * p1.first - p1.first * g -  ((g & 1) ? ((g - 1) / 2 * g) : (g / 2 * (g - 1))) * (a / g);
      /*if (val != val_true) {
        assert(false);
      }*/
      return val;
    } else {
      __int128_t a = p2.second - p1.second, b = p1.first - p2.first;
      if (a == 0) return 0;
      //__int128_t c = (__int128_t)p1.second * b + (__int128_t)p1.first * a;
      __int128_t g = gcd(a, b);
      //__int128_t val_true = T(-b, c, a, p2.second - 1) - T(-b, c, a, p1.second - 1) - p1.second * g -  ((g & 1) ? ((g - 1) / 2 * g) : (g / 2 * (g - 1))) * (a / g);
      __int128_t points = (__int128_t)(p2.second - p1.second) + g + p2.first + p1.first;
      __int128_t internal = ((__int128_t)(p2.first + p1.first) * (p2.second - p1.second) - points) / 2 + 1;
      __int128_t val2 = internal + p1.first - 1 + g;
      __int128_t val = T(-b, (__int128_t)p1.first * a, a, a - 1) + val2 * p1.second - p1.second * g -  ((g & 1) ? ((g - 1) / 2 * g) : (g / 2 * (g - 1))) * (a / g);
      /*if (val != val_true) {
        assert(false);
      }*/
      return val;
    }
  };
  std::pair<long long, long long> now = hull.back();
  std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
  std::for_each(hull.rbegin(), hull.rend(), [&](const std::pair<long long, long long> &p) -> void {
                  ret += section(now, p, true);
                  now = p;
                });
  //printf("time = %lf\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count() / 1000.0);

  std::swap(now.first, now.second);
  start_time = std::chrono::system_clock::now();
  std::for_each(hull.begin(), hull.end(), [&](const std::pair<long long, long long> &p) -> void {
                  ret += section(now, p);
                  now = p;
                });
  //printf("time = %lf\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count() / 1000.0);

  /*for (long long y = last.second - 1; y >= 2; --y) {
    std::pair<long long, long long> p = std::make_pair(N / y + 1, y);
    ret += section(now, p);
    now = p;
  }*/
  for (long long y = 3, right = N / 2, left; y <= last.second; ++y) {
    left = N / y;
    if (left + 1 <= last.first) {
      ret += ((__int128_t)right - last.first + 1) * ((__int128_t)right + last.first) / 2 * (y - 1);
      break;
    }
    ret += ((__int128_t)right - left) * ((__int128_t)right + left + 1) / 2 * (y - 1);
    right = left;
  }
  /*for (size_t i = 0; i + 1 < hull.size(); ++i) {
    const std::pair<long long, long long> &p1 = hull.at(i), &p2 = hull.at(i + 1);
    ret += section(p1, p2);
    if (i % 10000000 == 0) printf("%d finished, total = %d\n", (int)i, (int)hull.size());
  }*/
  __uint128_t temp = (N & 1) ? (((__uint128_t)N + 1) / 2 * N) : ((__uint128_t)N / 2 * (N + 1));
  ret -= temp;
  return ret;
}

void print(__uint128_t x) {
  if (x < 10) {printf("%d", (int)x); return; }
  print(x / 10);
  printf("%d", (int)(x % 10));
}

int main() {
  // 9223372036854775807LL
  /*for (unsigned long long n = 9223372036854775807LL; n <= 9223372036854775807LL; ++n) {
    __uint128_t correct = 0;
    unsigned long long v = sqrtl(n);
    for (unsigned long long beg : {1, 5}) for (unsigned long long i = beg; i <= v; i += 6) {
      for (unsigned long long M = n / i, j = i; j <= v; j *= 3, M /= 3) {
        for (unsigned long long L = M, k = j; k <= v; k <<= 1, L >>= 1) {
          correct += __uint128_t(L) * (L + 1 + 2 * k);
        }
      }
    }
    correct /= 2;
    correct = correct - __uint128_t(v) * v * (v + 1) / 2;
    __uint128_t temp = (n & 1) ? (((__uint128_t)n + 1) / 2 * n) : ((__uint128_t)n / 2 * (n + 1));
    correct -= temp;
    print(correct),printf("\n");
    __int128 my = S(n);
    print(my),printf("\n");
    assert(correct == my);
  }*/
  int testcases;
  scanf("%d", &testcases);
  for (int i = 0; i < testcases; ++i) {
    long long n;
    scanf("%lld", &n);
    print(S(n)),printf("\n");
  }
  return 0;
}
