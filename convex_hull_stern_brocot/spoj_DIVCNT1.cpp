#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("avx2")
#pragma GCC optimize("Os")

#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <utility>
#include <vector>

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

long long gcd(long long x, long long y) {
  return y ? gcd(y, x % y) : x;
}
__int128 S(long long N) {
  if (N < 12) {
    long long ret = 0;
    for (long long i = 1; i <= N; ++i)
      ret += N / i;
    return ret;
  }
  std::function<double(long long)> f = [&](long long x) { return N * 1.0 / x; };
  std::function<double(long long)> df = [&](long long x) {
    return -N * 1.0 / x / x;
  };
  std::function<double(long long, long long)> outside =
      [&](long long x, long long y) { return (__int128)x * y > (__int128)N; };

  std::vector<std::pair<long long, long long>> hull;
  long long x0 = sqrtl(N + 0.5), x1 = cbrtl(N) * cbrtl(N);
  convex_hull(N, x0 + 1, N / (x0 + 1) + 1, x1, f, df, outside, hull);
  for (long long y = hull.back().second - 1; y >= 2; --y)
    hull.emplace_back(N / y + 1, y);
  //for (const auto &i : hull)
  //  printf("(%lld, %lld), ", i.first, i.second);
  __int128 ret = 1 + (__int128)N - (N / 2 + 1);
  for (size_t i = 0; i + 1 < hull.size(); ++i) {
    const std::pair<long long, long long> &p1 = hull.at(i), &p2 = hull.at(i + 1);
    long long g = gcd(p2.first - p1.first, p1.second - p2.second);
    __int128 points = (__int128)(p2.first - p1.first) + g + p2.second + p1.second;
    __int128 internal = ((__int128)(p2.second + p1.second) * (p2.first - p1.first) - points) / 2 + 1;
    ret = ret + internal + p1.second - 1;
  }
  return ret * 2 + x0 * x0;
}

void print(__int128 x) {
  if (x < 10) {printf("%d", (int)x); return; }
  print(x / 10);
  printf("%d", (int)(x % 10));
}

int main() {
  /*for (int n = 1; n <= 100000; ++n) {
    long long correct = 0;
    for (int i = 1; i <= n; ++i)
      correct += n / i;
    long long my = S(n);
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
