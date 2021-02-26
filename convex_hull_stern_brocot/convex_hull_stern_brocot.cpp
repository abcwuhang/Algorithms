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
  if (!convex && inside(x1, y1) || convex && !inside(x1, y1))
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
            last = std::make_pair(x, y - 1);
          } else
            break;
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

std::pair<long long, long long>
brute_force_convex_hull(const long long &N, const long long &x1,
                        const long long &y1, const long long &x2,
                        const std::function<double(long long)> &f,
                        const std::function<double(long long)> &df,
                        const std::function<bool(long long, long long)> &inside,
                        std::vector<std::pair<long long, long long>> &ret) {
  assert(x1 < x2);

  std::function<long long(long long, long long)> gcd;
  gcd = [&gcd](long long x, long long y) -> long long {
    return y ? gcd(y, x % y) : x;
  };
  auto simplify = [&gcd](long long a,
                         long long b) -> std::pair<long long, long long> {
    long long g = gcd(a, b);
    return std::make_pair(a / g, b / g);
  };

  const bool convex = df(x1) < df(x1 + 1);
  if (!convex && inside(x1, y1) || convex && !inside(x1, y1))
    ret.emplace_back(std::make_pair(x1, y1));
  if (x1 + 1 == x2)
    return std::make_pair(0, 0);

  const int sign = (df(x1) < 0) ? -1 : 1;
  if (!convex) {
    std::vector<std::pair<long long, long long>> points;
    for (long long x = x1; x < x2; ++x) {
      long long y = f(x) + eps;
      if (x == x1)
        y = y1;
      while (!inside(x, y))
        --y;
      points.emplace_back(x, y);
    }
    long long x = x1, y = y1;
    if (sign == -1) {
      for (long long now = 0;;) {
        std::pair<long long, long long> min_slope;
        size_t next = -1;
        for (size_t i = now + 1; i < points.size(); ++i) {
          std::pair<long long, long long> slope =
              simplify(abs(points.at(i).second - points.at(now).second),
                       abs(points.at(i).first - points.at(now).first));
          if (next == -1 ||
              next != -1 &&
                  min_slope.first * slope.second >=
                      slope.first * min_slope.second) {
            next = i;
            min_slope = slope;
          }
        }
        if (next == -1)
          break;
        ret.emplace_back(points.at(next));
        now = next;
      }
    } else {
      for (long long now = 0;;) {
        std::pair<long long, long long> max_slope;
        size_t next = -1;
        for (size_t i = now + 1; i < points.size(); ++i) {
          std::pair<long long, long long> slope =
              simplify(abs(points.at(i).second - points.at(now).second),
                       abs(points.at(i).first - points.at(now).first));
          if (next == -1 ||
              next != -1 &&
                  max_slope.first * slope.second <=
                      slope.first * max_slope.second) {
            next = i;
            max_slope = slope;
          }
        }
        if (next == -1)
          break;
        ret.emplace_back(points.at(next));
        now = next;
      }
    }
  } else {
    const auto &outside = inside;
    std::vector<std::pair<long long, long long>> points;
    for (long long x = x1; x < x2; ++x) {
      long long y = f(x) + eps;
      if (x == x1)
        y = y1;
      while (outside(x, y))
        --y;
      points.emplace_back(x, y);
    }
    long long x = x1, y = y1;
    if (sign == -1) {
      for (long long now = 0;;) {
        std::pair<long long, long long> max_slope;
        size_t next = -1;
        for (size_t i = now + 1; i < points.size(); ++i) {
          std::pair<long long, long long> slope =
              simplify(abs(points.at(i).second - points.at(now).second),
                       abs(points.at(i).first - points.at(now).first));
          if (next == -1 ||
              next != -1 &&
                  max_slope.first * slope.second <=
                      slope.first * max_slope.second) {
            next = i;
            max_slope = slope;
          }
        }
        if (next == -1)
          break;
        ret.emplace_back(points.at(next));
        now = next;
      }
    } else {
      for (long long now = 0;;) {
        std::pair<long long, long long> min_slope;
        size_t next = -1;
        for (size_t i = now + 1; i < points.size(); ++i) {
          std::pair<long long, long long> slope =
              simplify(abs(points.at(i).second - points.at(now).second),
                       abs(points.at(i).first - points.at(now).first));
          if (next == -1 ||
              next != -1 &&
                  min_slope.first * slope.second >=
                      slope.first * min_slope.second) {
            next = i;
            min_slope = slope;
          }
        }
        if (next == -1)
          break;
        ret.emplace_back(points.at(next));
        now = next;
      }
    }
  }
  if ((int)ret.size() < 2)
    return std::make_pair(0, 0);
  return simplify(abs(ret.at(0).first - ret.at(1).first),
                  abs(ret.at(0).second - ret.at(1).second));
}

void test_concave_decrease(long long N, bool print_detail = false) {
  std::function<double(long long)> f = [&](long long x) {
    return sqrt(N - x * x);
  };
  std::function<double(long long)> df = [&](long long x) {
    return -x / sqrt(N - x * x);
  };
  std::function<double(long long, long long)> inside =
      [&](long long x, long long y) { return x * x + y * y <= N; };

  std::vector<std::pair<long long, long long>> hull;
  long long x0 = 1;
  std::pair<long long, long long> first_slope =
      convex_hull(N, x0, f(x0), sqrt(N), f, df, inside, hull);
  if (print_detail) {
    for (const auto &i : hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", first_slope.second, first_slope.first);
  }

  std::vector<std::pair<long long, long long>> correct_hull;
  std::pair<long long, long long> correct_first_slope = brute_force_convex_hull(
      N, x0, f(x0), sqrt(N), f, df, inside, correct_hull);
  if (print_detail) {
    for (const auto &i : correct_hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", correct_first_slope.second,
           correct_first_slope.first);
  }

  assert(hull == correct_hull);
  assert(first_slope == correct_first_slope);

  printf("pass test_concave_decrease %lld\n", N);
}

void test_concave_increase(long long N, bool print_detail = false) {
  std::function<double(long long)> f = [&](long long x) {
    return sqrt(N - pow(sqrt(N) - x, 2));
  };
  std::function<double(long long)> df = [&](long long x) {
    return (sqrt(N) - x) / sqrt(N - pow(sqrt(N) - x, 2));
  };
  std::function<double(long long, long long)> inside = [&](long long x,
                                                           long long y) {
    return (sqrt(N) - x) * (sqrt(N) - x) + y * y <= N - eps;
  };

  std::vector<std::pair<long long, long long>> hull;
  long long x0 = 1;
  std::pair<long long, long long> first_slope =
      convex_hull(N, x0, f(x0), sqrt(N), f, df, inside, hull);
  if (print_detail) {
    for (const auto &i : hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", first_slope.second, first_slope.first);
  }

  std::vector<std::pair<long long, long long>> correct_hull;
  std::pair<long long, long long> correct_first_slope = brute_force_convex_hull(
      N, x0, f(x0), sqrt(N), f, df, inside, correct_hull);
  if (print_detail) {
    for (const auto &i : correct_hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", correct_first_slope.second,
           correct_first_slope.first);
  }

  assert(hull == correct_hull);
  assert(first_slope == correct_first_slope);

  printf("pass test_concave_increase %lld\n", N);
}

void test_convex_decrease(long long N, bool eq, bool print_detail = false) {
  std::function<double(long long)> f = [&](long long x) { return N * 1.0 / x; };
  std::function<double(long long)> df = [&](long long x) {
    return -N * 1.0 / x / x;
  };
  std::function<double(long long, long long)> outside =
      [&](long long x, long long y) { return eq ? x * y >= N : x * y > N; };

  std::vector<std::pair<long long, long long>> hull;
  long long x0 = 1;
  std::pair<long long, long long> first_slope =
      convex_hull(N, x0, f(x0) - eq, N + 1, f, df, outside, hull);
  if (print_detail) {
    for (const auto &i : hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", first_slope.second, first_slope.first);
  }

  std::vector<std::pair<long long, long long>> correct_hull;
  std::pair<long long, long long> correct_first_slope = brute_force_convex_hull(
      N, x0, f(x0) - eq, N + 1, f, df, outside, correct_hull);
  if (print_detail) {
    for (const auto &i : correct_hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", correct_first_slope.second,
           correct_first_slope.first);
  }

  assert(hull == correct_hull);
  assert(first_slope == correct_first_slope);

  printf("pass test_convex_decrease %lld %d\n", N, (int)eq);
}

void test_convex_increase(long long N, bool eq, bool print_detail = false) {
  std::function<double(long long)> f = [&](long long x) {
    return N * 1.0 / (N + 1 - x);
  };
  std::function<double(long long)> df = [&](long long x) {
    return N * 1.0 / (N + 1 - x) / (N + 1 - x);
  };
  std::function<double(long long, long long)> outside = [&](long long x,
                                                            long long y) {
    return eq ? ((N + 1 - x) * y >= N) : ((N + 1 - x) * y > N);
  };

  std::vector<std::pair<long long, long long>> hull;
  long long x0 = 1;
  std::pair<long long, long long> first_slope =
      convex_hull(N, x0, f(x0) - eq, N + 1, f, df, outside, hull);
  if (print_detail) {
    for (const auto &i : hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", first_slope.second, first_slope.first);
  }

  std::vector<std::pair<long long, long long>> correct_hull;
  std::pair<long long, long long> correct_first_slope = brute_force_convex_hull(
      N, x0, f(x0) - eq, N + 1, f, df, outside, correct_hull);
  if (print_detail) {
    for (const auto &i : correct_hull)
      printf("(%lld, %lld), ", i.first, i.second);
    printf("\n%lld/%lld\n", correct_first_slope.second,
           correct_first_slope.first);
  }

  assert(hull == correct_hull);
  assert(first_slope == correct_first_slope);

  printf("pass test_convex_increase %lld %d\n", N, (int)eq);
}

int main() {
  test_concave_decrease(100);
  test_concave_decrease(101);
  test_concave_decrease(1e9 + 7);

  test_concave_increase(100);
  test_concave_increase(101);
  test_concave_increase(1e9 + 7);

  test_convex_decrease(100, false);
  test_convex_decrease(101, false);
  test_convex_decrease(1e5 + 7, false);

  test_convex_decrease(100, true);
  test_convex_decrease(101, true);
  test_convex_decrease(1e5 + 7, true);

  test_convex_increase(100, false);
  test_convex_increase(101, false);
  test_convex_increase(1e6 + 7, false);

  test_convex_increase(100, true);
  test_convex_increase(101, true);
  test_convex_increase(1e6 + 7, true);

  return 0;
}
