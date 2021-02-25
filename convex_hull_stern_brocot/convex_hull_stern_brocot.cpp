#include <cstdio>
#include <cmath>
#include <vector>
#include <cassert>
#include <utility>
#include <functional>

using namespace std;
using i64 = int64_t;

const double eps=1e-9;
std::pair<long long,long long>
    convex_hull(const long long &N,
        const long long &x1,const long long &y1,
        const long long &x2,
        const std::function<double(long long)> &f,
        const std::function<double(long long)> &df,
        const std::function<bool(long long,long long)> &inside,
        std::vector<std::pair<long long,long long>> &ret)
{
    // x1: from (inclusive)  x2: to (exclusive)
    // Note that for concave function, the inside function means the points on f should be included
    // but for convex function, the inside function means the outside
    // return:
    //   ret: the vertices of the convex hull under the function f
    //   return value: the slope of the convex hull (x,y) at x1

    // Assumptions:
    // df(x) is monotonic
    // df(x) has the same sign on x1<=x<x2
    // (x1, y1) should be on the convex hull (I don't check it)

    assert(x1<x2);
    const bool convex=df(x1)<df(x1+1);
    if (!convex && inside(x1,y1) || convex && !inside(x1,y1)) ret.emplace_back(std::make_pair(x1,y1));
    if (x1+1==x2) return std::make_pair(0,0);

    const int sign=(df(x1)<0)?-1:1;
    std::vector<std::pair<long long,long long>> stac;  // (x, y)
    std::pair<long long,long long> first_slope;
    if (!convex)
    {
        long long x=x1,y=y1;
        if (sign==-1)
        {
            stac.emplace_back(0,1);
            if ((long long)f(x1)==(long long)f(x1+1))
            {
                first_slope=std::make_pair(0,1);
                stac.emplace_back(1,0);
            }
            else
            {
                std::pair<long long,long long> left=std::make_pair(1,0),right=std::make_pair(0,1);
                while (true)
                {
                    std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                    if (inside(x+mid.first,y+sign*mid.second))
                    {
                        right=mid;
                        stac.emplace_back(mid);
                    }
                    else
                    {
                        if (x+mid.first>=x2) break;
                        if (sign*df(x+mid.first)*right.first>=right.second) break;
                        left=mid;
                    }
                }
                /*printf("my: ");
                for (const auto &i:stac)
                    printf("%lld/%lld, ",i.second,i.first);
                printf("\n");*/
                first_slope=stac.back();
            }
            while (true)
            {
                std::pair<long long,long long> left,right;
                right=stac.back();
                stac.pop_back();
                bool assigned=false;
                std::pair<long long,long long> last;
                while (inside(x+right.first,y+sign*right.second))
                {
                    x+=right.first;
                    y+=sign*right.second;
                    if (x<x2)
                    {
                        assigned=true;
                        last=std::make_pair(x,y);
                    }
                }
                if (assigned) ret.emplace_back(last);
                left=right;
                while (!stac.empty())
                {
                    right=stac.back();
                    if (inside(x+right.first,y+sign*right.second)) break;
                    stac.pop_back();
                    left=right;
                }
                if (stac.empty()) break;
                while (true)
                {
                    std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                    if (inside(x+mid.first,y+sign*mid.second))
                    {
                        right=mid;
                        stac.emplace_back(mid);
                    }
                    else
                    {
                        if (x+mid.first>=x2) break;
                        if (sign*df(x+mid.first)*right.first>=right.second) break;
                        left=mid;
                    }
                }
            }
            return first_slope;
        }
        else
        {
            stac.emplace_back(1,0);
            std::pair<long long,long long> left=std::make_pair(1,0),right=std::make_pair(0,1);
            while (true)
            {
                std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                if (inside(x+mid.first,y+sign*mid.second))
                {
                    left=mid;
                    stac.emplace_back(mid);
                }
                else
                {
                    if (x+mid.first>=x2) break;
                    if (sign*df(x+mid.first)*right.first>=right.second) break;
                    right=mid;
                }
            }
            /*printf("my: ");
            for (const auto &i:stac)
                printf("%lld/%lld, ",i.second,i.first);
            printf("\n");*/
            first_slope=stac.back();
            while (true)
            {
                std::pair<long long,long long> left,right;
                left=stac.back();
                stac.pop_back();
                bool assigned=false;
                std::pair<long long,long long> last;
                while (inside(x+left.first,y+sign*left.second))
                {
                    x+=left.first;
                    y+=sign*left.second;
                    if (x<x2)
                    {
                        assigned=true;
                        last=std::make_pair(x,y);
                    }
                }
                if (assigned) ret.emplace_back(last);
                right=left;
                while (!stac.empty())
                {
                    left=stac.back();
                    if (inside(x+left.first,y+sign*left.second)) break;
                    stac.pop_back();
                    right=left;
                }
                if (stac.empty()) break;
                while (true)
                {
                    std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                    if (inside(x+mid.first,y+sign*mid.second))
                    {
                        left=mid;
                        stac.emplace_back(mid);
                    }
                    else
                    {
                        if (x+mid.first>=x2) break;
                        double v=sign*df(x+mid.first);
                        if (sign*df(x+mid.first)*right.first>=right.second) break;
                        right=mid;
                    }
                }
            }
            return first_slope;
        }
    }
    else
    {
        const auto &outside=inside;
        if (sign==-1)
        {
            long long x=x1,y=y1+1;
            stac.emplace_back(1,0);
            std::pair<long long,long long> left=std::make_pair(0,1),right=std::make_pair(1,0);
            while (true)
            {
                std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                if (outside(x+mid.first,y+sign*mid.second)) // take care! It should be equivalent to f(x+mid.first)<=y+sign*mid.second
                {
                    right=mid;
                    stac.emplace_back(mid);
                }
                else
                {
                    if (x+mid.first>=x2) break;
                    if (sign*df(x+mid.first)*right.first<=right.second) break;
                    left=mid;
                }
            }
            printf("my: ");
            for (const auto &i:stac)
                printf("%lld/%lld, ",i.second,i.first);
            printf("\n");
            first_slope=stac.back();
            while (true)
            {
                std::pair<long long,long long> left,right;
                right=stac.back();
                stac.pop_back();
                bool assigned=false;
                std::pair<long long,long long> last;
                while (outside(x+right.first,y+sign*right.second))
                {
                    x+=right.first;
                    y+=sign*right.second;
                    if (x<x2)
                    {
                        assigned=true;
                        last=std::make_pair(x,y-1);
                    }
                    else break;
                }
                if (assigned) ret.emplace_back(last);
                left=right;
                while (!stac.empty())
                {
                    right=stac.back();
                    if (outside(x+right.first,y+sign*right.second)) break;
                    left=right;
                    stac.pop_back();
                }
                if (stac.empty()) break;
                while (true)
                {
                    std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                    if (outside(x+mid.first,y+sign*mid.second))
                    {
                        right=mid;
                        stac.emplace_back(mid);
                    }
                    else
                    {
                        if (x+mid.first>=x2) break;
                        if (sign*df(x+mid.first)*right.first<=right.second) break;
                        left=mid;
                    }
                }
            }
            return first_slope;
        }
        else
        {
            long long x=x1,y=ceil(f(x1));
            if (outside(x,y)) ++y;
            stac.emplace_back(0,1);
            long long y2=f(x1+1);
            while (outside(x1+1,y2)) y2++;
            if (y==y2)
            {
                first_slope=std::make_pair(1,0);
                stac.emplace_back(1,0);
            }
            else
            {
                std::pair<long long,long long> left=std::make_pair(0,1),right=std::make_pair(1,0);
                while (true)
                {
                    std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                    if (outside(x+mid.first,y+sign*mid.second)) // take care! It should be equivalent to f(x+mid.first)<=y+sign*mid.second
                    {
                        left=mid;
                        stac.emplace_back(mid);
                    }
                    else
                    {
                        if (x+mid.first>=x2) break;
                        if (sign*df(x+mid.first)*right.first<=right.second) break;
                        right=mid;
                    }
                }
                printf("my: ");
                for (const auto &i:stac)
                    printf("%lld/%lld, ",i.second,i.first);
                printf("\n");
                first_slope=stac.back();
            }
            while (true)
            {
                std::pair<long long,long long> left,right;
                left=stac.back();
                if (!left.first) break;
                stac.pop_back();
                bool assigned=false;
                std::pair<long long,long long> last;
                while (outside(x+left.first,y+sign*left.second))
                {
                    x+=left.first;
                    y+=sign*left.second;
                    if (x<x2)
                    {
                        assigned=true;
                        last=std::make_pair(x,y-1);
                    }
                    else break;
                }
                if (assigned) ret.emplace_back(last);
                right=left;
                while (!stac.empty())
                {
                    left=stac.back();
                    if (outside(x+left.first,y+sign*left.second)) break;
                    right=left;
                    stac.pop_back();
                }
                if (stac.empty()) break;
                while (true)
                {
                    std::pair<long long,long long> mid=std::make_pair(left.first+right.first,left.second+right.second);
                    if (outside(x+mid.first,y+sign*mid.second)) // take care! It should be equivalent to f(x+mid.first)<=y+sign*mid.second
                    {
                        left=mid;
                        stac.emplace_back(mid);
                    }
                    else
                    {
                        if (x+mid.first>=x2) break;
                        if (sign*df(x+mid.first)*right.first<=right.second) break;
                        right=mid;
                    }
                }
            }
            return first_slope;
        }
    }
}

std::pair<long long,long long>
    brute_force(const long long &N,
        const long long &x1,const long long &y1,
        const long long &x2,
        const std::function<bool(long long,long long)> &inside,
        const bool &include,std::vector<std::pair<long long,long long>> &ret)
{
    //for (const auto &i)
}

int main() {
  i64 N = 10;

    /*std::function<double(long long)> f = [&](long long x) {return sqrt(N-x*x);};
    std::function<double(long long)> df = [&](long long x) {return -x/sqrt(N-x*x));};
    std::function<double(long long,long long)> inside = [&] (i64 x, i64 y) {return x*x+y*y<=N;};*/

    /*std::function<double(long long)> f = [&](long long x) {return sqrt(N-pow(sqrt(N)-x,2));};
    std::function<double(long long)> df = [&](long long x) {return (sqrt(N)-x)/sqrt(N-pow(sqrt(N)-x,2));};
    std::function<double(long long,long long)> inside = [&] (i64 x, i64 y) {return (sqrt(N)-x)*(sqrt(N)-x) + y * y < N-eps;};*/

    /*std::function<double(long long)> f = [&](long long x) {return N*1.0/x;};
    std::function<double(long long)> df = [&](long long x) {return -N*1.0/x/x;};
    std::function<double(long long,long long)> outside = [&] (i64 x, i64 y) {return x*y>=N;};*/

    std::function<double(long long)> f = [&](long long x) {return N*1.0/(N+1-x);};
    std::function<double(long long)> df = [&](long long x) {return N*1.0/(N+1-x)/(N+1-x);};
    std::function<double(long long,long long)> outside = [&] (i64 x, i64 y) {return (N+1-x)*y>N;};

    std::vector<std::pair<long long,long long>> ret;
    long long x0=1;
    std::pair<long long,long long> first_slope=convex_hull(N,x0,f(x0),N+1,f,df,outside,ret);
    for (const auto &i:ret)
        printf("(%lld, %lld)\n",i.first,i.second);
    printf("%lld/%lld\n",first_slope.second,first_slope.first);
  return 0;
}
