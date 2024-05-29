#ifndef UTIL_H
#define UTIL_H

#include "log.h"

#ifdef WIN32
#undef min
#undef max
#endif

namespace backend {

using std::min;
using std::max;
using std::swap;

inline int sign(const float& x)
{
    return std::signbit(x) ? -1 : 1;
}

template<class T>
inline T sqr(const T& x)
{
    return x * x;
}

template<class T>
inline T cube(const T& x)
{
    return x * x * x;
}

template<class T>
inline T min(T a1, T a2, T a3)
{
    return min(a1, min(a2, a3));
}

template<class T>
inline T min(T a1, T a2, T a3, T a4)
{
    return min(min(a1, a2), min(a3, a4));
}

template<class T>
inline T second(T a1, T a2, T a3) {
    return a1 + a2 + a3 - max(a1, max(a2, a3)) - min(a1, min(a2, a3));
}

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5)
{
    return min(min(a1, a2), min(a3, a4), a5);
}

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return min(min(a1, a2), min(a3, a4), min(a5, a6));
}

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8)
{
    return min(min(a1, a2), min(a3, a4), min(a5, a6, a7, a8));
}

template<class T>
inline T max(T a1, T a2, T a3)
{
    return max(a1, max(a2, a3));
}

template<class T>
inline T max(T a1, T a2, T a3, T a4)
{
    return max(max(a1, a2), max(a3, a4));
}

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5)
{
    return max(max(a1, a2), max(a3, a4), a5);
}

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return max(max(a1, a2), max(a3, a4), max(a5, a6));
}

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8)
{
    return max(max(a1, a2), max(a3, a4), max(a5, a6, a7, a8));
}

template<class T>
inline void minmax(T a1, T a2, T& amin, T& amax)
{
    if (a1 < a2) {
        amin = a1;
        amax = a2;
    }
    else {
        amin = a2;
        amax = a1;
    }
}

template<class T>
inline void minmax(T a1, T a2, T a3, T& amin, T& amax)
{
    if (a1 < a2) {
        if (a1 < a3) {
            amin = a1;
            if (a2 < a3) amax = a3;
            else amax = a2;
        }
        else {
            amin = a3;
            if (a1 < a2) amax = a2;
            else amax = a1;
        }
    }
    else {
        if (a2 < a3) {
            amin = a2;
            if (a1 < a3) amax = a3;
            else amax = a1;
        }
        else {
            amin = a3;
            amax = a1;
        }
    }
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T& amin, T& amax)
{
    if (a1 < a2) {
        if (a3 < a4) {
            amin = min(a1, a3);
            amax = max(a2, a4);
        }
        else {
            amin = min(a1, a4);
            amax = max(a2, a3);
        }
    }
    else {
        if (a3 < a4) {
            amin = min(a2, a3);
            amax = max(a1, a4);
        }
        else {
            amin = min(a2, a4);
            amax = max(a1, a3);
        }
    }
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T a5, T& amin, T& amax)
{
    //@@@ the logic could be shortcircuited a lot!
    amin = min(a1, a2, a3, a4, a5);
    amax = max(a1, a2, a3, a4, a5);
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T a5, T a6, T& amin, T& amax)
{
    //@@@ the logic could be shortcircuited a lot!
    amin = min(a1, a2, a3, a4, a5, a6);
    amax = max(a1, a2, a3, a4, a5, a6);
}

template<class T>
inline void update_minmax(T a1, T& amin, T& amax)
{
    if (a1 < amin) amin = a1;
    else if (a1 > amax) amax = a1;
}

template<class T>
inline void sort(T& a, T& b, T& c)
{
    T temp;
    if (a < b) {
        if (a < c) {
            if (c < b) { // a<c<b
                temp = c; c = b; b = temp;
            } // else: a<b<c
        }
        else { // c<a<b
            temp = c; c = b; b = a; a = temp;
        }
    }
    else {
        if (b < c) {
            if (a < c) { //b<a<c
                temp = b; b = a; a = temp;
            }
            else { // b<c<a
                temp = b; b = c; c = a; a = temp;
            }
        }
        else { // c<b<a
            temp = c; c = a; a = temp;
        }
    }
}

template<class T>
inline T clamp(T a, T lower, T upper)
{
    if (a < lower) return lower;
    else if (a > upper) return upper;
    else return a;
}

// only makes sense with T=float or double
template<class T>
inline T smooth_step(T r)
{
    if (r < 0) return 0;
    else if (r > 1) return 1;
    return r * r * r * (10 + r * (-15 + r * 6));
}

// only makes sense with T=float or double
template<class T>
inline T smooth_step(T r, T r_lower, T r_upper, T value_lower, T value_upper)
{
    return value_lower + smooth_step((r - r_lower) / (r_upper - r_lower)) * (value_upper - value_lower);
}

// only makes sense with T=float or double
template<class T>
inline T ramp(T r)
{
    return smooth_step((r + 1) / 2) * 2 - 1;
}

#ifdef _MSC_VER

inline double remainder(double x, double y)
{
    return x - std::floor(x / y + 0.5) * y;
}
#endif

inline unsigned int round_up_to_power_of_two(unsigned int n)
{
    int exponent = 0;
    --n;
    while (n) {
        ++exponent;
        n >>= 1;
    }
    return 1 << exponent;
}

inline unsigned int round_down_to_power_of_two(unsigned int n)
{
    int exponent = 0;
    while (n > 1) {
        ++exponent;
        n >>= 1;
    }
    return 1 << exponent;
}

// Transforms even the sequence 0,1,2,3,... into reasonably good random numbers 
// Challenge: improve on this in speed and "randomness"!
// This seems to pass several statistical tests, and is a bijective map (of 32-bit unsigned ints)
inline unsigned int randhash(unsigned int seed)
{
    unsigned int i = (seed ^ 0xA3C59AC3u) * 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    return i;
}

// the inverse of randhash
inline unsigned int unhash(unsigned int h)
{
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= 0xA3C59AC3u;
    return h;
}

// returns repeatable stateless pseudo-random number in [0,1]
inline double randhashd(unsigned int seed)
{
    return randhash(seed) / (double)UINT_MAX;
}
inline float randhashf(unsigned int seed)
{
    return randhash(seed) / (float)UINT_MAX;
}

// returns repeatable stateless pseudo-random number in [a,b]
inline double randhashd(unsigned int seed, double a, double b)
{
    return (b - a) * randhash(seed) / (double)UINT_MAX + a;
}
inline float randhashf(unsigned int seed, float a, float b)
{
    return ((b - a) * randhash(seed) / (float)UINT_MAX + a);
}

inline int intlog2(int x)
{
    int exp = -1;
    while (x) {
        x >>= 1;
        ++exp;
    }
    return exp;
}

template<class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
{
    T s = std::floor(x);
    i = (int)s;
    if (i < i_low) {
        i = i_low;
        f = 0;
    }
    else if (i > i_high - 2) {
        i = i_high - 2;
        f = 1;
    }
    else
        f = (T)(x - s);
}

template<class S, class T>
inline S lerp(const S& value0, const S& value1, T f)
{
    return (1 - f) * value0 + f * value1;
}

template<class S, class T>
inline S bilerp(const S& v00, const S& v10,
    const S& v01, const S& v11,
    T fx, T fy)
{
    return lerp(lerp(v00, v10, fx),
        lerp(v01, v11, fx),
        fy);
}

template<class S, class T>
inline S trilerp(const S& v000, const S& v100,
    const S& v010, const S& v110,
    const S& v001, const S& v101,
    const S& v011, const S& v111,
    T fx, T fy, T fz)
{
    return lerp(bilerp(v000, v100, v010, v110, fx, fy),
        bilerp(v001, v101, v011, v111, fx, fy),
        fz);
}

template<class S, class T>
inline S quadlerp(const S& v0000, const S& v1000,
    const S& v0100, const S& v1100,
    const S& v0010, const S& v1010,
    const S& v0110, const S& v1110,
    const S& v0001, const S& v1001,
    const S& v0101, const S& v1101,
    const S& v0011, const S& v1011,
    const S& v0111, const S& v1111,
    T fx, T fy, T fz, T ft)
{
    return lerp(trilerp(v0000, v1000, v0100, v1100, v0010, v1010, v0110, v1110, fx, fy, fz),
        trilerp(v0001, v1001, v0101, v1101, v0011, v1011, v0111, v1111, fx, fy, fz),
        ft);
}

// f should be between 0 and 1, with f=0.5 corresponding to balanced weighting between w0 and w2
template<class T>
inline void quadratic_bspline_weights(T f, T& w0, T& w1, T& w2)
{
    w0 = T(0.5) * sqr(f - 1);
    w1 = T(0.75) - sqr(f - T(0.5));;
    w2 = T(0.5) * sqr(f);
}

// f should be between 0 and 1
template<class T>
inline void cubic_interp_weights(T f, T& wneg1, T& w0, T& w1, T& w2)
{
    T f2(f * f), f3(f2 * f);
    wneg1 = -T(1. / 3) * f + T(1. / 2) * f2 - T(1. / 6) * f3;
    w0 = 1 - f2 + T(1. / 2) * (f3 - f);
    w1 = f + T(1. / 2) * (f2 - f3);
    w2 = T(1. / 6) * (f3 - f);
}

template<class S, class T>
inline S cubic_interp(const S& value_neg1, const S& value0, const S& value1, const S& value2, T f)
{
    T wneg1, w0, w1, w2;
    cubic_interp_weights(f, wneg1, w0, w1, w2);
    return wneg1 * value_neg1 + w0 * value0 + w1 * value1 + w2 * value2;
}

template<class T>
void write_matlab(std::ostream& output, const std::vector<T>& a, const char* variable_name, bool column_vector = true, int significant_digits = 18)
{
    output << variable_name << "=[";
    std::streamsize old_precision = output.precision();
    output.precision(significant_digits);
    for (unsigned int i = 0; i < a.size(); ++i) {
        output << a[i] << " ";
    }
    output << "]";
    if (column_vector)
        output << "'";
    output << ";" << std::endl;
    output.precision(old_precision);
}

inline double interpolate_value(const double& i, const double& j, const MatrixXd& f) {
    // i, j >= 0 here    
    Assert(i >= 0 && j >= 0, "interpolate_value", "i and j should >= 0.");
    Assert(i < f.rows() && j < f.cols(), "interpolate_value", "out of index.");
    int i0 = static_cast<int>(i);
    int j0 = static_cast<int>(j);
    double a11 = (i - i0) * (j - j0);
    double a12 = (i0 + 1 - i) * (j - j0);
    double a21 = (i - i0) * (j0 + 1 - j);
    double a22 = (i0 + 1 - i) * (j0 + 1 - j);
    return f(i0, j0) * a22 + f(i0+1, j0) * a21 + f(i0, j0+1) * a12 + f(i0+1, j0+1) * a11;
}

inline double interpolate_value(const Vector2d& pos, const MatrixXd& f) {    
    return interpolate_value(pos(0), pos(1), f);
}

inline void interpolate_gradient(Vector2d& grad, const Vector2d& pos, const MatrixXd& f) {
    Assert(pos(0) >= 0 && pos(1) >= 0, "interpolate_value", "i and j should >= 0.");
    Assert(pos(0) < f.rows() && pos(1) < f.cols(), "interpolate_value", "out of index.");
    int i0 = static_cast<int>(pos(0));
    int j0 = static_cast<int>(pos(1));
    double grad_x1, grad_x2, grad_y1, grad_y2;
    grad_x1 = f(i0+1, j0) - f(i0, j0);
    grad_x2 = f(i0+1, j0+1) - f(i0, j0+1);
    grad_y1 = f(i0, j0+1) - f(i0, j0);
    grad_y2 = f(i0+1, j0+1) - f(i0+1, j0);
    grad(0) = grad_x1 * (j0+1-pos(1)) + grad_x2 * (pos(1)-j0);
    grad(1) = grad_y1 * (i0+1-pos(0)) + grad_y2 * (pos(0)-i0);
    return;
}

inline double computeFraction(const double& phi0, const double& phi1){
    if (phi0 <= 0 && phi1 <= 0) return 1;
    if (phi0 > 0 && phi1 > 0) return 0;
    if (phi0 <= 0) return -phi0/(phi1-phi0);
    if (phi1 <= 0) return -phi1/(phi0-phi1);
}

template<class T>
inline std::vector<T> union_phi(const std::vector<T>& phi1, const std::vector<T>& phi2) {
    if (phi1.size() != phi2.size()) {
        abort();
    }
    std::vector<T> ans = phi1;
    for (int i = 0; i < ans.size(); i++) {
        if (abs(phi1[i]) > abs(phi2[i])) {
            ans[i] = phi2[i];
        }
    }
    return ans;
}

class CleanupHelper {
public:
	CleanupHelper(const std::string& folderPath) : folderPath_(folderPath) {}

	~CleanupHelper() {
		try {
			for (const auto& entry : std::filesystem::directory_iterator(folderPath_)) {
				std::filesystem::remove(entry.path());
			}
		}
		catch (const std::filesystem::filesystem_error& err) {
			std::cerr << "Error during cleanup: " << err.what() << std::endl;
		}
	}

private:
	std::string folderPath_;
};

}

#endif