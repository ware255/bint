#pragma once

#include <cstdio>
#include <complex>
#include <vector>
#include <string>

using List = std::vector<int>;
using String = std::string;

using cd = std::complex<double>;
#define PI 3.14159265358979323846

class bint {
    int base = 10000;
    int carry, borrow;
    int sign;
    List z;

    int pow(int b, int e);
    void reverse(List& a);
    void fft(std::vector<cd>& a, bool invert);
    int cmpbint(const List& a, const List& b);
public:
    bint() : z({0}), sign(0) {}
    bint(const int& num) {
        if (num > 0 && num <= 9999) {
            sign = 1;
            z = {num};
        }
        else
            to_bint(std::to_string(num));
    }
    bint(const String& num) { to_bint(num); }
    bint(const bint& num) : z(num.z), sign(num.sign) {}
    void to_bint(String);
    friend String to_string(bint);

    bint& operator=(const int&);
    bint& operator=(const bint&);
    bint& operator=(const String&);

    bool operator<(const bint&);
    bool operator>(const bint&);
    bool operator<=(const bint&);
    bool operator>=(const bint&);
    bool operator==(const bint&);
    bool operator!=(const bint&);

    void add(const bint&);
    void sub(const bint&);
    void lmul(const int&);
    void mul(const bint&);
    void mul_fft(const bint&);
    void ldiv(const int&);
    void div(const bint&);

    void shr(const bint&);
    void shl(const bint&);

    bint operator+(const bint&) const;
    bint operator-(const bint&) const;
    bint operator*(const bint&) const;
    bint operator/(const bint&) const;
    bint operator%(const bint&) const;

    bint& operator++();
    bint& operator--();
    bint operator++(int);
    bint operator--(int);

    bint operator>>(const bint&) const;
    bint operator<<(const bint&) const;

    void doCarry(List&);
    void print();

    friend std::istream& operator>>(std::istream&, bint&);
    friend std::ostream& operator<<(std::ostream&, const bint&);
};

bint& bint::operator=(const int& num) {
    if (num > 0 && num <= 9999) {
        sign = 1;
        z = {num};
    }
    else
        to_bint(std::to_string(num));
    return *this;
}

bint& bint::operator=(const bint& num) {
    z = std::move(num.z);
    sign = num.sign;
    return *this;
}

bint& bint::operator=(const String& num) {
    to_bint(num);
    return *this;
}

int bint::pow(int b, int e) {
    int result = 1;
    while (e > 0) {
        if ((e & 1) == 1) {
            e -= 1;
            result *= b;
            if (e == 0) break;
        }
        e >>= 1;
        b *= b;
    }
    return result;
}

void bint::to_bint(String str) {
    if (str[0] == '-') {
        sign = -1;
        str.erase(str.begin());
    }
    else if (str == "0")
        sign = 0;
    else
        sign = 1;

    int N = (str.size() + 4 - 1) / 4, i, j;
    z.resize(N);
    for (int i = 0; i < z.size(); i++)
        z[i] = 0;

    while (str.size() & 3 != 0)
        str = '0' + str;

    for (i = 0; i < N; i++)
        for (j = 0; j < 4; j++)
            if ((i * 4 + j) < str.size())
                z[i] += (str[i * 4 + j] - '0') * bint::pow(10, 3 - j);
}

String to_string(bint num) {
    String value = "";
    for (int i = num.z.size() - 1; i > 0; --i) {
        value = std::to_string(num.z[i]) + value;
        while (value.length() & 3 != 0)
            value = '0' + value;
    }
    value = std::to_string(num.z[0]) + value;

    if (num.sign == -1)
        value = '-' + value;

    return value;
}

int bint::cmpbint(const List& a, const List& b) {
    int NA = a.size();
    int NB = b.size();
    if (NA > NB) return +1;
    if (NA < NB) return -1;
    for (int i = 0; i < NA; ++i) {
        if (a[i] > b[i]) return +1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

bool bint::operator<(const bint& num) {
    if (this->sign != num.sign)
        return this->sign < num.sign;

    if (this->sign == 0)
        return false;

    if (this->sign == 1)
        if (bint::cmpbint(this->z, num.z) < 0)
            return true;
        else
            return false;
    else
        if (bint::cmpbint(this->z, num.z) > 0)
            return true;
        else
            return false;
}

bool bint::operator>(const bint& num) {
    bint tmp(num);
    return tmp < *this;
}

bool bint::operator<=(const bint& num) {
    bint tmp(num);
    return !(tmp < *this);
}

bool bint::operator>=(const bint& num) {
    return !(*this < num);
}

bool bint::operator==(const bint& num) {
    if (this->sign != num.sign)
        return false;

    if (this->sign == 0 && num.sign == 0)
        return true;

    return bint::cmpbint(this->z, num.z) == 0;
}

bool bint::operator!=(const bint& num) {
    return !(*this == num);
}

void bint::add(const bint& num) {
    const List& a = this->z;
    const List& b = num.z;

    if (this->sign == 1 && num.sign == -1) {
        bint res;
        if (*this >= num) {
            sign = 1;
            res = *this;
            res.sub(num);
        }
        else {
            sign = -1;
            res = num;
            res.sub(*this);
        }
        z = std::move(res.z);
        return;
    }

    if (sign == -1 && num.sign == 1) {
        bint res;
        if (*this <= num) {
            sign = 1;
            res = num;
            res.sub(*this);
        }
        else {
            sign = -1;
            res = *this;
            res.sub(num);
        }
        z = std::move(res.z);
        return;
    }

    List res;
    int N = std::max(a.size(), b.size());
    int A, B;
    res.resize(N);
    carry = 0;

    for (int i = 0; i < N; ++i) {
        A = i < a.size() ? a[a.size() - i - 1] : 0;
        B = i < b.size() ? b[b.size() - i - 1] : 0;
        res[N - i - 1] = A + B + carry;

        if (res[N - i - 1] < base)
            carry = 0;
        else {
            res[N - i - 1] -= base;
            carry = 1;
        }
    }

    if (carry > 0)
        res.insert(res.begin(), carry);

    if (this->sign == -1 && num.sign == -1)
        sign = -1;
    else
        sign = 1;

    z = std::move(res);
}

bint bint::operator+(const bint& num) const {
    bint res = *this;
    res.add(num);
    return res;
}

void bint::sub(const bint& num) {
    List a = this->z;
    List b = num.z;

    if (this->sign == 1 && num.sign == -1) {
        sign = 1;
        bint temp(num), res(*this);
        temp.sign = 1;
        res.add(temp);
        z = std::move(res.z);
        return;
    }

    if (this->sign == -1 && num.sign == 1) {
        sign = -1;
        bint temp(*this);
        temp.sign = 1;
        temp.add(num);
        z = std::move(temp.z);
        return;
    }

    if (this->sign == 1 && num.sign == 1) {
        if (*this > num)
            sign = 1;
        else {
            sign = -1;
            a = num.z;
            b = this->z;
        }
    }

    if (this->sign == -1 && num.sign == -1) {
        if (*this > num)
            sign = -1;
        else {
            sign = 1;
            a = num.z;
            b = this->z;
        }
    }

    List res;
    int N = a.size() > b.size() ? a.size() : b.size();
    int A, B;
    res.resize(N);
    borrow = 0;

    for (int i = 0; i < N; ++i) {
        A = i < a.size() ? a[a.size() - i - 1] : 0;
        B = i < b.size() ? b[b.size() - i - 1] : 0;
        res[N - i - 1] = A - B - borrow;

        if (res[N - i - 1] >= 0)
            borrow = 0;
        else {
            res[N - i - 1] += base;
            borrow = 1;
        }
    }

    while (!res.empty() && res[0] == 0)
        res.erase(res.begin());

    if (res.empty()) {
        z = {0};
        sign = 0;
        return;
    }

    z = std::move(res);
}

bint bint::operator-(const bint& num) const {
    bint res = *this;
    res.sub(num);
    return res;
}

void bint::lmul(const int& num) {
    if (num == 0) {
        sign = 0;
        z = {0};
        return;
    }

    sign = 1;

    if (this->sign == -1 && num < 0)
        sign = 1;
    else if (this->sign == -1 || num < 0)
        sign = -1;

    const List& a = this->z;

    List res;
    res.resize(a.size());
    carry = 0;

    for (int i = a.size() - 1; i >=0; i--) {
        int w = a[i];
        res[i] = (w * num + carry) % base;
        carry = (w * num + carry) / base;
    }

    if (carry > 0)
        res.insert(res.begin(), carry);
    
    z = std::move(res);
}

void bint::mul(const bint& num) {
    const List& a = this->z;
    const List& b = num.z;

    sign = 1;

    if (this->sign == -1 && num.sign == -1)
        sign = 1;
    else if (this->sign == -1 || num.sign == -1)
        sign = -1;

    List res;
    res.resize(a.size() + b.size() - 1);

    int i, j;

    for(i = 0; i < res.size(); i++)
        res[i] = 0;

    for (j = 0; j < b.size(); j++)
        for (i = 0; i < a.size(); i++)
            res[j + i] += a[i] * b[j];

    doCarry(res);

    while (!res.empty() && res[0] == 0)
        res.erase(res.begin());

    if (res.empty()) {
        z = {0};
        sign = 0;
        return;
    }

    z = std::move(res);
}

void bint::fft(std::vector<cd>& a, bool invert) {
    int n = a.size(), i, j;

    for (i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (i = 0; i < n; i += len) {
            cd w(1);
            for (j = 0; j < (len >> 1); j++) {
                cd u = a[i+j], v = a[i+j+(len>>1)] * w;
                a[i+j] = u + v;
                a[i+j+(len>>1)] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert)
        for (cd & x : a) x /= n;
}

void bint::reverse(List& a) {
    for (int i = 0; i < a.size() >> 1; i++) {
        int temp = a[i];
        a[i] = a[a.size() - 1 - i];
        a[a.size() - 1 - i] = temp;
    }
}

void bint::mul_fft(const bint& num) {
    List x, y;
    int i;

    sign = 1;

    if (this->sign == -1 && num.sign == -1)
        sign = 1;
    else if (this->sign == -1 || num.sign == -1)
        sign = -1;

    String str_a = to_string(*this);
    for (i = 0; i < str_a.size(); ++i)
        x.emplace_back(str_a[str_a.size() - i - 1] - '0');

    String str_b = to_string(num);
    for (i = 0; i < str_b.size(); ++i)
        y.emplace_back(str_b[str_b.size() - i - 1] - '0');

    std::vector<cd> fa(x.begin(), x.end()), fb(y.begin(), y.end());

    int n = 1;
    while (n < x.size() + y.size()) 
        n <<= 1;

    fa.resize(n);
    fb.resize(n);

    bint::fft(fa, false);
    bint::fft(fb, false);

    for (i = 0; i < n; i++)
        fa[i] *= fb[i];

    bint::fft(fa, true);

    List result(n);
    for (i = 0; i < n; i++)
        result[i] = std::round(fa[i].real());

    carry = 0;
    for (i = 0; i < n; i++) {
        result[i] += carry;
        carry = std::floor(result[i] / 10);
        result[i] %= 10;
    }

    bint::reverse(result);

    while (!result.empty() && result[0] == 0)
        result.erase(result.begin());
    
    String str;
    for (i = 0; i < result.size(); i++)
        str += result[i] + '0';

    bint::to_bint(str);
}

bint bint::operator*(const bint& num) const {
    bint res = *this;

    if (bint(9999) >= num)
        res.lmul(num.z[0]);
    else if ((this->z.size() + num.z.size()) < 32)
        res.mul(num);
    else
        res.mul_fft(num);

    return res;
}

void bint::ldiv(const int& num) {
    if (num == 0)
        throw std::invalid_argument("Division by zero");

    sign = 1;

    if (this->sign == -1 && num < 0)
        sign = 1;
    else if (this->sign == -1 || num < 0)
        sign = -1;

    const List& a = this->z;

    List res;
    res.resize(a.size());
    int remainder = 0;

    for (int i = 0; i < a.size(); i++) {
        int w = a[i];
        res[i] = (w + remainder) / num;
        remainder = ((w + remainder) % num) * base;
    }

    while (!res.empty() && res[0] == 0)
        res.erase(res.begin());

    if (res.empty()) {
        z = {0};
        sign = 0;
        return;
    }

    z = std::move(res);
}

void bint::div(const bint& num) {
    if (num.z.empty() || (num.z.size() == 1 && num.z[0] == 0))
        throw std::invalid_argument("Division by zero");

    sign = 1;

    if (this->sign == -1 && num.sign == -1)
        sign = 1;
    else if (this->sign == -1 || num.sign == -1)
        sign = -1;

    if (*this < num) {
        z = {0};
        sign = 0;
        return;
    }

    bint rem(*this);
    bint m(num), s(1), quotient(0);

    while (rem > m) {
        m.lmul(2);
        s.lmul(2);
    }

    while (rem >= num) {
        while (rem < m) {
            m.ldiv(2);
            s.ldiv(2);
        }
        rem.sub(m);
        quotient.add(s);
    }

    z = std::move(quotient.z);
}

bint bint::operator/(const bint& num) const {
    bint res = *this;

    if (bint(9999) >= num)
        res.ldiv(num.z[0]);
    else
        res.div(num);

    return res;
}

bint bint::operator%(const bint& num) const {
    bint res, x, y;

    x = num;
    if (x == 2)
        return this->z[this->z.size() - 1] & 1;

    x = *this / num;
    y = x * num;
    res = *this - y;

    return res;
}

bint& bint::operator++() {
    *this = *this + 1;
    return *this;
}

bint& bint::operator--() {
    *this = *this - 1;
    return *this;
}

bint bint::operator++(int) {
    bint temp = *this;
    *this = *this + 1;
    return temp;
}

bint bint::operator--(int) {
    bint temp = *this;
    *this = *this - 1;
    return temp;
}

void bint::shr(const bint& num) {
    bint res = *this;
    for (bint i = 0; i < num; i++)
        res.ldiv(2);
    sign = res.sign;
    z = std::move(res.z);
}

bint bint::operator>>(const bint& num) const {
    bint res = *this;
    res.shr(num);
    return res;
}

void bint::shl(const bint& num) {
    bint res = *this;
    for (bint i = 0; i < num; i++)
        res.lmul(2);
    sign = res.sign;
    z = std::move(res.z);
}

bint bint::operator<<(const bint& num) const {
    bint res = *this;
    res.shl(num);
    return res;
}

void bint::doCarry(List& a) {
    carry = 0;
    for (int i = a.size() - 1; i >= 0; i--) {
        a[i] += carry;
        if (a[i] < 0)
            carry = -(-(a[i] + 1) / base + 1);
        else
            carry = a[i] / base;
        a[i] -= carry * base;
    }

    if (carry > 0)
        a.insert(a.begin(), carry);
}

void bint::print() {
    if (sign < 0)
        printf("-%d", z[0]);
    else
        printf("%d", z[0]);
    for (int i = 1; i < z.size(); i++)
        printf("%04d", z[i]);
    printf("\n");
}

std::istream& operator>>(std::istream& in, bint& num) {
    std::string input;
    in >> input;
    num = bint(input);
    return in;
}

std::ostream& operator<<(std::ostream& out, const bint& num) {
    String value = to_string(num);
    out << value;
    return out;
}
