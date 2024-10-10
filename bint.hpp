#pragma once

#include <cstdio>
#include <vector>
#include <string>

typedef long long int64;

using List = std::vector<int64>;
using String = std::string;

class bint {
    const int digit = 8;
    int base = bint::pow(10, digit);
    int64 carry, borrow;
    int sign;
    List z;

    int64 pow(int64, int64);
    int cmpbint(const List&, const List&);
public:
    bint();
    bint(const int64&);
    bint(const String& num) { to_bint(num); }
    bint(const bint&);
    void to_bint(String);
    String to_string() const;

    bint& operator=(const int64&);
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
    void mul_karatsuba(const bint&);
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

    //for debugging
    void print();

    friend std::istream& operator>>(std::istream&, bint&);
    friend std::ostream& operator<<(std::ostream&, const bint&);
};

bint::bint() {
    z = {}, sign = 0;
}

bint::bint(const int64& num) {
    if (num > 0 && num <= (base - 1)) {
        sign = 1;
        z = {num};
    }
    else
        to_bint(std::to_string(num));
}

bint::bint(const bint& num) {
    z = {num.z}, sign = num.sign;
}

bint& bint::operator=(const int64& num) {
    if (num > 0 && num <= (base - 1)) {
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

int64 bint::pow(int64 b, int64 e) {
    int64 result = 1;
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

    int N = (static_cast<int>(str.size()) + digit - 1) / digit, i, j;
    z.resize(N);

    while (str.size() % digit != 0)
        str = '0' + str;

    for (i = 0; i < N; i++)
        for (j = 0; j < digit; j++)
            if ((i * digit + j) < static_cast<int>(str.size()))
                z[i] += (str[i * digit + j] - '0') * bint::pow(10, digit - 1 - j);
}

String bint::to_string() const {
    String value = "";
    for (int i = static_cast<int>(this->z.size()) - 1; i > 0; --i) {
        value = std::to_string(this->z[i]) + value;
        while (value.length() % digit != 0)
            value = '0' + value;
    }
    value = std::to_string(this->z[0]) + value;

    if (this->sign == -1)
        value = '-' + value;

    return value;
}

int bint::cmpbint(const List& a, const List& b) {
    int NA = static_cast<int>(a.size());
    int NB = static_cast<int>(b.size());
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
    int N = std::max(this->z.size(), num.z.size());
    int64 A, B;
    res.resize(N);
    carry = 0;

    for (int i = 0; i < N; ++i) {
        A = i < static_cast<int>(this->z.size()) ? this->z[this->z.size() - i - 1] : 0;
        B = i < static_cast<int>(num.z.size()) ? num.z[num.z.size() - i - 1] : 0;
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

    List a = this->z;
    List b = num.z;

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
    int N = std::max(a.size(), b.size());
    int64 A, B;
    res.resize(N);
    borrow = 0;

    for (int i = 0; i < N; ++i) {
        A = i < static_cast<int>(a.size()) ? a[a.size() - i - 1] : 0;
        B = i < static_cast<int>(b.size()) ? b[b.size() - i - 1] : 0;
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

    List res;
    res.resize(this->z.size());
    carry = 0;

    for (int i = static_cast<int>(this->z.size()) - 1; i >=0; i--) {
        int64 w = this->z[i];
        res[i] = (w * num + carry) % base;
        carry = (w * num + carry) / base;
    }

    if (carry > 0)
        res.insert(res.begin(), carry);
    
    z = std::move(res);
}

void bint::mul(const bint& num) {
    sign = 1;

    if (this->sign == -1 && num.sign == -1)
        sign = 1;
    else if (this->sign == -1 || num.sign == -1)
        sign = -1;

    List res;
    res.resize(this->z.size() + num.z.size() - 1);

    int i, j;

    for (j = 0; j < static_cast<int>(num.z.size()); j++)
        for (i = 0; i < static_cast<int>(this->z.size()); i++)
            res[j + i] += this->z[i] * num.z[j];

    carry = 0;
    for (int i = static_cast<int>(res.size()) - 1; i >= 0; i--) {
        res[i] += carry;
        if (res[i] < 0)
            carry = -(-(res[i] + 1) / base + 1);
        else
            carry = res[i] / base;
        res[i] -= carry * base;
    }

    if (carry > 0)
        res.insert(res.begin(), carry);

    while (!res.empty() && res[0] == 0)
        res.erase(res.begin());

    if (res.empty()) {
        z = {0};
        sign = 0;
        return;
    }

    z = std::move(res);
}

void bint::mul_karatsuba(const bint& num) {
    sign = 1;

    if (this->sign == -1 && num.sign == -1)
        sign = 1;
    else if (this->sign == -1 || num.sign == -1)
        sign = -1;

    int length = std::max(this->z.size(), num.z.size());

    List lhs_padded = this->z;
    List rhs_padded = num.z;
    lhs_padded.resize(length, 0);
    rhs_padded.resize(length, 0);

    if (this->z.size() <= 64 && num.z.size() <= 64) {
        bint res = *this;
        res.mul(num);
        z = std::move(res.z);
        return;
    }

    int half = (length + 1) / 2;
    bint lhs0, lhs1, rhs0, rhs1;
    lhs0.z.assign(lhs_padded.begin(), lhs_padded.begin() + half);
    lhs1.z.assign(lhs_padded.begin() + half, lhs_padded.end());
    rhs0.z.assign(rhs_padded.begin(), rhs_padded.begin() + half);
    rhs1.z.assign(rhs_padded.begin() + half, rhs_padded.end());

    bint p0 = lhs0;
    p0.mul_karatsuba(rhs0);
    bint p1 = lhs1;
    p1.mul_karatsuba(rhs1);
    bint p2 = lhs0 + lhs1;
    p2.mul_karatsuba(rhs0 + rhs1);
    bint p3 = p2 - (p0 + p1);

    p0.z.resize(p0.z.size() + ((length - half) << 1), 0);
    p3.z.resize(p3.z.size() + (length - half), 0);

    bint r = (p0 + p1) + p3;

    while (r.z.size() > 1 && r.z[0] == 0)
        r.z.erase(r.z.begin());

    z = std::move(r.z);
}

bint bint::operator*(const bint& num) const {
    bint res = *this;

    if (bint((base - 1)) >= num)
        res.lmul(num.z[0]);
    else if ((res.z.size() + num.z.size() - 1) <= 128)
        res.mul(num);
    else
        res.mul_karatsuba(num);

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

    List res;
    res.resize(this->z.size());
    int64 remainder = 0;

    for (int i = 0; i < static_cast<int>(this->z.size()); i++) {
        int64 w = this->z[i];
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

    if (bint((base - 1)) >= num)
        res.ldiv(num.z[0]);
    else
        res.div(num);

    return res;
}

bint bint::operator%(const bint& num) const {
    bint res, x, y;

    if (bint(2) == num)
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

//for debugging
void bint::print() {
    if (sign < 0)
        printf("-%lld", z[0]);
    else
        printf("%lld", z[0]);
    for (int i = 1; i < static_cast<int>(z.size()); i++)
        printf("%08lld", z[i]);
    printf("\n");
}

std::istream& operator>>(std::istream& in, bint& num) {
    std::string input;
    in >> input;
    num = bint(input);
    return in;
}

std::ostream& operator<<(std::ostream& out, const bint& num) {
    String value = num.to_string();
    out << value;
    return out;
}
