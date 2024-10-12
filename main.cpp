//
// Example Code
//
#include <iostream>
#include <array>
#include <algorithm>
#include <functional>
#include <random>

#include "bint.hpp"

bint pow(bint b, int e) {
    bint result = 1;
    while (e > 0) {
        if ((e & 1) == 1) {
            e -= 1;
            result = result * b;
            if (e == 0) break;
        }
        e >>= 1;
        b = b * b;
    }
    return result;
}

bint GetRandNum(int length = 10) {
    std::array<
        std::seed_seq::result_type,
        std::mt19937::state_size
    > seed_data;

    std::random_device seed_gen;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(seed_gen));

    std::seed_seq seq(seed_data.begin(), seed_data.end());

    std::mt19937_64 engine(seq);

    bint n;
    for (int i = 0; i < length; i++)
        n = n + bint(engine() % 10) * pow(bint(10), i);
    return n;
}

bint gcd(bint x, bint y) {
    while (1) {
        if (y == 0) return x;
        x = x % y;
        if (x == 0) return y;
        y = y % x;
    }
}

const size_t k = 3;

int max(const int& a, const int& b) { return a >= b ? a : b; }

bint exp_sliding_window(const bint& x, bint y, const bint& mod) {
    std::vector<bint> table(1 << k, 0);
    table[1] = x;

    bint aa = (x * x) % mod;
    for (int i = 1; i < (1 << (k - 1)); ++i)
        table[(i << 1) + 1] = (table[(i << 1) - 1] * aa) % mod;

    std::string bin_exp;
    while (y > 0) {
        if (y % 2 == 1)
            bin_exp += '1';
        else
            bin_exp += '0';
        y = y / 2;
    }
    size_t bit_length = bin_exp.length();

    bint z = 1;

    int i = bit_length - 1;
    while (i >= 0) {
        if (bin_exp[i] == '0') {
            z = (z * z) % mod;
            i -= 1;
        }
        else {
            int s = max(i - k + 1, 0);

            while (s <= i && bin_exp[s] == '0')
                s += 1;

            for (int j = 0; j < i - s + 1; ++j)
                z = (z * z) % mod;

            std::string window_bits = bin_exp.substr(s, i - s + 1);
            std::reverse(window_bits.begin(), window_bits.end());
            int window = std::stoi(window_bits, nullptr, 2);
            z = (z * table[window]) % mod;
            i = s - 1;
        }
    }

    return z;
}

std::vector<bool> IsPrime;

void sieve(size_t max) {
    if (max + 1 > IsPrime.size())
        IsPrime.resize(max+1, true);

    IsPrime[0] = false;
    IsPrime[1] = false;

    for (size_t i = 2; i * i <= max; ++i)
        if (IsPrime[i])
            for (size_t j = 2; i * j <= max; ++j)
                IsPrime[i*j] = false;
}

bool Fermat_Test(bint num) {
    int a[4] = { 2, 3, 5, 7};

    for (int i = 0; i < 4; ++i) {
        if (gcd(num, a[i]) != 1)
            return false;
        if (exp_sliding_window(a[i], num - 1, num) != 1)
            return false;
    }

    return true;
}

bint prime_gen() {
    bint num = GetRandNum(32);

    if (num % 2 == 0)
        num = num + 1;

    constexpr int N = 1000;
    sieve(N);

    while (1) {
        for (int i = 3; i < N; i+=2) {
            if (IsPrime[i])
                if (num % i == 0) {
                    num = num + 2;
                    continue;
                }
        }

        if (Fermat_Test(num))
            break;

        num = num + 2;
    }

    return num;
}

bint factorial(unsigned int n) {
    bint res = 1;
    for (unsigned int i = 2; i <= n; i++)
        res = res * i;
    return res;
}

int main() {
    bint a = GetRandNum(32), b = GetRandNum(32);
    bint c = GetRandNum(32), d = 9999;
    
    std::cout << a << " + " << b << " = " << a + b << std::endl;
    std::cout << a << " - " << b << " = " << a - b << std::endl;
    std::cout << c << " * " << d << " = " << c * d << std::endl;
    std::cout << c << " / " << d << " = " << c / d << std::endl;
    std::cout << c << " % " << d << " = " << c % d << std::endl;

    std::cout << "4545! = " << factorial(4545) << std::endl;

    std::cout << "prime number = " << prime_gen() << std::endl;

    return 0;
}
