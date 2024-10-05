#include "bint.hpp"

int main() {
    bint a("90127890567834561234");
    bint b("56783456123490127890");
    bint c("127890567834561234");
    bint d("9999");
    bint res;

    a.print();
    b.print();
    c.print();
    d.print();

    puts("");

    res = a;
    res.add(b);
    res.print();

    res = a;
    res.sub(b);
    res.print();

    res = a;
    res.lmul(2);
    res.print();

    res = a;
    res.mul(b);
    res.print();

    res = a;
    res.mul_fft(b);
    res.print();

    res = c;
    res.ldiv(9999);
    res.print();

    res = c;
    res.div(d);
    res.print();
}