#include <iostream>
#include "bint.hpp"

int main() {
    bint a("90127890567834561234");
    bint b("56783456123490127890");
    bint c("127890567834561234");
    bint d("9999");

    std::cout << a << " + " << b << " = " << a + b << std::endl;
    std::cout << a << " - " << b << " = " << a - b << std::endl;
    std::cout << c << " * " << d << " = " << c * d << std::endl;
    std::cout << c << " / " << d << " = " << c / d << std::endl;
    std::cout << c << " % " << d << " = " << c % d << std::endl;

    return 0;
}
