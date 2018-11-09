#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <complex>

using namespace std;

const complex<double> I(0, 1); // imaginary unit
const double PI = 3.14159265359;


template<class T>
struct duplet {
    int ind;
    T val;

    duplet(int p, T v) {
        ind = p;
        val = v;
    }


};


template<class T>
struct sparse_vec {
    double tol = 1e-6;
    vector<duplet<T> > duplets;
    int len = 0;

    sparse_vec(int l) {
        len = l;
    }

    void append(int ind, T val) {
        if (abs(val) >= tol) {
            this->duplets.push_back(duplet(ind, val));
        }

    }

    void cleanup() {

        if (duplets.size() == 0) {
            return;
        }


        for (int i = 0; i < duplets.size(); ++i) {
            if (duplets[i].ind > len - 1 || abs(duplets[i].val) <= tol) {
                duplets.erase(duplets.begin() + i);
            }
        }

        sort(duplets.begin(), duplets.end(), [](duplet<T> x, duplet<T> y) { return x.ind < y.ind; });


        for (int j = 0; j < duplets.size() - 1; ++j) {
            if (duplets[j].ind == duplets[j + 1].ind) {
                duplets[j].val += duplets[j + 1].val;
                duplets.erase(duplets.begin() + j + 1);
                j--;

            }
        }

        for (int k = 0; k < duplets.size(); ++k) {
            if (abs(real(duplets[k].val)) <= tol) {

                duplets[k].val = complex<double>(0, imag(duplets[k].val));
            }
            if (abs(imag(duplets[k].val)) <= tol) {
                duplets[k].val = complex<double>(real(duplets[k].val), 0);
            }

        }


    }

    T get_val(int ind) const {
//        BinarySearch
        long links = 0;
        long rechts = duplets.size() - 1;

        while (links <= rechts) {
            long mitte = links + (rechts - links) / 2;

            if (duplets[mitte].ind == ind) {
                return duplets[mitte].val;
            } else if (duplets[mitte].ind > ind) {
                rechts = mitte - 1;
            } else {
                links = mitte + 1;
            }

        }


        return 0;
    }

    static sparse_vec swap(const sparse_vec &a) {
        sparse_vec vec(a.len);
        for (int i = 0; i < a.duplets.size(); ++i) {
            double realteil = real(a.duplets[i].val);
            double imag2 = imag(a.duplets[i].val);
            int ind = a.duplets[i].ind;
            vec.append(ind, complex<double>(imag2, realteil));
        }
        vec.cleanup();
    }

    static sparse_vec cwise_mult(const sparse_vec &a, const sparse_vec &b) {
        sparse_vec out(max(a.len, b.len));

        int i = 0;
        int j = 0;
        while (i < a.len && j < b.len) {
            if (a.duplets[i].ind == b.duplets[j].ind) {
                out.append(a.duplets[i].ind, a.duplets[i].val * b.duplets[j].val);
                i++;
                j++;
            } else if (a.duplets[i].ind > b.duplets[j].ind) {
                j++;
            } else {
                i++;
            }
        }


        return out;
    }

    static sparse_vec conv(const sparse_vec &a, const sparse_vec &b) {
//        Ist schneller möglich...
        int size = a.len + b.len - 1;
        sparse_vec out(size);

        for (int k = 0; k < size; ++k) {
            T x = 0;
            for (int i = 0; i < k; ++i) {
                x += a.get_val(i) * b.get_val(k - i);
            }

            out.append(k, x);
        }
        out.cleanup();

        return out;
    }

    static sparse_vec fft(const sparse_vec &x) {
        int n = x.len;
        sparse_vec tot(n);

        for (int i = 0; i < n; ++i) {

            T a = 0;

            for (int j = 0; j < n; ++j) {
                auto num = (2 * PI * i * j) / n;
                auto real = cos(num);
                auto imag = sin(num);
                a += x.get_val(j) * (real - I * imag);


            }


            tot.append(i, a);
        }
        tot.cleanup();
        return tot;
    }

    static sparse_vec ifft(const sparse_vec &x) {
        double n = x.len;
        sparse_vec vec = swap(x);
        fft(vec);
        return swap(vec);

//
//        sparse_vec vec = x;
//
//        int left = 1;
//        int right = n;
//
//
//        for (int i = 1; i < n; ++i) {
//            vec.duplets[i].ind = n - vec.duplets[i].ind;
//        }
//        vec.cleanup();
//
//        vec = fft(vec);
//
//        for (int j = 0; j < n; ++j) {
//            vec.duplets[j].val /= 4;
//        }
//
//
//        return vec;
    }

    static sparse_vec conv_fft(sparse_vec a, sparse_vec b) {
        int size = a.len + b.len - 1;
        a.len = size;
        b.len = size;

        auto fftA = fft(a);
        print(fftA);
        auto fftB = fft(b);
        print(fftB);

        auto product = cwise_mult(fftA, fftB);
//        auto product = fftA * fftB;
        print(product);

        auto ifftAB = ifft(product);


        return ifftAB;
    }

    std::string to_string() const {
        std::stringstream ss;
        for (auto p : this->duplets) {
            ss << "(" << p.ind << "," << p.val << "),";
        }
        ss << "\n";
        std::string out = ss.str();
        return out;
    }


};

void print(sparse_vec<complex<double> > &x) {

    for (int i = 0; i < x.duplets.size(); ++i) {

        cout << x.duplets[i].ind << ": " << x.duplets[i].val.real() << ", " << x.duplets[i].val.imag() << endl;

    }
    cout << endl;
}

/***** TESTING ******/

int main() {
    sparse_vec<complex<double> > example(4);
    example.append(0, complex<double>(1, 0));
    example.append(1, complex<double>(2, -1));
    example.append(2, complex<double>(0, -1));
    example.append(3, complex<double>(-1, +2));
    example.cleanup();

    print(example);

//    sparse_vec<complex<double> > example2(4);
//    example2.append(0, complex<double>(1, 0));
//    example2.append(1, complex<double>(2, -1));
//    example2.append(2, complex<double>(0, -1));
//    example2.append(3, complex<double>(-1, +2));
//    print(example2);
//    auto result = sparse_vec<complex<double> >::conv_fft(example, example2);
//    print(result);


    sparse_vec<complex<double> > x(5);
    x.append(4, complex<double>(0, 4));
    x.append(0, complex<double>(8.2, 0));
    x.append(1, complex<double>(1, -2));
    x.append(3, complex<double>(-3, 4.66));
    x.cleanup();


    sparse_vec<complex<double> > y(4);
    y.append(1, complex<double>(5, 0));
    y.append(2, complex<double>(1.21, -4));
    y.append(3, complex<double>(4, 2.4));
    y.cleanup();

    auto m = sparse_vec<complex<double> >::cwise_mult(x, y);
    m.cleanup();
    cout << "TESTS. Correct componentwise multiplication between x and y: (1,(5,-10)),(3,(-23.184,11.44)),\n";
    cout << "cwise_mult(x,y) = " << m.to_string();

    auto c = sparse_vec<complex<double> >::conv(x, y);
    c.cleanup();
    cout
            << "Correct exact discrete convolution between x and y: (1,(41,0)),(2,(14.922,-42.8)),(3,(26.01,13.26)),(4,(-6.2,17.7)),(5,(15.01,37.6386)),(6,(-7.184,16.28)),(7,(-9.6,16)),\n";
    cout << "conv(x,y) = " << c.to_string();
    auto cf = sparse_vec<complex<double> >::conv_fft(x, y);
    cf.cleanup();
    cout << "conv_fft(x,y) = " << cf.to_string();
}


