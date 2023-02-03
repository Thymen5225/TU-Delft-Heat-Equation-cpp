#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

template<typename T>
class Vector {
public:
    Vector() : length(0), data(nullptr) {}

    Vector(const Vector &other) : length(other.length), data(new T[length]) {
        for (int i = 0; i < length; i++) {
            data[i] = other.data[i];
        }
    }

    Vector(Vector &&other) : length(other.length), data(other.data) {
        other.length = 0;
        other.data = nullptr;
    }

    Vector(int len) : length(len), data(new T[len]) {}

    Vector(std::initializer_list<T> list) : length(list.size()), data(new T[length]) {
        int i = 0;
        for (auto &item: list) {
            data[i++] = item;
        }
    }

    ~Vector() {
        delete[] data;
    }

    Vector &operator=(const Vector &other) {
        if (this != &other) {
            delete[] data;

            length = other.length;
            data = new T[length];
            for (int i = 0; i < length; i++) {
                data[i] = other.data[i];
            }
        }
        return *this;
    }

    Vector &operator=(Vector &&other) {
        if (this != &other) {
            delete[] data;

            length = other.length;
            data = other.data;
            other.length = 0;
            other.data = nullptr;
        }
        return *this;
    }

    T &operator[](int i) {
        return data[i];
    }

    const T &operator[](int i) const {
        return data[i];
    }

    Vector operator+(const Vector &other) const {
        if (length != other.length) {
            throw std::invalid_argument("Vectors have different lengths");
        }

        Vector result(length);
        for (int i = 0; i < length; i++) {
            result[i] = data[i] + other.data[i];
        }
        return result;
    }

    Vector operator-(const Vector &other) const {
        if (length != other.length) {
            throw std::invalid_argument("Vectors have different lengths");
        }

        Vector result(length);
        for (int i = 0; i < length; i++) {
            result[i] = data[i] - other.data[i];
        }
        return result;
    }

    template<typename S>
    Vector operator*(const S &scalar) const {
        Vector result(length);
        for (int i = 0; i < length; i++) {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    template<typename S>
    friend Vector operator*(const S &scalar, const Vector &vector) {
        return vector * scalar;
    }

    int len() const {
        return length;
    }


private:
    int length;
    T *data;
};

template<typename T, typename U>
typename std::common_type<T, U>::type
dot(const Vector<T> &lhs, const Vector<U> &rhs) {
    if (lhs.len() != rhs.len()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    typename std::common_type<T, U>::type result = 0;
    for (int i = 0; i < lhs.len(); i++) {
        result += lhs[i] * rhs[i];
    }

    return result;
}

template<typename T>
class Matrix {
public:
    Matrix(int rows, int cols) : num_rows(rows), num_cols(cols) {}

    ~Matrix() {}

    T &operator[](const std::pair<int, int> &ij) {
        return data[ij];
    }

    const T &operator()(const std::pair<int, int> &ij) const {
        auto it = data.find(ij);
        if (it == data.end())
            throw std::out_of_range("Matrix entry not present");
        return it->second;
    }

    int get_rows() const {
        return num_rows;
    }

    int get_cols() const {
        return num_cols;
    }

    bool has_entry(const std::pair<int, int> &ij) const {
        auto it = data.find(ij);
        if (it == data.end())
            return false;
        return true;
    }


private:
    std::map<std::pair<int, int>, T> data;
    int num_rows;
    int num_cols;
};

template<typename T, typename U>
Vector<typename std::common_type<T, U>::type>
operator*(const Matrix<T> &lhs,
          const Vector<U> &rhs) {
    int rows = lhs.get_rows();
    int cols = lhs.get_cols();
    int rhs_len = rhs.len();

    if (cols != rhs_len) {
        throw std::invalid_argument("Matrix and Vector dimensions are not compatible.");
    }

    using result_type = typename std::common_type<T, U>::type;
    Vector<result_type> result(rows);
    for (int i = 0; i < rows; i++) {
        result_type sum = 0;
        for (int j = 0; j < cols; j++) {
            double aa = lhs({i, j});
            sum += aa * rhs[j];
        }
        result[i] = sum;
    }
    return result;
}

template<typename T>
int cg(const Matrix<T> &A, const Vector<T> &b, Vector<T> &x, T tol = (T) 1e-8, int maxiter = 100) {
    Vector<T> p = b - A * x;
    Vector<T> r = p;

    Vector<T> p_new = p;
    Vector<T> r_new = r;
    Vector<T> x_new = x;

    int k = 0;

    for (k = 0; k < maxiter; k++){
        auto alpha = dot(r, r) / dot(p, p);
        x_new = x + alpha * p;
        r_new = r - alpha * A * p;

        if (dot(r_new, r_new) < tol * tol)
            break;

        auto beta = dot(r_new, r_new) / dot(r, r);
        p_new = r_new + beta * p;

        p = p_new;
        r = r_new;
        x = x_new;
    }

    if (k < maxiter)
        return k;
    if (k >= maxiter)
        return -1;
}

template <int n, typename T>
class Heat
{
public:
    Heat(T alpha, int m, T dt) : alpha(alpha), m(m), dt(dt) {
        M(m, m);
    }

private:
    T alpha;
    int m;
    T dt;
    Matrix<T> M;
};


int main(int argc, char *argv[]) {
    Vector<double> V1(5);
    Vector<double> V2 = {1.0, 2.0, 3.0};
    Vector<double> V3 = {3.0, 2.0, 1.0};
    Vector<double> V4 = V2 + V3;

    auto a = dot(V2, V3);

    std::cout << "Vector Testing..." << std::endl;
    std::cout << V4[2] << std::endl;
    std::cout << V2.len() << std::endl;
    std::cout << a << std::endl;

    Matrix<double> M1(3, 3);
    M1[{0, 0}] = 1.0; // set value at row 0, column 0 to 1.0
    M1[{1, 2}] = 2.0; // set value at row 1, column 2 to 2.0
    bool tt = M1.has_entry({0, 0});
    std::cout << tt << std::endl;

//
//    Vector<double> la(M1.get_cols());
//    for (int i = 0; i < M1.get_cols(); i++){
//        if (M1.has_entry({0, i}))
//            la[i] = M1({0,i});
//    }
//
//    std::cout << la[0] << std::endl;
//    std::cout << la[1] << std::endl;
//    std::cout << la[2] << std::endl;

//    Vector<double> V6 = M1 * V2;
//
//    std::cout << "Matrix testing..." << std::endl;
//    std::cout << M1({0,0}) << std::endl; // prints 1.0
//    std::cout << M1({1,2}) << std::endl; // throws an exception
//    std::cout << V6[3] << std::endl;
//


    return 0;
}






//    template<typename T1, typename U>
//    Vector<typename std::common_type<T1,U>::type>
//    operator*(const Matrix<T1>& lhs,
//              const Vector<U>& rhs)
//    {
//        int rows = lhs.getRows();
//        int cols = lhs.getCols();
//        int rhs_len = rhs.len();
//
//        if (cols != rhs_len) {
//            throw std::invalid_argument("Matrix and Vector dimensions are not compatible.");
//        }
//
//        using result_type = typename std::common_type<T1,U>::type;
//        Vector<result_type> result(rows);
//        for (int i = 0; i < rows; i++) {
//            result_type sum = 0;
//            for (int j = 0; j < cols; j++) {
//                sum += lhs({i, j}) * rhs[j];
//            }
//            result[i] = sum;
//        }
//        return result;
//    }


//    template<typename U>
//    friend Vector<typename std::common_type<T,U>::type>
//    operator*(const Matrix<T>& lhs, const Vector<U>& rhs)
//    {
//        if (lhs.num_cols != rhs.length)
//            throw std::invalid_argument("Matrix and Vector dimensions are not compatible");
//
//        Vector<typename std::common_type<T,U>::type> result(lhs.num_rows);
//        for (const auto& entry : lhs.data)
//        {
//            int i = entry.first.first;
//            int j = entry.first.second;
//            auto k = entry.second;
//
//            Vector<typename std::common_type<T,U>::type> new_vec(lhs.num_cols);
//            for (int t = 0; t < lhs.num_cols; t++){
//                if (lhs.has_entry({i, t}))
//                    new_vec[i] = lhs({i, t});
//            }
//
//            result[i] += k * new_vec[j];
//        }
//
//        return result;
//    }
