// File: P1.cpp
// Date: 4/17/2021
// Author: Brian Hilnbrand

#include <iostream>
#include <iomanip>
#include <vector>
#include <bits/stdc++.h>
#include <math.h>
#include <cmath>

// Calculate greatest common denominator
long long gcd(long long a, long long b) {
   if (b == 0)
   return a;
   return gcd(b, a % b);
}

// Number class contains both fraction and double representation of number
class Number {
    // long long num, den;
    double value;
    // bool isDouble;
public:
    long long num, den;
    bool isDouble;

    // Construct fraction
    Number(long long n, long long d) {
        num = n;
        den = d;
        value = toDouble();
        isDouble = false;
    }

    // Construct double
    Number(double v) {
        value = v;
        isDouble = true;
    }

    // Reduce fraction
    void reduce() {
        if (num == 0) {
            den = 1;
        }
        else if(num == den) {
            num = den = 1;
        }
        else {
            long long gcdVal = gcd(num, den);
            num /= gcdVal;
            den /= gcdVal;
        }
    }

    // Use double instead of fraction. This will convert the fraction to a double.
    void useDouble() {
        if (!isDouble) {
            value = toDouble();
            isDouble = true;
            num = 1;
            den = 1;
        }
    }

    // Overloaded multiply operator to another Number object
    Number operator*(Number const& rhs) {
        Number result = *this;
        if (result.isDouble || rhs.isDouble) {
            result.useDouble();
            result.value *= rhs.toDouble();
        }
        else {
            result.num *= rhs.num;
            result.den *= rhs.den;
            result.reduce();
            result.value = result.toDouble();
        }
        return result;
    }

    // Overloaded multiply operator to a double
    Number operator*(double val) {
        Number result = *this;
        if (result.isDouble) {
            result.value *= val;
        }
        else {
            result.num *= val;
            result.reduce();
            result.value = result.toDouble();
        }
        return result;
    }

    // Overloaded divide operator. If both numbers are fractions,
    // multiply this fraction by inverse of the other fraction
    Number operator/(Number const& rhs) {
        Number result = *this;
        if (result.isDouble || rhs.isDouble) {
            result.useDouble();
            result.value /= rhs.toDouble();
        }
        else {
            result = result * Number(rhs.den, rhs.num);
        }
        return result;
    }

    // Overloaded addition operator
    Number operator+(Number const& rhs) {
        Number result = *this;
        if (result.isDouble || rhs.isDouble) {
            result.useDouble();
            result.value += rhs.toDouble();
        }
        else {
            long long lcd = result.den * rhs.den;
            result.num *= rhs.den;
            result.num += rhs.num * result.den;
            result.den = lcd;
            result.reduce();
            result.value = result.toDouble();
        }
        return result;
    }

    // Overloaded subtraction operator
    Number operator-(Number const& rhs) {
        Number result = *this;
        if (result.isDouble || rhs.isDouble) {
            result.useDouble();
            result.value -= rhs.toDouble();
        }
        else {
            long long lcd = result.den * rhs.den;
            result.num *= rhs.den;
            result.num -= rhs.num * result.den;
            result.den = lcd;
            result.reduce();
            result.value = result.toDouble();
        }
        return result;
    }

    // Return double representation of this number
    double toDouble() const {
        if (isDouble) {
            return value;
        }
        double num_d = static_cast<double>(num);
        double den_d = static_cast<double>(den);
        return num_d/den_d;
    }

    // Return the delta (absolute distance) from this number to the other number
    double delta(Number& other) {
        return std::abs(toDouble() - other.toDouble());
    }

    // Return true if this number is equivalent to 0
    bool isZero() {
        if (isDouble) {
            return value == 0.0;
        }
        return num == 0;
    }

    // Return true if this number is equivalent to 1
    bool isOne() {
        if (isDouble) {
            return value == 1.0;
        }
        return num == den;
    }

    // Return true if this number is greater than the given double value
    bool operator>(double val) {
        return toDouble() > val;
    }

    // Overloaded output operator
    friend std::ostream& operator<<(std::ostream& os, const Number& f);
};

// Implement the overloaded output operator with the given number.
// If the number is a fraction, print as a fraction.
// If the number is a decimal, print as a decimal.
std::ostream& operator<<(std::ostream& os, const Number& f)
{
    if (f.isDouble) {
        os << f.value;
    }
    else if (f.num == 0) {
        os << 0;
    }
    else if(f.num == f.den) {
        os << 1;
    }
    else {
        long long gcdVal = gcd(f.num, f.den);
        os << f.num / gcdVal << '/' << f.den / gcdVal;
    }
    return os;
}

// Define Matrix as a vector of vector of Numnbers
typedef std::vector< std::vector<Number> > Matrix;

// Print a matrix. Decimals will be printed with 15 decimal place precision.
void printMatrix(Matrix &matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        std::cout << "  ";
        for (int j = 0; j < matrix[i].size(); ++j) {
            std::cout << std::setprecision(15) << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

// Scale a matrix by the given Number
Matrix scaleMatrix(Matrix &matrix, Number scale) {
    Matrix ret;
    for (int i = 0; i < matrix.size(); ++i) {
        std::vector<Number> temp;
        for (int j = 0; j < matrix[i].size(); ++j) {
            temp.push_back(matrix[i][j] * scale);
        }
        ret.push_back(temp);
    }
    return ret;
}

// Create a new matrix with size [rows, cols] with each element initial value as initial_value.
Matrix createNewMatrix(int rows, int cols, Number initial_value) {
    std::vector<Number> row;
    Matrix mat;
    for (int i = 0; i < cols; ++i) {
        row.push_back(initial_value);
    }
    for (int i = 0; i < rows; ++i) {
        mat.push_back(row);
    }
    return mat;
}

// Add two matricies
Matrix matAdd(Matrix &a, Matrix &b) {
    Matrix sum = createNewMatrix(a.size(), a[0].size(), Number(1,1));
    for(int i = 0; i < a.size(); ++i) {
        for(int j = 0; j < a[0].size(); ++j) { 
            sum[i][j] = a[i][j] + b[i][j];
        }
    }
    return sum;
}

// Multiply two matricies
Matrix matMul(Matrix &a, Matrix &b) {
    const int rows1 = a.size();
    const int cols1 = a[0].size();
    const int cols2 = b[0].size();
    std::vector<Number> m;
    Matrix mult = createNewMatrix(rows1, cols2, Number(0,1));

    // Multiplying matrix a and b and storing in array mult.
    for(int i = 0; i < rows1; ++i) {
        for(int j = 0; j < cols2; ++j) {
            for(int k = 0; k < cols1; ++k)
            {
                mult[i][j] = (a[i][k] * b[k][j]) + mult[i][j];
            }
        }
    }
    return mult;
}

// Return true if the given string is a number
// (either whole number with no decimal places, or decimal number with only 1 '.')
bool is_number(const std::string& s)
{
    bool isValidNumber = (s.find_first_not_of("0123456789.") == std::string::npos);
    isValidNumber &= (s.find_first_of(".") == s.find_last_of("."));
    return isValidNumber;
}

// Convert all values in the given matrix to doubles
void convertMatrixToDouble(Matrix &matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            matrix[i][j].useDouble();
        }
    }
}

// Read matrix from stdin with the following rules/functions:
// 1. Number of rows must match number of columns in the first inputted row.
// 2. All values must be a valid number.
// 3. Ignore all comments (from '#' to the end of the line)
// 4. Ignore eroneous whitespace
// 5. Exit program if matrix has bad formatting
Matrix readMatrix() {
    std::string input;
    Matrix matrix;
    int matrix_size = 0;
    bool initialized = false;
    bool useDouble = false;
    bool notNumeric = false;
    int row_count = 0;

    do {
        bool found_comment = false;
        std::string word;
        int word_count = 0;
        std::vector<Number> temp_row;

        std::getline (std::cin,input);
        if (input == "") {
            break;
        }
        std::istringstream ss(input);
        
        // Read each row column by column
        while (!found_comment && ss >> word) {
            // Remove comment
            int comment_start = word.find('#');
            if (comment_start != std::string::npos) {
                if (comment_start == 0) {
                    found_comment = true;
                }
                word = word.substr(0,comment_start);
            }
            
            // If value contains values
            if (!found_comment) {
                int del = word.find('/');

                // If value is a fraction
                if (del != std::string::npos) {
                    double num, den;
                    std::istringstream num_ss(word.substr(0,del));
                    std::istringstream den_ss(word.substr(del+1));
                    if (!is_number(num_ss.str()) || !is_number(den_ss.str())) {
                        std::cout << "BAD FORMATTING, QUITTING.\n";
                        exit(0);
                    }
                    num_ss >> num;
                    den_ss >> den;
                    temp_row.push_back(Number(num,den));
                }

                // If value is not a fraction
                else {
                    
                    // If not a valid number
                    if (!is_number(word)) {
                        std::cout << "BAD FORMATTING, QUITTING.\n";
                        exit(0);
                    }
                    double val;
                    std::istringstream(word) >> val;

                    // If value is a decimal
                    if (word.find('.') != std::string::npos) {
                        temp_row.push_back(Number(val));
                        useDouble = true;
                    }

                    // If value is a whole number (not a fraction but no decimal places, i.e. '1')
                    // add as fraction (1/1)
                    else {
                        temp_row.push_back(Number(val,1));
                    }
                }
                // If matrix has not been initialized yet, set expected matrix size
                if (!initialized) {
                    ++matrix_size;
                }
                ++word_count;
            }
        }

        // If good row, push to matrix
        if (!temp_row.empty()) {
            matrix.push_back(temp_row);
            initialized = true;
            ++row_count;

            // If matrix is not initialized or row count is not as expected, there is bad formatting, quit.
            if (row_count > matrix_size || word_count != matrix_size) {
                std::cout << "BAD FORMATTING, QUITTING.\n";
                exit(0);
            }
        }
    } while (input != "");

    if (!initialized) {
        std::cout << "BAD FORMATTING, QUITTING.\n";
        exit(0);
    }
    
    for (int c = 0; c < matrix_size; ++c) {
        Number col_total = Number(0,1);
        for (int r = 0; r < matrix_size; ++r) {
            col_total = matrix[r][c] + col_total;
        }
        if (col_total.toDouble() != 1 && col_total.toDouble() != 0) {
            std::cout << "BAD FORMATTING, QUITTING.\n";
            exit(0);
        }
    }

    // If there is a single decimal, convert entire matrix to double
    if (useDouble) {
        convertMatrixToDouble(matrix);
    }

    return matrix;
}

// Remove sink nodes
std::vector<int> reduceMatrix(Matrix &inputMatrix) {
    std::vector<int> nodes_to_remove;
    bool reduced = false;
    for (int col = 0; col < inputMatrix[0].size(); ++col) {
        int zero_count = 0;
        for (int row = 0; row < inputMatrix.size(); ++row) {
            if (inputMatrix[row][col].isZero()) {
                zero_count++;
            }
            else if (inputMatrix[row][col].isOne() && row == col) {
                nodes_to_remove.push_back(col);
                break;
            }
        }
        if (zero_count == inputMatrix[col].size()) {
            nodes_to_remove.push_back(col);
        }
    }

    if (nodes_to_remove.size() > 0) {
        // Remove node columns and rows
        // Only remove first node (node with the lowest index). The rest will be removed recursively.
        int node_to_remove = nodes_to_remove[0];
        for (int row = 0; row < inputMatrix.size(); ++row) {
            inputMatrix[row].erase(inputMatrix[row].begin()+node_to_remove);
        }
        inputMatrix.erase(inputMatrix.begin()+node_to_remove);
        reduced = true;

        nodes_to_remove.clear();
        nodes_to_remove.push_back(node_to_remove);

        // Recursively call reduceMatrix on the current matrix and collect returned removed nodes indices
        if (reduced) {
            std::vector<int> rec_ret = reduceMatrix(inputMatrix);
            for (int i = 0; i < rec_ret.size(); ++i) {
                if (rec_ret[i] >= node_to_remove) {
                    rec_ret[i] += 1;
                }
            }
            nodes_to_remove.insert(nodes_to_remove.end(), rec_ret.begin(), rec_ret.end());
        }

        // Redistribute weights
        for (int col = 0; col < inputMatrix[0].size(); ++col) {
            Number total = Number(0,1);
            for (int row = 0; row < inputMatrix.size(); ++row) {
                total = inputMatrix[row][col] + total;
            }
            for (int row = 0; row < inputMatrix.size(); ++row) {
                inputMatrix[row][col] = inputMatrix[row][col] / total;
            }
        }
    }

    return nodes_to_remove;
}

// Complete one iteration with random walk.
// B is given.
// Keep track of iteration count. If > 5 iterations, convert matrix to double.
bool iterateWithRandomWalk(Matrix &matrix, Matrix &rank, Number B, int iterationCount) {

    // If > 5 iterations, convert matrix to double
    if (iterationCount == 6) {
        convertMatrixToDouble(matrix);
    }

    // Scale matrix by B
    Matrix scaledMatrix = scaleMatrix(matrix, B);

    // Create random walk matrix with initial value (1/rows * (1 - B))
    double rows = matrix.size();
    Matrix randomWalkMatrix = createNewMatrix(rows, rows, Number(1,rows) * (Number(1,1)-B));

    // Add scaled matrix to scaled random walk matrix
    Matrix sum = matAdd(scaledMatrix, randomWalkMatrix);

    // Multiply resulting matrix with rank matrix to create a new rank
    Matrix newrank = matMul(sum, rank);

    // Check if any value has changed more then epsilon. If so, allow another iteration.
    double e = 0.000000000001;
    bool cont = false;
    for (int i = 0; i < newrank.size(); ++i) {
        if (rank[i][0].delta(newrank[i][0]) > e) {
            cont = true;
            break;
        }
    }
    rank = newrank;

    return cont;
}

// Run PageRank algorithm via stdin
int main() {

    // Read matrix from stdin
    Matrix matrix = readMatrix();

    // Reduce matrix
    std::vector<int> removed_nodes = reduceMatrix(matrix);

    // Print input matrix size as well as those nodes that have been removed
    sort(removed_nodes.begin(), removed_nodes.end());
    std::cout << "Input matrix A with dimension " << matrix.size() << "x" << matrix[0].size();
    if (removed_nodes.size() > 0) {
        std::cout << " by removing node ";
        for (int i = 0; i < removed_nodes.size(); ++i) {
            std::cout << removed_nodes[i] + 1;
            if (i < removed_nodes.size()-1) {
                std::cout << ", ";
            }
        }
    }
    
    // Print input matrix
    std::cout << ":\n";
    printMatrix(matrix);

    // Run PageRank iterations. Print rank after every iteration.
    // Stop when delta < epsilon. Stop at 100000 iterations to avoid never-ending program.
    int iterationCount = 0;
    Matrix rank = createNewMatrix(matrix.size(), 1, Number(1, matrix.size()));
    std::cout << "Output matrix C with dimension " << rank.size() << "x" << rank[0].size() << ":\n";
    while (iterateWithRandomWalk(matrix, rank, Number(7,8), ++iterationCount)) {
        std::cout << "--Iteration " << iterationCount << ":" << std::endl;
        printMatrix(rank);
        if (iterationCount > 100000) {
            break;
        }
    }

    // Print final rank
    std::cout << "--Iteration " << iterationCount << ":" << std::endl;
    printMatrix(rank);

    return 0;
}

