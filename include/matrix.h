#pragma once

#include "array.h"
#include "vec2.h"

struct Matrix {
    int rows;
    int cols;
    TArray<double> data;

    Matrix(int r, int c, double value = 0.0) : rows(r), cols(c) {
        
        data.resize(rows*cols);
        
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                
                data[i*cols + j] = value;
            }
        }
        
    }

    ~Matrix() {
        data.clear();
    }

    void set_value(int r, int c, double value) {
        if (r >= 0 && r < rows && c >= 0 && c < cols) {
            data[r* cols + c] = value;
        }
    }

    double& operator()(int r, int c) {

        assert(r >= 0 && r < rows && c >= 0 && c < cols); 
        return data[r*cols + c];
        
    }

    void reset(int N = 0, double value = 0.0) {
        if (N == 0)
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    data[i*cols + j] = value;
                }
            }
        else
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < N; j++) {
                    data[i*cols + j] = value;
                }
            }
    }

    double get_value(int r, int c) const {
        if (r >= 0 && r < rows && c >= 0 && c < cols) {
            return data[r*cols + c];
        }
        return 0.0;
    }



    double operator()(int r, int c) const {
        if (r >= 0 && r < rows && c >= 0 && c < cols) {
            return data[r*cols + c];
        }
    }

    void mvp(const Vec2d &v, Vec2d &result) const {
        result = {0.0, 0.0};
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += data[i*cols + j] * v[j];
            }
        }
        return;
    }

    void mvp(const double *v, double *result) const {
        for (int i = 0; i < rows; i++) {
            result[i] = 0.0;
            for (int j = 0; j < cols; j++) {
                result[i] += data[i*cols + j] * v[j];
            }
        }
        return;
    }
};