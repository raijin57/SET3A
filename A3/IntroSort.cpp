//
// Created by arsen on 15.11.2025.
//
#include <iostream>
#include <random>
#include <vector>

void InsertionSort(std::vector<int> &a, int left, int right) {
    for (int i = left + 1; i <= right; ++i) {
        int key = a[i];
        int j = i - 1;
        while (j >= left && a[j] > key) {
            a[j + 1] = a[j];
            --j;
        }
        a[j + 1] = key;
    }
}
void SinkDown(std::vector<int>& A, size_t i, size_t n) {
    while (2 * i + 1 < n) {
        size_t j = 2 * i + 1;
        if (j + 1 < n && A[j] < A[j + 1]) {
            ++j;
        }
        if (A[i] >= A[j]) {
            break;
        }
        std::swap(A[i], A[j]);
        i = j;
    }
}
void Heapify(std::vector<int>& A, size_t i) {
    SinkDown(A, i, A.size());
}
void BuildMaxHeap(std::vector<int>& A) {
    if (A.empty()) {
        return;
    }
    for (int i = static_cast<int>(A.size() / 2) - 1; i >= 0; --i) {
        Heapify(A, i);
    }
}
void HeapSort(std::vector<int>& A, int left, int right) {
    if (right <= left) {
        return;
    }
    std::vector tmp(A.begin() + left, A.begin() + right + 1);
    BuildMaxHeap(tmp);
    size_t n = tmp.size() - 1;
    while (n > 0) {
        std::swap(tmp[0], tmp[n]);
        SinkDown(tmp, 0, n);
        --n;
    }
    std::copy(tmp.begin(), tmp.end(), A.begin() + left);
}

int Partition(std::vector<int>& A, int left, int right) {
    thread_local std::mt19937 rng(42);
    std::uniform_int_distribution dist(left, right);
    const int random = dist(rng);
    std::swap(A[random], A[right]);
    const int pivot = A[right];
    int i = left - 1;
    for (int j = left; j <= right - 1; ++j) {
        if (A[j] <= pivot) {
            ++i;
            std::swap(A[i], A[j]);
        }
    }
    std::swap(A[i + 1], A[right]);
    return i + 1;
}

void IntroSort(std::vector<int>& A, int left = 0, int right = -1, int depth = -1) {
    if (right == -1) {
        right = A.size() - 1;
    }
    if (left >= right) {
        return;
    }
    int size = right - left + 1;
    if (depth == -1) {
        depth = 2 * std::floor(std::log2(size));
    }
    if (size <= 16) {
        InsertionSort(A, left, right);
        return;
    }
    if (depth == 0) {
        HeapSort(A, left, right);
        return;
    }
    int pivot = Partition(A, left, right);
    IntroSort(A, left, pivot - 1, depth - 1);
    IntroSort(A, pivot + 1, right, depth - 1);
}

int main() {
    int n;
    std::cin >> n;
    std::vector<int> A(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> A[i];
    }
    IntroSort(A);
    for (int x : A) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;
    return 0;
}