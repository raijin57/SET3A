//
// Created by arsen on 15.11.2025.
//
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <cmath>

class ArrayGenerator {
public:
    enum class Type { RANDOM, REVERSE_SORTED, ALMOST_SORTED };
    explicit ArrayGenerator(size_t maxSize = 100000, int maxValue = 6000, double swaps = 0.001, unsigned int seed = 42) :
        maxSize_(maxSize), maxValue_(maxValue), swaps_(swaps) {
        rng_.seed(seed);
        GenerateRandomBase();
        GenerateSortedBase();
        GenerateReverseSortedBase();
        GenerateAlmostSortedBase();
    }
    [[nodiscard]] std::vector<int> GetSubarray(Type type, size_t size) const {
        switch (type) {
            case Type::RANDOM:
                return std::vector(randomBase_.begin(), randomBase_.begin() + size);
            case Type::REVERSE_SORTED:
                return std::vector(reverseSortedBase_.begin(), reverseSortedBase_.begin() + size);
            case Type::ALMOST_SORTED:
                return std::vector(almostSortedBase_.begin(), almostSortedBase_.begin() + size);
            default:
                return {};
        }
    }

    [[nodiscard]] size_t max_size() const { return maxSize_; }
    [[nodiscard]] int max_value() const { return maxValue_; }

private:
    size_t maxSize_;
    int maxValue_;
    double swaps_;
    std::mt19937 rng_;
    std::vector<int> randomBase_;
    std::vector<int> sortedBase_;
    std::vector<int> reverseSortedBase_;
    std::vector<int> almostSortedBase_;
    void GenerateRandomBase() {
        randomBase_.resize(maxSize_);
        std::uniform_int_distribution dist(0, maxValue_);
        for (size_t i = 0; i < maxSize_; ++i)
            randomBase_[i] = dist(rng_);
    }
    void GenerateSortedBase() {
        sortedBase_.resize(maxSize_);
        std::uniform_int_distribution dist(0, maxValue_);
        for (size_t i = 0; i < maxSize_; ++i)
            sortedBase_[i] = dist(rng_);
        std::sort(sortedBase_.begin(), sortedBase_.end());
    }
    void GenerateReverseSortedBase() {
        reverseSortedBase_ = sortedBase_;
        std::reverse(reverseSortedBase_.begin(), reverseSortedBase_.end());
    }
    void GenerateAlmostSortedBase() {
        almostSortedBase_ = sortedBase_;
        size_t num_pairs = std::max<size_t>(1, static_cast<int>(maxSize_ * swaps_));
        std::uniform_int_distribution<size_t> pos_dist(0, maxSize_ - 1);
        for (size_t k = 0; k < num_pairs; ++k) {
            size_t i = pos_dist(rng_);
            size_t j = pos_dist(rng_);
            if (i != j) {
                std::swap(almostSortedBase_[i], almostSortedBase_[j]);
            }
        }
    }
};

int64_t MedianVec(std::vector<int64_t> v) {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n == 0) {
        return 0;
    }
    if (n % 2 == 1) {
        return v[n / 2];
    }
    return (v[n / 2 - 1] + v[n / 2]) / 2;
}
double MeanVec(const std::vector<int64_t> &v) {
    if (v.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (auto x: v) {
        sum += static_cast<double>(x);
    }
    return sum / static_cast<double>(v.size());
}
double StddevVec(const std::vector<int64_t> &v, double mean) {
    if (v.size() < 2) {
        return 0.0;
    }
    double s = 0.0;
    for (auto x: v) {
        double d = static_cast<double>(x) - mean;
        s += d * d;
    }
    return std::sqrt(s / static_cast<double>(v.size() - 1));
}

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
void QuickSort(std::vector<int>& A, int left, int right) {
    if (left < right) {
        int pivot = Partition(A, left, right);
        QuickSort(A, left, pivot - 1);
        QuickSort(A, pivot + 1, right);
    }
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

class SortTester {
public:
    SortTester(ArrayGenerator &arrayGen, size_t startN, size_t maxN, size_t stepN, int repeats) :
        arrayGen_(arrayGen), startN_(startN), maxN_(maxN), stepN_(stepN), repeats_(repeats) {}
    void RunQuick(const std::string &outFilename) const {
        std::ofstream ofs(outFilename, std::ios::out | std::ios::trunc);
        ofs << "type,n,min_us,max_us,mean_us,median_us,stddev_us\n";
        ofs.flush();
        std::vector types = {ArrayGenerator::Type::RANDOM, ArrayGenerator::Type::REVERSE_SORTED,
                             ArrayGenerator::Type::ALMOST_SORTED};
        size_t jobCnt = 0;
        std::cout << "Quick Sort started\n";
        for (auto curType: types) {
            for (size_t n = startN_; n <= maxN_; n += stepN_) {
                ++jobCnt;
                std::vector<int64_t> timesUs;
                timesUs.reserve(repeats_);
                for (int rep = 0; rep < repeats_; ++rep) {
                    std::vector<int> base = arrayGen_.GetSubarray(curType, n);
                    if (rep == 0) {
                        auto warm = base;
                        QuickSort(warm, 0, warm.size() - 1);
                    }
                    auto test = base;
                    auto t1 = std::chrono::high_resolution_clock::now();
                    QuickSort(test, 0, test.size() - 1);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto durUs = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
                    timesUs.push_back(durUs);
                }
                int64_t minUs = *std::min_element(timesUs.begin(), timesUs.end());
                int64_t maxUs = *std::max_element(timesUs.begin(), timesUs.end());
                double meanUs = MeanVec(timesUs);
                int64_t medianUs = MedianVec(timesUs);
                double stddevUs = StddevVec(timesUs, meanUs);
                ofs << static_cast<int>(curType) << ',' << n << ',' << minUs << ',' << maxUs << ',' << std::fixed
                    << std::setprecision(3) << meanUs << ',' << medianUs << ',' << std::fixed << std::setprecision(3)
                    << stddevUs << '\n';
                if (jobCnt % 10 == 0) {
                    ofs.flush();
                }
            }
        }
        std::cout << "Quick Sort finished -> " << outFilename << '\n';
        ofs.close();
    }
    void RunIntro(const std::string &outFilename) const {
        std::ofstream ofs(outFilename, std::ios::out | std::ios::trunc);
        ofs << "type,n,min_us,max_us,mean_us,median_us,stddev_us\n";
        ofs.flush();
        std::vector types = {ArrayGenerator::Type::RANDOM, ArrayGenerator::Type::REVERSE_SORTED,
                             ArrayGenerator::Type::ALMOST_SORTED};
        size_t jobCnt = 0;
        std::cout << "Intro Sort started\n";
        for (auto curType: types) {
            for (size_t n = startN_; n <= maxN_; n += stepN_) {
                ++jobCnt;
                std::vector<int64_t> timesUs;
                timesUs.reserve(repeats_);
                for (int rep = 0; rep < repeats_; ++rep) {
                    std::vector<int> base = arrayGen_.GetSubarray(curType, n);
                    if (rep == 0) {
                        auto warm = base;
                        IntroSort(warm);
                    }
                    auto test = base;
                    auto t1 = std::chrono::high_resolution_clock::now();
                    IntroSort(test);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto durUs = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
                    timesUs.push_back(static_cast<int64_t>(durUs));
                }
                int64_t minUs = *std::min_element(timesUs.begin(), timesUs.end());
                int64_t maxUs = *std::max_element(timesUs.begin(), timesUs.end());
                double meanUs = MeanVec(timesUs);
                int64_t medianUs = MedianVec(timesUs);
                double stddevUs = StddevVec(timesUs, meanUs);
                ofs << static_cast<int>(curType) << ',' << n << ',' << minUs << ',' << maxUs << ','
                    << std::fixed << std::setprecision(3) << meanUs << ',' << medianUs << ',' << std::fixed
                    << std::setprecision(3) << stddevUs << '\n';
                if (jobCnt % 10 == 0) {
                    ofs.flush();
                }
            }
        }
        std::cout << "Intro ended" << '\n';
        ofs.close();
    }

private:
    ArrayGenerator &arrayGen_;
    size_t startN_;
    size_t maxN_;
    size_t stepN_;
    int repeats_;
};

int main() {
    constexpr size_t maxN = 30000;
    constexpr size_t startN = 500;
    constexpr size_t stepN = 100;
    constexpr int runs = 9;
    constexpr int maxValue = 6000;
    constexpr double swapsFraction = 0.001;
    constexpr unsigned int seed = 42;
    ArrayGenerator arrayGen(maxN, maxValue, swapsFraction, seed);
    SortTester tester(arrayGen, startN, maxN, stepN, runs);
    const std::string quickOutput = "/home/arsen/CLionProjects/SET3/A3/quick_standard_results.csv";
    tester.RunQuick(quickOutput);
    const std::string introOut = "/home/arsen/CLionProjects/SET3/A3/intro_results.csv";
    tester.RunIntro(introOut);
    std::cout << "Все эксперименты завершены.\n";
    return 0;
}