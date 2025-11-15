//
// Created by arsen on 15.11.2025.
//
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

class ArrayGenerator {
public:
    enum class Type { RANDOM, REVERSE_SORTED, ALMOST_SORTED };
    ArrayGenerator(size_t maxSize = 100000, int maxValue = 6000, double swaps = 0.001, unsigned int seed = 42) :
        maxSize_(maxSize), maxValue_(maxValue), swaps_(swaps) {
        rng_.seed(seed);
        GenerateRandomBase();
        GenerateSortedBase();
        GenerateReverseSortedBase();
        GenerateAlmostSortedBase();
    }
    std::vector<int> GetSubarray(Type type, size_t size) const {
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
    size_t max_size() const { return maxSize_; }
    int max_value() const { return maxValue_; }

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

void MergeRange(std::vector<int> &a, int left, int mid, int right, std::vector<int> &tmpBuffer) {
    int i = left;
    int j = mid + 1;
    int k = left;
    while (i <= mid && j <= right) {
        if (a[i] <= a[j]) {
            tmpBuffer[k++] = a[i++];
        } else {
            tmpBuffer[k++] = a[j++];
        }
    }
    while (i <= mid) {
        tmpBuffer[k++] = a[i++];
    }
    while (j <= right) {
        tmpBuffer[k++] = a[j++];
    }
    for (int idx = left; idx <= right; ++idx) {
        a[idx] = tmpBuffer[idx];
    }
}
void MergeSortRec(std::vector<int> &a, int left, int right, std::vector<int> &tmpBuffer) {
    if (left >= right) {
        return;
    }
    int mid = left + (right - left) / 2;
    MergeSortRec(a, left, mid, tmpBuffer);
    MergeSortRec(a, mid + 1, right, tmpBuffer);
    MergeRange(a, left, mid, right, tmpBuffer);
}
void MergeSort(std::vector<int> &a) {
    if (a.empty()) {
        return;
    }
    std::vector<int> tmpBuffer(a.size());
    MergeSortRec(a, 0, static_cast<int>(a.size()) - 1, tmpBuffer);
}

void InsertionSortRangeHybrid(std::vector<int> &a, int left, int right) {
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
void MergeSortHybridRec(std::vector<int> &a, int left, int right, std::vector<int> &tmpBuffer, int threshold) {
    int len = right - left + 1;
    if (len <= threshold) {
        InsertionSortRangeHybrid(a, left, right);
        return;
    }
    int mid = left + (right - left) / 2;
    MergeSortHybridRec(a, left, mid, tmpBuffer, threshold);
    MergeSortHybridRec(a, mid + 1, right, tmpBuffer, threshold);
    MergeRange(a, left, mid, right, tmpBuffer);
}
void MergeSortHybrid(std::vector<int> &a, int threshold) {
    if (a.empty()) {
        return;
    }
    std::vector<int> tmpBuffer(a.size());
    MergeSortHybridRec(a, 0, static_cast<int>(a.size()) - 1, tmpBuffer, threshold);
}

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

class SortTester {
public:
    SortTester(ArrayGenerator &arrayGen, size_t startN, size_t maxN, size_t stepN, int repeats) :
        arrayGen_(arrayGen), startN_(startN), maxN_(maxN), stepN_(stepN), repeats_(repeats) {}
    void RunMerge(const std::string &outFilename) {
        std::ofstream ofs(outFilename, std::ios::out | std::ios::trunc);
        ofs << "type,n,min_us,max_us,mean_us,median_us,stddev_us\n";
        ofs.flush();
        std::vector types = {ArrayGenerator::Type::RANDOM, ArrayGenerator::Type::REVERSE_SORTED,
                             ArrayGenerator::Type::ALMOST_SORTED};
        size_t jobCnt = 0;
        std::cout << "Merge Sort started -> " << outFilename << '\n';
        for (auto curType: types) {
            for (size_t n = startN_; n <= maxN_; n += stepN_) {
                ++jobCnt;
                std::vector<int64_t> timesUs;
                timesUs.reserve(repeats_);
                for (int rep = 0; rep < repeats_; ++rep) {
                    std::vector<int> base = arrayGen_.GetSubarray(curType, n);
                    if (rep == 0) {
                        auto warm = base;
                        MergeSort(warm);
                    }
                    auto test = base;
                    auto t1 = std::chrono::high_resolution_clock::now();
                    MergeSort(test);
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
        std::cout << "Merge Sort finished -> " << outFilename << '\n';
        ofs.close();
    }
    void RunHybrid(const std::string &outFilename, int threshold) {
        std::ofstream ofs(outFilename, std::ios::out | std::ios::trunc);
        ofs << "type,n,min_us,max_us,mean_us,median_us,stddev_us\n";
        ofs.flush();
        std::vector types = {ArrayGenerator::Type::RANDOM, ArrayGenerator::Type::REVERSE_SORTED,
                             ArrayGenerator::Type::ALMOST_SORTED};
        size_t jobCnt = 0;
        std::cout << "Merge+Insertion started (threshold = " << threshold << ") -> " << outFilename << '\n';
        for (auto curType: types) {
            for (size_t n = startN_; n <= maxN_; n += stepN_) {
                ++jobCnt;
                std::vector<int64_t> timesUs;
                timesUs.reserve(repeats_);
                for (int rep = 0; rep < repeats_; ++rep) {
                    std::vector<int> base = arrayGen_.GetSubarray(curType, n);
                    if (rep == 0) {
                        auto warm = base;
                        MergeSortHybrid(warm, threshold);
                    }
                    auto test = base;
                    auto t1 = std::chrono::high_resolution_clock::now();
                    MergeSortHybrid(test, threshold);
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
        std::cout << "Merge+Insertion ended" << '\n';
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
    try {
        constexpr size_t maxN = 30000;
        constexpr size_t startN = 500;
        constexpr size_t stepN = 100;
        constexpr int runs = 9;
        constexpr int maxValue = 6000;
        constexpr double swapsFraction = 0.001;
        constexpr unsigned int seed = 42;
        ArrayGenerator arrayGen(maxN, maxValue, swapsFraction, seed);
        SortTester tester(arrayGen, startN, maxN, stepN, runs);
        const std::string standardOut = "/home/arsen/CLionProjects/SET3/A2/merge_standard_results.csv";
        tester.RunMerge(standardOut);
        const std::string hybridOut = "/home/arsen/CLionProjects/SET3/A2/merge_insertion_results.csv";
        constexpr int threshold = 15;
        tester.RunHybrid(hybridOut, threshold);
        std::cout << "Все эксперименты завершены.\n";
        return 0;
    } catch (...) {
        return 101;
    }
}
