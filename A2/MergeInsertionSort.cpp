//
// Created by arsen on 15.11.2025.
//
#include <vector>
#include <iostream>
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

int main() {
    int n;
    std::cin >> n;
    std::vector<int> a(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> a[i];
    }
    MergeSortHybrid(a, 15);
    for (int i = 0; i < a.size(); ++i) {
        std::cout << a[i] << " ";
    }
}