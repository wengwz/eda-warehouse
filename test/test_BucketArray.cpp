#include <iostream>
#include <vector>
#include <cassert>

#include "bucket_array.h"
using namespace partition;
int main() {

    // test case
    int max_id = 32;
    std::vector<int> gains = {-1, 2, 8, 9, 20, 22, 99, 100, -8, 9, 11, 2, 8};
    std::vector<int> sorted_gains = {100, 99, 22, 20, 11, 9, 9, 8, 8, 2, 2, -1, -8};

    PriorityBucketArray<int, int> bucket_array(max_id);

    // test insert, pop and getters
    std::cout << "Test insert, pop, and getters" << std::endl;

    std::cout << "Test insert and get_size functions" << std::endl;
    for (int i = 0; i < gains.size(); i++) {
        bucket_array.insert(i, gains[i]);
    }
    assert(bucket_array.get_size() == gains.size());

    std::cout << "Test get_gain_of function" << std::endl;
    for (int i = 0; i < gains.size(); i++) {
        int gain = bucket_array.get_value_of(i);
        assert(gain == gains[i]);
    }

    std::cout << "Test get_max_gain_id and pop function" << std::endl;
    for (int i = 0; i < sorted_gains.size(); i++) {
        int id = bucket_array.get_top_id();
        int gain = bucket_array.get_value_of(id);
        assert(gain == sorted_gains[i]);
        bucket_array.pop();
    }

    std::cout << "Test empty function" << std::endl;
    assert(bucket_array.empty());

    // test update
    std::vector<int> new_gains = {-1, 2, 8, 9, 20, 19, 99, 8, -8, 9, 11, 2, 8};
    std::vector<int> new_sorted_gains = {99, 20, 19, 11, 9, 9, 8, 8, 8, 2, 2, -1, -8};
    std::cout << "Test update and remove function of bucket array" << std::endl;
    for (int i = 0; i < gains.size(); i++) {
        bucket_array.insert(i, gains[i]);
    }

    std::cout << "Test update function" << std::endl;
    for (int i = 0; i < new_gains.size(); i++) {
        bucket_array.update(i, new_gains[i]);
    }
    
    for (int i = 0; i < new_sorted_gains.size(); i++) {
        int id = bucket_array.get_top_id();
        int gain = bucket_array.get_value_of(id);
        assert(gain == new_sorted_gains[i]);
        bucket_array.remove(id);
    }

    bucket_array.insert(1, 5);
    std::cout << "All tests passed" << std::endl;

    return 0;
}