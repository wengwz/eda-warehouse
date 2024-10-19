#ifndef BUCKET_ARRAY_H
#define BUCKET_ARRAY_H

#include <unordered_map>
#include <queue>

namespace partition {

template<typename ValueType, typename IndexType, typename Compare = std::less<ValueType>>
class PriorityBucketArray {
private:
struct bucket_elem
{
    ValueType value;
    IndexType id;
    bucket_elem* pre_elem = nullptr;
    bucket_elem* nxt_elem = nullptr;
};

public:
    PriorityBucketArray(IndexType max_elem_count) : max_elem_count(max_elem_count), elem_count(0) {
        id2bucket_elem.resize(max_elem_count, nullptr);
    };

    
    ~PriorityBucketArray() {
        for (IndexType id = 0; id < max_elem_count; id++) {
            if (id2bucket_elem[id] != nullptr) {
                delete id2bucket_elem[id];
            }
        }
        bucket_array.clear();
    }


    void insert(IndexType id, ValueType value) {
        assert(id >= 0 && id < max_elem_count);

        if (id2bucket_elem[id] != nullptr) {
            throw std::runtime_error("Element already exists in bucket array");
        }

        bucket_elem* new_elem = new bucket_elem();
        new_elem->value = value;
        new_elem->id = id;

        //
        id2bucket_elem[id] = new_elem;
        //
        if (bucket_array.count(value) == 0) {
            bucket_array[value] = new_elem;
            value_queue.push(value);
        } else {
            bucket_elem* head_elem = bucket_array[value];
            if (head_elem != nullptr) {
                head_elem->pre_elem = new_elem;
            }
            new_elem->nxt_elem = head_elem;
            bucket_array[value] = new_elem;
        }
        //
        elem_count++;
    }

    bool remove(IndexType id) {
        assert(id < max_elem_count);

        if (id2bucket_elem[id] == nullptr) {
            return false; // return false if element does not exist
        }

        bucket_elem* elem = id2bucket_elem[id];

        if (elem->pre_elem == nullptr) {
            bucket_array[elem->value] = elem->nxt_elem;
        } else {
            elem->pre_elem->nxt_elem = elem->nxt_elem;
        }

        if (elem->nxt_elem != nullptr) {
            elem->nxt_elem->pre_elem = elem->pre_elem;
        }

        delete id2bucket_elem[id];
        id2bucket_elem[id] = nullptr;
        elem_count--;
        return true;
    }

    bool update(IndexType id, ValueType new_val) {
        assert(id >= 0 && id < max_elem_count);

        if (id2bucket_elem[id] == nullptr) {
            return false; // return false if element does not exist
        }

        if (id2bucket_elem[id]->value == new_val) {
            return true; // return true if value does not change
        }
        
        remove(id);
        insert(id, new_val);

        return true;
    }

    IndexType pop() {
        IndexType id = get_top_id();
        remove(id);
        return id;
    }

    // void resize(IndexType max_elem_count) {
    //     this->max_elem_count = max_elem_count;
    //     id2bucket_elem.resize(max_elem_count, nullptr);
    // }

// getters
    bool has_elem(IndexType id) const {
        assert(id >=0 && id < max_elem_count);
        return id2bucket_elem[id] != nullptr;
    }
    
    bool empty() const {
        return elem_count == 0;
    }

    IndexType get_top_id() {
        if (elem_count == 0) {
            return -1;
        }

        ValueType max_val = value_queue.top();

        while (bucket_array.at(max_val) == nullptr) {
            value_queue.pop();
            bucket_array.erase(max_val);

            max_val = value_queue.top();
        }

        return bucket_array[max_val]->id;
    }

    IndexType get_size() const {
        return elem_count;
    }

    ValueType get_value_of(IndexType id) const {
        assert(id >= 0 && id < max_elem_count);

        if (id2bucket_elem[id] == nullptr) {
            throw std::runtime_error("Element does not exist in bucket array");
        }
        return id2bucket_elem[id]->value;
    }

private:
    IndexType max_elem_count;
    IndexType elem_count;
    std::unordered_map<ValueType, bucket_elem*> bucket_array;
    std::vector<bucket_elem*> id2bucket_elem;

    std::priority_queue<ValueType, std::vector<ValueType>, Compare> value_queue;

};

template<typename ValueType, typename IndexType>
using MaxPriorityBucketArray = PriorityBucketArray<ValueType, IndexType, std::less<ValueType>>;

template<typename ValueType, typename IndexType>
using MinPriorityBucketArray = PriorityBucketArray<ValueType, IndexType, std::greater<ValueType>>;

} // namespace partitioner

#endif
