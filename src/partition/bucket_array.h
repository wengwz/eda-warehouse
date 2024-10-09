#ifndef BUCKET_ARRAY_H
#define BUCKET_ARRAY_H

#include <unordered_map>
#include <queue>

namespace fm_partitioner {

template<typename GainType, typename IndexType>
class GainBucketArray {
private:
struct bucket_elem
{
    GainType gain;
    IndexType id;
    bucket_elem* pre_elem = nullptr;
    bucket_elem* nxt_elem = nullptr;
};

public:
    GainBucketArray(IndexType max_elem_count) : max_elem_count(max_elem_count), elem_count(0) {
        id2bucket_elem.resize(max_elem_count, nullptr);
    };

    GainBucketArray() = delete;
    
    ~GainBucketArray() {
        for (IndexType i = 0; i < id2bucket_elem.size(); i++) {
            if (id2bucket_elem[i] != nullptr) {
                delete id2bucket_elem[i];
            }
        }
        bucket_array.clear();
    }


    void insert(IndexType id, GainType gain) {
        assert(id < max_elem_count);

        if (id2bucket_elem[id] != nullptr) {
            throw std::runtime_error("Element already exists in bucket array");
        }

        bucket_elem* new_elem = new bucket_elem();
        new_elem->gain = gain;
        new_elem->id = id;

        //
        id2bucket_elem[id] = new_elem;
        //
        if (bucket_array.count(gain) == 0) {
            bucket_array[gain] = new_elem;
            gain_queue.push(gain);
        } else {
            bucket_elem* head_elem = bucket_array[gain];
            if (head_elem != nullptr) {
                head_elem->pre_elem = new_elem;
            }
            new_elem->nxt_elem = head_elem;
            bucket_array[gain] = new_elem;
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
            bucket_array[elem->gain] = elem->nxt_elem;
        } else {
            elem->pre_elem->nxt_elem = elem->nxt_elem;
        }

        if (elem->nxt_elem != nullptr) {
            elem->nxt_elem->pre_elem = elem->pre_elem;
        }

        id2bucket_elem[id] = nullptr;
        delete elem;
        elem_count--;
        return true;
    }

    bool update(IndexType id, GainType new_gain) {
        assert(id < max_elem_count);

        if (id2bucket_elem[id] == nullptr) {
            return false; // return false if element does not exist
        }

        remove(id);
        insert(id, new_gain);

        return true;
    }

    IndexType pop() {
        IndexType id = get_max_gain_id();
        remove(id);
        return id;
    }

// getters
    bool has_elem(IndexType id) const {
        assert(id < max_elem_count);
        return id2bucket_elem[id] != nullptr;
    }
    
    bool empty() const {
        return elem_count == 0;
    }

    IndexType get_max_gain_id() {
        if (elem_count == 0) {
            return -1;
        }

        GainType max_gain = gain_queue.top();

        while (bucket_array.at(max_gain) == nullptr) {
            gain_queue.pop();
            bucket_array.erase(max_gain);

            max_gain = gain_queue.top();
        }

        return bucket_array[max_gain]->id;
    }

    IndexType get_size() const {
        return elem_count;
    }

    GainType get_gain_of(IndexType id) const {
        assert(id < max_elem_count);

        if (id2bucket_elem[id] == nullptr) {
            throw std::runtime_error("Element does not exist in bucket array");
        }
        return id2bucket_elem[id]->gain;
    }

private:
    IndexType max_elem_count;
    IndexType elem_count;
    std::unordered_map<GainType, bucket_elem*> bucket_array;
    std::vector<bucket_elem*> id2bucket_elem;
    std::priority_queue<GainType> gain_queue;

};
} // namespace fm_partitioner
#endif
