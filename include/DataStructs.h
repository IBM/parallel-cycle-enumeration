#ifndef _DATA_STRUCTS_H_
#define _DATA_STRUCTS_H_

#include <list>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <tbb/combinable.h>
#include <tbb/spin_mutex.h>

#include "Macros.h"

using namespace std;
using namespace tbb;

/// Concurrent counter for collecting statistics
template<typename T = unsigned long>
class TConcurrentCounter {
public:
    TConcurrentCounter() {}

    void increment(T inc = 1) {
        bool exists = true;
        T& local_cnt = pt_counter.local(exists);
        if(exists) local_cnt += inc; else local_cnt = inc;
    }

    void decrement() {
        bool exists = true;
        T& local_cnt = pt_counter.local(exists);
        if(exists) local_cnt--; else local_cnt = 0;
    }

    T getResult() {
        if(!combined) {
            combined = true;
            result = 0;
            pt_counter.combine_each([&](T num) {
                result += num;
            });
        }
        return result;
    }

private:
    tbb::combinable<T> pt_counter;

    bool combined = false;
    T result = 0;
};

typedef TConcurrentCounter<unsigned long> ConcurrentCounter;

/// Wrapper for dynamic array
template<typename T>
class VectorPath {
public:
    VectorPath(int s) : vect(new T [s]) {}
    ~VectorPath() {
        if(vect) delete [] vect;
    }

    void push_back(T x) {
        ++it;
        vect[it] = x;
    }

    T back() {
        if(it == -1) return T();
        return vect[it];
    }

    void pop_back() {
        if(it >= 0) it--;
    }

    int size() { return it+1; }
private:
    T* vect = NULL;
    int it = -1;
};

/// Wrapper for std::unordered_set
class HashSet {
public:
    HashSet() : elems() {};
    HashSet(int s) : elems() {};
    HashSet(const HashSet& hs) : elems(hs.elems) {};

    void insert(int el) { elems.insert(el); }

    void remove(int el) {
        auto it = elems.find(el);
        if(it != elems.end()) elems.erase(it);
    }

    bool exists(int el) {
        if(elems.find(el) != elems.end()) return true;
        return false;
    }
    void include(const HashSet& other) {
        for(auto el : other.elems)
            insert(el);
    }

    int size() { return elems.size(); }
    void clear() { elems.clear(); }

    template <typename TF>
    void for_each(TF&& f){ for(auto el : elems) f(el); }
    std::unordered_set<int>::iterator begin() { return elems.begin(); }
    std::unordered_set<int>::iterator end() { return elems.end(); }
    std::unordered_set<int>::iterator erase(std::unordered_set<int>::iterator it) { return elems.erase(it); }
private:
    friend class HashSetStack;
    unordered_set<int> elems;
};

class HashSetStack {
private:
    typedef tbb::spin_mutex HashSetMutexType;
    HashSetMutexType HashSetMutex;

public:
    HashSetStack(bool conc = false) : curLevel(0), elems(), concurrent(conc) {};
    HashSetStack(int s, bool conc = false) : curLevel(0), elems(), concurrent(conc){};
    HashSetStack(const HashSetStack& hs) :
        curLevel(hs.curLevel), elems(hs.elems), concurrent(hs.concurrent)  {};

    HashSetStack* clone() {
        if(concurrent) HashSetMutex.lock();
        HashSetStack* ret =  new HashSetStack(*this);
        if(concurrent) HashSetMutex.unlock();
        return ret;
    }

    HashSetStack* clone(int lvl) {
        HashSetStack* ret = new HashSetStack();
        ret->curLevel = lvl;
        ret->concurrent = concurrent;

        if(concurrent) HashSetMutex.lock();
        for(auto it = elems.begin(); it != elems.end(); ++it) {
            if(it->second <= lvl) ret->elems.insert({it->first, it->second});
        }
        if(concurrent) HashSetMutex.unlock();

        return ret;
    }

    void reserve(int s) {}

    void incrementLevel() {
        if(concurrent) HashSetMutex.lock();
        curLevel++;
        if(concurrent) HashSetMutex.unlock();
    }

    void decrementLevel() {
        if(concurrent) HashSetMutex.lock();
        curLevel--;
        for(auto it = elems.begin(); it != elems.end(); ) {
            if(it->second > curLevel) it = elems.erase(it);
            else it++;
        }
        if(concurrent) HashSetMutex.unlock();
    }

    void setLevel(int lvl) {
        if(concurrent) HashSetMutex.lock();
        if(lvl < curLevel) {
            curLevel = lvl;
            for(auto it = elems.begin(); it != elems.end(); ) {
                if(it->second > curLevel) it = elems.erase(it);
                else it++;
            }
        }
        if(concurrent) HashSetMutex.unlock();
    }

    void insert(int el) {
        if(concurrent) HashSetMutex.lock();
        if(elems.find(el) == elems.end()) elems[el] = curLevel;
        if(concurrent) HashSetMutex.unlock();
    }

    void remove(int el) {
        auto it = elems.find(el);
        if(it != elems.end()) elems.erase(it);
    }

    bool exists(int el) {
        if(elems.find(el) != elems.end()) return true;
        return false;
    }

    void include(const HashSet& other) {
        for(auto el : other.elems) {
            insert(el);
        }
    }

    template <typename TF>
    void for_each(TF&& f){
        for(auto el : elems) {
            f(el);
        }
    }

    void exclude(const HashSet& other) {
        for(auto el : other.elems) {
            remove(el);
        }
    }

    int size() { return elems.size(); }

    void clear() {
        elems.clear();
    }

    void copy(const HashSet& other) {
        elems.clear();

        for(auto el : other.elems) {
            insert(el);
        }
    }

private:
    int curLevel = 0;
    unordered_map<int, int> elems;
    bool concurrent;

};

/// Wrapper for std::unordered_map
class HashMap {
public:
    HashMap() : elems() {};
    HashMap(int s) : elems() {};
    HashMap(const HashMap& hs) : elems(hs.elems) {};

    void insert(int el, int num) { if(num != 0) elems[el] = num; else elems.erase(el); }
    void erase(int el) { elems.erase(el); }

    bool exists(int el) {
        auto it = elems.find(el);
        if(it != elems.end()) return !!it->second;
        return false;
    }

    bool exists(int el, int ts) {
        auto it = elems.find(el);
        if(it == elems.end()) return false;
        if(it->second == -1 || ts >= it->second) return true;
        return false;
    }

    int at(int el) {
        auto it = elems.find(el);
        if(it == elems.end()) return 0;
        return it->second;
    }

    template <typename TF>
    void for_each(TF&& f){ for(auto el : elems) f(el.first); }
    int size() { return elems.size(); }
private:
    friend class HashMapStack;
    unordered_map<int, int> elems;
};

/// Hash map with a stack
class HashMapStack {
private:
    struct StackElem { int cltime; int level; StackElem(int ct = 0, int l = 0) : cltime(ct), level(l) {}; };
    typedef tbb::spin_mutex HashMapMutexType;
    HashMapMutexType HashMapMutex;

public:
    HashMapStack(bool conc = false) : curLevel(0), elems(), concurrent(conc) {};
    HashMapStack(const HashMapStack& hs) : curLevel(hs.curLevel), elems(hs.elems), concurrent(hs.concurrent)  {};

    HashMapStack* clone() {
        if(concurrent) HashMapMutex.lock();
        HashMapStack* ret =  new HashMapStack(*this);
        if(concurrent) HashMapMutex.unlock();
        return ret;
    }

    HashMapStack* clone(int lvl) {
        HashMapStack* ret = new HashMapStack();
        ret->curLevel = lvl;
        ret->concurrent = concurrent;

        if(concurrent) HashMapMutex.lock();
        for(auto it = elems.begin(); it != elems.end(); ++it) {
            auto& vect = it->second;
            for(int ind = vect.size() - 1; ind >= 0; ind--) {
                if(vect[ind].level <= lvl) {
                    ret->elems.insert({it->first, vector<StackElem>(1, vect[ind])});
                    break;
                }
            }
        }
        if(concurrent) HashMapMutex.unlock();
        return ret;
    }

    void incrementLevel() {
        if(concurrent) HashMapMutex.lock();
        curLevel++;
        if(concurrent) HashMapMutex.unlock();
    }

    void decrementLevel() {
        if(concurrent) HashMapMutex.lock();
        curLevel--;
        for(auto it = elems.begin(); it != elems.end(); ) {
            if(it->second.back().level > curLevel) it->second.pop_back();
            if(it->second.empty()) it = elems.erase(it); else ++it;
        }
        if(concurrent) HashMapMutex.unlock();
    }

    void setLevel(int lvl) {
        if(concurrent) HashMapMutex.lock();
        if(lvl < curLevel) {
            curLevel = lvl;
            for(auto it = elems.begin(); it != elems.end(); ) {
                while(!it->second.empty() && it->second.back().level > curLevel) it->second.pop_back();
                if(it->second.empty()) it = elems.erase(it); else ++it;
            }
        }
        if(concurrent) HashMapMutex.unlock();
    }

    void insert(int el, int num) {
        if(concurrent) HashMapMutex.lock();
        if(num != 0) {
            auto it = elems.find(el);
            if(it == elems.end()) elems[el].push_back(StackElem(num, curLevel));
            else if(it->second.back().level < curLevel)  it->second.push_back(StackElem(num, curLevel));
            else if(it->second.back().level == curLevel) { auto& last = it->second.back(); last.cltime = num; }
        }
        else elems.erase(el);
        if(concurrent) HashMapMutex.unlock();
    }

    bool exists(int el) {
        auto it = elems.find(el);
        if(it != elems.end()) return !!it->second.back().cltime;
        return false;
    }

    bool exists(int el, int ts) {
        auto it = elems.find(el);
        if(it == elems.end()) return false;
        int closeTime = it->second.back().cltime;
        if(closeTime == -1 || ts >= closeTime) return true;
        return false;
    }

    void include(const HashMap& other) { for(auto el : other.elems) insert(el.first, el.second); }

    int at(int el) {
        auto it = elems.find(el);
        if(it == elems.end()) return 0;
        return it->second.back().cltime;
    }

    int size() { return elems.size(); }

private:
    int curLevel = 0;
    unordered_map<int, vector<StackElem>> elems;
    bool concurrent;
};

/// Concurrent List
template<typename T>
class ConcurrentList {
public:
    ConcurrentList(bool conc = false) : elems(), concurrent(conc) {}
    ConcurrentList(const ConcurrentList& cl) : elems(cl.elems), concurrent(cl.concurrent) {}
    ConcurrentList(const ConcurrentList& cl, int len) : elems(cl.elems.begin(), cl.elems.begin()+len), concurrent(cl.concurrent) {}

    ConcurrentList* clone() {
        if(concurrent) CListMutex.lock();
        ConcurrentList* ret = new ConcurrentList(*this);
        if(concurrent) CListMutex.unlock();
        return ret;
    }

    ConcurrentList* clone(int len) {
        if(concurrent) CListMutex.lock();
        ConcurrentList* ret = new ConcurrentList(*this, len);
        if(concurrent) CListMutex.unlock();
        return ret;
    }

    void push_back(T x) {
        if(concurrent) CListMutex.lock();
        elems.push_back(x);
        if(concurrent) CListMutex.unlock();
    }

    T front() { return elems.front(); }
    T back() { return elems.back(); }

    void pop_back() {
        if(concurrent) CListMutex.lock();
        elems.pop_back();
        if(concurrent) CListMutex.unlock();
    }

    void pop_back_until(int sz = 1) {
        while(elems.size() > sz) {
            if(concurrent) CListMutex.lock();
            elems.pop_back();
            if(concurrent) CListMutex.unlock();
        }
    }

    int size() { return elems.size(); }

    template <typename TF>
    void for_each(TF&& f){ for(auto el : elems) f(el); }
    T& at(int idx){ return elems[idx]; }

    typename std::vector<T>::iterator begin() { return elems.begin(); }
    typename std::vector<T>::iterator end() { return elems.end(); }
    typename std::vector<T>::reverse_iterator rbegin() { return elems.rbegin(); }
    typename std::vector<T>::reverse_iterator rend() { return elems.rend(); }
private:
    bool concurrent = false;
    vector<T> elems;

    typedef tbb::spin_mutex CListMutexType;
    CListMutexType CListMutex;
};

#endif//_DATA_STRUCTS_H_