#include <cstddef>
#include <cstdio>
#include <memory>
#include <sys/mman.h>
#include <type_traits>
#include <vector>
#define ALLOC_LEN 0x40000000
template<typename T>
struct allocator_state {
    T* start_addr;
    T* storage_end;
    T* allocated_end;
    /*#ifndef NDEBUG
    std::vector<T*> last_allocate_ptr;
    #endif*/
    allocator_state() {
        start_addr=(T*) mmap(0,ALLOC_LEN,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
        storage_end=(T*) ((char*) start_addr+ALLOC_LEN);
        allocated_end=start_addr;
    }
    ~allocator_state() {
        munmap(start_addr, ALLOC_LEN);
    }
};
template<typename T>
struct stack_allocator {
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::false_type propagate_on_container_move_assignment;
    typedef std::false_type is_always_equal;
    template<typename U>
    struct rebind {
        typedef stack_allocator<U> other;
    };
  private:
    allocator_state<T>& state;
  public:
    stack_allocator(allocator_state<T>& state):state(state) {}
    T* allocate(size_t n) {
        auto to_return=state.allocated_end;
        state.allocated_end+=n;
        /*#ifndef NDEBUG
                    state.last_allocate_ptr.push_back(to_return);
        #endif*/
        return to_return;
    }
    void deallocate(T* p,size_t n) {
        /*#ifndef NDEBUG
                    assert(p==state.last_allocate_ptr.back());
                    assert(p+n==state.allocated_end);
                    state.last_allocate_ptr.pop_back();
        #endif*/
        state.allocated_end-=n;
    }
    bool empty()const {
        return(state.start_addr==state.allocated_end);
    }
};
