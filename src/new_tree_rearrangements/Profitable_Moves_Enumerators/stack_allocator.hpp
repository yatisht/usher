#include <cstddef>
#include <sys/mman.h>
#include <type_traits>
#include <vector>
#define ALLOC_LEN 0x4000
template<typename T>
struct stack_allocator{
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::false_type propagate_on_container_move_assignment;
    typedef std::false_type rebind;
    typedef std::false_type is_always_equal;
    private:
        T* start_addr;
        T* storage_end;
        T* allocated_end;
        #ifndef NDEBUG
        std::vector<T*> last_allocate_ptr;
        #endif
    public:
        stack_allocator(){
            start_addr=(T*) mmap(0,ALLOC_LEN,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
            storage_end=(char*) start_addr+ALLOC_LEN;
            allocated_end=start_addr;
        }
        ~stack_allocator(){
            munmap(start_addr, ALLOC_LEN);
        }
        T* allocate(size_t n){
            auto to_return=allocated_end;
            allocated_end+=n;
#ifndef NDEBUG
            last_allocated_ptr.push_back(to_return);
#endif
            return to_return;
        }
        void deallocate(T* p,size_t n){
#ifndef NDEBUG
            assert(p==last_allocate_ptr.back());
            assert(p+n==allocated_end);
#endif
            allocated_end-=n;
        }

};