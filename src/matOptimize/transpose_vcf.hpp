#include <cstdint>
struct Pos_Mut{
    int position;
    uint8_t mut;
    bool operator<(const Pos_Mut& other)const{
        return position<other.position;
    }
};