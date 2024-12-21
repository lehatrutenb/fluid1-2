#include <vector>

template<class T> // primary template
void sort(std::vector<T>& v) {

}
template<>        // specialization for T = int
void sort(std::vector<int>&) {
    
}
int main()
{
}