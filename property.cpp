#include<functional>
#include<iostream>
template<typename T, bool RO=false, bool WO=false>
class Property {
public:
    template<typename ...Args>
    Property(Args&&... args):t(std::forward<Args>(args)...){
        getter = std::bind(&Property::get,this);
        setter = std::bind(&Property::set,this,std::placeholders::_1);
    }

    operator const T& () const {
        static_assert(!WO, "Cannot access write-only property.");
        return getter();
    }
    const T& operator = (const T& other) {
        static_assert(!RO, "Cannot set read-only property.");
        return setter(other);
    }
    bool operator == (const T& other) const {
        // Static cast makes sure our getter operator is called, so we could use overrides if those are in place
        return static_cast<const T&>(*this) == other;
    }

    // Assign getter and setter to these properties
    std::function<const T&()> getter;
    std::function<const T&(const T&)> setter;

    // Use this to always get without overrides, useful for use with overriding implementations
    const T& get() const {
        return t;
    }
    // Use this to always set without overrides, useful for use with overriding implementations
    const T& set(const T& other) {
        return t = other;
    }

private:
    T t;
};

class Test {
public:
    Test(int _v=3):
        val(_v) {
        val.setter = std::bind(&Test::set, this, std::placeholders::_1);
        val.getter = std::bind(&Test::get, this);
    }
    Property<int> val;
    
    const int& get() const {
        std::cout << "get\n";
        return val.get();
    }

    const int& set(const int& v) {
        std::cout << "set\n";
        return val.set(v);
    }
};

int main() {
    Test t(10);
    std::cout << t.val << std::endl;
    t.val = 14;
    return 0;
}
