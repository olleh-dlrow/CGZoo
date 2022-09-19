/*
c++ scan_line.cpp -std=c++17 -o scan_line && ./scan_line

扫描线填充算法
https://www.zhihu.com/search?q=扫描线算法&utm_content=search_suggestion&type=content
*/

#include<cstdio>
#include<cmath>
#include<vector>
struct Color24 {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

struct Vec2 {
    int x;
    int y;
};

struct Polygon {
    std::vector<Vec2> vertices;
};

class Image {
public:
    Image(int _w, int _h):w(_w), h(_h) { 
        colors = new Color24[w*h];
        for(int i = 0; i < w*h; i++) {
            colors[i] = {255, 255, 255};
        }
    }
    void fill_color(int x, int y, const Color24& color) {
        int index = (h - 1 - y) * w + x;
        if(x >= 0 && y >= 0 && x < w && y < h)
            colors[index] = color;
    }
    void write(const char* filename) {
        FILE* f = fopen(filename, "w");
        fprintf(f, "P3\n%d %d\n255\n", w, h);
        for (int i = 0; i < w*h; i++) 
            fprintf(f, "%d %d %d ", colors[i].r, colors[i].g, colors[i].b);
        fclose(f);  
    }
    ~Image() { delete[] colors; }
    Color24* colors;
    int w;  int h;
};

/*
 IMPLEMENTATION
*/
struct Node {
    // 这个不是较小的x，而是较小的y对应的x，即下端x
    double x_min;
    // 斜率的倒数，与x轴平行的线条被舍弃，因此不存在无穷大的情况
    double delta_x;
    // 较大的y
    int y_max;
    Node* next = nullptr;

    bool operator<(const Node& rhs) const {
        if(x_min != rhs.x_min)return x_min < rhs.x_min;
        return delta_x < rhs.delta_x;
    }

    bool operator==(const Node& rhs) const {
        return (x_min == rhs.x_min) && (delta_x == rhs.delta_x);
    }

    bool operator>=(const Node& rhs) const {
        return !(*this < rhs);
    }
};

/*
ET表
记录线条数据，方便扫描
其中，线条中较小的y（先被扫描到的点）作为表索引
表中记录线条信息，包括：
当前x
x的变化率
终点y

每个索引记录一个有序链表，排序方式为：
先按照x排序，再按照x的变化率排序
这样能够保证线条是从左往右绘制的
*/
class EdgeTable {
public:
    EdgeTable(int size, Polygon& _poly) {
        table = std::vector<Node*>(size, 0);
        poly = &_poly;
        build();
    }

    void build() {
        int sz = poly->vertices.size();
        for(int i = 0; i < sz; i++) {
            Vec2& p1 = poly->vertices[i];
            Vec2& p2 = poly->vertices[(i+1)%sz];
            int y_min = std::min(p1.y, p2.y);
            // 与x轴平行的边不计入
            if(p1.y != p2.y) {
                Node* node = new Node {
                    .x_min = (double)((p1.y < p2.y) ? p1.x : p2.x),
                    .delta_x = (double)(p1.x - p2.x) / (double)(p1.y - p2.y),
                    .y_max = std::max(p1.y, p2.y),
                    .next = nullptr
                };
                insert(y_min, node);
            }
        }
    }

    void insert(int y_min, Node* node) {
        if(table[y_min] == nullptr) {
            table[y_min] = node;
            node->next = nullptr;
            return;
        }
        if(*table[y_min] >= *node) {
            node->next = table[y_min];
            table[y_min] = node;
            return;
        }
        Node* p = table[y_min];
        Node* pre = p;
        while(p != nullptr && *p < *node) {
            pre = p;
            p = p->next;
        }
        node->next = p;
        pre->next = node;
    }

    ~EdgeTable() {
        for (int i = 0; i < table.size(); i++)
        {
            Node* p = table[i];
            while(p) {
                Node* tmp = p;
                p = p->next;
                delete tmp;
            }
        }
        
    }

    std::vector<Node*> table;
    Polygon* poly;
};

/*
AET表，记录当前要画的线条信息
结点类型和ET表的结点类型一样
*/
class ActiveEdgeTable {
public:
    void insert(Node* node) {
        if(head == nullptr) {
            head = node;
            node->next = nullptr;
            return;
        }
        if(*head >= *node) {
            node->next = head;
            head = node;
            return;
        }
        Node* p = head;
        Node* pre = p;
        while(p != nullptr && *p < *node) {
            pre = p;
            p = p->next;
        }
        node->next = p;
        pre->next = node;
    }
    ~ActiveEdgeTable() {
        Node* p = head;
        while(p) {
            Node* tmp = p;
            p = p->next;
            delete tmp;
        }
    }

    Node* head = nullptr;
};

/*
扫描线算法
自底向上扫描，根据ET表的数据更新AET表，使用AET表画横线
步骤：
1. 遍历y = [y_min, y_max]
2. 如果ET[y]非空，将其加入到AET表中，注意都是有序插入的
3. 遍历AET表结点，判断结点的y是否等于当前y，如果是，说明该结点对应的线条已经画完，
可以将其删除
4. 如果AET表不为空，两两画线
TIPS：如何保证AET表此时一定是偶数个结点？
归纳法：假设前一个y是偶数结点，当遇到顶点时，有2种情况：
前一个线条被删除，后一个被添加，总结点数不变；
顶点对应的两个线条都是新添加的；
而初始时，顶点数为0
因此假设成立
5. 更新x，即x = x + delta_x

*/
class ScanLineAlgorithm {
public:
    ScanLineAlgorithm(Image& _img, Color24 col) {
        img = &_img;
        color = col;
    }

    void draw(Polygon& poly) {
        EdgeTable et(img->h, poly);
        ActiveEdgeTable aet;
        
        int sz = et.table.size();
        for(int i = 0; i < sz; i++) {
            // et非空时，将其全部插入aet表
            if(et.table[i] != nullptr) {
                Node* p = et.table[i];
                while(p != nullptr) {
                    et.table[i] = p->next;
                    aet.insert(p);
                    p = et.table[i];
                }
            }
            // aet表非空时，检查y_max是否等于当前y, 如果是，则删除
            if(aet.head != nullptr) {
                Node* p = aet.head;
                Node* pre = p;
                while(p != nullptr) {
                    if(p->y_max == i) {
                        if(p == aet.head) {
                            Node* tmp = p;
                            aet.head = p->next;
                            delete tmp;
                            p = aet.head;
                            pre = p;
                        } else {
                            Node* tmp = p;
                            p = p->next;
                            delete tmp;
                            pre->next = p;
                        }
                    } else {
                        pre = p;
                        p = p->next;
                    }
                }
            }
            // 两两上色
            if(aet.head != nullptr) {
                Node* p = aet.head;
                while(p != nullptr) {
                    int x1 = std::round(p->x_min);
                    int x2 = std::round(p->next->x_min);
                    for(int x = x1; x <= x2; x++) {
                        img->fill_color(x, i, color);
                    }
                    p = p->next->next;
                }
                // 使用delta_x迭代
                p = aet.head;
                while(p != nullptr) {
                    p->x_min += p->delta_x;
                    p = p->next;
                }
            }
        }
    }

    Image* img;
    Color24 color;
};

Polygon build_polygon() {
    Polygon poly;
    poly.vertices = {
        {550, 900},{300, 512}, {650, 512}, {450, 200}, {750, 570}, {400, 570}
        //{200, 200},{500, 100},{1000, 300},{1000, 800},{500, 500},{200, 700}
    };
    return poly;
}


int main() {
    Image img(1024, 1024);
    Polygon poly = build_polygon();
    ScanLineAlgorithm algo(img, {255, 0, 0});
    algo.draw(poly);
    img.write("scan_line.ppm");
    return 0;
}
