#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <time.h>
#include <sys/time.h>

using namespace std;

FILE* netFile;
FILE* solutionFile;

int Grid_H_numb, Grid_V_numb;
int H_capacity, V_capacity;
int NumbNets;

struct pin {
    int x;
    int y;
};

struct path {
    int x;
    int y;
    path(int _x, int _y) : x(_x), y(_y) {}
};

struct net {
    int name;
    int numbOfNodesInBox;
    vector<pin> pins;
    vector<path> paths;
};

struct weightG {
    int original_index;
    int adjacent_id;
    int weight;
    weightG(int index, int id, int w) : original_index(index), adjacent_id(id), weight(w) {}
};

net* nets;
vector<net> Net_vec;
vector<vector<weightG>> WeightStateGraph;

int** map;
int** mazeMap;

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)) return 0;
    
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time() {
    return (double)clock() / CLOCKS_PER_SEC;
}

// --------------------- Paring Files ---------------------
bool parseFileName(int argc, char** argv) {
    if (argc != 3 || !strstr(argv[1], ".modified.txt") || !strstr(argv[2], ".result")) {
        cerr << "Commend is Wrong!!!\n";
        return 0;
    }
    
    if(!(netFile = fopen(argv[1], "r"))) {
        cerr << " Can not open .modified.txt file\n";
        return 0;
    }
    
    if(!(solutionFile = fopen(argv[2], "r"))) {
        cerr << " Can not open .result file\n";
        return 0;
    }
    
    return 1;
}

void parseNets(FILE* file) {
    fscanf(file, "grid %d %d\n", &Grid_H_numb, &Grid_V_numb);
    fscanf(file, "vertical capacity %d\n", &V_capacity);
    fscanf(file, "horizontal capacity %d\n", &H_capacity);
    fscanf(file, "num net %d\n", &NumbNets);
    
    cout << "grid H: " << Grid_H_numb << " grid V: " << Grid_V_numb << endl;
    cout << "V capacity: " << V_capacity << endl;
    cout << "H capacity: " << H_capacity << endl;
    cout << "NumbNets: " << NumbNets << endl;
    
    nets = new net[NumbNets];
    map = new int*[Grid_V_numb];
    mazeMap = new int*[Grid_V_numb];
    
    for (int i = 0; i < Grid_V_numb; i++) {
        map[i] = new int[Grid_H_numb];
        mazeMap[i] = new int[Grid_H_numb];
    }
    
    for (int i = 0; i < Grid_V_numb; i++) {
        for (int j = 0; j < Grid_H_numb ; j++) {
            map[i][j] = 0;
            mazeMap[i][j] = -1;
        }
    }
    
    for (int i = 0; i < NumbNets; i++) {
        pin tmp1, tmp2;
        fscanf(file, "net%d %d 2\n", &nets[i].name, &nets[i].name);
        fscanf(file, "%d %d\n", &tmp1.x, &tmp1.y);
        fscanf(file, "%d %d\n", &tmp2.x, &tmp2.y);
        nets[i].pins.push_back(tmp1);
        nets[i].pins.push_back(tmp2);
        map[tmp1.y][tmp1.x] = 1;
        map[tmp2.y][tmp2.x] = 1;
    }
}

void writeSolution(FILE* file) {
    for (int i = 0; i < NumbNets; i++) {
        fprintf(file, "net%d %d\n", Net_vec[i].name, Net_vec[i].name);
        for (int j = 0; j < Net_vec[i].paths.size() - 1; j++) {
            fprintf(file, "(%d, %d, 1)-(%d, %d, 1)\n", Net_vec[i].paths[j].x, Net_vec[i].paths[j].y, Net_vec[i].paths[j + 1].x, Net_vec[i].paths[j + 1].y);
        }
        fprintf(file, "!\n");
    }
}

// --------------------- Global Routing ---------------------
void countNodesNumb() {
    for (int i = 0; i < NumbNets; i++) {
        int minX = min(nets[i].pins[0].x, nets[i].pins[1].x);
        int minY = min(nets[i].pins[0].y, nets[i].pins[1].y);
        int maxX = max(nets[i].pins[0].x, nets[i].pins[1].x);
        int maxY = max(nets[i].pins[0].y, nets[i].pins[1].y);
        int block_num = 0;
        
        for (int row = minY; row <= maxY; row++) {
            for (int col = minX; col <= maxX; col++) {
                if (!(col == minX && row == minY)) {
                    if (!(col == maxX && row == maxY)) {
                        if (map[row][col] == 1) block_num++;
                    }
                }
            }
        }
        
        nets[i].numbOfNodesInBox = block_num;
    }
}

vector<net> createVector() {
    vector<net> colletions;
    
    for (int i = 0; i < NumbNets; i++) {
        net tmp;
        tmp.name = nets[i].name;
        tmp.numbOfNodesInBox = nets[i].numbOfNodesInBox;
        tmp.pins = nets[i].pins;
        colletions.push_back(tmp);
    }
    
    return colletions;
}

void quickSort1(int L, int R) {
    int i, j, mid;
    float piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = Net_vec[mid].numbOfNodesInBox;

    while (i < R || j > L) {
        while (Net_vec[i].numbOfNodesInBox < piv) i++;
        while (Net_vec[j].numbOfNodesInBox > piv) j--;

        if (i <= j) {
            swap(Net_vec[i], Net_vec[j]);
            i++;
            j--;
        } else {
            if (i < R) quickSort1(i, R);
            if (j > L) quickSort1(L, j);
            return;
        }
    }
}

void quickSort2(int L, int R) {
    int i, j, mid;
    float piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = Net_vec[mid].name;

    while (i < R || j > L) {
        while (Net_vec[i].name < piv) i++;
        while (Net_vec[j].name > piv) j--;

        if (i <= j) {
            swap(Net_vec[i], Net_vec[j]);
            i++;
            j--;
        } else {
            if (i < R) quickSort2(i, R);
            if (j > L) quickSort2(L, j);
            return;
        }
    }
}

void initialWeightGraph() {
    for (int i = 0; i < Grid_V_numb; i++) {
        for (int j = 0; j < Grid_H_numb; j++) {
            vector<weightG> tmp;
            if (j - 1 >= 0) tmp.push_back(weightG(nets[i * Grid_H_numb + j].name, i * Grid_H_numb + j - 1, 1));  // left
            if (j + 1 < Grid_H_numb) tmp.push_back(weightG(nets[i * Grid_H_numb + j].name, i * Grid_H_numb + j + 1, 1));   // right
            if (i - 1 >= 0) tmp.push_back(weightG(nets[i * Grid_H_numb + j].name, (i - 1) * Grid_H_numb + j, 1)); // top
            if(i + 1 < Grid_V_numb) tmp.push_back(weightG(nets[i * Grid_H_numb + j].name, (i + 1) * Grid_H_numb + j, 1));   // below
            WeightStateGraph.push_back(tmp);
        }
    }
}

void updateEdgeWeight(int start, int target) {
    int index, index0;
    
    for (index = 0; index < WeightStateGraph[start].size(); index++) {
        if (WeightStateGraph[start][index].adjacent_id == target) {
            WeightStateGraph[start][index].weight++;
            break;
        } else if (index == WeightStateGraph[start].size() - 1 && WeightStateGraph[start][index].adjacent_id != target) {
            cout << "cannot find corresponding neighbor!!!----1" << endl;
        }
    }
    
    for (index0 = 0; index0 < WeightStateGraph[target].size(); index0++) {
        if (WeightStateGraph[target][index0].adjacent_id == start) {
            WeightStateGraph[target][index0].weight++;
            break;
        } else if (index0 == WeightStateGraph[target].size() - 1 && WeightStateGraph[target][index0].adjacent_id != start) {
            cout << "cannot find corresponding neighbor!!!----2" << endl;
        }
    }
}

void updateWeight(int id) {
    for (int i = 0; i < Net_vec[id].paths.size() - 1; i++) {
        if (Net_vec[id].paths[i].x == Net_vec[id].paths[i + 1].x) {
            int x_cor = Net_vec[id].paths[i].x;
            
            for (int j = Net_vec[id].paths[i].y; j < Net_vec[id].paths[i + 1].y - 1; j++) {
                int start = x_cor + j * Grid_H_numb;
                int target = x_cor + (j + 1) * Grid_H_numb;
                
                updateEdgeWeight(start, target);
            }
        } else {
            int y_cor = Net_vec[id].paths[i].y;
            
            for (int j = Net_vec[id].paths[i].x; j < Net_vec[id].paths[i + 1].x - 1; j++) {
                int start = j + y_cor * Grid_H_numb;
                int target = j + 1 + y_cor * Grid_H_numb;
                
                updateEdgeWeight(start, target);
            }
        }
    }
}

void updateMazeMap(int start, int target) {
    int min_x = min(start % Grid_H_numb, target % Grid_H_numb);
    int max_x = max(start % Grid_H_numb, target % Grid_H_numb);
    int min_y = min(start / Grid_H_numb, target / Grid_H_numb);
    int max_y = max(start / Grid_H_numb, target / Grid_H_numb);
    int minX = (min_x - 3 >= 0) ? (min_x - 3) : 0;
    int maxX = (max_x + 3 < Grid_H_numb) ? (max_x + 3) : (Grid_H_numb - 1);
    int minY = (min_y - 3 >= 0) ? (min_y - 3) : 0;
    int maxY = (max_y + 3 < Grid_V_numb) ? (max_y + 3) : (Grid_V_numb - 1);
    
    for (int i = 0; i < Grid_V_numb; i++) {
        for (int j = 0; j < Grid_H_numb ; j++) {
            mazeMap[i][j] = -1;
        }
    }
        
    for (int row = minY; row <= maxY; row++) {
        for (int col = minX; col <= maxX; col++) {
            int tmp = (abs(col - start % Grid_H_numb) + abs(row - start / Grid_H_numb)) % 3;
            mazeMap[row][col] = (tmp == 0) ? 3 : tmp;
        }
    }
}

void findPath(int start, int target, int id) {
    stack<int> corner_id;
    bool reachTarget = false;
    int backTracePoint = target;
    int cur_cornet_id = target;
    updateMazeMap(start, target);
    
    corner_id.push(target);
    
    while (!reachTarget) {
        int cur_id = backTracePoint;
        int cur_index = mazeMap[backTracePoint / Grid_H_numb][backTracePoint % Grid_H_numb];
        int next_maze_index = (cur_index == 1) ? 3 : (cur_index - 1);
        int next_id = 0;
        int minW = 10000;
                
        for (int i = 0; i < WeightStateGraph[cur_id].size(); i++) {
            int tmp_x = WeightStateGraph[cur_id][i].adjacent_id % Grid_H_numb;
            int tmp_y = WeightStateGraph[cur_id][i].adjacent_id / Grid_H_numb;
            
            if (mazeMap[tmp_y][tmp_x] == next_maze_index && WeightStateGraph[cur_id][i].weight < minW) {
                next_id = WeightStateGraph[cur_id][i].adjacent_id;
                minW = WeightStateGraph[cur_id][i].weight;
                backTracePoint = next_id;
            }
        }
                
        updateEdgeWeight(cur_id, next_id);
        
        int cur_x = cur_cornet_id % Grid_H_numb;
        int cur_y = cur_cornet_id / Grid_H_numb;
        int next_x = backTracePoint % Grid_H_numb;
        int next_y = backTracePoint / Grid_H_numb;
        
        if (cur_x != next_x && cur_y != next_y) {
            corner_id.push(cur_id);
            cur_cornet_id = cur_id;
        }
        
        if (next_id == start) reachTarget = true;
    }
    
    corner_id.push(start);
    
    while (!corner_id.empty()) {
        int tmp_id = corner_id.top(); corner_id.pop();
        Net_vec[id].paths.push_back(path(tmp_id % Grid_H_numb, tmp_id / Grid_H_numb));
    }
}

void Global_Routing() {
    countNodesNumb();
    initialWeightGraph();
    Net_vec = createVector();
    quickSort1(0, NumbNets - 1);
    
    for (int i = 0; i < Net_vec.size(); i++) {
        cout << i << endl;
        if (Net_vec[i].numbOfNodesInBox == 0) {    // Two pins on same row or column
            if (Net_vec[i].pins[0].x == Net_vec[i].pins[1].x || Net_vec[i].pins[0].y == Net_vec[i].pins[1].y) {
                Net_vec[i].paths.push_back(path(Net_vec[i].pins[0].x, Net_vec[i].pins[0].y));
                Net_vec[i].paths.push_back(path(Net_vec[i].pins[1].x, Net_vec[i].pins[1].y));
            } else {                            // Two pins has L relationship
                Net_vec[i].paths.push_back(path(Net_vec[i].pins[0].x, Net_vec[i].pins[0].y));
                Net_vec[i].paths.push_back(path(Net_vec[i].pins[1].x, Net_vec[i].pins[0].y));
                Net_vec[i].paths.push_back(path(Net_vec[i].pins[1].x, Net_vec[i].pins[1].y));
            }
            updateWeight(i);
        } else {
            findPath(Net_vec[i].pins[0].x + Net_vec[i].pins[0].y * Grid_H_numb, Net_vec[i].pins[1].x + Net_vec[i].pins[1].y * Grid_H_numb, i);
        }
    }

    quickSort2(0, NumbNets - 1);
}

// --------------------- Main ---------------------
int main(int argc, char * argv[]) {
    if (!parseFileName(argc, argv)) return 0;
    
    parseNets(netFile);
    fclose(netFile);
    
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    
    cout << "start global_routing" << endl;
    Global_Routing();
    
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    
    cout << "\nTime(s) : Abacus : Wall Time : " << wall1 - wall0 << ", CPU Time : " << cpu1  - cpu0 << std::endl;
    
    writeSolution(solutionFile);
    fclose(solutionFile);
    return 0;
}
