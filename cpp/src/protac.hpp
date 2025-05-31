#pragma once
#include <vector>
#include <array>
#include <memory>
#include <tuple>
#include <random>
#include <future>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>  
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include "protein.hpp"
#include "getgriden.hpp"
#include "minimize_h.hpp"
#include "align.hpp"
#include "parameters.hpp"


class Protac {
public:
    // 构造函数
    Protac(const GRID& grid_anchor, const GRID& grid_flex, int processes = 1);
    
    // 内部结构体定义
    struct Solution {
        std::vector<double> dihedrals;
        double energy;
        std::vector<double> parameters;  // 参数向量
        bool operator<(const Solution& other) const {
            return energy < other.energy;
        }
    };    
    

    void init(RDKit::ROMol* protac,
              RDKit::ROMol* w_anch,
              RDKit::ROMol* w_flex,
              const std::string& fpro_flex,
              bool verbose = false);
    // 采样相关函数
    Solution sample_single(bool verbose = false);
    void sample(int ntotal = 100, int nsolu = 100, bool verbose = false);
    void search(bool verbose = false);
    double e_intra(const RDKit::ROMol* mol) const;
    double score(const std::vector<double>& dihe);
    
    // 新增score_only功能
    double score_only(bool verbose = false);
    double score_only(const std::vector<double>& dihe, bool verbose = false);
    
    // 公有访问器方法
    const std::shared_ptr<RDKit::ROMol>& getProtac() const { return protac_; }
    const std::vector<std::array<int, 4>>& getRotatableDihedrals() const { return rot_dihe_; }
    std::vector<Solution>& getSolutions() { return solutions_; }
    void clearSolutions() { solutions_.clear(); }
    void addSolution(const Solution& solution) { solutions_.push_back(solution); }
    
    void output(RDKit::SDWriter& w, 
                std::ostream& fpro, 
                int nKeep,
                bool fpro_w = false,
                double rmsd_cutoff = 1.0);
    void printProtacInfo();

private:
    // 成员变量
    GRID grid_anchor_;
    GRID grid_flex_;
    int processes_;

    std::shared_ptr<RDKit::ROMol> protac_;
    std::vector<std::array<int, 4>> rot_dihe_;
    double E_intra_ref_;
    std::vector<Solution> solutions_;

    
    Protein protein_;

    Coords coord_subs_var; //coord_subs_var = moveToOrigin(coords_flex_conf);
    
    std::vector<int> idx_;
    std::vector<IQHb> q_flex_;
    std::vector<IQHb> q_anchor_;
    std::array<double, 3> translation_;
    std::vector<int> flex_warhead_atoms_;  // 存储弹头原子的索引，用于RMSD计算
    
    // 新增：保存anchor warhead对齐信息
    std::vector<int> anchor_warhead_atoms_;  // 保存anchor warhead原子的索引
    
    struct VdwParam {
        int atom1_idx;
        int atom2_idx;
        double param1;
        double param1_1;
        double param2;
        double param3;
        double param4;
        double param5;
    };
    
    struct DihedralParam {
        int atom1_idx;
        int atom2_idx;
        int atom3_idx;
        int atom4_idx;
        double v1;
        double v2;
        double v3;
        std::vector<double> lowEnerDiheAngles;
    };

    std::vector<VdwParam> list_vdw_;
    std::vector<DihedralParam> list_dihe_;

    void list(const std::vector<int>& warheads, 
              const std::vector<std::pair<int, int>>& rbond,
              bool print_info = false);
    
    // 对齐相关的私有函数
    void alignProtacToFlexWarhead(RDKit::ROMol* w_flex, bool verbose = false);
    void alignProtacToAnchorWarhead(RDKit::ROMol* w_anch, bool verbose = false);
    
    // 可旋转键和二面角查找函数
    void findRotatableDihedrals(const std::vector<int>& linker, bool verbose = false);
    
    void calHeavyAtomsCharge(RDKit::ROMol& mol);
    // 电荷和氢键类型计算函数
    void calculateQAnchor(const std::vector<int>& hb_donors, const std::vector<int>& hb_acceptors, bool verbose = false);
    void calculateQFlex(const std::vector<int>& hb_donors, const std::vector<int>& hb_acceptors, bool verbose = false);
    
    // 氢键供受体识别函数
    std::pair<std::vector<int>, std::vector<int>> findHB_DA(bool verbose = false);
    
    // 坐标系统设置函数
    void setupCoordinateSystem(RDKit::ROMol* w_flex, const std::string& fpro_flex, bool verbose = false);
    
    // Linker原子识别函数
    std::vector<int> findLinkerAtoms(bool verbose = false);
    
    // 能量组分详细输出的辅助函数
    void printEnergyComponents(const RDKit::ROMol& mol, 
                              const std::vector<double>& dihedrals, 
                              double total_energy, 
                              bool is_custom_dihedrals);
    
    // 扫描单个二面角的低能构象
    std::vector<double> scanTorsion(double v1, double v2, double v3);
    std::vector<double> listTorsion(double grid_step);
    
    // 随机数生成器
    std::mt19937 rng_;
    std::uniform_real_distribution<double> angle_dist_;
    std::uniform_real_distribution<double> uniform_dist_;

    // Powell最小化的辅助函数
    Solution powell_minimize(const std::vector<double>& initial_guess, RDKit::ROMol* mol_copy, double tol = 0.01);
    Solution search_single(const Solution& initial_solution);

    //
    thread_local static std::unique_ptr<RDKit::ROMol> thread_local_mol_;
    double thread_safe_score(const std::vector<double>& dihe, RDKit::ROMol* mol_copy, bool print_info = false); // 只声明

    static std::mutex output_mutex_;

    // Helper methods for conformer clustering
    double computeRMSD(const RDKit::ROMol& mol, int confId1, int confId2);
    std::vector<std::vector<int>> clusterConformers(const RDKit::ROMol& mol, double rmsdThreshold);
    

};

class ThreadPool {
public:
    ThreadPool(size_t num_threads);
    ~ThreadPool();

    template<class F>
    auto enqueue(F&& f) -> std::future<typename std::result_of<F()>::type>;

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

class ProgressBar {
public:
    ProgressBar(int total) : total_(total), current_(0), lastProgress_(0) {}

    void update(int completed) {
        std::lock_guard<std::mutex> lock(mutex_); // 确保线程安全
        current_ = completed; // 更新当前完成的任务数
        int progress = static_cast<int>(static_cast<double>(current_) / total_ * 100 + 0.5); // 计算进度百分比，加0.5进行四舍五入
        progress = std::min(progress, 100); // 确保百分比不会超过100
        int barWidth = 50; // 进度条的宽度

        if (progress != lastProgress_) { // 只有当进度变化时才更新进度条
            // 计算已填充的部分
            int filledWidth = progress * barWidth / 100;
            std::cout << "\r[" << std::string(filledWidth, '=') 
                      << std::string(barWidth - filledWidth, ' ') // 填充剩余部分
                      << "] " << progress << "%"; // 只显示一次百分比
            std::cout.flush();
            lastProgress_ = progress;
        }
    }

    void finish() {
        std::lock_guard<std::mutex> lock(mutex_); // 确保线程安全
        // 完成时确保输出100%且没有多余的百分号
        std::cout << "\r[" << std::string(50, '=') 
                  << "] 100%" << std::endl;
    }

private:
    int total_;
    int current_;
    int lastProgress_; // 上次显示的进度
    std::mutex mutex_; // 互斥锁
};
