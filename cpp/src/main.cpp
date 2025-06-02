#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <limits>
#include <cstdio>  // for std::remove
#include <sstream>
#include <memory>  // for std::unique_ptr and std::make_unique
#include <iomanip>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <RDGeneral/Invariant.h>
#include <Geometry/point.h>

#include "ternify.hpp"
#include "protac.hpp"
#include "grid.hpp"


void read_parameters(const std::string& filename, Parameters& params) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open parameter file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        // 去除行首行尾的空白字符
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        
        // 跳过空行
        if (line.empty()) {
            continue;
        }
        
        // 跳过以#开头的注释行
        if (line[0] == '#') {
            continue;
        }
        
        // 处理行内注释：找到#的位置，截取#之前的内容
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
            // 再次去除可能的尾部空白字符
            line.erase(line.find_last_not_of(" \t") + 1);
        }
        
        // 跳过处理注释后变为空的行
        if (line.empty()) {
            continue;
        }
        
        std::istringstream iss(line);
        std::string key, value;
        
        if (std::getline(iss, key, ':') && std::getline(iss, value)) {
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            if (key == "Interface") {
                std::istringstream value_stream(value);
                std::string number;
                std::vector<double> numbers;
                
                // 解析所有数字
                while (std::getline(value_stream, number, ',')) {
                    number.erase(0, number.find_first_not_of(" \t"));
                    number.erase(number.find_last_not_of(" \t") + 1);
                    if (!number.empty()) {
                        numbers.push_back(std::stod(number));
                    }
                }

                // 检查是否有足够的数字
                if (numbers.size() >= 6) {
                    for (int i = 0; i < 3; ++i) {
                        params.interface[i][0] = numbers[i * 2];
                        params.interface[i][1] = numbers[i * 2 + 1];
                    }
                } else {
                    throw std::runtime_error("Interface parameter requires 6 numbers");
                }
            }
            else if (key == "PROTACs") params.protacs_file = value;
            else if (key == "Warhead_anchor") params.warhead_anchor_file = value;
            else if (key == "Warhead_flex") params.warhead_flex_file = value;
            else if (key == "Protein_anchor") params.protein_anchor_file = value;
            else if (key == "Protein_flex") params.protein_flex_file = value;
            else if (key == "Output_protac") params.output_protac_file = value;
            else if (key == "Output_protein") params.output_protein_file = value;
            else if (key == "RMSD_cutoff") params.output_rmsd_cutoff = std::stod(value);
            else if (key == "N_ini") params.n_ini = std::stoi(value);
            else if (key == "N_search") params.n_search = std::stoi(value);
            else if (key == "N_keep") params.n_keep = std::stoi(value);
            else if (key == "N_processes") params.n_processes = std::stoi(value);
            else if (key == "Verbose") params.verbose = std::stoi(value);
            else if (key == "Score_only") params.score_only = (std::stoi(value) != 0);
        }
    }
}

void run_ternify(const Parameters& params) {
    // 读取warhead分子
    RDKit::SDMolSupplier w_anch_supplier(params.warhead_anchor_file,true,true,true); //读取的时候，分子并没有较好地检查，应该是bug
    RDKit::SDMolSupplier w_flex_supplier(params.warhead_flex_file,true,true, true);
    
    if (w_anch_supplier.atEnd() || w_flex_supplier.atEnd()) {
        throw std::runtime_error("Failed to read warhead molecules");
    }
    
    std::unique_ptr<RDKit::ROMol> w_anch(w_anch_supplier.next());
    std::unique_ptr<RDKit::ROMol> w_flex(w_flex_supplier.next());

    if (w_anch && w_flex) {
            std::cout << "READ: Number of atoms in WANCHOR: " << w_anch->getNumAtoms() << std::endl;
            std::cout << "READ: Number of bonds in WANCHOR: " << w_anch->getNumBonds() << std::endl;
            std::cout << "READ: Number of atoms in WFLEX: " << w_flex->getNumAtoms() << std::endl;
            std::cout << "READ: Number of bonds in WFLEX: " << w_flex->getNumBonds() << std::endl;
    }
    else {
        throw std::runtime_error("Failed to read warhead molecules");
    }
    
    // 检查是否需要输出蛋白质坐标
    bool write_protein = !params.output_protein_file.empty() && 
                        params.output_protein_file != "None" && 
                        params.output_protein_file != "none";
    
    // 在程序开始时删除已存在的输出文件（重建）
    std::cout << "\nPreparing output files..." << std::endl;
    
    // 删除PROTAC输出文件（如果存在）
    if (std::remove(params.output_protac_file.c_str()) == 0) {
        std::cout << "Removed existing PROTAC output file: " << params.output_protac_file << std::endl;
    }
    
    // 删除蛋白质输出文件（如果存在且需要输出蛋白质）
    if (write_protein) {
        if (std::remove(params.output_protein_file.c_str()) == 0) {
            std::cout << "Removed existing protein output file: " << params.output_protein_file << std::endl;
        }
        std::cout << "PROTAC output: " << params.output_protac_file << std::endl;
        std::cout << "Protein output: " << params.output_protein_file << std::endl;
    } else {
        std::cout << "PROTAC output: " << params.output_protac_file << std::endl;
        std::cout << "Protein output disabled (Output_protein = None or not specified)" << std::endl;
    }

    // 创建网格
    std::cout << "\nCreating grid for anchor_pro..." << std::endl;
    std::unique_ptr<GRID> grid_anchor;
    bool grid_anchor_initialized = false;
    try {
        grid_anchor = std::make_unique<GRID>(Grid(params.protein_anchor_file, params.interface, params.n_processes));
        grid_anchor_initialized = true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create grid for anchor_pro: " << e.what() << std::endl;
    }
    // 检查网格是否已初始化
    if (!grid_anchor_initialized) {
        std::cerr << "Warning: grid_anchor was not initialized." << std::endl;
    }
    
    // 计算柔性蛋白网格的位置
    std::cout << "Creating grid for flex_pro..." << std::endl;
    const RDKit::Conformer& conf_flex = w_flex->getConformer();
    VolRegion site_flex;
    
    for (int i = 0; i < 3; ++i) {
        double min_val = std::numeric_limits<double>::max();
        double max_val = std::numeric_limits<double>::lowest();
        
        for (unsigned int j = 0; j < conf_flex.getNumAtoms(); ++j) {
            const RDGeom::Point3D& pos = conf_flex.getAtomPos(j);
            double coord = (i == 0) ? pos.x : ((i == 1) ? pos.y : pos.z);
            min_val = std::min(min_val, coord);
            max_val = std::max(max_val, coord);
        }
        
        site_flex[i][0] = min_val - 5.0;  // 使用数组索引而不是 push_back
        site_flex[i][1] = max_val + 5.0;    
    }

    std::unique_ptr<GRID> grid_flex;
    bool grid_flex_initialized = false;
    try {
        grid_flex = std::make_unique<GRID>(Grid(params.protein_flex_file, site_flex, params.n_processes));
        grid_flex_initialized = true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create grid for flex_pro: " << e.what() << std::endl;
    }


    if (!grid_flex_initialized) {
        std::cerr << "Warning: grid_flex was not initialized." << std::endl;
    }
    
    // 处理PROTAC分子
    std::cout << "Reading protac file..." << std::endl;
    RDKit::SDMolSupplier protacs_supplier(params.protacs_file,true,true,false);
    size_t total_protacs = protacs_supplier.length();
    std::cout << "There is(are) " << protacs_supplier.length() << " protac molecule(s) to be processed..." << std::endl;
    
    // 创建输出文件流（在循环外创建，用于append模式）
    std::unique_ptr<std::ofstream> protac_file_stream;
    std::unique_ptr<RDKit::SDWriter> protac_writer;
    std::unique_ptr<std::ofstream> protein_file_stream;
    
    std::unique_ptr<RDKit::ROMol> mol;
    for (size_t protac_index = 0; protac_index < total_protacs; ++protac_index) {
        mol = std::unique_ptr<RDKit::ROMol>(protacs_supplier[protac_index]);
        // Check if the molecule was successfully read and matches the substructure
        if (mol) {
            std::cout << "\n=== Processing PROTAC molecule " << (protac_index + 1) << " of " << total_protacs << " ===" << std::endl;
            std::cout << "READ: Number of atoms in protac: " << mol->getNumAtoms() << std::endl;
            std::cout << "READ: Number of bonds in protac: " << mol->getNumBonds() << std::endl;
            
            if (RDKit::SubstructMatch(*mol, *w_anch).empty() || 
                RDKit::SubstructMatch(*mol, *w_flex).empty() ){
                std::cout << "WARNING: Skipping molecule " << protac_index + 1 << ": this protac does not match substructures with anch/flex." << std::endl;
                continue;
            }
            
            
            // 检查网格是否都已成功初始化
            if (!grid_anchor_initialized || !grid_flex_initialized) {
                std::cerr << "Error: Cannot process PROTAC " << (protac_index + 1) 
                         << " because grids are not properly initialized." << std::endl;
                continue;
            }
            
            // 为第一个PROTAC分子创建输出文件，后续分子使用append模式
            if (protac_index == 0 || !protac_writer) {
                // 第一个分子：创建新文件
                protac_file_stream.reset(new std::ofstream(params.output_protac_file));
                protac_writer.reset(new RDKit::SDWriter(protac_file_stream.get()));
                
                if (write_protein) {
                    protein_file_stream.reset(new std::ofstream(params.output_protein_file));
                }
            } else {
                // 后续分子：关闭当前writer，以append模式重新打开
                protac_writer.reset();  // 关闭当前writer
                protac_file_stream.reset();  // 关闭当前文件流
                
                // 以append模式打开文件
                protac_file_stream.reset(new std::ofstream(params.output_protac_file, std::ios::app));
                protac_writer.reset(new RDKit::SDWriter(protac_file_stream.get()));
                
                if (write_protein && protein_file_stream) {
                    protein_file_stream.reset();  // 关闭当前writer
                    protein_file_stream.reset(new std::ofstream(params.output_protein_file, std::ios::app));  // append模式
                }
            }
            
            Protac PROTac(*grid_anchor, *grid_flex, params.n_processes);
            std::cout << "Initializing protac..." << std::endl;
            // Initialize Protac object with the current molecule and other parameters
            PROTac.init(mol.get(), w_anch.get(), w_flex.get(),
                    params.protein_flex_file, params.verbose > 0);
            // Print Protac information based on verbose setting
            if (params.verbose > 0) {
                std::cout << "Verbose mode enabled (level: " << params.verbose << ")" << std::endl;
                PROTac.printProtacInfo();
            }

            if (params.score_only) {
                std::cout << "Score Only Mode Enabled..." << std::endl;
                double energy = PROTac.score_only(params.verbose > 0);
                std::cout << "Final energy score: " << std::fixed << std::setprecision(3) << energy << " kcal/mol" << std::endl;
                
                // 在score_only模式下，我们不需要采样和搜索，但可以输出当前构象
                std::cout << "Saving current conformation to output..." << std::endl;
                // 创建一个简单的solution用于输出
                Protac::Solution current_solution;
                const RDKit::Conformer& conf = PROTac.getProtac()->getConformer();
                const auto& rot_dihe = PROTac.getRotatableDihedrals();
                for (size_t i = 0; i < rot_dihe.size(); ++i) {
                    const auto& atoms = rot_dihe[i];
                    double angle = MolTransforms::getDihedralDeg(conf, atoms[0], atoms[1], atoms[2], atoms[3]);
                    current_solution.dihedrals.push_back(angle);
                }
                current_solution.energy = energy;
                
                // 清空solutions_并添加当前解
                PROTac.clearSolutions();
                PROTac.addSolution(current_solution);
            } else {
                std::cout << "Sampling and Searching..." << std::endl;
                PROTac.sample(params.n_ini, params.n_search, params.verbose > 0);
            }
        
            std::cout << "Writing output..." << std::endl;
            // Write the output to the specified writer
            PROTac.output(*protac_writer, protein_file_stream ? *protein_file_stream : std::cout, params.n_keep, write_protein, params.output_rmsd_cutoff);            
        } 
        else {
            std::cout << "Warning: Molecule " << protac_index + 1 << " is invalid or could not be read." << std::endl;
        }
        
    }

}

int main(int argc, char* argv[]) {
    try {
        // 命令行参数解析
        std::string param_file = "tcs.inp";
        if (argc > 1) {
            for (int i = 1; i < argc; i++) {
                std::string arg = argv[i];
                if (arg == "-p" && i + 1 < argc) {
                    param_file = argv[++i];
                }
            }
        }

        // 打印版本信息
        std::cout << "+---------------------------------------------------------------+\n"
                  << "|TERNIFY: Efficient Sampling of PROTAC-Induced Ternary Complexes|\n"
                  << "|                FOR ACADEMIC USAGE  ONLY                       |\n"
                  << "|---------------------------------------------------------------|\n"
                  << "|Hongtao Zhao, PhD                                              |\n"
                  << "|Ximing XU, PhD [C++] xuximing@ouc.edu.cn                       |\n" 
                  << "|Version: 2025.06.02                                            |\n" 
                  << "+---------------------------------------------------------------+\n"
                  << std::endl;

        // 读取参数并运行
        Parameters params;
        std::cout << "Reading parameters from " << param_file << std::endl;
        read_parameters(param_file, params);
        run_ternify(params);

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}